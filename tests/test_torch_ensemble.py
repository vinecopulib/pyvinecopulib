"""Many vines, one cascade: the stacked ensemble is a reordering.

`BatchedVineEnsemble` stacks M vines' pair copulas so a tree level is one
call over all of them rather than M. Every pair sees exactly the values it
would have seen alone, and no term is regrouped, so the reference is the
per-vine *batched* cascade and the gate is a last bit rather than a
tolerance in any meaningful sense: `_LAST_BIT` is 1e-14 relative, about
forty-five times the double-precision epsilon and two orders tighter than
the batched-vs-unbatched gates below.

It is not exactly zero, and the reason is worth recording. On x86-64 it
*is* zero -- measured on cpu and cuda over 120 shape configurations, and
that is what pinned the vine-major slot ordering, since edge-major
ordering differs in 4 of those on cpu and 10 on cuda. On arm64 macOS and
on Windows a handful of elements differ by one unit in the last place.
That is not the density product's reduction order: `rosenblatt` carries no
reduction at all on the cached path -- gathers, `lerp`, `addcmul`, a
division -- and it moves too. What differs is the *leading extent* of the
tensors those elementwise kernels run on, `M * (d - t - 1)` rows against
`d - t - 1`, and a kernel is free to pick a different vectorized path,
with different fused-multiply-add contraction, for a different size.
Batching more rows into one call is the whole point of the class, so this
is a property of the approach and not of a mistake in it.

The comparison against the *unbatched* cascade is a genuine tolerance, for
two reasons that predate this class: `_pdf` accumulates one product across
all edges of all trees where the batched path takes one product per level,
and `TorchBicop._prep` trims each pair's inputs where the batched path
clamps them. Those are the published 1e-12 / 1e-13, checked here too so a
regression cannot hide behind the tighter gate.
"""

import pickle
from typing import Any

import numpy as np
import pytest

torch = pytest.importorskip("torch")

import pyvinecopulib as pv  # noqa: E402
from pyvinecopulib.core import (  # noqa: E402
  NonSimplifiedContext,
  VinecopLike,
)
from pyvinecopulib.torch import (  # noqa: E402
  BatchedVineEnsemble,
  FitControlsTorchBicop,
  FitControlsTorchVinecop,
  TorchVinecop,
)
from pyvinecopulib.torch._batched import BatchedTreeLevel  # noqa: E402

from .helpers import assert_on_device  # noqa: E402

_FAMILIES = [pv.families.tll, pv.families.indep]

#: Cross-path gate: the stacked cascade regroups nothing, so anything
#: above a last bit is a wiring bug rather than arithmetic. See the module
#: docstring for why this is not exactly zero off x86-64.
_LAST_BIT = {"rtol": 1e-14, "atol": 1e-16}


def _simulate(d: int, n: int, seed: int = 0) -> np.ndarray:
  rng = np.random.default_rng(seed)
  base = rng.standard_normal(size=(n, 1))
  noise = rng.standard_normal(size=(n, d))
  return pv.to_pseudo_obs(0.6 * base + 0.4 * noise)


def _eval_grid(n: int, d: int, seed: int = 0) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return rng.uniform(0.02, 0.98, size=(n, d))


def _fit_vines(
  m: int,
  d: int,
  n: int = 400,
  seed: int = 0,
  trunc_lvl: int | None = None,
  **kwargs,
) -> list[TorchVinecop]:
  """M vines on one dataset, each on a different sampled structure.

  Different structures are the point: they give every vine its own
  ``order`` and its own wiring, which is what the stacked levels have to
  keep straight.
  """
  u = _simulate(d, n, seed)
  controls = (
    pv.FitControlsVinecop(family_set=_FAMILIES, num_threads=1)
    if trunc_lvl is None
    else pv.FitControlsVinecop(
      family_set=_FAMILIES, num_threads=1, trunc_lvl=trunc_lvl
    )
  )
  return [
    TorchVinecop.from_vinecop(
      pv.Vinecop.from_data(
        u,
        structure=pv.RVineStructure.sample(d, seeds=[k + 1]),
        controls=controls,
      ),
      **kwargs,
    )
    for k in range(m)
  ]


def _torch_vine(d: int, seed: int = 3, k: int = 1, **control_kwargs):
  """One vine fitted natively in torch, for the mismatch cases."""
  u = torch.from_numpy(_simulate(d, 400, seed))
  return TorchVinecop.from_data(
    u,
    pv.RVineStructure.sample(d, seeds=[k]),
    controls=FitControlsTorchVinecop(**control_kwargs),
  )


def _loop(vines, method: str, u):
  return torch.stack([getattr(v, method)(u, batched=True) for v in vines], 0)


# --- parity against the per-vine loop ---------------------------------- #


@pytest.mark.parametrize("method", ("pdf", "rosenblatt"))
@pytest.mark.parametrize("m", (1, 2, 7))
@pytest.mark.parametrize("n", (1, 2, 37, 300))
def test_matches_the_per_vine_loop(method: str, m: int, n: int) -> None:
  """Agrees with the loop to a last bit, at every shape including n = 1.

  `n = 1` is not padding. The density's per-level product is the one
  operation whose shape changes, from `(N_t, n)` reduced over dim 0 to
  `(M, N_t, n)` reduced over dim 1, and a trailing axis of length one is
  where the two could pick different reduction geometries. Vine-major
  slots keep them the same; edge-major ordering was measured to differ at
  n = 1 and nowhere else, which is what pinned the layout. Windows then
  found a 1-ULP discrepancy at n = 1 anyway, from a different cause -- see
  the module docstring -- so the shape stays in the sweep.
  """
  vines = _fit_vines(m, d=6, seed=11)
  ens = BatchedVineEnsemble(vines)
  u = torch.from_numpy(_eval_grid(n, 6, seed=99))
  got = getattr(ens, method)(u)
  assert got.shape == ((m, n) if method == "pdf" else (m, n, 6))
  torch.testing.assert_close(got, _loop(vines, method, u), **_LAST_BIT)


@pytest.mark.parametrize(
  ("method", "tol"), (("pdf", 1e-12), ("rosenblatt", 1e-13))
)
def test_matches_the_unbatched_cascade(method: str, tol: float) -> None:
  """The published batched-vs-unbatched gates carry over to the stack."""
  vines = _fit_vines(4, d=6, seed=12)
  ens = BatchedVineEnsemble(vines)
  u = torch.from_numpy(_eval_grid(200, 6, seed=98))
  ref = torch.stack([getattr(v, method)(u, batched=False) for v in vines], 0)
  torch.testing.assert_close(getattr(ens, method)(u), ref, atol=tol, rtol=tol)


@pytest.mark.parametrize("chunk", (1, 3, 5, 20))
def test_chunking_does_not_change_the_result(chunk: int) -> None:
  """A chunk is a scheduling unit, including the ragged last one."""
  vines = _fit_vines(11, d=5, seed=13)
  ens = BatchedVineEnsemble(vines, chunk_size=chunk)
  u = torch.from_numpy(_eval_grid(64, 5, seed=97))
  assert len(ens.chunks) == -(-11 // min(chunk, 11))
  for method in ("pdf", "rosenblatt"):
    torch.testing.assert_close(
      getattr(ens, method)(u), _loop(vines, method, u), **_LAST_BIT
    )


@pytest.mark.parametrize("d", (2, 3))
def test_the_smallest_vines(d: int) -> None:
  """d = 2 is one pair and one tree level -- every loop runs once.

  The degenerate end of the cascade: a level whose stacked pair count is
  the vine count itself, a density product over a single edge, and for
  d = 2 an ``out_idx`` that is a permutation of two rows. Nothing here is
  special-cased, so this is guarding that nothing needs to be.
  """
  vines = _fit_vines(3, d=d, seed=30 + d)
  ens = BatchedVineEnsemble(vines)
  assert ens.trunc_lvl == d - 1
  u = torch.from_numpy(_eval_grid(17, d, seed=87))
  for method in ("pdf", "rosenblatt"):
    torch.testing.assert_close(
      getattr(ens, method)(u), _loop(vines, method, u), **_LAST_BIT
    )


def test_a_uniformly_truncated_set() -> None:
  """Truncation shortens the cascade; it does not change its shape."""
  vines = _fit_vines(4, d=7, seed=14, trunc_lvl=2)
  assert [v.trunc_lvl for v in vines] == [2] * 4
  ens = BatchedVineEnsemble(vines)
  assert ens.trunc_lvl == 2
  u = torch.from_numpy(_eval_grid(50, 7, seed=96))
  for method in ("pdf", "rosenblatt"):
    torch.testing.assert_close(
      getattr(ens, method)(u), _loop(vines, method, u), **_LAST_BIT
    )


def test_a_set_mixing_indep_and_tll_pairs() -> None:
  """An independence pair has a 2x2 sentinel grid; the level still stacks.

  The substitute density the bake puts in its place has to be the same one
  every vine in the set gets, or two vines' independence pairs would not
  agree with each other.
  """
  rng = np.random.RandomState(83)
  base = rng.normal(size=800)
  x = np.column_stack(
    [
      0.95 * base + 0.05 * rng.normal(size=800),
      0.95 * base + 0.05 * rng.normal(size=800),
      rng.normal(size=800),
      rng.normal(size=800),
      0.9 * base + 0.1 * rng.normal(size=800),
    ]
  )
  u = pv.to_pseudo_obs(x)
  controls = pv.FitControlsVinecop(family_set=_FAMILIES, num_threads=1)
  vines, families = [], set()
  for k in range(4):
    cop = pv.Vinecop.from_data(
      u, structure=pv.RVineStructure.sample(5, seeds=[k + 1]), controls=controls
    )
    families |= {
      cop.get_pair_copula(t, e).family
      for t in range(cop.trunc_lvl)
      for e in range(cop.dim - 1 - t)
    }
    vines.append(TorchVinecop.from_vinecop(cop))
  assert {pv.BicopFamily.tll, pv.BicopFamily.indep} <= families

  ens = BatchedVineEnsemble(vines)
  u_t = torch.from_numpy(u[:64])
  for method in ("pdf", "rosenblatt"):
    torch.testing.assert_close(
      getattr(ens, method)(u_t), _loop(vines, method, u_t), **_LAST_BIT
    )


# --- what it refuses --------------------------------------------------- #


def test_an_empty_set_is_refused() -> None:
  with pytest.raises(ValueError, match="at least one vine"):
    BatchedVineEnsemble([])


def test_a_non_vine_is_refused() -> None:
  """The check earns its keep on callers the type checker never sees."""
  # Annotated `list[Any]` because the point is to pass what the signature
  # forbids, and a type checker is right to object to the literal form.
  not_vines: list[Any] = [_torch_vine(5), 3]
  with pytest.raises(ValueError, match=r"vines\[1\] has type int"):
    BatchedVineEnsemble(not_vines)


def test_a_mismatched_dimension_is_refused() -> None:
  with pytest.raises(ValueError, match=r"vines\[1\] has dimension 6"):
    BatchedVineEnsemble([_torch_vine(5), _torch_vine(6)])


def test_a_mismatched_truncation_is_refused() -> None:
  """The message has to name the workaround: group by `trunc_lvl`."""
  u = _simulate(6, 400, 15)
  vines = [
    TorchVinecop.from_vinecop(
      pv.Vinecop.from_data(
        u,
        controls=pv.FitControlsVinecop(
          family_set=_FAMILIES, num_threads=1, trunc_lvl=tl
        ),
      )
    )
    for tl in (5, 2)
  ]
  with pytest.raises(ValueError, match="group the set by trunc_lvl"):
    BatchedVineEnsemble(vines)


def test_a_mismatched_grid_size_is_refused() -> None:
  with pytest.raises(ValueError, match="interpolation grid"):
    BatchedVineEnsemble(
      [
        _torch_vine(5, k=1),
        _torch_vine(5, k=2, bicop_controls=FitControlsTorchBicop(grid_size=20)),
      ]
    )


def test_a_mismatched_grid_type_is_refused() -> None:
  """The one a shape comparison cannot catch.

  A normal and a linear grid of the same `grid_size` have the same shape
  and different coordinates, so the single-vine precondition -- which
  compares shapes, every pair there having come out of one fit -- passes
  them. Stacked, they would be silently evaluated against one of the two.
  """
  with pytest.raises(ValueError, match="interpolation grid"):
    BatchedVineEnsemble(
      [
        _torch_vine(5, k=1),
        _torch_vine(
          5, k=2, bicop_controls=FitControlsTorchBicop(grid_type="linear")
        ),
      ]
    )


def test_a_mismatched_cache_setting_is_refused() -> None:
  """Degrading the set to the uncached path would break its own parity.

  The cached and on-the-fly integrals agree to rounding, not exactly, so a
  level that dropped to quadrature because *one* vine lacked tables would
  stop being bit-identical for all the others.
  """
  with pytest.raises(ValueError, match="cached integrals at tree level"):
    BatchedVineEnsemble(
      [_torch_vine(5, k=1), _torch_vine(5, k=2, cache_integrals=False)]
    )


def test_a_mismatched_dtype_is_refused() -> None:
  with pytest.raises(ValueError, match="float32"):
    BatchedVineEnsemble(
      [_torch_vine(5, k=1), _torch_vine(5, k=2, dtype=torch.float32)]
    )


def test_a_discrete_vine_is_refused() -> None:
  u = _simulate(3, 300, 9)
  u[:, 0] = np.ceil(u[:, 0] * 4) / 4
  wide = np.column_stack([u, np.maximum(u[:, 0] - 0.25, 0.0)])
  vine = TorchVinecop.from_data(
    torch.from_numpy(wide),
    pv.RVineStructure.sample(3, seeds=[1]),
    var_types=["d", "c", "c"],
  )
  with pytest.raises(ValueError, match="discrete variables"):
    BatchedVineEnsemble([vine])


def test_a_non_simplified_vine_is_refused() -> None:
  vine = _torch_vine(5)
  conditional = TorchVinecop(
    vine.pair_copulas, vine.structure, context=NonSimplifiedContext()
  )
  with pytest.raises(ValueError, match="non-simplified"):
    BatchedVineEnsemble([conditional])


def test_a_grad_tracking_grid_is_refused() -> None:
  """A snapshot cannot see a later `requires_grad_`, so it declines to try."""
  vine = _torch_vine(5)
  vine._pair_module(0, 0).interp_grid.values.requires_grad_(True)
  with pytest.raises(ValueError, match="tracks gradients"):
    BatchedVineEnsemble([vine])


def test_a_bad_chunk_size_is_refused() -> None:
  with pytest.raises(ValueError, match="chunk_size must be >= 1"):
    BatchedVineEnsemble(_fit_vines(2, d=5, seed=16), chunk_size=0)


def test_a_wrong_width_input_is_refused() -> None:
  ens = BatchedVineEnsemble(_fit_vines(2, d=5, seed=17))
  with pytest.raises(ValueError, match=r"\(n, 5\)"):
    ens.pdf(torch.rand(10, 6, dtype=torch.float64))


# --- torch plumbing ---------------------------------------------------- #


def test_it_is_not_a_vinecoplike() -> None:
  """Deliberately: its `pdf` is `(M, n)`, which no vine consumer expects.

  `VinecopLike` is `runtime_checkable`, so a class carrying the method
  names would pass `isinstance` and could be handed to `Vinedist` or the
  sklearn estimators, which would read the stacked density as one vine's.
  """
  ens = BatchedVineEnsemble(_fit_vines(2, d=5, seed=18))
  assert not isinstance(ens, VinecopLike)


def test_the_stacked_state_round_trips_through_state_dict() -> None:
  vines = _fit_vines(3, d=5, seed=19)
  ens = BatchedVineEnsemble(vines)
  u = torch.from_numpy(_eval_grid(40, 5, seed=95))
  want = ens.pdf(u)

  fresh = BatchedVineEnsemble(vines)
  fresh.load_state_dict(ens.state_dict(), strict=True)
  torch.testing.assert_close(fresh.pdf(u), want, atol=0.0, rtol=0.0)


def test_evaluation_leaves_the_state_dict_alone() -> None:
  """No lazily-derived cache leaks in, as it would on a vine."""
  ens = BatchedVineEnsemble(_fit_vines(2, d=5, seed=20))
  before = set(ens.state_dict())
  ens.pdf(torch.from_numpy(_eval_grid(16, 5, seed=94)))
  ens.rosenblatt(torch.from_numpy(_eval_grid(16, 5, seed=94)))
  assert set(ens.state_dict()) == before


def test_it_pickles_after_an_evaluation() -> None:
  ens = BatchedVineEnsemble(_fit_vines(3, d=5, seed=21))
  u = torch.from_numpy(_eval_grid(32, 5, seed=93))
  want = ens.pdf(u)
  again = pickle.loads(pickle.dumps(ens))
  torch.testing.assert_close(again.pdf(u), want, atol=0.0, rtol=0.0)


def test_it_pickles_after_nn_module_compile() -> None:
  """`nn.Module.compile()` installs a callable that has to be dropped.

  Distinct from :func:`test_it_pickles_after_an_evaluation`: that one
  covers the ensemble's own ``_compiled`` cache, while this covers
  ``_compiled_call_impl``, which torch installs on the module itself.
  Copying ``__dict__`` instead of delegating to ``nn.Module.__getstate__``
  keeps it, and pickling then raises ``PicklingError``.
  """
  vines = _fit_vines(2, d=5, seed=27)
  ens = BatchedVineEnsemble(vines)
  u = torch.from_numpy(_eval_grid(20, 5, seed=88))
  want = ens.pdf(u)
  ens.compile()
  again = pickle.loads(pickle.dumps(ens))
  torch.testing.assert_close(again.pdf(u), want, atol=0.0, rtol=0.0)


def test_every_buffer_follows_to_device(device: str) -> None:
  """No `extra=`: unlike a vine's bake, the stacked state is registered."""
  ens = BatchedVineEnsemble(_fit_vines(3, d=5, seed=22)).to(device)
  u = torch.from_numpy(_eval_grid(32, 5, seed=92)).to(device)
  out = (ens.pdf(u), ens.rosenblatt(u))
  assert_on_device(ens, device, out)


def test_a_device_move_drops_the_compiled_cache() -> None:
  ens = BatchedVineEnsemble(_fit_vines(2, d=5, seed=23))
  ens.compile_cascades = True
  ens._compiled["pdf"] = object()
  ens.to(torch.float64)
  assert ens._compiled == {}


def test_gradients_flow_to_the_input() -> None:
  vines = _fit_vines(3, d=5, seed=24)
  ens = BatchedVineEnsemble(vines)
  u = torch.from_numpy(_eval_grid(48, 5, seed=91)).requires_grad_(True)
  (grad,) = torch.autograd.grad(ens.pdf(u).sum(), u)
  assert torch.isfinite(grad).all() and grad.abs().sum() > 0

  u_ref = u.detach().clone().requires_grad_(True)
  (want,) = torch.autograd.grad(
    torch.stack([v.pdf(u_ref, batched=True) for v in vines], 0).sum(), u_ref
  )
  torch.testing.assert_close(grad, want, rtol=1e-10, atol=1e-12)


def test_the_number_of_stacked_calls_does_not_scale_with_m(
  monkeypatch: pytest.MonkeyPatch,
) -> None:
  """The whole point: one call per tree level, whatever M is.

  Counted at the stacked-level call rather than in kernel launches, so it
  is a deterministic assertion instead of a profiler reading.
  """
  calls = []
  original = BatchedTreeLevel.pdf_h1_h2

  def counted(self, *args, **kwargs):
    calls.append(1)
    return original(self, *args, **kwargs)

  monkeypatch.setattr(BatchedTreeLevel, "pdf_h1_h2", counted)
  u = torch.from_numpy(_eval_grid(32, 6, seed=90))
  counts = []
  for m in (1, 8):
    ens = BatchedVineEnsemble(_fit_vines(m, d=6, seed=25), chunk_size=m)
    calls.clear()
    ens.pdf(u)
    counts.append(len(calls))
  assert counts[0] == counts[1] == 5  # d - 1 tree levels


@pytest.mark.compile
@pytest.mark.parametrize("method", ("pdf", "rosenblatt"))
def test_compiled_matches_eager(method: str) -> None:
  """Inductor fuses the cascade; it does not change what it computes."""
  vines = _fit_vines(4, d=5, seed=26)
  ens = BatchedVineEnsemble(vines)
  u = torch.from_numpy(_eval_grid(64, 5, seed=89))
  eager = getattr(ens, method)(u)
  ens.compile_cascades = True
  torch.testing.assert_close(
    getattr(ens, method)(u), eager, atol=1e-12, rtol=1e-12
  )
