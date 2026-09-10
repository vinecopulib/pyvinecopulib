"""The torch cascade on CUDA: same numbers, same placement, no host detours.

The C++-parity suite in ``test_torch_tll_bicop`` / ``test_torch_vinecop`` pins
torch against the compiled library on the cpu, in float64. This file pins
cuda against *that*, so the compiled-library agreement carries over to the
device without running the whole parity suite twice.

``float32`` is compared only to our own ``float64`` answer, never to the
compiled library: it is a lower precision, so the float64 tolerances do not
apply to it and loosening them would weaken the check that does.
"""

from typing import Any

import numpy as np
import pytest

torch = pytest.importorskip("torch")

import pyvinecopulib as pv  # noqa: E402
from pyvinecopulib.torch import (  # noqa: E402
  FitControlsTorchVinecop,
  TorchTllBicop,
  TorchKde1d,
  TorchVinecop,
  TorchVinedist,
)
from tests.helpers import assert_on_device, count_transfers  # noqa: E402

# cuda vs cpu at the same dtype. Tight enough that the cpu-vs-C++ 1e-10
# tolerances carry over to cuda by the triangle inequality.
DEVICE_TOL = {torch.float64: 1e-12, torch.float32: 2e-5}
# float32 against our own float64 result.
F32_VS_F64_TOL = 5e-5

_EVAL_OPS = ("pdf", "rosenblatt", "inverse_rosenblatt")


def _u(d: int, n: int, seed: int) -> np.ndarray:
  rng = np.random.default_rng(seed)
  a = rng.standard_normal((d, d))
  s = a @ a.T + d * np.eye(d)
  sd = np.sqrt(np.diag(s))
  return pv.to_pseudo_obs(
    rng.multivariate_normal(np.zeros(d), s / np.outer(sd, sd), size=n)
  )


@pytest.fixture
def cpp_vine() -> pv.Vinecop:
  """A fitted reference vine.

  Function-scoped like the rest of the suite: a module-scoped compiled
  object outlives pytest's own teardown and nanobind reports it as a leak
  at interpreter shutdown.
  """
  return pv.Vinecop.from_data(
    _u(5, 600, 7),
    controls=pv.FitControlsVinecop(family_set=[pv.BicopFamily.tll]),
  )


@pytest.fixture
def u_eval() -> np.ndarray:
  return _u(5, 400, 8)


def _np(t: "torch.Tensor") -> np.ndarray:
  return t.detach().cpu().numpy()


@pytest.mark.parametrize("op", _EVAL_OPS)
@pytest.mark.parametrize("cache", [False, True])
def test_eval_matches_cpu(
  device: str, op: str, cache: bool, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """Every evaluation op returns the cpu answer, on the device."""
  ref = TorchVinecop.from_vinecop(
    cpp_vine, cache_integrals=cache, device=torch.device("cpu")
  )
  got = TorchVinecop.from_vinecop(
    cpp_vine, cache_integrals=cache, device=torch.device(device)
  )
  a = getattr(ref, op)(torch.as_tensor(u_eval))
  b = getattr(got, op)(torch.as_tensor(u_eval, device=device))
  assert b.device.type == torch.device(device).type
  np.testing.assert_allclose(
    _np(b),
    _np(a),
    rtol=DEVICE_TOL[torch.float64],
    atol=DEVICE_TOL[torch.float64],
  )


@pytest.mark.parametrize(
  ("op", "tol"),
  [("pdf", 1e-12), ("rosenblatt", 1e-12), ("inverse_rosenblatt", 0.0)],
)
def test_batched_matches_unbatched_on_device(
  device: str, op: str, tol: float, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """The batched fast path agrees with the per-edge cascade on the device.

  ``inverse_rosenblatt`` is held to equality rather than a tolerance: its
  batched path groups cells of the dependency graph without changing what any
  one of them computes, so it is a reordering and not a second implementation.
  """
  vine = TorchVinecop.from_vinecop(cpp_vine, device=torch.device(device))
  ut = torch.as_tensor(u_eval, device=device)
  np.testing.assert_allclose(
    _np(getattr(vine, op)(ut, batched=True)),
    _np(getattr(vine, op)(ut, batched=False)),
    rtol=tol,
    atol=tol / 10.0,
  )


def test_batched_sample_matches_unbatched_on_device(
  device: str, cpp_vine: pv.Vinecop
) -> None:
  """Seeded draws are identical either way, since the inverse is a reordering."""
  vine = TorchVinecop.from_vinecop(cpp_vine, device=torch.device(device))
  np.testing.assert_array_equal(
    _np(vine.sample(300, seeds=[5, 6, 7], batched=True)),
    _np(vine.sample(300, seeds=[5, 6, 7], batched=False)),
  )


def test_float32_tracks_float64(
  device: str, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """float32 reproduces the float64 answer to float32 precision.

  Pinned against our own float64 result rather than the compiled library:
  the C++ tolerances describe float64 and do not transfer.
  """
  f64 = TorchVinecop.from_vinecop(cpp_vine, device=torch.device(device))
  f32 = TorchVinecop.from_vinecop(
    cpp_vine, device=torch.device(device), dtype=torch.float32
  )
  a = _np(f64.pdf(torch.as_tensor(u_eval, device=device)))
  b = _np(f32.pdf(torch.as_tensor(u_eval, device=device, dtype=torch.float32)))
  assert np.isfinite(b).all()
  np.testing.assert_allclose(b, a, rtol=F32_VS_F64_TOL, atol=F32_VS_F64_TOL)


def test_float32_keeps_arguments_inside_the_unit_square(device: str) -> None:
  """The trim survives float32.

  ``1 - 1e-10`` rounds to exactly ``1.0`` in float32, so a bound written as
  that literal would admit the value it exists to exclude -- and ``ndtri(1)``
  is an infinity for whatever refits on the result.
  """
  u = pv.to_pseudo_obs(
    np.random.default_rng(3).multivariate_normal(
      [0, 0], [[1, 0.9], [0.9, 1]], size=1500
    )
  )
  cpp = pv.Bicop.from_data(
    u, controls=pv.FitControlsBicop(family_set=[pv.BicopFamily.tll])
  )
  bc = TorchTllBicop.from_bicop(
    cpp, device=torch.device(device), dtype=torch.float32
  )
  edge = torch.tensor(
    [[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0], [1.0, 1.0], [0.0, 0.0]],
    dtype=torch.float32,
    device=device,
  )
  for name in ("hfunc1", "hfunc2", "cdf"):
    out = getattr(bc, name)(edge)
    assert (out > 0.0).all(), f"{name} returned 0"
    assert (out < 1.0).all(), f"{name} returned 1"
    assert torch.isfinite(torch.special.ndtri(out)).all()


def test_float16_keeps_arguments_inside_the_unit_square(
  device: str, cpp_vine: pv.Vinecop
) -> None:
  """The lower trim bound remains representable in half precision."""
  bc = TorchTllBicop.from_bicop(
    cpp_vine.get_pair_copula(0, 0),
    device=torch.device(device),
    dtype=torch.float16,
  )
  edge = torch.tensor(
    [[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0]],
    dtype=torch.float16,
    device=device,
  )
  for name in ("hfunc1", "hfunc2", "cdf"):
    out = getattr(bc, name)(edge)
    assert (out > 0.0).all(), f"{name} returned 0"
    assert (out < 1.0).all(), f"{name} returned 1"


def test_every_buffer_follows_to_device(
  device: str, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """``.to(device)`` moves the whole object, batched cache included."""
  vine = TorchVinecop.from_vinecop(cpp_vine)
  vine.to(device)
  ut = torch.as_tensor(u_eval, device=device)
  out = vine.pdf(ut, batched=True)
  assert_on_device(vine, device, out, extra=(vine._batched,))


@pytest.mark.parametrize("op", _EVAL_OPS + ("sample",))
def test_evaluation_does_not_round_trip_through_the_host(
  device: str, op: str, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """An evaluation call must not move data back to the host.

  ``fit`` / ``select`` legitimately do -- Kendall's tau goes through the
  compiled ``wdm`` -- which is why this is scoped to evaluation. The first
  call is untimed: it builds the batched cache, which reads the structure
  from the compiled extension.
  """
  vine = TorchVinecop.from_vinecop(cpp_vine, device=torch.device(device))
  ut = torch.as_tensor(u_eval, device=device)
  call = (
    (lambda: vine.sample(64, seeds=[1]))
    if op == "sample"
    else (lambda: getattr(vine, op)(ut))
  )
  call()
  with count_transfers(device) as c:
    call()
  c.assert_no_d2h(f"TorchVinecop.{op}")


def test_batched_fit_peak_memory_stays_bounded(device: str) -> None:
  """A wide level on a long sample does not scale the footprint with either.

  The batched fit's peak lives in the kernel evaluation's temporaries, which
  a fixed grid block grows with both `n` and the level width: this vine held
  about 1.8 GiB that way, on a card most users have 8 of. Sizing the block
  from the two instead holds the peak near the budget, which is what makes a
  `d = 20` fit on real-sized data something a laptop can run at all.
  """
  if torch.device(device).type != "cuda":
    pytest.skip("peak allocation is only observable on cuda")
  from pyvinecopulib.torch._bicop_fit_tll import _KDE_MEM_BUDGET_BYTES

  d, n = 20, 8000
  u_np = _u(d, n, 7)
  structure = pv.Vinecop.from_data(
    u_np,
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.tll], num_threads=1, trunc_lvl=20
    ),
  ).structure
  u = torch.as_tensor(u_np, device=device)
  controls = FitControlsTorchVinecop(
    device=torch.device(device), batched_fit=True
  )
  TorchVinecop.from_data(u, structure=structure, controls=controls)  # warm
  torch.cuda.empty_cache()
  torch.cuda.reset_peak_memory_stats()
  TorchVinecop.from_data(u, structure=structure, controls=controls)
  peak = torch.cuda.max_memory_allocated()
  # Generous against the budget, the data and the fitted pairs sitting
  # outside it, and far under the ~1.8 GiB a fixed block took.
  assert peak < 3 * _KDE_MEM_BUDGET_BYTES, f"peak {peak / 2**20:.0f} MiB"


def test_fit_and_select_run_on_device(device: str) -> None:
  """Fitting and structure selection work with device-resident data."""
  u = _u(4, 600, 11)
  ut = torch.as_tensor(u, device=device)
  ctl = FitControlsTorchVinecop(device=torch.device(device))
  fixed = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.BicopFamily.tll])
  )
  vine = TorchVinecop.from_data(ut, structure=fixed.structure, controls=ctl)
  assert vine.pdf(ut).device.type == torch.device(device).type
  selected = TorchVinecop.from_data(ut, controls=ctl)
  assert selected.pdf(ut).device.type == torch.device(device).type


@pytest.mark.parametrize("threshold", [0.3, 0.5])
def test_thresholded_pairs_land_where_the_data_is(
  device: str, threshold: float
) -> None:
  """A thresholded edge follows the data, not `controls.device`.

  A thresholded edge is not fitted, so its pair is constructed rather than
  derived from `u`, and it has to be placed explicitly. Taking the
  placement from `controls.device` gets it wrong whenever the caller left
  that `None` and let the data choose -- which is the documented way to pass
  an already-resident tensor.

  Both thresholds are here because the failure has two faces. At 0.3 the
  level is mixed, so fitted pairs sit on the data's device and thresholded
  ones on the cpu, and evaluation raises. At 0.5 every pair on this fixture
  is thresholded, nothing raises, and the whole vine quietly lands on the
  cpu with its data on the accelerator.
  """
  d, n = 6, 400
  u = _u(d, n, 7)
  ut = torch.as_tensor(u, device=device)
  # `controls.device` left None: the data carries the placement.
  vine = TorchVinecop.from_data(
    ut, controls=FitControlsTorchVinecop(trunc_lvl=20, threshold=threshold)
  )
  want = torch.device(device).type
  rows: Any = vine.pair_copulas
  for t in range(vine.trunc_lvl):
    for e in range(d - t - 1):
      got = rows[t][e].interp_grid.values.device.type
      assert got == want, f"pair ({t}, {e}) on {got}, expected {want}"
  assert vine.pdf(ut).device.type == want


@pytest.mark.parametrize("var_types", [["d", "c", "c"], ["c", "d", "d"]])
def test_discrete_vine_matches_cpu(device: str, var_types: list[str]) -> None:
  """A vine with atoms evaluates identically on either device.

  The discrete cascade differences the distribution function over an atom's
  width, which amplifies any absolute error by ``~4/(w1 w2)`` -- so this is
  the path where a device-dependent rounding difference would show first.
  """
  rng = np.random.default_rng(19)
  d = len(var_types)
  x = rng.multivariate_normal(np.zeros(d), np.eye(d) * 0.5 + 0.5, size=800)
  cols, lims = [], []
  for j, t in enumerate(var_types):
    if t == "d":
      k = np.floor(4.0 * pv.to_pseudo_obs(x[:, [j]])[:, 0])
      cols.append((k + 1.0) / 4.0)
      lims.append(k / 4.0)
    else:
      cols.append(pv.to_pseudo_obs(x[:, [j]])[:, 0])
  u = np.column_stack(cols + lims)
  cpp = pv.Vinecop.from_data(
    u,
    var_types=var_types,
    controls=pv.FitControlsVinecop(family_set=[pv.BicopFamily.tll]),
  )
  ref = TorchVinecop.from_vinecop(cpp, device=torch.device("cpu"))
  got = TorchVinecop.from_vinecop(cpp, device=torch.device(device))
  a = _np(ref.pdf(torch.as_tensor(u)))
  b = _np(got.pdf(torch.as_tensor(u, device=device)))
  assert b.shape == a.shape
  np.testing.assert_allclose(b, a, rtol=1e-12, atol=1e-12)


def test_vinedist_matches_cpu(device: str) -> None:
  """The data-scale distribution agrees across devices, margins included."""
  rng = np.random.default_rng(23)
  d = 3
  y = rng.multivariate_normal(np.zeros(d), np.eye(d) * 0.4 + 0.6, size=600)
  cpp = pv.Vinecop.from_data(
    pv.to_pseudo_obs(y),
    controls=pv.FitControlsVinecop(family_set=[pv.BicopFamily.tll]),
  )

  def build(dev: str) -> TorchVinedist:
    cop = TorchVinecop.from_vinecop(cpp, device=torch.device(dev))
    margins = [
      TorchKde1d.from_kde1d(pv.core.Kde1d().fit(y[:, j])).to(dev)
      for j in range(d)
    ]
    return TorchVinedist(cop, margins)

  ref, got = build("cpu"), build(device)
  yt = torch.as_tensor(y)
  a = _np(ref.logpdf(yt))
  b = _np(got.logpdf(torch.as_tensor(y, device=device)))
  assert np.isfinite(b).all()
  np.testing.assert_allclose(b, a, rtol=1e-11, atol=1e-11)
  assert_on_device(got, device)


@pytest.mark.compile
def test_compiled_output_survives_the_next_call(
  device: str, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """A returned density is the caller's, not a buffer the next call reuses.

  On CUDA the compiled cascade replays as one graph, whose result lands in a
  buffer the next replay overwrites -- so holding two results, as dividing one
  density by another does, has to keep working.
  """
  vine = TorchVinecop.from_vinecop(cpp_vine, device=torch.device(device))
  a = torch.as_tensor(u_eval, device=device)
  b = torch.as_tensor(_u(5, 400, 9), device=device)
  eager_a, eager_b = vine.pdf(a).clone(), vine.pdf(b).clone()
  vine.compile_cascades = True
  held = vine.pdf(a)
  later = vine.pdf(b)
  np.testing.assert_allclose(_np(held), _np(eager_a), rtol=1e-12, atol=1e-13)
  np.testing.assert_allclose(_np(later), _np(eager_b), rtol=1e-12, atol=1e-13)


@pytest.mark.parametrize("d", [5, 9])
def test_batched_fit_matches_the_per_edge_fit_on_device(
  device: str, d: int
) -> None:
  """`batched_fit` is a schedule, not a model, on either device.

  One tolerance for both devices, because the mechanism is not a device
  asymmetry: stacking a level changes how many elements the bandwidth
  search's `pow` sees, and torch selects an elementwise kernel by element
  count. A cpu that vectorizes at a lower lane count than another therefore
  diverges where the other does not, which is why this cannot be pinned at
  zero on the strength of one machine agreeing.
  """
  u_fit = _u(d, 1200, 600 + d)
  structure = pv.Vinecop.from_data(
    u_fit,
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.tll], num_threads=1, trunc_lvl=20
    ),
  ).structure
  u_t = torch.as_tensor(u_fit, device=device)
  fits = {
    flag: TorchVinecop.from_data(
      u_t,
      structure=structure,
      controls=FitControlsTorchVinecop(
        device=torch.device(device), batched_fit=flag
      ),
    )
    for flag in (False, True)
  }
  u_eval = torch.as_tensor(_u(d, 300, 610 + d), device=device)
  np.testing.assert_allclose(
    _np(fits[True].pdf(u_eval)),
    _np(fits[False].pdf(u_eval)),
    rtol=1e-9,
    atol=1e-11,
  )


def test_a_declared_placement_serves_a_host_that_is_not_a_module(
  device: str,
) -> None:
  """The mixin resolves a declaration, not only registered tensors.

  A pair copula that is not an ``nn.Module`` -- backend
  estimators, a device handle and Python scalars, no tensor -- registers
  nothing for ``reference_tensor`` to find. Before this it reached
  ``self.parameters()`` and raised ``AttributeError``; the array-API inference
  it falls back to instead returned its argument untouched, which is how a
  host ``x`` met a device ``u`` inside a ``column_stack`` downstream.
  """
  from pyvinecopulib.core import BicopBase
  from pyvinecopulib.torch import TensorPlacementMixin

  class _Declared(TensorPlacementMixin, BicopBase[torch.Tensor]):
    supports_covariates = True

    def __init__(self) -> None:
      self.device = torch.device(device)
      self.dtype = torch.float64

    def pdf(self, u: Any, *, x: Any = None) -> Any:
      u = self._prep_args(u)
      assert isinstance(u, torch.Tensor) and u.device.type == device
      assert x is None or (
        isinstance(x, torch.Tensor) and x.device.type == device
      )
      return torch.ones(u.shape[0], dtype=u.dtype, device=u.device)

    def hfunc1(self, u: Any, *, x: Any = None) -> Any:
      return self._prep_args(u)[:, 1]

    def hfunc2(self, u: Any, *, x: Any = None) -> Any:
      return self._prep_args(u)[:, 0]

  pair = _Declared()
  u_np = np.array([[0.3, 0.5], [0.7, 0.2]])

  placed = pair._prep(u_np)
  assert placed.dtype is torch.float64 and placed.device.type == device
  # The dtype is normalized rather than inherited: `torch.tensor([0.3])` is
  # float32, and a downstream test had been pinning that as correct.
  assert pair._prep(u_np.astype(np.float32)).dtype is torch.float64
  # Both arguments arrive placed, so they can meet in one expression.
  assert float(pair.loglik(u_np, x=np.array([[1.0], [2.0]]))) == 0.0


def test_a_member_named_parameters_does_not_decide_the_placement(
  device: str,
) -> None:
  """Being a module selects the registered-tensor step, not the member name.

  A foreign estimator may carry a ``parameters`` of its own -- a dict of
  hyperparameters is the obvious one -- and it says nothing about where the
  numerics run. Reading it as though it were ``nn.Module.parameters`` raised
  ``TypeError`` out of ``_prep`` rather than falling back to the declaration
  the host had made, which is a refusal in place of an answer that was
  available.
  """
  from pyvinecopulib.core import BicopBase
  from pyvinecopulib.torch import TensorPlacementMixin

  class _Coincidental(TensorPlacementMixin, BicopBase[torch.Tensor]):
    # Not callable, and not tensors: exactly what a wrapped estimator's own
    # hyperparameter record looks like.
    parameters = {"n_estimators": 400, "depth": 6}

    def __init__(self) -> None:
      self.device = torch.device(device)
      self.dtype = torch.float64

    def pdf(self, u: Any, *, x: Any = None) -> Any:
      return torch.ones(u.shape[0], dtype=u.dtype, device=u.device)

    def hfunc1(self, u: Any, *, x: Any = None) -> Any:
      return self._prep_args(u)[:, 1]

    def hfunc2(self, u: Any, *, x: Any = None) -> Any:
      return self._prep_args(u)[:, 0]

  placed = _Coincidental()._prep(np.array([[0.3, 0.5]], dtype=np.float32))
  assert placed.dtype is torch.float64 and placed.device.type == device


def test_a_host_with_parameters_but_no_buffers_still_resolves(
  device: str,
) -> None:
  """Each member is read on its own, so a partial one is not fatal.

  ``parameters()`` and ``buffers()`` used to be read as one expression, and
  both were evaluated before either was consumed -- so a host exposing only
  the first raised ``AttributeError`` even when that first one held exactly
  the tensor being looked for.
  """
  from pyvinecopulib.torch import reference_tensor

  ref = torch.zeros(1, dtype=torch.float32, device=device)

  class _HalfModule:
    def parameters(self) -> Any:
      return iter([ref])

  # Read directly: `reference_tensor` is exported, so it is reachable with an
  # object that is no module at all.
  assert reference_tensor(_HalfModule()) is ref

  class _BuffersOnly:
    def buffers(self) -> Any:
      return iter([ref])

  assert reference_tensor(_BuffersOnly()) is ref
  # An object exposing neither member answers `None` rather than raising, and
  # so does one whose member of that name is something else entirely.
  assert reference_tensor(object()) is None

  class _Coincidental:
    parameters = {"depth": 6}

  assert reference_tensor(_Coincidental()) is None

  # An integer tensor is no answer here: every caller has a floating default.
  class _IntOnly:
    def parameters(self) -> Any:
      return iter([torch.zeros(1, dtype=torch.int64, device=device)])

  assert reference_tensor(_IntOnly()) is None


def test_a_declared_host_keeps_a_tracked_gradient(device: str) -> None:
  """The declared fallback goes through ``as_tensor``, which preserves grad.

  ``place`` cannot promise this -- it converts through the namespace's own
  ``asarray``, whose answer for a tracked tensor changed between torch 2.11
  and 2.13 -- so the guarantee belongs to this route and is worth pinning.
  """
  from pyvinecopulib.torch import TensorPlacementMixin

  on = torch.device(device)

  class _Declared(TensorPlacementMixin):
    device = on
    dtype = torch.float64

  tracked = torch.ones(
    2, dtype=torch.float64, device=device, requires_grad=True
  )
  placed = _Declared()._prep(tracked)
  assert placed.requires_grad
  assert placed.dtype is torch.float64 and placed.device.type == device


def test_an_undeclared_host_gets_the_documented_default() -> None:
  """No registered tensor and no declaration is the third step, not a raise."""
  from pyvinecopulib.core import BicopBase
  from pyvinecopulib.torch import TensorPlacementMixin

  class _Bare(TensorPlacementMixin, BicopBase[torch.Tensor]):
    def pdf(self, u: Any, *, x: Any = None) -> Any:
      return torch.ones(u.shape[0], dtype=u.dtype)

    def hfunc1(self, u: Any, *, x: Any = None) -> Any:
      return u[:, 1]

    def hfunc2(self, u: Any, *, x: Any = None) -> Any:
      return u[:, 0]

  placed = _Bare()._prep(np.array([[0.3, 0.5]], dtype=np.float32))
  assert placed.dtype is torch.float64 and placed.device.type == "cpu"


def test_fit_and_select_place_their_data(device: str) -> None:
  """The two re-estimators placed nothing, unlike every other entry point.

  `TorchVinecop.fit(numpy_u)` on a CUDA vine raised `can't convert cuda:0
  device type tensor to numpy`: the engines allocate their per-tree scratch in
  the data's namespace, so a vine whose pairs answer elsewhere then assigns
  across namespaces. `from_data` was unaffected -- it places from the controls.
  """
  u_np = pv.to_pseudo_obs(np.random.default_rng(0).normal(size=(150, 3)))
  controls = FitControlsTorchVinecop(
    device=torch.device(device), dtype=torch.float64
  )
  vine = TorchVinecop.from_data(
    torch.as_tensor(u_np, dtype=torch.float64, device=device), controls
  )
  assert vine.fit(u_np, controls) is vine
  assert vine.select(u_np, controls) is vine
  out = vine.pdf(torch.as_tensor(u_np[:4], dtype=torch.float64, device=device))
  assert out.shape == (4,) and out.device.type == device
