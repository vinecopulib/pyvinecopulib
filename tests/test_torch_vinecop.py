"""Tests for the optional ``pyvinecopulib.torch.TorchVinecop`` class.

Skipped when PyTorch isn't installed. Compares the torch ``TorchVinecop``
against the C++ ``pv.Vinecop`` (with TLL pair copulas) for ``pdf`` /
``rosenblatt`` / ``inverse_rosenblatt``; verifies the round-trip
``inverse_rosenblatt(rosenblatt(u)) ≈ u``; spot-checks an
independent vine; and confirms the truncation and discrete paths.
"""

from __future__ import annotations

import pickle
from typing import Any, Callable

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import DiscretePair

torch = pytest.importorskip("torch")

from pyvinecopulib.torch import (  # noqa: E402
  FitControlsTorchBicop,
  FitControlsTorchVinecop,
  TorchBicop,
  TorchVinecop,
)

_TLL_BICOP = pv.FitControlsBicop(family_set=[pv.families.tll])
_TLL_CONTROLS = pv.FitControlsVinecop(
  family_set=[pv.families.tll], num_threads=1
)


def _fit_tll_vine(u: np.ndarray, **extra) -> pv.Vinecop:
  ctl = pv.FitControlsVinecop(
    family_set=[pv.families.tll], num_threads=1, **extra
  )
  return pv.Vinecop.from_data(u, controls=ctl)


def _eval_grid(n: int, d: int, seed: int = 0) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return rng.uniform(0.02, 0.98, size=(n, d))


def _simulate(d: int, n: int, seed: int = 0) -> np.ndarray:
  # Use a non-trivial structure: Gaussian pair-copula vine via Clayton
  # samples per pair would require building a parametric Vinecop. Easier:
  # take pseudo-observations from a multivariate normal with banded
  # correlation, which gives a smoothly varying density per pair.
  rng = np.random.default_rng(seed)
  base = rng.standard_normal(size=(n, 1))
  noise = rng.standard_normal(size=(n, d))
  x = 0.6 * base + 0.4 * noise
  return pv.to_pseudo_obs(x)


def test_pdf_matches_pvvinecop() -> None:
  u_fit = _simulate(d=5, n=2000, seed=1)
  cop_tll = _fit_tll_vine(u_fit)

  # Pin cache=False — this test verifies vine-level parity with the C++
  # on-the-fly cascade at 1e-10. The default cache_integrals=True uses
  # bilinear interp; its precision floor is pinned separately in
  # test_cache_integrals_precision_floor.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(500, d=5, seed=11)

  out_torch = bc.pdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_rosenblatt_matches_pvvinecop() -> None:
  u_fit = _simulate(d=5, n=2000, seed=2)
  cop_tll = _fit_tll_vine(u_fit)

  # Pin cache=False — see test_pdf_matches_pvvinecop.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(500, d=5, seed=12)

  out_torch = bc.rosenblatt(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.rosenblatt(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize(
  "method,cache,max_abs,mean_abs",
  [
    # cache=False is bisection-precise at every step (forward cascades are
    # exact trapezoidal). Vine-level error stays at machine precision for
    # pdf/rosenblatt and at bisection convergence for inverse_rosenblatt.
    ("pdf", False, 1e-11, 1e-13),
    ("rosenblatt", False, 1e-10, 1e-13),
    ("inverse_rosenblatt", False, 1e-8, 1e-9),
    # cache=True reconstructs the same integrals in closed form from the
    # prefix tables, so it differs from cache=False only in summation
    # order: the same bounds hold. `inverse_rosenblatt` shares one code
    # path across both modes, so its two rows are bit-identical.
    ("pdf", True, 1e-11, 1e-13),
    ("rosenblatt", True, 1e-10, 1e-13),
    ("inverse_rosenblatt", True, 1e-8, 1e-9),
  ],
)
def test_cache_integrals_precision_floor(
  method: str, cache: bool, max_abs: float, mean_abs: float
) -> None:
  """Pin the precision floor for each cache mode at a moderately-deep vine.

  Both modes preserve precision through the cascade: ``cache_integrals=False``
  integrates and inverts per call, ``True`` reads the closed-form value off the
  prefix tables, and the two differ only in the order the same terms are
  summed. So the pins are the same in each column, and a regression in the
  cached path shows up as a bound violation rather than as a widened tolerance.
  """
  d = 10
  u_fit = _simulate(d=d, n=2000, seed=1)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=cache)
  u_eval = _eval_grid(500, d=d, seed=2)
  out_cpp = getattr(cop_tll, method)(u_eval)
  out_torch = getattr(bc, method)(torch.from_numpy(u_eval)).numpy()
  diff = np.abs(out_torch - out_cpp)
  assert diff.max() < max_abs, f"{method}/{cache}: max {diff.max():.2e}"
  assert diff.mean() < mean_abs, f"{method}/{cache}: mean {diff.mean():.2e}"


def test_inverse_rosenblatt_matches_pvvinecop() -> None:
  u_fit = _simulate(d=5, n=2000, seed=3)
  cop_tll = _fit_tll_vine(u_fit)

  # Pin cache=False — C++ parity at 1e-9 via the bisection-precise
  # on-the-fly cascade. The cached path's precision floor is checked
  # by test_cache_integrals_precision_floor.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  # Inverse-rosenblatt takes independent uniforms; sample fresh.
  w = _eval_grid(400, d=5, seed=13)

  out_torch = bc.inverse_rosenblatt(torch.from_numpy(w)).numpy()
  out_cpp = cop_tll.inverse_rosenblatt(w)
  # Bisection at n_iter=35 gives bracket-width accuracy 0.5**35 ≈ 6e-11.
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-9, rtol=1e-9)


def test_inverse_rosenblatt_roundtrip() -> None:
  u_fit = _simulate(d=5, n=2000, seed=4)
  cop_tll = _fit_tll_vine(u_fit)

  # Pin cache=False so the per-pair forward/inverse cascade is
  # bisection-precise; the round-trip identity requires the on-the-fly
  # path to hit the 99-percentile tolerance below.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  rng = np.random.default_rng(14)
  u_eval = rng.uniform(0.1, 0.9, size=(400, 5))
  u_t = torch.from_numpy(u_eval)
  back = bc.inverse_rosenblatt(bc.rosenblatt(u_t)).numpy()
  diff = np.abs(back - u_eval)
  # ITP makes per-call precision ~1e-15. The cascade's irreducible noise
  # is the [TRIM_LO, TRIM_HI] = [1e-10, 1-1e-10] clamp inside hfunc /
  # hinv, which creates plateaus the solver cannot resolve at the
  # boundaries. The bulk of u_eval roundtrips at machine precision; a
  # tiny fraction (~0.1%) of points whose forward pass saturates against
  # the clamp don't, so we measure the 99th percentile rather than the
  # max.
  assert np.percentile(diff, 99) < 1e-10, (
    f"99th percentile diff = {np.percentile(diff, 99):.2e}; "
    f"max = {diff.max():.2e}"
  )


@pytest.mark.parametrize("batched", [False, True])
def test_independent_vine_short_circuits(batched: bool) -> None:
  d = 4
  structure = pv.RVineStructure.from_order([1, 2, 3, 4])
  # Build a vine with all-independent pair copulas.
  pair_copulas = [
    [pv.Bicop(family=pv.families.indep) for _ in range(d - 1 - t)]
    for t in range(d - 1)
  ]
  cop = pv.Vinecop.from_structure(
    structure=structure, pair_copulas=pair_copulas
  )

  # cache_integrals doesn't matter for the all-independence short-circuit
  # (every pair short-circuits to pdf=1, hfunc=u, etc.), but we pin
  # cache=False so the 1e-10 / 1e-12 tolerances below hold under either
  # default.
  bc = TorchVinecop.from_vinecop(cop, cache_integrals=False)
  u = _eval_grid(200, d=d, seed=21)
  u_t = torch.from_numpy(u)

  np.testing.assert_allclose(
    bc.pdf(u_t, batched=batched).numpy(), np.ones(200), atol=1e-12, rtol=0
  )
  np.testing.assert_allclose(
    bc.rosenblatt(u_t, batched=batched).numpy(), u, atol=1e-10, rtol=1e-10
  )
  np.testing.assert_allclose(
    bc.inverse_rosenblatt(u_t, batched=batched).numpy(),
    u,
    atol=1e-10,
    rtol=1e-10,
  )


@pytest.mark.parametrize("batched", [False, True])
def test_truncated_vine_pdf(batched: bool) -> None:
  u_fit = _simulate(d=6, n=2000, seed=5)
  cop_tll = _fit_tll_vine(u_fit, trunc_lvl=2)
  assert cop_tll.structure.trunc_lvl == 2

  # cache=False for 1e-10 C++ parity.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  assert bc.trunc_lvl == 2

  u_eval = _eval_grid(300, d=6, seed=15)
  np.testing.assert_allclose(
    bc.pdf(torch.from_numpy(u_eval), batched=batched).numpy(),
    cop_tll.pdf(u_eval),
    atol=1e-10,
    rtol=1e-10,
  )


def test_from_vinecop_rejects_unsupported_family() -> None:
  d = 3
  structure = pv.RVineStructure.from_order([1, 2, 3])
  pair_copulas = [
    [
      pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
      for _ in range(d - 1 - t)
    ]
    for t in range(d - 1)
  ]
  cop = pv.Vinecop.from_structure(
    structure=structure, pair_copulas=pair_copulas
  )
  with pytest.raises(ValueError, match="tll and indep"):
    TorchVinecop.from_vinecop(cop)


def test_pdf_shape_validation() -> None:
  u_fit = _simulate(d=4, n=1000, seed=6)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)

  with pytest.raises(ValueError, match=r"shape \(n, 4\)"):
    bc.pdf(torch.zeros(10, 5, dtype=torch.float64))


# --------------------------------------------------------------------- #
# Batched cascade                                                        #
# --------------------------------------------------------------------- #


@pytest.mark.parametrize("cache", [False, True])
def test_pdf_batched_matches_unbatched(cache: bool) -> None:
  u_fit = _simulate(d=10, n=800, seed=300)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=cache)
  u_t = torch.from_numpy(_eval_grid(500, d=10, seed=310))
  ref = bc.pdf(u_t, batched=False).numpy()
  got = bc.pdf(u_t, batched=True).numpy()
  np.testing.assert_allclose(got, ref, atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("cache", [False, True])
def test_rosenblatt_batched_matches_unbatched(cache: bool) -> None:
  u_fit = _simulate(d=10, n=800, seed=301)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=cache)
  u_t = torch.from_numpy(_eval_grid(500, d=10, seed=311))
  ref = bc.rosenblatt(u_t, batched=False).numpy()
  got = bc.rosenblatt(u_t, batched=True).numpy()
  np.testing.assert_allclose(got, ref, atol=1e-13, rtol=1e-13)


def test_inverse_rosenblatt_batched_matches_sequential() -> None:
  """The wave-scheduled inverse is a reordering, so it agrees bit-for-bit.

  The inverse cascade's dependency graph is genuinely 2-D -- a cell depends on
  one at a *different* tree level -- so the per-tree-level wavefront that works
  for pdf / rosenblatt does not apply, and the batched path groups cells by
  longest-path level of that graph instead. Every cell still consumes exactly
  the values it consumed sequentially, only sooner, so equality is exact rather
  than approximate.
  """
  u_fit = _simulate(d=6, n=600, seed=302)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  w_t = torch.from_numpy(_eval_grid(300, d=6, seed=312))
  out_seq = bc.inverse_rosenblatt(w_t, batched=False)
  out_bat = bc.inverse_rosenblatt(w_t, batched=True)
  torch.testing.assert_close(out_bat, out_seq, atol=0.0, rtol=0.0)


def test_batched_matches_cpp_pdf() -> None:
  """End-to-end: torch batched pdf vs pv.Vinecop on the same fit."""
  u_fit = _simulate(d=8, n=1000, seed=320)
  cop_tll = _fit_tll_vine(u_fit)
  # cache=False for 1e-10 C++ parity.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(400, d=8, seed=330)
  out_torch = bc.pdf(torch.from_numpy(u_eval), batched=True).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_batched_matches_cpp_rosenblatt() -> None:
  u_fit = _simulate(d=8, n=1000, seed=321)
  cop_tll = _fit_tll_vine(u_fit)
  # cache=False for 1e-10 C++ parity.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(400, d=8, seed=331)
  out_torch = bc.rosenblatt(torch.from_numpy(u_eval), batched=True).numpy()
  out_cpp = cop_tll.rosenblatt(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("d", [5, 10])
def test_from_data_matches_from_vinecop(d: int) -> None:
  """``TorchVinecop.from_data(u, structure)`` should agree with
  ``pv.Vinecop.from_data(u, controls=..., structure=structure)`` to machine
  precision. Note: this is the *fixed-structure* C++ path. Freely-selected
  C++ Vinecop.from_data produces a different pair-copula orientation per
  edge (the graph-traversal order during selection differs from the order
  used when the structure is passed in), so the comparison target is the
  fixed-structure refit.
  """
  u_fit = _simulate(d=d, n=1500, seed=400 + d)
  # First, get a structure (any TLL fit will do).
  structure = _fit_tll_vine(u_fit).structure
  # Now refit C++ with that structure passed in.
  cop_cpp_fixed = pv.Vinecop.from_data(
    u_fit, controls=_TLL_CONTROLS, structure=structure
  )

  # Pin cache=False so the per-pair cascade is bisection-precise — needed
  # to match at 1e-9 (the cached path round-trips only to ~1e-3).
  bc_cpp = TorchVinecop.from_vinecop(cop_cpp_fixed, cache_integrals=False)
  bc_torch = TorchVinecop.from_data(
    torch.from_numpy(u_fit),
    structure,
    controls=FitControlsTorchVinecop(cache_integrals=False),
  )

  u_eval = torch.from_numpy(_eval_grid(400, d=d, seed=410 + d))
  np.testing.assert_allclose(
    bc_torch.pdf(u_eval).numpy(),
    bc_cpp.pdf(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )
  np.testing.assert_allclose(
    bc_torch.rosenblatt(u_eval).numpy(),
    bc_cpp.rosenblatt(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )
  np.testing.assert_allclose(
    bc_torch.inverse_rosenblatt(u_eval).numpy(),
    bc_cpp.inverse_rosenblatt(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )


def test_from_data_auto_selects_structure() -> None:
  """``from_data(structure=None)`` selects the R-vine natively in torch
  (``VinecopBase.select``, TLL, ``trunc_lvl=20``) and reuses the
  selection-time pairs — an exact port of ``pv.Vinecop.from_data``, so the
  selected structure (identical matrix encoding), the reoriented pairs, and
  the resulting density must all match."""
  u_fit = _simulate(d=5, n=1500, seed=777)
  cop_cpp = pv.Vinecop.from_data(
    u_fit,
    structure=None,
    var_types=[],
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
    ),
  )
  got = TorchVinecop.from_data(
    torch.from_numpy(u_fit),
    controls=FitControlsTorchVinecop(cache_integrals=False),
  )
  # Identical structure, down to the matrix encoding.
  np.testing.assert_array_equal(
    np.asarray(got.structure.matrix), np.asarray(cop_cpp.structure.matrix)
  )
  # Identical density: same selection, same reused pairs (torch TLL matches
  # the ``Bicop`` TLL fit to machine precision), same flips.
  reference = TorchVinecop.from_vinecop(cop_cpp, cache_integrals=False)
  u_eval = torch.from_numpy(_eval_grid(400, d=5, seed=778))
  np.testing.assert_allclose(
    got.pdf(u_eval).numpy(),
    reference.pdf(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )


def test_batched_to_device_invalidates() -> None:
  """``.to()`` should drop the lazily-built BatchedVine so the next
  batched call rebuilds it on the new device."""
  u_fit = _simulate(d=5, n=500, seed=322)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  u_t = torch.from_numpy(_eval_grid(50, d=5, seed=332))
  bc.pdf(u_t, batched=True)  # build the BatchedVine
  assert bc._batched is not None
  bc.to(torch.float32)
  assert bc._batched is None
  # And the re-build on the new dtype produces a finite output.
  out = bc.pdf(u_t.to(torch.float32), batched=True)
  assert torch.isfinite(out).all()


# --------------------------------------------------------------------------- #
# User-facing input validation for TorchVinecop                                #
# --------------------------------------------------------------------------- #


def test_from_data_rejects_non_2d_u() -> None:
  """``TorchVinecop.from_data`` rejects 1-D / 3-D inputs."""
  structure = pv.RVineStructure.sample(5, seeds=[1])
  with pytest.raises(ValueError, match="must be 2-D"):
    TorchVinecop.from_data(torch.zeros(100, dtype=torch.float64), structure)
  with pytest.raises(ValueError, match="must be 2-D"):
    TorchVinecop.from_data(
      torch.zeros(100, 5, 2, dtype=torch.float64), structure
    )


def test_from_data_rejects_structure_dim_mismatch() -> None:
  """``TorchVinecop.from_data`` rejects ``structure.dim != u.shape[1]``."""
  structure = pv.RVineStructure.sample(5, seeds=[1])  # d=5
  with pytest.raises(ValueError, match="does not match"):
    TorchVinecop.from_data(torch.zeros(100, 7, dtype=torch.float64), structure)


@pytest.mark.parametrize("op", ["pdf", "rosenblatt", "inverse_rosenblatt"])
def test_eval_rejects_wrong_input_shape(op: str) -> None:
  """Every eval entry point requires ``(n, d)`` and raises otherwise."""
  d = 5
  u_fit = _simulate(d=d, n=400, seed=520)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  fn = getattr(bc, op)
  with pytest.raises(ValueError):
    fn(torch.zeros(100, dtype=torch.float64))  # 1-D
  with pytest.raises(ValueError):
    fn(torch.zeros(100, d + 1, dtype=torch.float64))  # wrong second dim


# ---------------------------------------------------------------------------
# API alignment with pv.Vinecop (from_structure, sample, cdf, num_threads)
# ---------------------------------------------------------------------------


def test_from_structure_pair_copulas_matches_init():
  """from_structure(structure, pair_copulas) builds the same vine as
  TorchVinecop(pair_copulas, structure)."""
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (300, 3))
  cop = _fit_tll_vine(U)
  pcs = [[TorchBicop.from_bicop(b) for b in row] for row in cop.pair_copulas]
  via_init = TorchVinecop(pair_copulas=pcs, structure=cop.structure)
  via_factory = TorchVinecop.from_structure(
    structure=cop.structure, pair_copulas=pcs
  )

  u_t = torch.from_numpy(U[:10])
  assert torch.allclose(via_init.pdf(u_t), via_factory.pdf(u_t))


def test_from_structure_empty_is_independence_vine():
  """No pair_copulas -> every edge is an independence TorchBicop;
  pdf is identically 1."""
  s = pv.RVineStructure.sample(4, seeds=[1, 2, 3, 4, 5])
  tv = TorchVinecop.from_structure(structure=s)
  u = torch.rand(50, 4, dtype=torch.float64)
  assert torch.allclose(
    tv.pdf(u), torch.ones(50, dtype=torch.float64), atol=1e-12
  )


def test_from_structure_matrix_path():
  """Matrix kwarg routes through RVineStructure.from_matrix."""
  s = pv.RVineStructure.sample(3, seeds=[1, 2, 3, 4, 5])
  mat = np.asarray(s.matrix, dtype=np.uint64)
  tv = TorchVinecop.from_structure(matrix=mat)
  assert tv.d == 3
  assert tv.trunc_lvl == 2


def test_from_structure_rejects_both_or_neither():
  s = pv.RVineStructure.sample(3, seeds=[1, 2, 3, 4, 5])
  with pytest.raises(ValueError, match="exactly one"):
    TorchVinecop.from_structure()
  with pytest.raises(ValueError, match="exactly one"):
    TorchVinecop.from_structure(
      structure=s, matrix=np.asarray(s.matrix, dtype=np.uint64)
    )


def test_simulate_seeded_reproducible():
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (300, 3))
  tv = TorchVinecop.from_vinecop(_fit_tll_vine(U))

  a = tv.sample(n=200, seeds=[42])
  b = tv.sample(n=200, seeds=[42])
  assert torch.allclose(a, b)
  c = tv.sample(n=200, seeds=[43])
  assert not torch.allclose(a, c)


def test_simulate_qrng_runs_and_returns_uniforms():
  """qrng=True draws Halton/Sobol via pv.utils.sample_uniform, then
  pushes through inverse_rosenblatt. We check basic shape + range."""
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (300, 3))
  tv = TorchVinecop.from_vinecop(_fit_tll_vine(U))

  s = tv.sample(n=500, qrng=True, seeds=[1, 2, 3])
  assert s.shape == (500, 3)
  assert (s > 0).all() and (s < 1).all()


def test_cdf_matches_cpp_within_mc_error():
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (300, 3))
  cop = _fit_tll_vine(U)
  tv = TorchVinecop.from_vinecop(cop)
  query = torch.from_numpy(U[:10])

  cdf_torch = tv.cdf(query, N=20000, qrng=True, seeds=[1, 2, 3]).numpy()
  cdf_cpp = np.asarray(cop.cdf(U[:10], N=20000, seeds=[1, 2, 3]))
  # Both are MC estimates with N=20000; agree to ~3 sig figs.
  np.testing.assert_allclose(cdf_torch, cdf_cpp, atol=2e-2)


def test_pdf_accepts_num_threads_no_op():
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (200, 3))
  tv = TorchVinecop.from_vinecop(_fit_tll_vine(U))
  u_t = torch.from_numpy(U[:10])
  p1 = tv.pdf(u_t, num_threads=1)
  p4 = tv.pdf(u_t, num_threads=4)
  assert torch.allclose(p1, p4)


def test_rosenblatt_inverse_rosenblatt_accept_num_threads():
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (200, 3))
  tv = TorchVinecop.from_vinecop(_fit_tll_vine(U))
  u_t = torch.from_numpy(U[:5])
  # Just confirm the kwarg is accepted; behavior is identical.
  r1 = tv.rosenblatt(u_t, num_threads=1)
  r4 = tv.rosenblatt(u_t, num_threads=4)
  assert torch.allclose(r1, r4)
  ir1 = tv.inverse_rosenblatt(u_t, num_threads=1)
  ir4 = tv.inverse_rosenblatt(u_t, num_threads=4)
  assert torch.allclose(ir1, ir4)


def test_pdf_autograd_through_grid_param() -> None:
  """Grad flows through ``TorchVinecop.pdf`` to a ``TorchBicop`` grid value.

  Guards the Stage-2a cascade extraction onto ``VinecopBase``: ``pdf`` /
  ``rosenblatt`` must stay autograd-capable (only inverse / sample / cdf are
  wrapped in ``no_grad``), and the cascade must not detach pair outputs.
  """
  u_fit = _simulate(d=4, n=800, seed=700)
  cop = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop, cache_integrals=False)
  # The stored module, not what `_get_pair_copula` hands the cascade: that is
  # the grid the parameters live on, and it is the same object here (this vine
  # is continuous, so nothing is wrapped).
  pair = bc._pair_module(0, 0)
  pair.interp_grid.values.requires_grad_(True)
  u = torch.from_numpy(_eval_grid(64, d=4, seed=701))
  out = bc.pdf(u, batched=False)
  assert out.requires_grad
  out.sum().backward()
  grad = pair.interp_grid.values.grad
  assert grad is not None
  assert torch.isfinite(grad).all()
  assert grad.abs().sum() > 0


def test_conditional_cdf_raises() -> None:
  """``cdf`` with external covariates ``x`` is unsupported (raises)."""
  u_fit = _simulate(d=4, n=400, seed=710)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  u = torch.from_numpy(_eval_grid(5, d=4, seed=711))
  x = torch.zeros(5, 1, dtype=torch.float64)
  with pytest.raises(NotImplementedError, match="Conditional cdf"):
    bc.cdf(u, x=x)


def test_sample_conditional_requires_row_per_sample() -> None:
  """``sample(n, x=...)`` requires exactly one covariate row per sample."""
  u_fit = _simulate(d=3, n=400, seed=712)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  x = torch.zeros(4, 1, dtype=torch.float64)
  with pytest.raises(ValueError, match="one covariate row per sample"):
    bc.sample(10, x=x)


def test_batched_cache_stays_out_of_the_state_dict() -> None:
  # The batched state is a memo derived from the pair copulas, not part of
  # the model. Registering it as a child module would put derived buffers in
  # every checkpoint taken after a batched call, and `load_state_dict` on a
  # fresh model would then reject them as unexpected keys.
  rng = np.random.default_rng(0)
  d = 4
  u = pv.to_pseudo_obs(
    rng.standard_normal((400, d)) @ rng.standard_normal((d, d))
  )
  cpp = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.tll])
  )
  vine = TorchVinecop.from_vinecop(cpp)
  data = torch.as_tensor(u)

  keys_before = set(vine.state_dict())
  vine.pdf(data, batched=True)
  vine.rosenblatt(data, batched=True)
  assert set(vine.state_dict()) == keys_before

  fresh = TorchVinecop.from_vinecop(cpp)
  fresh.load_state_dict(vine.state_dict(), strict=True)
  torch.testing.assert_close(fresh.pdf(data), vine.pdf(data))


# --- gradients ---------------------------------------------------------------- #


def _central_fd(fn, q, eps: float = 1e-6):
  """Central finite differences of a scalar-per-row function in every column."""
  fd = torch.zeros_like(q)
  with torch.no_grad():
    for j in range(q.shape[1]):
      hi, lo = q.clone(), q.clone()
      hi[:, j] += eps
      lo[:, j] -= eps
      fd[:, j] = (fn(hi) - fn(lo)) / (2 * eps)
  return fd


@pytest.mark.parametrize("cache_integrals", [True, False])
@pytest.mark.parametrize("batched", [False, True])
def test_pdf_gradient_matches_finite_differences(
  cache_integrals: bool, batched: bool
) -> None:
  """`d pdf / d u` must be right on all four (cache, batched) paths.

  It was not: with `cache_integrals=True` the cached `hfunc1` / `hfunc2` were
  evaluated under `torch.no_grad()`, so the cascade's conditioning arguments were
  detached and only the final density factors contributed. The error was 111% of
  the largest finite-difference entry — and only on the non-batched path, because
  the batched level calls the interpolation primitives directly and was never
  wrapped. So
  the two paths silently disagreed on gradients while agreeing on values.

  Query points are drawn away from the grid nodes: the interpolant is piecewise
  bilinear, so a finite difference straddling a knot compares two different
  polynomials.
  """
  cop = _fit_tll_vine(_simulate(d=5, n=400, seed=4))
  bc = TorchVinecop.from_vinecop(cop, cache_integrals=cache_integrals)
  q = torch.from_numpy(
    np.random.default_rng(23).uniform(0.15, 0.85, size=(8, 5))
  )

  x = q.clone().requires_grad_(True)
  (grad,) = torch.autograd.grad(bc.pdf(x, batched=batched).sum(), x)
  fd = _central_fd(lambda uu: bc.pdf(uu, batched=batched), q)

  assert grad is not None
  rel = (grad - fd).abs().max().item() / fd.abs().max().item()
  assert rel < 1e-6, f"cache={cache_integrals} batched={batched}: {rel:.3e}"


def test_cached_and_uncached_gradients_agree_in_direction() -> None:
  """Bounded, not equal: with a cache the surrogate *is* the model.

  The cached members interpolate a table of the integrals bilinearly, so their
  derivative is the surrogate's, not the exact integral's. That difference is
  the same order as the documented value gap, so pin the agreement as a
  direction rather than a value.
  """
  cop = _fit_tll_vine(_simulate(d=4, n=400, seed=5))
  q = torch.from_numpy(
    np.random.default_rng(24).uniform(0.15, 0.85, size=(16, 4))
  )
  grads = []
  for cache in (True, False):
    bc = TorchVinecop.from_vinecop(cop, cache_integrals=cache)
    x = q.clone().requires_grad_(True)
    (g,) = torch.autograd.grad(bc.pdf(x, batched=False).sum(), x)
    grads.append(g.flatten())

  cosine = torch.nn.functional.cosine_similarity(grads[0], grads[1], dim=0)
  assert cosine.item() > 0.99, cosine.item()


def test_the_batched_cache_re_bakes_when_grad_tracking_changes() -> None:
  """`requires_grad_` after a batched call must not leave a stale cache.

  The bake copies each pair's grid into a stacked tensor and `requires_grad_`
  mutates a flag in place, so flipping it afterwards used to leave the copy
  behind. It did not raise where it was read: the batched cascade returned a
  value detached from the grid, and it surfaced as torch's generic "does not
  require grad" out of `backward()`. Both orderings must now give the same
  gradient, and it must equal the non-batched one.
  """
  fitted = _fit_tll_vine(_simulate(d=3, n=600, seed=900))
  u = torch.from_numpy(_eval_grid(48, d=3, seed=901))

  def lift():
    return TorchVinecop.from_vinecop(fitted, cache_integrals=False)

  cop = lift()
  # Bake first -- the ordering a caller hits by evaluating before deciding to
  # optimize -- then start tracking the grid.
  baked = cop.pdf(u, batched=True)
  assert not baked.requires_grad
  values = cop._get_pair_copula(0, 0).interp_grid.values
  values.requires_grad_(True)

  (g_batched,) = torch.autograd.grad(cop.pdf(u, batched=True).sum(), values)
  (g_plain,) = torch.autograd.grad(cop.pdf(u, batched=False).sum(), values)
  assert float(g_batched.abs().sum()) > 0.0
  torch.testing.assert_close(g_batched, g_plain, rtol=1e-10, atol=1e-12)

  # The other ordering agrees, so neither is the special case.
  other = lift()
  other_values = other._get_pair_copula(0, 0).interp_grid.values
  other_values.requires_grad_(True)
  (g_first,) = torch.autograd.grad(
    other.pdf(u, batched=True).sum(), other_values
  )
  torch.testing.assert_close(g_batched, g_first, rtol=1e-10, atol=1e-12)


#: Cumulative atom masses of a Binomial(4, 0.5) count variable: its top atom
#: reaches ``F(x) = 1`` and its bottom ``F(x^-) = 0``, both places where the
#: cascade's clamp and a pair copula's own trimming have to agree.
_COUNT_CDF = np.cumsum([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0


def _discrete_data(
  var_types: list[str], n: int = 1500, seed: int = 3
) -> np.ndarray:
  """An ``(n, d + k)`` sample: counts at ``"d"``, continuous ranks at ``"c"``."""
  rng = np.random.default_rng(seed)
  d = len(var_types)
  base = rng.standard_normal((n, 1))
  z = 0.7 * base + 0.7 * rng.standard_normal((n, d))
  values, limits = [], []
  for j, t in enumerate(var_types):
    p = pv.to_pseudo_obs(z[:, [j]]).ravel()
    if t == "d":
      level = np.searchsorted(_COUNT_CDF[:-1], p)
      values.append(_COUNT_CDF[level])
      limits.append(np.where(level > 0, _COUNT_CDF[level - 1], 0.0))
    else:
      values.append(p)
  return np.column_stack(values + limits)


def _discrete_vinecop(var_types: list[str], u: np.ndarray) -> pv.Vinecop:
  return pv.Vinecop.from_data(u, var_types=var_types, controls=_TLL_CONTROLS)


@pytest.mark.parametrize(
  "var_types", [["d", "c", "c"], ["c", "d", "d"], ["d", "d", "d"]]
)
def test_from_vinecop_matches_discrete_vinecop(var_types: list[str]) -> None:
  # A discrete C++ vine lifted into torch. The stored grids are continuous and
  # `DiscretePair` supplies the mixed-discrete surface, so what is compared is
  # two difference quotients over the same grid: measured 1.0e-14 / 2.4e-14 /
  # 5.6e-14 across the three type patterns. `cache_integrals=None` resolves to
  # `False` here, which is what makes that floor reachable -- differencing the
  # bilinearly interpolated cached cdf gives 38% instead.
  u = _discrete_data(var_types)
  cop = _discrete_vinecop(var_types, u)
  bc = TorchVinecop.from_vinecop(cop)
  assert bc.var_types == var_types
  u_t = torch.from_numpy(u)
  np.testing.assert_allclose(
    bc.pdf(u_t).numpy(), cop.pdf(u), atol=1e-12, rtol=1e-12
  )
  np.testing.assert_allclose(
    bc.rosenblatt(u_t, randomize_discrete=False).numpy(),
    cop.rosenblatt(u, randomize_discrete=False),
    atol=1e-12,
    rtol=1e-12,
  )


def _first_discrete_edge(
  vine: TorchVinecop,
) -> tuple[int, int, tuple[str, str]]:
  """The first edge whose variable types include an atom.

  Which edge that is depends on the selected structure, so it has to be looked
  up rather than assumed: a `["d", "c", "c"]` vine may well pair its two
  continuous variables at ``(0, 0)``.
  """
  for tree in range(vine.trunc_lvl):
    for edge in range(vine.d - 1 - tree):
      types = vine.pair_var_types(tree, edge)
      if "d" in types:
        return tree, edge, types
  raise AssertionError("no discrete edge in a vine with discrete variables")


def _pair_edge(types: tuple[str, str], seed: int) -> torch.Tensor:
  """A pair-copula input for an edge with the given variable types.

  The pair contract is the expanded four-column layout `[u1, u2, u1^-, u2^-]`
  whenever either argument has atoms, and a continuous argument's left limit is
  its own value. Atoms are `0.125` wide, what eight equiprobable levels give,
  so the quotients are the ones a fitted count vine evaluates.
  """
  rng = np.random.default_rng(seed)
  width = 0.125
  vals = rng.uniform(0.2 + width, 0.95, size=(300, 2))
  limits = [
    vals[:, j] - width if ty == "d" else vals[:, j]
    for j, ty in enumerate(types)
  ]
  return torch.from_numpy(np.column_stack([vals[:, 0], vals[:, 1], *limits]))


@pytest.mark.parametrize(
  "var_types", [["d", "c", "c"], ["c", "d", "d"], ["d", "d", "d"]], ids=str
)
def test_the_integral_cache_is_allowed_on_a_discrete_vine(
  var_types: list[str],
) -> None:
  """A discrete vine takes the cache like any other, and agrees without it.

  It used to be refused: the cache was a bilinear interpolation of tabulated
  integrals, and a discrete edge reads its density from *differences* over an
  atom's width -- 38% maximum relative error on a ``("d","d")`` density. The
  tables are exact now, so both modes difference the same quantity and
  differ only in the order its terms are summed. The pin is that agreement
  rather than a refusal.
  """
  u = _discrete_data(var_types, n=400, seed=6)
  cop = _discrete_vinecop(var_types, u)
  u_t = torch.from_numpy(u)

  cached = TorchVinecop.from_vinecop(cop, cache_integrals=True)
  plain = TorchVinecop.from_vinecop(cop, cache_integrals=False)
  assert cached._pair_module(0, 0)._cache_integrals
  assert not plain._pair_module(0, 0)._cache_integrals
  # `None` no longer resolves differently for a discrete vine.
  assert TorchVinecop.from_vinecop(cop)._pair_module(0, 0)._cache_integrals

  # Tier one, and the tight one: at the pair level the two modes are the same
  # formula summed differently, so this is where a cache regression shows -- the
  # bilinear cache this replaced read 38% here.
  tree, edge_ix, types = _first_discrete_edge(cached)
  edge = _pair_edge(types, seed=61)
  for method in ("pdf", "hfunc1", "hfunc2"):
    torch.testing.assert_close(
      getattr(DiscretePair(cached._pair_module(tree, edge_ix), types), method)(
        edge
      ),
      getattr(DiscretePair(plain._pair_module(tree, edge_ix), types), method)(
        edge
      ),
      rtol=1e-12,
      atol=1e-13,
    )

    # Tier two: the vine. Looser than the pair's because the discrete cascade
    # amplifies -- tree 1's atom widths are themselves differences of tree 0's
    # h-functions -- but only to 1.6e-13 relative, measured, since both modes
    # difference the same exact distribution function and only the summation
    # order separates them.
    #
    # `rosenblatt` needs its randomization pinned. It jitters the discrete
    # columns by default, so two calls on the *same* vine disagree by O(10) --
    # comparing two vines without pinning it measures the draw, not the cache.
    torch.testing.assert_close(
      cached.pdf(u_t), plain.pdf(u_t), rtol=1e-11, atol=1e-10
    )
    torch.testing.assert_close(
      cached.rosenblatt(u_t, randomize_discrete=False),
      plain.rosenblatt(u_t, randomize_discrete=False),
      rtol=1e-11,
      atol=1e-10,
    )


def test_from_structure_keeps_discrete_var_types() -> None:
  s = pv.RVineStructure.sample(3, seeds=[1, 2, 3, 4, 5])
  vine = TorchVinecop.from_structure(structure=s, var_types=["c", "d", "c"])
  assert vine.var_types == ["c", "d", "c"]
  # An independence vine has density 1 whatever the variable types are, and the
  # discrete edges still consume the four-column input.
  u = _discrete_data(["c", "d", "c"], n=64, seed=5)
  np.testing.assert_allclose(
    vine.pdf(torch.from_numpy(u)).numpy(), np.ones(64), atol=1e-12
  )


@pytest.mark.parametrize("var_types", [["d", "c"], ["c", "d"], ["d", "d"]])
def test_the_discrete_pair_fit_reproduces_the_compiled_grid(
  var_types: list[str],
) -> None:
  """The fitted grid itself, not the vine that hosts it.

  ``TllBicop::fit`` selects the bandwidth from randomly-jittered ranks and then
  fits on the *latent sample* drawn with it, so a discrete fit needs both steps.
  ``TorchBicop.from_data`` reuses the compiled ``find_latent_sample`` for the
  second, since it is a stochastic iterative reconstruction over a spatial index
  and a reimplementation would put this parity at its mercy rather than at the
  same code's. Measured 1.5e-13 / 1.4e-13 / 5.5e-14.
  """
  u = _discrete_data(var_types, n=2000, seed=5)
  # The pair contract is the expanded four-column layout: a continuous column
  # is its own left limit.
  k, limits = 0, []
  for j, ty in enumerate(var_types):
    if ty == "d":
      limits.append(u[:, 2 + k])
      k += 1
    else:
      limits.append(u[:, j])
  u4 = np.column_stack([u[:, 0], u[:, 1], *limits])

  ref = pv.Bicop.from_data(u, controls=_TLL_BICOP, var_types=var_types)
  bc = TorchBicop.from_data(
    torch.from_numpy(u4), cache_integrals=False, var_types=var_types
  )
  np.testing.assert_allclose(
    bc.interp_grid.values.numpy(),
    np.asarray(ref.parameters),
    rtol=1e-12,
    atol=1e-12,
  )


@pytest.mark.parametrize("var_types", [["d", "c", "c"], ["d", "d", "d"]])
def test_from_data_matches_discrete_vinecop(var_types: list[str]) -> None:
  # End to end: the torch TLL fit on data with atoms, against the compiled vine
  # fitted on the same data with the same structure. Measured 2.0e-14 / 5.6e-14.
  u = _discrete_data(var_types)
  cop = _discrete_vinecop(var_types, u)
  bc = TorchVinecop.from_data(
    torch.from_numpy(u),
    var_types=var_types,
    controls=FitControlsTorchVinecop(),
  )
  assert np.array_equal(
    np.asarray(bc.structure.matrix), np.asarray(cop.structure.matrix)
  )
  u_t = torch.from_numpy(u)
  np.testing.assert_allclose(
    bc.pdf(u_t).numpy(), cop.pdf(u), atol=1e-12, rtol=1e-12
  )
  np.testing.assert_allclose(
    bc.rosenblatt(u_t, randomize_discrete=False).numpy(),
    cop.rosenblatt(u, randomize_discrete=False),
    atol=1e-12,
    rtol=1e-12,
  )


def test_discrete_vine_declines_the_batched_path() -> None:
  # The batched level carries no distribution-function grid, which a discrete
  # edge's h-functions are difference quotients of, so an explicit
  # `batched=True` must fall back rather than evaluate something else.
  var_types = ["d", "c", "c"]
  u = _discrete_data(var_types, n=400, seed=8)
  bc = TorchVinecop.from_vinecop(_discrete_vinecop(var_types, u))
  u_t = torch.from_numpy(u)
  np.testing.assert_array_equal(
    bc.pdf(u_t, batched=True).numpy(), bc.pdf(u_t, batched=False).numpy()
  )


def test_discrete_vine_round_trips_through_pickle() -> None:
  # `var_types` is plain Python state on the vine, and the per-edge type table
  # is derived from it at bind time — so both have to survive a round trip, or a
  # reloaded vine would evaluate a continuous density on the same data.
  var_types = ["d", "c", "c"]
  u = _discrete_data(var_types, n=400, seed=9)
  bc = TorchVinecop.from_vinecop(_discrete_vinecop(var_types, u))
  clone = pickle.loads(pickle.dumps(bc))
  assert clone.var_types == var_types
  assert clone.pair_var_types(0, 0) == bc.pair_var_types(0, 0)
  u_t = torch.from_numpy(u)
  np.testing.assert_array_equal(clone.pdf(u_t).numpy(), bc.pdf(u_t).numpy())


def test_compile_flag_defaults_off_and_survives_a_round_trip() -> None:
  """`compile_cascades` is off unless asked for, and is a cache, not state."""
  u_fit = _simulate(d=4, n=400, seed=71)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  assert bc.compile_cascades is False
  assert TorchVinecop.from_data(
    torch.from_numpy(u_fit), controls=FitControlsTorchVinecop(compile=True)
  ).compile_cascades

  bc.compile_cascades = True
  # The state is everything but the cache -- `__getstate__` cannot lean on an
  # inherited one, which `nn.Module` only grew in torch 2.1 and `object` in
  # Python 3.11, below both of this project's floors.
  state = bc.__getstate__()
  assert state["_compiled"] == {}
  assert {"_parameters", "_buffers", "_modules"} <= set(state)
  # Nothing is compiled until a cascade runs, so this round-trips either way.
  again = pickle.loads(pickle.dumps(bc))
  assert again.compile_cascades is True
  torch.testing.assert_close(
    again.pdf(torch.from_numpy(_eval_grid(64, d=4, seed=72)), batched=False),
    bc.pdf(torch.from_numpy(_eval_grid(64, d=4, seed=72)), batched=False),
    atol=0.0,
    rtol=0.0,
  )


@pytest.mark.compile
@pytest.mark.parametrize("op", ("pdf", "rosenblatt", "inverse_rosenblatt"))
def test_compiled_cascades_match_eager(op: str) -> None:
  """Inductor fuses the cascade; it does not change what the cascade computes."""
  u_fit = _simulate(d=4, n=400, seed=73)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  u_t = torch.from_numpy(_eval_grid(128, d=4, seed=74))
  eager = getattr(bc, op)(u_t, batched=True)
  bc.compile_cascades = True
  torch.testing.assert_close(
    getattr(bc, op)(u_t, batched=True), eager, atol=1e-12, rtol=1e-12
  )


def test_pickles_after_an_evaluation() -> None:
  """Evaluating first must not make a vine unpicklable.

  The cascades memoize their array namespace -- a module, which no pickle can
  carry -- so the state has to drop it. Pickling before the first call cannot
  see this, which is what every other pickle test here does.
  """
  cop_tll = _fit_tll_vine(_simulate(d=4, n=400, seed=81))
  bc = TorchVinecop.from_vinecop(cop_tll)
  u_t = torch.from_numpy(_eval_grid(64, d=4, seed=82))
  bc.pdf(u_t)
  again = pickle.loads(pickle.dumps(bc))
  torch.testing.assert_close(again.pdf(u_t), bc.pdf(u_t), atol=0.0, rtol=0.0)


@pytest.mark.parametrize("op", ("pdf", "rosenblatt", "inverse_rosenblatt"))
def test_batched_handles_a_vine_mixing_indep_and_tll(op: str) -> None:
  """An independence pair has its own 2x2 grid; the level still stacks.

  `TorchBicop` gives an independence copula a two-point sentinel grid and no
  prefix tables, because none of its own methods read either. Stacking a tree
  level reads both, so the bake substitutes an independence density built on
  the shared grid -- otherwise a vine whose fit chose `indep` anywhere raises
  from `torch.stack`, which is most vines.
  """
  rng = np.random.RandomState(83)
  base = rng.normal(size=800)
  x = np.column_stack(
    [
      0.95 * base + 0.05 * rng.normal(size=800),
      0.95 * base + 0.05 * rng.normal(size=800),
      rng.normal(size=800),
      rng.normal(size=800),
    ]
  )
  u = pv.to_pseudo_obs(x)
  cop = pv.Vinecop.from_data(
    u,
    controls=pv.FitControlsVinecop(
      family_set=[pv.BicopFamily.tll, pv.BicopFamily.indep]
    ),
  )
  fams = {
    cop.get_pair_copula(t, e).family
    for t in range(cop.trunc_lvl)
    for e in range(cop.dim - 1 - t)
  }
  assert {pv.BicopFamily.tll, pv.BicopFamily.indep} <= fams

  bc = TorchVinecop.from_vinecop(cop)
  u_t = torch.from_numpy(u[:32])
  torch.testing.assert_close(
    getattr(bc, op)(u_t, batched=True),
    getattr(bc, op)(u_t, batched=False),
    atol=0.0,
    rtol=0.0,
  )


@pytest.mark.parametrize("first", ["sample", "cdf", "inverse_rosenblatt"])
def test_a_no_grad_cascade_does_not_detach_the_bake(first: str) -> None:
  """`sample` / `cdf` / the inverse evaluate under `no_grad`; `pdf` still fits.

  Those three bake the stacked grids inside `torch.no_grad()`, which copies
  them detached while the grids themselves still track grad -- so the flags
  the re-bake watches do not change and the detached copy would be kept.
  """
  cop_tll = _fit_tll_vine(_simulate(d=4, n=400, seed=84))
  bc = TorchVinecop.from_vinecop(cop_tll)
  values = bc._pair_module(0, 0).interp_grid.values
  values.requires_grad_(True)
  u_t = torch.from_numpy(_eval_grid(32, d=4, seed=85))
  if first == "sample":
    bc.sample(8, batched=True)
  elif first == "cdf":
    bc.cdf(u_t, N=64, batched=True)
  else:
    bc.inverse_rosenblatt(u_t, batched=True)
  grad = torch.autograd.grad(bc.pdf(u_t, batched=True).sum(), values)[0]
  assert torch.isfinite(grad).all()
  assert grad.abs().sum() > 0


@pytest.mark.parametrize("value", [0.0, 1.0])
def test_batched_inverse_matches_at_the_unit_square_boundary(
  value: float,
) -> None:
  """The waves agree with the per-edge cascade at the trimmed boundary too.

  Random interior points exercise neither clamp, and the two paths trim in
  different places -- the vine trims its input, each pair trims again.
  """
  cop_tll = _fit_tll_vine(_simulate(d=5, n=600, seed=86))
  bc = TorchVinecop.from_vinecop(cop_tll)
  u_t = torch.full((16, 5), value, dtype=torch.float64)
  torch.testing.assert_close(
    bc.inverse_rosenblatt(u_t, batched=True),
    bc.inverse_rosenblatt(u_t, batched=False),
    atol=0.0,
    rtol=0.0,
  )


@pytest.mark.parametrize("d", [5, 9])
def test_batched_fit_matches_the_per_edge_fit(d: int) -> None:
  """A level fitted in one call gives the vine a per-edge fit would.

  `batched_fit` is a scheduling choice, so it must not move the model --
  but it does move the last bits, on every device rather than only on CUDA.
  Stacking a level changes how many elements the bandwidth search's `pow`
  sees, and torch picks an elementwise kernel by element count. What a
  schedule must not touch is whether a lane's companions change its answer,
  and that is pinned exactly in `test_torch_bicop`.
  """
  u_fit = _simulate(d=d, n=1200, seed=500 + d)
  structure = _fit_tll_vine(u_fit).structure
  u_t = torch.from_numpy(u_fit)
  fits = {
    flag: TorchVinecop.from_data(
      u_t,
      structure,
      controls=FitControlsTorchVinecop(batched_fit=flag),
    )
    for flag in (False, True)
  }
  u_eval = torch.from_numpy(_eval_grid(200, d=d, seed=510 + d))
  for op in ("pdf", "rosenblatt"):
    torch.testing.assert_close(
      getattr(fits[True], op)(u_eval),
      getattr(fits[False], op)(u_eval),
      atol=1e-11,
      rtol=1e-9,
    )


def test_batched_fit_defaults_to_off_on_cpu() -> None:
  """The default follows the device, as the evaluation cascade's does.

  Batching amortizes kernel launches, and cpu has none worth amortizing --
  measured slower there at one thread. Pinned because a default that
  silently flipped would change every cpu user's fit timing.
  """
  assert FitControlsTorchVinecop().batched_fit is None
  u_fit = _simulate(d=4, n=400, seed=61)
  structure = _fit_tll_vine(u_fit).structure
  seen: list[int] = []
  real = TorchBicop.from_data_batched

  def spy(u: torch.Tensor, *args: Any, **kwargs: Any) -> Any:
    seen.append(int(u.shape[0]))
    return real(u, *args, **kwargs)

  # The resolution is only observable through whether the level fitter runs,
  # so watch that rather than the field that was left unset.
  with pytest.MonkeyPatch.context() as mp:
    mp.setattr(TorchBicop, "from_data_batched", spy)
    TorchVinecop.from_data(torch.from_numpy(u_fit), structure)
    assert seen == [], "cpu resolved batched_fit on"
    TorchVinecop.from_data(
      torch.from_numpy(u_fit),
      structure,
      controls=FitControlsTorchVinecop(batched_fit=True),
    )
    assert seen == [3, 2, 1], f"explicit batched_fit=True not honored: {seen}"


def test_batched_fit_falls_back_for_a_discrete_level() -> None:
  """A discrete edge cannot stack, so its level stays per-edge.

  Asking for `batched_fit=True` on a discrete vine must fit it, not raise:
  the hook is an optimization, and a level it cannot serve is simply not
  handed to it.

  Exact rather than toleranced because a discrete level cannot stack at all,
  so `batched_fit=True` never reaches the batched fitter here and the two
  runs are the same arithmetic.
  """
  var_types = ["d", "c", "c"]
  wide = _discrete_data(var_types, n=600, seed=17)
  structure = pv.RVineStructure.sample(3, seeds=[2])
  fits = {
    flag: TorchVinecop.from_data(
      torch.from_numpy(wide),
      structure,
      controls=FitControlsTorchVinecop(batched_fit=flag),
      var_types=var_types,
    )
    for flag in (False, True)
  }
  u_eval = torch.from_numpy(wide[:64])
  torch.testing.assert_close(
    fits[True].pdf(u_eval), fits[False].pdf(u_eval), atol=0.0, rtol=0.0
  )


@pytest.mark.parametrize("d", [5, 9])
def test_batched_fit_selects_the_same_structure(d: int) -> None:
  """Selection is unmoved by how its pairs were scheduled.

  `structure=None` is the default path, so this is the common case rather
  than a corner.

  The matrix is compared exactly on every device even though the fit itself
  is not bit-identical anywhere, and what licenses that is rank
  quantization rather than the fit: the tree criterion is a function of
  ranks alone, so perturbing an h-value by ~1e-15 can only change a tau by
  swapping two neighbors, and their spacing at this `n` is ~1e-3. The same
  file already asserts exact matrix equality across a ~1e-11 pair-fit
  difference, four orders more perturbation than scheduling applies. Only
  the pdf comparison carries a tolerance.
  """
  u_fit = _simulate(d=d, n=1200, seed=700 + d)
  u_t = torch.from_numpy(u_fit)
  fits = {
    flag: TorchVinecop.from_data(
      u_t, controls=FitControlsTorchVinecop(batched_fit=flag)
    )
    for flag in (False, True)
  }
  np.testing.assert_array_equal(
    np.asarray(fits[True].structure.matrix),
    np.asarray(fits[False].structure.matrix),
  )
  # And still the structure the compiled selector picks.
  np.testing.assert_array_equal(
    np.asarray(fits[True].structure.matrix),
    np.asarray(_fit_tll_vine(u_fit, trunc_lvl=20).structure.matrix),
  )
  u_eval = torch.from_numpy(_eval_grid(200, d=d, seed=710 + d))
  torch.testing.assert_close(
    fits[True].pdf(u_eval), fits[False].pdf(u_eval), atol=1e-11, rtol=1e-9
  )


@pytest.mark.parametrize(
  ("make_structure", "label"),
  [
    (lambda d: pv.CVineStructure(list(range(1, d + 1))), "c-vine"),
    (lambda d: pv.CVineStructure(list(range(d, 0, -1))), "c-vine reversed"),
    (lambda d: pv.DVineStructure(list(range(1, d + 1))), "d-vine"),
    (lambda d: pv.RVineStructure.sample(d, seeds=[11]), "r-vine sampled"),
    (lambda d: pv.RVineStructure.sample(d, seeds=[29]), "r-vine sampled 2"),
  ],
)
def test_fit_matches_cpp_on_a_given_structure(
  make_structure: Callable[[int], Any], label: str
) -> None:
  """The gathered level indexes correctly on structures selection never picks.

  `VinecopBase.fit` gathers a whole tree level's inputs before fitting any of
  it, which is licensed by every edge of a tree reading only columns
  finalized by earlier trees. That holds -- checked over every C-vine and
  D-vine order permutation up to d = 5 and 660 sampled R-vines up to d = 12,
  19888 slots, no violation -- so what is at risk is not the claim but this
  port's indexing of it, on a structure shape the suite does not otherwise
  produce.

  Comparing `batched_fit=True` against `False` cannot catch that: both read
  the same gathered `inputs`, so a mis-gathered level is mis-gathered
  identically for both and every such test still passes. The independent
  check is C++, which traverses per node rather than per column, so the two
  implementations of the same indexing disagree where one is wrong.

  Every other fixed-structure parity test uses a Dissmann-selected structure
  from one generator. These five are chosen rather than generated, and the
  C-vine and D-vine are *complementary* rather than redundant -- measured by
  injecting three mis-indexings into the gather and seeing which cases the
  1e-9 gate then catches:

  =================  ==========  ==========  ==========
  case               shifted     wrong h     swapped
  =================  ==========  ==========  ==========
  c-vine (both)      caught      **missed**  caught
  d-vine             **missed**  caught      caught
  r-vine (2 seeds)   caught      caught      caught
  =================  ==========  ==========  ==========

  So a C-vine, whose tree 0 is a hub, cannot see a wrong h-function, and a
  D-vine, whose tree 0 is a path, cannot see a column shifted by one. Drop
  either as duplicative and the suite acquires a blind spot; the sampled
  R-vines catch all three but are not a substitute, being one shape each.
  """
  del label
  d = 6
  u_fit = _simulate(d=d, n=700, seed=920)
  structure = make_structure(d)
  cpp = pv.Vinecop.from_data(u_fit, controls=_TLL_CONTROLS, structure=structure)
  fitted = TorchVinecop.from_data(
    torch.from_numpy(u_fit),
    structure,
    controls=FitControlsTorchVinecop(cache_integrals=False),
  )
  lifted = TorchVinecop.from_vinecop(cpp, cache_integrals=False)
  u_eval = torch.from_numpy(_eval_grid(200, d=d, seed=921))
  np.testing.assert_allclose(
    fitted.pdf(u_eval).numpy(),
    lifted.pdf(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )


def test_mixed_grids_refuse_the_batched_path() -> None:
  """Pairs that disagree on their grid fall back rather than stack.

  Every stacked pair is read against one level's knots, so pairs must agree
  on the grid itself and not merely on its size: two grids of one size and
  different spacing interpolate to different functions. Agreeing on the
  shape alone let a level be evaluated on the wrong one, silently -- the two
  cascades disagreed by 2.0 in absolute density on the vine below.
  """
  rng = np.random.default_rng(3)
  d, n = 3, 400
  u = pv.to_pseudo_obs(
    0.6 * rng.standard_normal((n, 1)) + 0.5 * rng.standard_normal((n, d))
  )

  def pair(a: int, b: int, grid_type: str) -> TorchBicop:
    return TorchBicop.from_data(
      np.column_stack([u[:, a], u[:, b]]),
      FitControlsTorchBicop(grid_type=grid_type, grid_size=30),
    )

  # Same grid_size, so the shapes match and only the spacing differs.
  vine = TorchVinecop(
    pair_copulas=[
      [pair(0, 1, "normal"), pair(1, 2, "linear")],
      [pair(0, 2, "normal")],
    ],
    structure=pv.DVineStructure([1, 2, 3]),
  )
  u_eval = torch.from_numpy(_eval_grid(64, d=d, seed=4))
  torch.testing.assert_close(
    vine.pdf(u_eval, batched=True),
    vine.pdf(u_eval, batched=False),
    atol=0.0,
    rtol=0.0,
  )


def test_load_state_dict_drops_the_stacked_bake() -> None:
  """Loading new grids must not leave the batched path on the old ones.

  The bake and the compiled cascades copy the grids rather than viewing
  them, so a load that replaces the grids leaves both answering from the
  density they were built with. Evaluating first is the point: it is what
  builds the bake that the load then has to invalidate.
  """
  d, n = 4, 400
  strong = _simulate(d=d, n=n, seed=5)
  weak = pv.to_pseudo_obs(np.random.default_rng(6).standard_normal((n, d)))
  structure = _fit_tll_vine(strong).structure
  controls = FitControlsTorchVinecop(cache_integrals=False)
  vine = TorchVinecop.from_data(
    torch.from_numpy(strong), structure, controls=controls
  )
  other = TorchVinecop.from_data(
    torch.from_numpy(weak), structure, controls=controls
  )
  u_eval = torch.from_numpy(_eval_grid(64, d=d, seed=7))
  vine.pdf(u_eval, batched=True)
  vine.load_state_dict(other.state_dict())
  torch.testing.assert_close(
    vine.pdf(u_eval, batched=True),
    other.pdf(u_eval, batched=True),
    atol=0.0,
    rtol=0.0,
  )


def test_cxi_criterion_thresholds_nothing_at_zero() -> None:
  """A negative criterion must not fall below a threshold of zero.

  `calculate_criterion` returns `std::fabs(w)` for *every* branch, `cxi`
  included. This port applied the absolute value inside the branches and the
  `cxi` one -- `max(xi12, xi21)`, mirroring `pairwise_cxi` -- did not, so a
  criterion could be negative. Chatterjee's xi routinely is on weak
  dependence: symmetrized, it came out negative in 102 of 400 independent
  samples.

  That was harmless while the criterion only set a spanning-tree weight, and
  stopped being harmless once a criterion below the threshold also decides
  whether an edge is fitted at all. With the default `threshold=0.0`, a
  negative criterion is below it, so edges silently became independence: on
  independent uniforms this thresholded something in 16 of 18 fixtures, up to
  7 pairs in one d = 6 vine.

  This fixture is the worst of those. Independent data is the point -- it is
  where xi has no signal and sits astride zero.
  """
  rng = np.random.default_rng(1000 + 17 * 2 + 6)
  u = rng.uniform(size=(400, 6))
  controls = FitControlsTorchVinecop(
    trunc_lvl=20, tree_criterion="cxi", threshold=0.0
  )
  fitted = TorchVinecop.from_data(torch.as_tensor(u), controls=controls)
  # A thresholded edge becomes the grid that is exactly one everywhere; a
  # fitted one never is, however weak the dependence.
  # `nn.Module.__getattr__` returns a union no checker can narrow.
  rows: Any = fitted.pair_copulas
  thresholded = [
    (t, e)
    for t in range(fitted.trunc_lvl)
    for e in range(6 - t - 1)
    if bool((rows[t][e].interp_grid.values == 1.0).all())
  ]
  assert thresholded == [], f"thresholded at threshold=0.0: {thresholded}"

  cpp = pv.Vinecop.from_data(
    u,
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.tll],
      num_threads=1,
      trunc_lvl=20,
      tree_criterion="cxi",
      threshold=0.0,
    ),
  )
  np.testing.assert_array_equal(
    np.asarray(fitted.structure.matrix), np.asarray(cpp.structure.matrix)
  )
  u_eval = rng.uniform(0.05, 0.95, size=(120, 6))
  np.testing.assert_allclose(
    fitted.pdf(torch.as_tensor(u_eval)).numpy(),
    cpp.pdf(u_eval),
    rtol=1e-9,
    atol=1e-11,
  )


@pytest.mark.parametrize("threshold", [0.0, 0.3, 0.5, 0.95])
def test_threshold_matches_pvvinecop(threshold: float) -> None:
  """A thresholded edge holds independence, as the compiled selector's does.

  `threshold` decides two things upstream: which edges the spanning tree
  prefers, and which of the survivors are fitted at all -- an edge whose
  criterion falls below it keeps a default-constructed pair, which is the
  independence copula. Only the first was ported, so a threshold that bites
  used to give a different model rather than a differently-weighted one.

  The four values are chosen against this fixture's candidate taus, which
  run from 0.47 to 0.52: 0.0 thresholds nothing and is the path every other
  test takes, 0.3 and 0.5 straddle the candidates, and 0.95 is above all of
  them so every pair is independence and both vines are exactly `pdf = 1`.
  """
  u = _simulate(d=6, n=800, seed=81)
  cpp = pv.Vinecop.from_data(
    u,
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.tll],
      num_threads=1,
      trunc_lvl=20,
      threshold=threshold,
    ),
  )
  fitted = TorchVinecop.from_data(
    torch.from_numpy(u),
    controls=FitControlsTorchVinecop(trunc_lvl=20, threshold=threshold),
  )
  np.testing.assert_array_equal(
    np.asarray(fitted.structure.matrix), np.asarray(cpp.structure.matrix)
  )
  u_eval = torch.from_numpy(_eval_grid(120, d=6, seed=82))
  np.testing.assert_allclose(
    fitted.pdf(u_eval).numpy(),
    cpp.pdf(u_eval.numpy()),
    rtol=1e-9,
    atol=1e-11,
  )


def test_batched_fit_runs_one_call_per_tree(
  monkeypatch: pytest.MonkeyPatch,
) -> None:
  """The batched path is taken, once per tree level, with every edge in it.

  Every other batched-fit test compares the batched result against the
  per-edge one, which a batched path that silently fell back would pass
  too: substituting a plain per-lane loop for `from_data_batched` leaves
  all of them green at exactly zero difference. So nothing else here
  separates a working batched fit from a dead one.

  This does, and it compares no floating-point value -- so it is immune to
  the batch-extent arithmetic that put a tolerance on the others. The lane
  counts also say the levels arrive whole, which is what
  `test_batched_selection_survives_a_shortened_cascade` claims and cannot
  see.
  """
  d = 6
  u_fit = _simulate(d=d, n=400, seed=91)
  structure = _fit_tll_vine(u_fit).structure
  seen: list[int] = []
  real = TorchBicop.from_data_batched

  def spy(u: torch.Tensor, *args: Any, **kwargs: Any) -> Any:
    seen.append(int(u.shape[0]))
    return real(u, *args, **kwargs)

  monkeypatch.setattr(TorchBicop, "from_data_batched", spy)
  TorchVinecop.from_data(
    torch.from_numpy(u_fit),
    structure,
    controls=FitControlsTorchVinecop(batched_fit=True),
  )
  assert seen == list(range(d - 1, 0, -1))

  # `structure=None` routes through `select`, which installs the hook
  # separately -- dropping it there left every select-path test green.
  seen.clear()
  TorchVinecop.from_data(
    torch.from_numpy(u_fit),
    controls=FitControlsTorchVinecop(batched_fit=True),
  )
  assert seen == list(range(d - 1, 0, -1)), f"select did not batch: {seen}"


@pytest.mark.parametrize(
  ("kwargs", "why"),
  [
    ({"threshold": 0.5}, "a threshold the candidates straddle"),
    ({"threshold": 0.95}, "a threshold above every candidate"),
    ({"trunc_lvl": 2}, "truncation shortens the cascade"),
  ],
)
def test_batched_selection_survives_a_shortened_cascade(
  kwargs: dict, why: str
) -> None:
  """The batched selector holds where a level could come out ragged.

  Both knobs change how many trees get built, and thresholding changes the
  edge weights that decide which edges survive -- so both are places a
  level could arrive with a different number of edges than the unbatched
  path saw, or with none. That the levels line up is what
  `test_batched_fit_runs_one_call_per_tree` pins; this checks the model they
  produce.

  Both thresholds are here because they exercise different things. At 0.95
  no candidate tau reaches the threshold for this fixture, so every MST
  weight is exactly 1.0 and the structure is dependence-blind -- which
  makes its matrix comparison hold for reasons unrelated to scheduling. At
  0.5 the candidates straddle it, so the weights carry dependence and the
  comparison has something to say.
  """
  del why
  u_fit = _simulate(d=6, n=800, seed=81)
  u_t = torch.from_numpy(u_fit)
  fits = {
    flag: TorchVinecop.from_data(
      u_t, controls=FitControlsTorchVinecop(batched_fit=flag, **kwargs)
    )
    for flag in (False, True)
  }
  assert fits[True].trunc_lvl == fits[False].trunc_lvl
  np.testing.assert_array_equal(
    np.asarray(fits[True].structure.matrix),
    np.asarray(fits[False].structure.matrix),
  )
  u_eval = torch.from_numpy(_eval_grid(120, d=6, seed=82))
  torch.testing.assert_close(
    fits[True].pdf(u_eval), fits[False].pdf(u_eval), atol=1e-11, rtol=1e-9
  )
