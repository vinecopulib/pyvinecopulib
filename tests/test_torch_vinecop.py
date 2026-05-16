"""Tests for the optional ``pyvinecopulib.torch.TorchVinecop`` class.

Skipped when PyTorch isn't installed. Compares the torch ``TorchVinecop``
against the C++ ``pv.Vinecop`` (with TLL pair copulas) for ``pdf`` /
``rosenblatt`` / ``inverse_rosenblatt``; verifies the round-trip
``inverse_rosenblatt(rosenblatt(u)) ≈ u``; spot-checks an
independent vine; and confirms truncation and discrete-rejection paths.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyvinecopulib as pv

torch = pytest.importorskip("torch")

from pyvinecopulib.torch import TorchVinecop  # noqa: E402

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

  bc = TorchVinecop.from_vinecop(cop_tll)
  u_eval = _eval_grid(500, d=5, seed=11)

  out_torch = bc.pdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_rosenblatt_matches_pvvinecop() -> None:
  u_fit = _simulate(d=5, n=2000, seed=2)
  cop_tll = _fit_tll_vine(u_fit)

  bc = TorchVinecop.from_vinecop(cop_tll)
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
    # cache=True replaces every integration / inversion with a bilinear
    # interp on a 30x30 grid — the per-pair ~1e-3 mean error compounds
    # through the d-1 cascade levels. pdf is unbounded so the *max* error
    # in absolute terms can be O(1); rosenblatt / inverse_rosenblatt are
    # u-space and bounded so their max stays around the bicop-level cap.
    # These bounds are loose — they pin the floor, not optimize against it.
    ("pdf", True, 1e1, 1e-1),
    ("rosenblatt", True, 1.0, 5e-2),
    ("inverse_rosenblatt", True, 5e-2, 5e-3),
  ],
)
def test_cache_integrals_precision_floor(
  method: str, cache: bool, max_abs: float, mean_abs: float
) -> None:
  """Pin the precision floor for each cache mode at a moderately-deep vine.

  ``cache_integrals=False`` is the right choice for likelihood / sampling
  applications: every per-pair integration and inversion is exact (modulo
  bisection convergence), so the cascade preserves precision. Setting
  ``cache=True`` is the right choice when ~1e-3 u-space error is
  acceptable and per-call speed matters — typically pdf-only workloads
  where the bilinear-interp gap is small compared to downstream noise.
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

  bc = TorchVinecop.from_vinecop(cop_tll)
  # Inverse-rosenblatt takes independent uniforms; sample fresh.
  w = _eval_grid(400, d=5, seed=13)

  out_torch = bc.inverse_rosenblatt(torch.from_numpy(w)).numpy()
  out_cpp = cop_tll.inverse_rosenblatt(w)
  # Bisection at n_iter=35 gives bracket-width accuracy 0.5**35 ≈ 6e-11.
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-9, rtol=1e-9)


def test_inverse_rosenblatt_roundtrip() -> None:
  u_fit = _simulate(d=5, n=2000, seed=4)
  cop_tll = _fit_tll_vine(u_fit)

  bc = TorchVinecop.from_vinecop(cop_tll)
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

  bc = TorchVinecop.from_vinecop(cop)
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

  bc = TorchVinecop.from_vinecop(cop_tll)
  assert bc.trunc_lvl == 2

  u_eval = _eval_grid(300, d=6, seed=15)
  np.testing.assert_allclose(
    bc.pdf(torch.from_numpy(u_eval), batched=batched).numpy(),
    cop_tll.pdf(u_eval),
    atol=1e-10,
    rtol=1e-10,
  )


def test_from_vinecop_rejects_discrete_var_types() -> None:
  d = 3
  structure = pv.RVineStructure.from_order([1, 2, 3])
  pair_copulas = [
    [pv.Bicop(family=pv.families.indep) for _ in range(d - 1 - t)]
    for t in range(d - 1)
  ]
  cop = pv.Vinecop.from_structure(
    structure=structure,
    pair_copulas=pair_copulas,
    var_types=["c", "d", "c"],
  )
  with pytest.raises(ValueError, match="continuous-only"):
    TorchVinecop.from_vinecop(cop)


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
# `impl="lazy"` equivalence: dict-based bookkeeping must give identical
# results to the legacy impl (same math, different bookkeeping).
# --------------------------------------------------------------------- #


def test_pdf_lazy_matches_legacy() -> None:
  u_fit = _simulate(d=5, n=2000, seed=101)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  u_eval = _eval_grid(500, d=5, seed=111)
  u_t = torch.from_numpy(u_eval)
  out_legacy = bc.pdf(u_t, impl="legacy").numpy()
  out_lazy = bc.pdf(u_t, impl="lazy").numpy()
  np.testing.assert_allclose(out_lazy, out_legacy, atol=1e-12, rtol=1e-12)


def test_rosenblatt_lazy_matches_legacy() -> None:
  u_fit = _simulate(d=5, n=2000, seed=102)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  u_eval = _eval_grid(500, d=5, seed=112)
  u_t = torch.from_numpy(u_eval)
  out_legacy = bc.rosenblatt(u_t, impl="legacy").numpy()
  out_lazy = bc.rosenblatt(u_t, impl="lazy").numpy()
  np.testing.assert_allclose(out_lazy, out_legacy, atol=1e-12, rtol=1e-12)


def test_inverse_rosenblatt_lazy_matches_legacy() -> None:
  u_fit = _simulate(d=5, n=2000, seed=103)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  rng = np.random.default_rng(113)
  w = rng.uniform(0.05, 0.95, size=(400, 5))
  w_t = torch.from_numpy(w)
  out_legacy = bc.inverse_rosenblatt(w_t, impl="legacy").numpy()
  out_lazy = bc.inverse_rosenblatt(w_t, impl="lazy").numpy()
  # Both paths invoke the same bisection root-finder, just in different
  # iteration orders; agree to machine precision modulo a few ULPs.
  np.testing.assert_allclose(out_lazy, out_legacy, atol=1e-12, rtol=1e-12)


def test_invalid_impl_raises() -> None:
  u_fit = _simulate(d=3, n=500, seed=104)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  u_t = torch.from_numpy(_eval_grid(20, d=3, seed=114))
  for fn in (bc.pdf, bc.rosenblatt, bc.inverse_rosenblatt):
    with pytest.raises(ValueError, match="impl must be"):
      fn(u_t, impl="bogus")


# --------------------------------------------------------------------- #
# Batched cascade                                                        #
# --------------------------------------------------------------------- #


@pytest.mark.parametrize("impl", ["legacy", "lazy"])
@pytest.mark.parametrize("cache", [False, True])
def test_pdf_batched_matches_legacy(impl: str, cache: bool) -> None:
  u_fit = _simulate(d=10, n=800, seed=300)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=cache)
  u_t = torch.from_numpy(_eval_grid(500, d=10, seed=310))
  ref = bc.pdf(u_t, impl=impl, batched=False).numpy()
  got = bc.pdf(u_t, impl=impl, batched=True).numpy()
  np.testing.assert_allclose(got, ref, atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("impl", ["legacy", "lazy"])
@pytest.mark.parametrize("cache", [False, True])
def test_rosenblatt_batched_matches_legacy(impl: str, cache: bool) -> None:
  u_fit = _simulate(d=10, n=800, seed=301)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=cache)
  u_t = torch.from_numpy(_eval_grid(500, d=10, seed=311))
  ref = bc.rosenblatt(u_t, impl=impl, batched=False).numpy()
  got = bc.rosenblatt(u_t, impl=impl, batched=True).numpy()
  np.testing.assert_allclose(got, ref, atol=1e-13, rtol=1e-13)


@pytest.mark.parametrize("impl", ["legacy", "lazy"])
def test_inverse_rosenblatt_batched_falls_back(impl: str) -> None:
  """``batched=True`` on inverse_rosenblatt routes to legacy in v1
  (cross-tree dependencies block a clean tree-level batched traversal).
  Output must be bit-equal to the non-batched path."""
  u_fit = _simulate(d=6, n=600, seed=302)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  w_t = torch.from_numpy(_eval_grid(300, d=6, seed=312))
  ref = bc.inverse_rosenblatt(w_t, impl=impl, batched=False).numpy()
  got = bc.inverse_rosenblatt(w_t, impl=impl, batched=True).numpy()
  np.testing.assert_array_equal(got, ref)


def test_batched_matches_cpp_pdf() -> None:
  """End-to-end: torch batched pdf vs pv.Vinecop on the same fit."""
  u_fit = _simulate(d=8, n=1000, seed=320)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  u_eval = _eval_grid(400, d=8, seed=330)
  out_torch = bc.pdf(torch.from_numpy(u_eval), batched=True).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_batched_matches_cpp_rosenblatt() -> None:
  u_fit = _simulate(d=8, n=1000, seed=321)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  u_eval = _eval_grid(400, d=8, seed=331)
  out_torch = bc.rosenblatt(torch.from_numpy(u_eval), batched=True).numpy()
  out_cpp = cop_tll.rosenblatt(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


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
