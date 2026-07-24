"""Tests for the optional ``pyvinecopulib.torch`` submodule.

Skipped when PyTorch isn't installed. Compares the torch ``TorchBicop``
against the C++ ``pv.Bicop`` (TLL family) for ``pdf`` / ``cdf`` / ``hfunc``
on the same fitted interpolation grid and verifies ``hinv`` round-trips.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyvinecopulib as pv

torch = pytest.importorskip("torch")

from pyvinecopulib.torch import FitControlsTorchBicop, TorchBicop  # noqa: E402

_TLL_CONTROLS = pv.FitControlsBicop(family_set=[pv.families.tll], num_threads=1)


def _fit_tll(u: np.ndarray) -> pv.Bicop:
  return pv.Bicop.from_data(u, controls=_TLL_CONTROLS)


def _eval_grid(n: int, seed: int = 0) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return rng.uniform(0.02, 0.98, size=(n, 2))


def test_pdf_matches_pvbicop() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  # Pin cache=False: this test verifies parity with the C++ on-the-fly
  # integration math at 1e-10. The default cache_integrals=True uses
  # bilinear-interp on the cache, which intentionally introduces a
  # ~1e-3 IAE gap and is covered by test_cached_*_smoke below.
  bc = TorchBicop.from_bicop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(500, seed=11)

  out_torch = bc.pdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_cdf_matches_pvbicop() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  # Pin cache=False — C++ parity at 1e-10; see test_pdf_matches_pvbicop.
  bc = TorchBicop.from_bicop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(500, seed=12)

  out_torch = bc.cdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.cdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_hfunc_matches_pvbicop() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  # Pin cache=False — C++ parity at 1e-10; see test_pdf_matches_pvbicop.
  bc = TorchBicop.from_bicop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(500, seed=13)
  u_t = torch.from_numpy(u_eval)

  np.testing.assert_allclose(
    bc.hfunc1(u_t).numpy(), cop_tll.hfunc1(u_eval), atol=1e-10, rtol=1e-10
  )
  np.testing.assert_allclose(
    bc.hfunc2(u_t).numpy(), cop_tll.hfunc2(u_eval), atol=1e-10, rtol=1e-10
  )


def test_hinv_roundtrip() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.7]]))
  u_fit = cop.simulate(2000, seeds=[4, 5, 6])
  cop_tll = _fit_tll(u_fit)

  # Pin cache=False: the closed-form inversion is the exact inverse of the
  # on-the-fly h-function, so the round-trip holds to machine precision.
  # The cached path round-trips only to ~1e-3 — covered by
  # test_cached_hinv_speedup.
  bc = TorchBicop.from_bicop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(400, seed=21)
  u_t = torch.from_numpy(u_eval)

  u2 = bc.hinv1(u_t).unsqueeze(-1)
  back = bc.hfunc1(torch.cat([u_t[:, 0:1], u2], dim=-1)).numpy()
  np.testing.assert_allclose(back, u_eval[:, 1], atol=1e-12, rtol=1e-12)

  u1 = bc.hinv2(u_t).unsqueeze(-1)
  back = bc.hfunc2(torch.cat([u1, u_t[:, 1:2]], dim=-1)).numpy()
  np.testing.assert_allclose(back, u_eval[:, 0], atol=1e-12, rtol=1e-12)


def test_hinv_closed_form_matches_cpp() -> None:
  """The no-cache torch h-inverse ports the C++ closed-form conditional
  quantile (vinecopulib#691), so torch and C++ agree to machine precision
  on the same fitted grid."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.7]]))
  u_fit = cop.simulate(2000, seeds=[4, 5, 6])
  cop_tll = _fit_tll(u_fit)

  bc = TorchBicop.from_bicop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(400, seed=22)
  u_t = torch.from_numpy(u_eval)

  np.testing.assert_allclose(
    bc.hinv1(u_t).numpy(), cop_tll.hinv1(u_eval), atol=1e-12, rtol=1e-12
  )
  np.testing.assert_allclose(
    bc.hinv2(u_t).numpy(), cop_tll.hinv2(u_eval), atol=1e-12, rtol=1e-12
  )


def test_inverse_integrate_1d() -> None:
  """Unit test of the closed-form conditional quantile on the grid itself:
  exact inverse of ``integrate_1d``, NaN propagation, shape validation."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.simulate(1500, seeds=[7, 8, 9])
  grid = TorchBicop.from_bicop(
    _fit_tll(u_fit), cache_integrals=False
  ).interp_grid

  u_eval = _eval_grid(300, seed=23)
  u_t = torch.from_numpy(u_eval)
  for cond_var in (1, 2):
    x = grid.inverse_integrate_1d(u_t, cond_var)
    assert ((x >= 0) & (x <= 1)).all()
    if cond_var == 1:
      u_back = torch.stack([u_t[:, 0], x], dim=-1)
      p = u_eval[:, 1]
    else:
      u_back = torch.stack([x, u_t[:, 1]], dim=-1)
      p = u_eval[:, 0]
    np.testing.assert_allclose(
      grid.integrate_1d(u_back, cond_var).numpy(), p, atol=1e-12, rtol=1e-12
    )

  # NaN rows return NaN (mirroring the C++ binaryExpr_or_nan wrapper).
  u_nan = u_t.clone()
  u_nan[0, 0] = torch.nan
  u_nan[1, 1] = torch.nan
  out = grid.inverse_integrate_1d(u_nan, 1)
  assert out[:2].isnan().all() and out[2:].isfinite().all()

  with pytest.raises(ValueError, match="shape"):
    grid.inverse_integrate_1d(u_t[:, 0], 1)


def test_from_bicop_rejects_rotated() -> None:
  """TLL pair-copulas in pyvinecopulib always have rotation=0; the
  TorchBicop wrapper enforces this at construction."""

  class _FakeCop:
    family = _fit_tll(_eval_grid(100, seed=0)).family
    rotation = 90
    parameters = np.eye(2)

  with pytest.raises(ValueError, match="rotation"):
    TorchBicop.from_bicop(_FakeCop())  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# TorchBicop.from_data — pure-torch TLL fit                                    #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("n", [500, 2000])
@pytest.mark.parametrize("rho", [0.3, 0.6, 0.9])
def test_from_data_matches_cpp(n: int, rho: float) -> None:
  """The pure-torch TLL constant fit produces the same density grid as
  ``pv.Bicop.from_data`` to machine precision after the standard
  ``normalize_margins(3)`` round in :class:`InterpolationGrid2D`.
  """
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[rho]]))
  u_np = cop.simulate(n, seeds=[1, 2, 3])
  cop_cpp = pv.Bicop.from_data(
    u_np,
    controls=pv.FitControlsBicop(family_set=[pv.families.tll], num_threads=1),
  )
  bc_torch = TorchBicop.from_data(u_np)
  np.testing.assert_allclose(
    bc_torch.interp_grid.values.numpy(),
    cop_cpp.parameters,
    atol=1e-11,
    rtol=1e-11,
  )


def test_from_data_evaluates_consistently() -> None:
  """Fit on data, evaluate pdf/cdf/hfunc/hinv — sanity checks: finite,
  in-range, hinv round-trips. Doesn't pin to a reference; the
  matches-vs-cpp test covers the fit accuracy.

  Pin cache_integrals=False so the hinv round-trip below is the exact
  closed-form inverse; the cached path round-trips only to ~1e-3."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_np = cop.simulate(1000, seeds=[11, 22, 33])
  bc = TorchBicop.from_data(u_np, cache_integrals=False)

  u_eval = _eval_grid(300, seed=99)
  u_t = torch.from_numpy(u_eval)

  pdf = bc.pdf(u_t)
  assert pdf.isfinite().all() and (pdf > 0).all()
  cdf = bc.cdf(u_t)
  assert ((cdf > 0) & (cdf < 1)).all()
  h1 = bc.hfunc1(u_t)
  h2 = bc.hfunc2(u_t)
  assert ((h1 >= 0) & (h1 <= 1)).all()
  assert ((h2 >= 0) & (h2 <= 1)).all()

  # hinv round-trip
  u2 = bc.hinv1(u_t).unsqueeze(-1)
  back = bc.hfunc1(torch.cat([u_t[:, 0:1], u2], dim=-1)).numpy()
  np.testing.assert_allclose(back, u_eval[:, 1], atol=1e-9, rtol=1e-9)


def test_cached_integrals_smoke() -> None:
  """Sanity check that ``cache_integrals=True`` produces finite outputs in
  range. We don't pin it to the trapezoidal path — the two paths agree only
  at grid nodes and the off-node gap depends on the local curvature.
  """
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.simulate(2000, seeds=[10, 11, 12])
  cop_tll = _fit_tll(u_fit)

  bc_cache = TorchBicop.from_bicop(cop_tll, cache_integrals=True)
  # All five caches must be populated for a non-indep pair.
  assert bc_cache._cdf_cache is not None
  assert bc_cache._hfunc1_cache is not None
  assert bc_cache._hfunc2_cache is not None
  assert bc_cache._hinv1_cache is not None
  assert bc_cache._hinv2_cache is not None

  u_t = torch.from_numpy(_eval_grid(200, seed=42))
  for fn in (
    bc_cache.cdf,
    bc_cache.hfunc1,
    bc_cache.hfunc2,
    bc_cache.hinv1,
    bc_cache.hinv2,
  ):
    out = fn(u_t)
    assert torch.isfinite(out).all()
    assert (out >= 0.0).all() and (out <= 1.0).all()


def test_cached_hinv_speedup() -> None:
  """``cache_integrals=True`` should make ``hinv1`` / ``hinv2`` a single
  bilinear interp instead of 35 iters of bisection. Verify the cached
  path matches the bisection path to within bilinear-interp precision
  (~1e-2 max, the usual gap between the cached integral and the
  trapezoidal recomputation at off-node points).
  """
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  bc_bisect = TorchBicop.from_bicop(cop_tll, cache_integrals=False)
  bc_cached = TorchBicop.from_bicop(cop_tll, cache_integrals=True)

  u_t = torch.from_numpy(_eval_grid(500, seed=77))
  for which in ("hinv1", "hinv2"):
    out_bisect = getattr(bc_bisect, which)(u_t).numpy()
    out_cached = getattr(bc_cached, which)(u_t).numpy()
    diff = np.abs(out_bisect - out_cached)
    # Bilinear-interp gap between the cached hinv and the bisection-on-
    # cached-h path (the latter converges to ~1e-9 of h's cache; the
    # former lives off-grid).
    assert diff.max() < 2e-2, f"{which}: max diff {diff.max():.3e}"
    assert diff.mean() < 3e-3, f"{which}: mean diff {diff.mean():.3e}"


def test_independent_bicop() -> None:
  bc = TorchBicop()  # default: independence
  u = torch.tensor([[0.2, 0.7], [0.5, 0.5], [0.9, 0.1]], dtype=torch.float64)
  assert torch.allclose(bc.pdf(u), torch.ones(3, dtype=torch.float64))
  assert torch.allclose(bc.cdf(u), u[:, 0] * u[:, 1])
  assert torch.allclose(bc.hfunc1(u), u[:, 1])
  assert torch.allclose(bc.hfunc2(u), u[:, 0])


def test_simulate_smoke() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.simulate(2000, seeds=[20, 21, 22])
  cop_tll = _fit_tll(u_fit)
  bc = TorchBicop.from_bicop(cop_tll)

  samples = bc.simulate(n=1000, qrng=False, seeds=[0])
  assert samples.shape == (1000, 2)
  assert torch.isfinite(samples).all()
  assert ((samples > 0.0) & (samples < 1.0)).all()

  samples_sobol = bc.simulate(n=500, qrng=True, seeds=[0])
  assert samples_sobol.shape == (500, 2)
  assert torch.isfinite(samples_sobol).all()
  assert ((samples_sobol > 0.0) & (samples_sobol < 1.0)).all()


def test_from_bicop_rejects_non_kernel_family() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  with pytest.raises(ValueError, match="kernel-family"):
    TorchBicop.from_bicop(cop)


# --------------------------------------------------------------------------- #
# TorchBicop linear-grid path                                                  #
# --------------------------------------------------------------------------- #


def test_linear_grid_roundtrip_and_range() -> None:
  """``grid_type='linear'`` produces a usable fit: hinv round-trips and
  all outputs are in-range. The linear-vs-normal disagreement on pdf is
  small (~5e-2 max) since both fit the same KDE in z-space and only the
  storage grid differs.
  """
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.simulate(2000, seeds=[10, 11, 12])
  # cache=False so the hinv round-trip below holds at 1e-9 (the cached
  # path round-trips only to ~1e-3 — that's covered by
  # test_cached_hinv_speedup).
  bc_lin = TorchBicop.from_data(
    u_fit, FitControlsTorchBicop(grid_type="linear"), cache_integrals=False
  )
  bc_nrm = TorchBicop.from_data(
    u_fit, FitControlsTorchBicop(grid_type="normal"), cache_integrals=False
  )

  assert bc_lin.interp_grid._is_linear is True
  assert bc_nrm.interp_grid._is_linear is False
  # Linear grid points are equispaced after the [0,1] endpoint clamp.
  gp = bc_lin.interp_grid.grid_points.numpy()
  d = np.diff(gp[1:-1])
  assert np.allclose(d, d[0], atol=1e-12), "interior spacing must be uniform"

  u_t = torch.from_numpy(_eval_grid(300, seed=99))
  pdf = bc_lin.pdf(u_t)
  cdf = bc_lin.cdf(u_t)
  h1 = bc_lin.hfunc1(u_t)
  h2 = bc_lin.hfunc2(u_t)
  assert (pdf > 0).all() and torch.isfinite(pdf).all()
  assert ((cdf > 0) & (cdf < 1)).all()
  assert ((h1 >= 0) & (h1 <= 1)).all()
  assert ((h2 >= 0) & (h2 <= 1)).all()

  # hinv round-trip
  u2 = bc_lin.hinv1(u_t).unsqueeze(-1)
  back = bc_lin.hfunc1(torch.cat([u_t[:, 0:1], u2], dim=-1)).numpy()
  np.testing.assert_allclose(back, u_t[:, 1].numpy(), atol=1e-9, rtol=1e-9)


def test_linear_grid_cell_index_matches_searchsorted() -> None:
  """The O(1) floor-based cell index for the linear grid must agree with
  the generic searchsorted-based cell index — same indices at every input.
  """
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.4]]))
  u_fit = cop.simulate(500, seeds=[1, 2, 3])
  bc_lin = TorchBicop.from_data(
    u_fit, FitControlsTorchBicop(grid_type="linear")
  )
  grid = bc_lin.interp_grid

  u = torch.linspace(0.0, 1.0, 1001, dtype=torch.float64)
  # O(1) floor path
  fast_idx = grid._cell_index(u)
  # Force the searchsorted path for the same grid
  grid._is_linear = False
  slow_idx = grid._cell_index(u)
  grid._is_linear = True
  assert torch.equal(fast_idx, slow_idx)


def test_linear_grid_cached_integrals_consistent() -> None:
  """``cache_integrals=True`` must build all five caches and produce in-
  range outputs on the linear grid as well as on the normal grid."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(1500, seeds=[7, 8, 9])
  bc = TorchBicop.from_data(
    u_fit, FitControlsTorchBicop(grid_type="linear"), cache_integrals=True
  )
  assert bc._cdf_cache is not None
  assert bc._hfunc1_cache is not None
  assert bc._hfunc2_cache is not None
  assert bc._hinv1_cache is not None
  assert bc._hinv2_cache is not None

  u_t = torch.from_numpy(_eval_grid(300, seed=21))
  for fn in (bc.cdf, bc.hfunc1, bc.hfunc2, bc.hinv1, bc.hinv2):
    out = fn(u_t)
    assert torch.isfinite(out).all()
    assert ((out >= 0.0) & (out <= 1.0)).all()


def test_linear_rejects_invalid_grid_type() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.4]]))
  u_fit = cop.simulate(200, seeds=[1, 2, 3])
  with pytest.raises(ValueError, match="grid_type"):
    TorchBicop.from_data(u_fit, FitControlsTorchBicop(grid_type="quadratic"))


# --------------------------------------------------------------------------- #
# User-facing input validation                                                #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
  "op", ["pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"]
)
def test_eval_rejects_wrong_input_shape(op: str) -> None:
  """Every eval entry point requires (n, 2) and raises ValueError otherwise."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.simulate(500, seeds=[1, 2, 3])
  bc = TorchBicop.from_data(u_fit)
  fn = getattr(bc, op)
  with pytest.raises(ValueError, match=r"shape \(n, 2\)"):
    fn(torch.zeros(100, dtype=torch.float64))  # 1-D
  with pytest.raises(ValueError, match=r"shape \(n, 2\)"):
    fn(torch.zeros(100, 3, dtype=torch.float64))  # wrong second dim


def test_from_data_rejects_wrong_input_shape() -> None:
  """`from_data` rejects inputs that aren't ``(n, 2)``."""
  with pytest.raises(ValueError, match="must have shape"):
    TorchBicop.from_data(torch.zeros(100, dtype=torch.float64))
  with pytest.raises(ValueError, match="must have shape"):
    TorchBicop.from_data(torch.zeros(100, 3, dtype=torch.float64))


def test_from_data_rejects_bad_args() -> None:
  """`from_data` rejects grid_size < 2 and mult <= 0."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.simulate(200, seeds=[1, 2, 3])
  with pytest.raises(ValueError, match="grid_size"):
    TorchBicop.from_data(u_fit, FitControlsTorchBicop(grid_size=1))
  with pytest.raises(ValueError, match="mult"):
    TorchBicop.from_data(u_fit, FitControlsTorchBicop(mult=0.0))


def test_from_data_rejects_unknown_method() -> None:
  """`FitControlsTorchBicop` rejects unknown ``method`` values up-front."""
  with pytest.raises(ValueError, match="unknown method"):
    FitControlsTorchBicop(method="bogus")


def test_simulate_rejects_nonpositive_n() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.simulate(200, seeds=[1, 2, 3])
  bc = TorchBicop.from_data(u_fit)
  with pytest.raises(ValueError, match="must be > 0"):
    bc.simulate(n=0)
  with pytest.raises(ValueError, match="must be > 0"):
    bc.simulate(n=-5)


# ---------------------------------------------------------------------------
# InterpolationGrid2D.normalize_margins(tol=...) early-stop
# ---------------------------------------------------------------------------


def test_normalize_margins_tol_converges() -> None:
  """``tol`` early-stops once margins are uniform to the requested
  precision. The default ``tol=None`` path keeps fixed-budget semantics
  for C++ TLL parity."""
  from pyvinecopulib.torch._interp import InterpolationGrid2D

  # Deliberately skewed initial density on a 16x16 linear grid.
  m = 16
  grid = torch.linspace(0.0, 1.0, m, dtype=torch.float64)
  rng = np.random.default_rng(0)
  values = torch.from_numpy(rng.uniform(0.5, 1.5, size=(m, m)))

  # times=50 with tol=1e-12 should converge to marginals within tol.
  ig = InterpolationGrid2D(grid, values, norm_times=50, norm_tol=1e-12)
  dgrid = ig.grid_points[1:] - ig.grid_points[:-1]
  row_int = 0.5 * ((ig.values[:, :-1] + ig.values[:, 1:]) * dgrid).sum(-1)
  col_int = 0.5 * (
    (ig.values[:-1, :] + ig.values[1:, :]) * dgrid.unsqueeze(-1)
  ).sum(0)
  assert (row_int - 1.0).abs().max().item() < 1e-10
  assert (col_int - 1.0).abs().max().item() < 1e-10


def test_normalize_margins_default_tol_is_fixed_budget() -> None:
  """``tol=None`` (default) preserves the fixed-budget loop — same number
  of divides every time, matching the C++ TLL pipeline's
  ``normalize_margins(3)`` byte-for-byte. We verify this indirectly by
  checking that the standard ``from_data`` path still matches C++ to
  machine precision (covered by ``test_from_data_matches_cpp``); here we
  just spot-check that the constructor accepts the default."""
  from pyvinecopulib.torch._interp import InterpolationGrid2D

  grid = torch.linspace(0.0, 1.0, 8, dtype=torch.float64)
  values = torch.ones(8, 8, dtype=torch.float64)
  # Should not raise; tol is optional.
  InterpolationGrid2D(grid, values, norm_times=3)


def test_flip_swaps_arguments() -> None:
  # Rotated Clayton data give an asymmetric TLL fit, so the flip is a genuine
  # change (not a no-op). flip() must give c'(u1, u2) = c(u2, u1) with the
  # h-functions swapped, exactly (it transposes the interpolation grid), and
  # leave the original unchanged.
  cop = pv.Bicop.from_family(
    family=pv.families.clayton, rotation=90, parameters=np.array([[3.0]])
  )
  u_fit = cop.simulate(4000, seeds=[1, 2, 3])
  tb = TorchBicop.from_data(
    u_fit, FitControlsTorchBicop(), cache_integrals=False
  )
  before = tb.pdf(torch.tensor([[0.3, 0.7]], dtype=torch.float64)).clone()
  flipped = tb.flip()
  u = torch.as_tensor(_eval_grid(300, seed=7), dtype=torch.float64)
  swapped = u[:, [1, 0]]
  np.testing.assert_allclose(
    flipped.pdf(u).numpy(), tb.pdf(swapped).numpy(), atol=1e-12, rtol=1e-12
  )
  np.testing.assert_allclose(
    flipped.hfunc1(u).numpy(),
    tb.hfunc2(swapped).numpy(),
    atol=1e-12,
    rtol=1e-12,
  )
  np.testing.assert_allclose(
    flipped.hfunc2(u).numpy(),
    tb.hfunc1(swapped).numpy(),
    atol=1e-12,
    rtol=1e-12,
  )
  # Genuinely different from the unflipped copula on an asymmetric fit, and
  # the original is left unchanged.
  assert not np.allclose(flipped.pdf(u).numpy(), tb.pdf(u).numpy())
  np.testing.assert_array_equal(
    tb.pdf(torch.tensor([[0.3, 0.7]], dtype=torch.float64)).numpy(),
    before.numpy(),
  )
