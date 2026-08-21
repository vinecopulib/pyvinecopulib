"""Tests for the optional ``pyvinecopulib.torch`` submodule.

Skipped when PyTorch isn't installed. Compares the torch ``TorchBicop``
against the C++ ``pv.Bicop`` (TLL family) for ``pdf`` / ``cdf`` / ``hfunc``
on the same fitted interpolation grid and verifies ``hinv`` round-trips.
"""

from __future__ import annotations

from fractions import Fraction

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
  u_fit = cop.sample(2000, seeds=[1, 2, 3])
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
  u_fit = cop.sample(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  # Pin cache=False — C++ parity at 1e-10; see test_pdf_matches_pvbicop.
  bc = TorchBicop.from_bicop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(500, seed=12)

  out_torch = bc.cdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.cdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_hfunc_matches_pvbicop() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.sample(2000, seeds=[1, 2, 3])
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
  u_fit = cop.sample(2000, seeds=[4, 5, 6])
  cop_tll = _fit_tll(u_fit)

  # Pin cache=False: the closed-form inversion is the exact inverse of the
  # on-the-fly h-function, so the round-trip holds to machine precision. The
  # cached path inverts the same quadratic but reads `hfunc1` off the prefix
  # tables, so it round-trips to summation-order noise rather than exactly.
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
  u_fit = cop.sample(2000, seeds=[4, 5, 6])
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
  u_fit = cop.sample(1500, seeds=[7, 8, 9])
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
  margin normalization in :class:`InterpolationGrid2D`.
  """
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[rho]]))
  u_np = cop.sample(n, seeds=[1, 2, 3])
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
  u_np = cop.sample(1000, seeds=[11, 22, 33])
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
  u_fit = cop.sample(2000, seeds=[10, 11, 12])
  cop_tll = _fit_tll(u_fit)

  bc_cache = TorchBicop.from_bicop(cop_tll, cache_integrals=True)
  # All three prefix tables must be populated for a non-indep pair.
  assert bc_cache._sy is not None
  assert bc_cache._sx is not None
  assert bc_cache._prefix is not None

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


def test_cached_and_uncached_hinv_agree_exactly() -> None:
  """``hinv1`` / ``hinv2`` do not read the prefix tables, in either mode.

  Locating the bracketing cell of the inverse needs the conditional cumulative
  along the whole free axis, which is O(m) to assemble — so there is no O(1)
  exact lookup to cache, and both modes run the same closed-form inversion of
  the same quadratic. The two are therefore bit-identical, not merely close.
  """
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.sample(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  bc_bisect = TorchBicop.from_bicop(cop_tll, cache_integrals=False)
  bc_cached = TorchBicop.from_bicop(cop_tll, cache_integrals=True)

  u_t = torch.from_numpy(_eval_grid(500, seed=77))
  for which in ("hinv1", "hinv2"):
    out_bisect = getattr(bc_bisect, which)(u_t)
    out_cached = getattr(bc_cached, which)(u_t)
    torch.testing.assert_close(out_cached, out_bisect, atol=0.0, rtol=0.0)


def test_independent_bicop() -> None:
  bc = TorchBicop()  # default: independence
  u = torch.tensor([[0.2, 0.7], [0.5, 0.5], [0.9, 0.1]], dtype=torch.float64)
  assert torch.allclose(bc.pdf(u), torch.ones(3, dtype=torch.float64))
  assert torch.allclose(bc.cdf(u), u[:, 0] * u[:, 1])
  assert torch.allclose(bc.hfunc1(u), u[:, 1])
  assert torch.allclose(bc.hfunc2(u), u[:, 0])


def test_simulate_smoke() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.sample(2000, seeds=[20, 21, 22])
  cop_tll = _fit_tll(u_fit)
  bc = TorchBicop.from_bicop(cop_tll)

  samples = bc.sample(n=1000, qrng=False, seeds=[0])
  assert samples.shape == (1000, 2)
  assert torch.isfinite(samples).all()
  assert ((samples > 0.0) & (samples < 1.0)).all()

  samples_sobol = bc.sample(n=500, qrng=True, seeds=[0])
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
  u_fit = cop.sample(2000, seeds=[10, 11, 12])
  # cache=False so the hinv round-trip below is the exact inverse of the
  # h-function it is composed with, holding at 1e-9.
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
  u_fit = cop.sample(500, seeds=[1, 2, 3])
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
  """``cache_integrals=True`` must build all three prefix tables and produce
  in-range outputs on the linear grid as well as on the normal grid."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.sample(1500, seeds=[7, 8, 9])
  bc = TorchBicop.from_data(
    u_fit, FitControlsTorchBicop(grid_type="linear"), cache_integrals=True
  )
  assert bc._sy is not None
  assert bc._sx is not None
  assert bc._prefix is not None

  u_t = torch.from_numpy(_eval_grid(300, seed=21))
  for fn in (bc.cdf, bc.hfunc1, bc.hfunc2, bc.hinv1, bc.hinv2):
    out = fn(u_t)
    assert torch.isfinite(out).all()
    assert ((out >= 0.0) & (out <= 1.0)).all()


def test_linear_rejects_invalid_grid_type() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.4]]))
  u_fit = cop.sample(200, seeds=[1, 2, 3])
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
  u_fit = cop.sample(500, seeds=[1, 2, 3])
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
  u_fit = cop.sample(200, seeds=[1, 2, 3])
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
  u_fit = cop.sample(200, seeds=[1, 2, 3])
  bc = TorchBicop.from_data(u_fit)
  with pytest.raises(ValueError, match="must be > 0"):
    bc.sample(n=0)
  with pytest.raises(ValueError, match="must be > 0"):
    bc.sample(n=-5)


# ---------------------------------------------------------------------------
# InterpolationGrid2D.normalize_margins
# ---------------------------------------------------------------------------


def test_normalize_margins_balances_both_margins() -> None:
  """Both margins end up uniform, and equally so.

  The old scheme was three fixed row-then-column sweeps, which left the second
  margin exact and dumped the whole residual on the first -- up to 3.3e-2 at
  strong dependence, so a fitted grid was not a copula density in one direction.
  Averaging the two sweep orders splits the residual, which is what makes the
  balance assertion meaningful rather than tautological.
  """
  from pyvinecopulib.torch._interp import InterpolationGrid2D, _trap_weights

  m = 16
  grid = torch.linspace(0.0, 1.0, m, dtype=torch.float64)
  rng = np.random.default_rng(0)
  values = torch.from_numpy(rng.uniform(0.5, 1.5, size=(m, m)))

  ig = InterpolationGrid2D(grid, values)
  w = _trap_weights(ig.grid_points)
  r = (ig.values @ w - 1.0).abs().max().item()
  c = (ig.values.t().contiguous() @ w - 1.0).abs().max().item()
  assert max(r, c) < 1e-10
  # Neither margin is allowed to carry the whole residual.
  assert max(r, c) < 10 * max(min(r, c), 1e-18)


def test_normalize_margins_commutes_with_transposition() -> None:
  """`normalize(V).T == normalize(V.T)`, exactly.

  This is what makes `flip` correct: it transposes an already-normalized grid
  without renormalizing, so if the normalization were not equivariant then
  `fit(a, b).flip()` and `fit(b, a)` would be different models. Under the old
  three-sweep scheme they differed by 2.7e-4.
  """
  from pyvinecopulib.torch._interp import InterpolationGrid2D

  m = 16
  grid = torch.linspace(0.0, 1.0, m, dtype=torch.float64)
  rng = np.random.default_rng(1)
  values = torch.from_numpy(rng.uniform(0.2, 2.0, size=(m, m)))

  direct = InterpolationGrid2D(grid, values).values
  swapped = InterpolationGrid2D(grid, values.t().contiguous()).values
  torch.testing.assert_close(direct, swapped.t(), rtol=1e-13, atol=1e-15)


def test_normalize_margins_leaves_a_normalized_grid_alone() -> None:
  """An already-uniform grid costs one margin check and no scaling."""
  from pyvinecopulib.torch._interp import InterpolationGrid2D

  grid = torch.linspace(0.0, 1.0, 8, dtype=torch.float64)
  values = torch.ones(8, 8, dtype=torch.float64)
  ig = InterpolationGrid2D(grid, values)
  torch.testing.assert_close(ig.values, values, rtol=0.0, atol=0.0)


def test_flip_swaps_arguments() -> None:
  # Rotated Clayton data give an asymmetric TLL fit, so the flip is a genuine
  # change (not a no-op). flip() must give c'(u1, u2) = c(u2, u1) with the
  # h-functions swapped, exactly (it transposes the interpolation grid), and
  # leave the original unchanged.
  cop = pv.Bicop.from_family(
    family=pv.families.clayton, rotation=90, parameters=np.array([[3.0]])
  )
  u_fit = cop.sample(4000, seeds=[1, 2, 3])
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


# --- gradients ---------------------------------------------------------------- #


def test_hinv_is_differentiable() -> None:
  """Both hinv paths carry a gradient, and it is the right one.

  Neither is a root-finder: cached is one bilinear lookup, uncached is the
  closed-form quadratic root of the interpolated conditional cdf. Both were
  nonetheless wrapped in `torch.no_grad()`, which detached them for no reason.
  """
  gauss = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  cop = _fit_tll(gauss.sample(2000, seeds=[6, 7, 8]))
  q = torch.from_numpy(
    np.random.default_rng(25).uniform(0.15, 0.85, size=(6, 2))
  )

  for cache in (True, False):
    bc = TorchBicop.from_bicop(cop, cache_integrals=cache)
    for name in ("hinv1", "hinv2"):
      method = getattr(bc, name)
      x = q.clone().requires_grad_(True)
      out = method(x)
      assert out.requires_grad, f"{name}, cache={cache}"
      (grad,) = torch.autograd.grad(out.sum(), x)

      fd = torch.zeros_like(q)
      with torch.no_grad():
        for j in range(2):
          hi, lo = q.clone(), q.clone()
          hi[:, j] += 1e-6
          lo[:, j] -= 1e-6
          fd[:, j] = (method(hi) - method(lo)) / 2e-6
      rel = (grad - fd).abs().max().item() / fd.abs().max().item()
      assert rel < 1e-6, f"{name}, cache={cache}: {rel:.3e}"


def test_cached_integrals_carry_a_grid_gradient() -> None:
  """The cached ``cdf`` is differentiable in the grid, and agrees with the
  on-the-fly one.

  The prefix tables are cumulative trapezoids of ``values``, so the closed-form
  ``cdf`` read off them is an exact function of the grid -- not the bilinear
  surrogate the previous cache interposed, which carried no such gradient at
  all. Both directions matter: turning ``requires_grad`` on **after**
  construction has to reach the tables (they are buffers, so a stale detached
  copy would silently zero the gradient), and the two modes must agree.
  """
  gauss = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  cop = _fit_tll(gauss.sample(2000, seeds=[7, 8, 9]))
  q = torch.from_numpy(
    np.random.default_rng(26).uniform(0.15, 0.85, size=(8, 2))
  )

  grads = {}
  for cache in (True, False):
    bc = TorchBicop.from_bicop(cop, cache_integrals=cache)
    bc.interp_grid.values.requires_grad_(True)
    assert bc.pdf(q).requires_grad
    out = bc.cdf(q)
    assert out.requires_grad, f"cache={cache}: cdf lost its grid gradient"
    (grad,) = torch.autograd.grad(out.sum(), bc.interp_grid.values)
    assert torch.isfinite(grad).all()
    assert grad.abs().sum().item() > 0.0
    grads[cache] = grad

  # Both are exact integrals of the same bilinear interpolant -- `integrate_2d`
  # is piecewise-exact in each direction, the tables are cumulative trapezoids
  # of the same pieces -- so they differ only in the order the terms are summed.
  torch.testing.assert_close(grads[True], grads[False], atol=1e-10, rtol=1e-7)


def _exact_rect_prob(
  grid: list[Fraction], values: list[list[Fraction]], rect: tuple[Fraction, ...]
) -> Fraction:
  """The four-corner difference of the renormalized distribution function.

  Exact rational arithmetic throughout, so it is the reference the float routes
  are measured against. ``C(u1, u2) = M(u1, u2) * u2 / M(1, u2)``, the object
  ``integrate_2d`` computes.
  """
  m = len(grid)

  def weights(lo: Fraction, hi: Fraction) -> list[Fraction]:
    w = [Fraction(0)] * m

    def add(k: int, s0: Fraction, s1: Fraction) -> None:
      h = grid[k + 1] - grid[k]
      q = (s1 * s1 - s0 * s0) / 2
      w[k] += h * (s1 - s0 - q)
      w[k + 1] += h * q

    def cell(x: Fraction) -> int:
      k = 0
      while k < m - 2 and grid[k + 1] <= x:
        k += 1
      return k

    ka, kb = cell(lo), cell(hi)
    if ka == kb:
      h = grid[ka + 1] - grid[ka]
      add(ka, (lo - grid[ka]) / h, (hi - grid[ka]) / h)
      return w
    add(ka, (lo - grid[ka]) / (grid[ka + 1] - grid[ka]), Fraction(1))
    for k in range(ka + 1, kb):
      h = (grid[k + 1] - grid[k]) / 2
      w[k] += h
      w[k + 1] += h
    add(kb, Fraction(0), (hi - grid[kb]) / (grid[kb + 1] - grid[kb]))
    return w

  def mass(a1: Fraction, b1: Fraction, a2: Fraction, b2: Fraction) -> Fraction:
    wx, wy = weights(a1, b1), weights(a2, b2)
    return sum(
      (
        wx[i] * wy[j] * values[i][j]
        for i in range(m)
        if wx[i]
        for j in range(m)
        if wy[j]
      ),
      Fraction(0),
    )

  def cdf(u1: Fraction, u2: Fraction) -> Fraction:
    total = mass(Fraction(0), Fraction(1), Fraction(0), u2)
    return (
      Fraction(0)
      if total == 0
      else mass(Fraction(0), u1, Fraction(0), u2) * u2 / total
    )

  a1, b1, a2, b2 = rect
  return cdf(b1, b2) - cdf(a1, b2) - cdf(b1, a2) + cdf(a1, a2)


# Both bounds move with the width, which is the point: differencing four
# distribution values amplifies any absolute error by ~4 / (w1 w2), where
# `rect_mass` amplifies by 1 / w2 alone -- one power instead of two, because
# only its `lam` difference cancels and that multiplies a term of order w1.
# Measuring both in one run is what shows the gap is the construction rather
# than the test data: the two agree at 3/8 and are 3000x apart at 1/8192.
# Dyadic endpoints keep `float(x) == x`, so the comparison is against the
# algorithm and not against how the query points round.
@pytest.mark.parametrize(
  "width,tol_rect,tol_diff",
  [
    # One decade of tolerance per decade of width for `rect_mass`, two for the
    # difference -- which is the finding, stated as the shape of the table.
    (Fraction(3, 8), 3e-15, 3e-15),
    (Fraction(1, 16), 3e-14, 1e-13),
    (Fraction(1, 64), 1e-13, 3e-12),
    (Fraction(1, 1024), 2e-12, 1e-9),
    (Fraction(1, 8192), 1e-11, 3e-8),
  ],
)
def test_rect_mass_is_exact_at_every_atom_width(
  width: Fraction, tol_rect: float, tol_diff: float
) -> None:
  """``rect_mass`` against exact rational truth, with the differencing route
  measured beside it.
  """
  from pyvinecopulib.torch._interp import InterpolationGrid2D

  m = 30
  grid_f = [Fraction(i, m - 1) for i in range(m)]
  raw = np.random.default_rng(4).integers(1, 4000, size=(m, m))
  values_f = [[Fraction(int(raw[i, j])) for j in range(m)] for i in range(m)]

  gp = torch.tensor([float(x) for x in grid_f], dtype=torch.float64)
  vals = torch.tensor(raw, dtype=torch.float64)
  grid = InterpolationGrid2D(gp, vals, norm_maxiter=0, is_linear=True)

  rng = np.random.default_rng(5)
  lim = int((1 - width) * 512)
  rects = []
  for _ in range(40):
    a1 = Fraction(int(rng.integers(1, lim)), 512)
    a2 = Fraction(int(rng.integers(1, lim)), 512)
    rects.append((a1, a1 + width, a2, a2 + width))
  truth = np.array(
    [float(_exact_rect_prob(grid_f, values_f, r)) for r in rects]
  )

  a1t, b1t, a2t, b2t = (
    torch.tensor([float(r[k]) for r in rects], dtype=torch.float64)
    for k in range(4)
  )
  got = grid.rect_mass(a1t, b1t, a2t, b2t).numpy()

  def cdf(x1: torch.Tensor, x2: torch.Tensor) -> np.ndarray:
    return grid.integrate_2d(torch.stack([x1, x2], dim=-1)).numpy()

  diff4 = cdf(b1t, b2t) - cdf(a1t, b2t) - cdf(b1t, a2t) + cdf(a1t, a2t)
  e_rect = np.max(np.abs(got - truth) / truth)
  e_diff = np.max(np.abs(diff4 - truth) / truth)
  assert e_rect < tol_rect, f"rect_mass: {e_rect:.3e}"
  # The differencing route must stay inside its own, looser bound; if the two
  # columns ever coincided, the parametrization would be measuring nothing.
  assert e_diff < tol_diff, f"difference: {e_diff:.3e}"


def test_rect_mass_reproduces_the_cdf_on_a_corner_rectangle() -> None:
  """``rect_mass`` is the ``cdf``, not a different integral of the same grid.

  On a rectangle anchored at the origin the four-corner difference collapses to
  one ``cdf`` value, so the two must agree with no width-dependent slack.
  Pinning it keeps them from drifting into two definitions.
  """
  gauss = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  cop = _fit_tll(gauss.sample(2000, seeds=[31, 32, 33]))
  bc = TorchBicop.from_bicop(cop)
  u = torch.from_numpy(_eval_grid(200, seed=34))
  zero = torch.zeros_like(u[:, 0])
  torch.testing.assert_close(
    bc.rect_mass(zero, u[:, 0], zero, u[:, 1]),
    bc.cdf(u),
    rtol=1e-12,
    atol=1e-12,
  )


def test_values_must_be_nonnegative() -> None:
  """A density grid with a negative node is refused, not silently integrated.

  The exact tables and ``rect_mass`` both rely on nonnegativity for their
  no-cancellation bound, so it is a precondition rather than a nicety.
  """
  from pyvinecopulib.torch._interp import InterpolationGrid2D

  gp = torch.linspace(0.0, 1.0, 8, dtype=torch.float64)
  vals = torch.ones(8, 8, dtype=torch.float64)
  vals[3, 4] = -1e-12
  with pytest.raises(ValueError, match="nonnegative"):
    InterpolationGrid2D(gp, vals, norm_maxiter=0, is_linear=True)
