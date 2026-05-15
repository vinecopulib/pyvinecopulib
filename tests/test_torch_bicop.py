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

from pyvinecopulib.torch import TorchBicop  # noqa: E402

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

  bc = TorchBicop.from_bicop(cop_tll)
  u_eval = _eval_grid(500, seed=11)

  out_torch = bc.pdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_cdf_matches_pvbicop() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  bc = TorchBicop.from_bicop(cop_tll)
  u_eval = _eval_grid(500, seed=12)

  out_torch = bc.cdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.cdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_hfunc_matches_pvbicop() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  bc = TorchBicop.from_bicop(cop_tll)
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

  bc = TorchBicop.from_bicop(cop_tll)
  u_eval = _eval_grid(400, seed=21)
  u_t = torch.from_numpy(u_eval)

  u2 = bc.hinv1(u_t).unsqueeze(-1)
  back = bc.hfunc1(torch.cat([u_t[:, 0:1], u2], dim=-1)).numpy()
  np.testing.assert_allclose(back, u_eval[:, 1], atol=1e-9, rtol=1e-9)

  u1 = bc.hinv2(u_t).unsqueeze(-1)
  back = bc.hfunc2(torch.cat([u1, u_t[:, 1:2]], dim=-1)).numpy()
  np.testing.assert_allclose(back, u_eval[:, 0], atol=1e-9, rtol=1e-9)


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
  matches-vs-cpp test covers the fit accuracy."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_np = cop.simulate(1000, seeds=[11, 22, 33])
  bc = TorchBicop.from_data(u_np)

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
  u_t = torch.from_numpy(_eval_grid(200, seed=42))
  for fn in (bc_cache.cdf, bc_cache.hfunc1, bc_cache.hfunc2):
    out = fn(u_t)
    assert torch.isfinite(out).all()
    assert (out >= 0.0).all() and (out <= 1.0).all()


def test_independent_bicop() -> None:
  bc = TorchBicop()  # default: independence
  u = torch.tensor([[0.2, 0.7], [0.5, 0.5], [0.9, 0.1]], dtype=torch.float64)
  assert torch.allclose(bc.pdf(u), torch.ones(3, dtype=torch.float64))
  assert torch.allclose(bc.cdf(u), u[:, 0] * u[:, 1])
  assert torch.allclose(bc.hfunc1(u), u[:, 1])
  assert torch.allclose(bc.hfunc2(u), u[:, 0])


def test_sample_smoke() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  u_fit = cop.simulate(2000, seeds=[20, 21, 22])
  cop_tll = _fit_tll(u_fit)
  bc = TorchBicop.from_bicop(cop_tll)

  samples = bc.sample(num_sample=1000, seed=0, is_sobol=False)
  assert samples.shape == (1000, 2)
  assert torch.isfinite(samples).all()
  assert ((samples > 0.0) & (samples < 1.0)).all()

  samples_sobol = bc.sample(num_sample=500, seed=0, is_sobol=True)
  assert samples_sobol.shape == (500, 2)
  assert torch.isfinite(samples_sobol).all()
  assert ((samples_sobol > 0.0) & (samples_sobol < 1.0)).all()


def test_from_bicop_rejects_non_kernel_family() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  with pytest.raises(ValueError, match="kernel-family"):
    TorchBicop.from_bicop(cop)
