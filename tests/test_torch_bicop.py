"""Tests for the optional ``pyvinecopulib.torch`` submodule.

Skipped when PyTorch isn't installed. Compares the torch ``TorchBicop``
against the C++ ``pv.Bicop`` (TLL family) for ``pdf`` / ``cdf`` / ``hfunc``
on the same fitted interpolation grid; verifies ``hinv`` round-trips; and
spot-checks the rotation handling on an asymmetric Clayton sample.
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


def test_pdf_matches_pvbicop_rotation_0() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  bc = TorchBicop.from_bicop(cop_tll)
  u_eval = _eval_grid(500, seed=11)

  out_torch = bc.pdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_cdf_matches_pvbicop_rotation_0() -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[1, 2, 3])
  cop_tll = _fit_tll(u_fit)

  bc = TorchBicop.from_bicop(cop_tll)
  u_eval = _eval_grid(500, seed=12)

  out_torch = bc.cdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.cdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_hfunc_matches_pvbicop_rotation_0() -> None:
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
  np.testing.assert_allclose(back, u_eval[:, 1], atol=1e-5, rtol=1e-5)

  u1 = bc.hinv2(u_t).unsqueeze(-1)
  back = bc.hfunc2(torch.cat([u1, u_t[:, 1:2]], dim=-1)).numpy()
  np.testing.assert_allclose(back, u_eval[:, 0], atol=1e-5, rtol=1e-5)


def _rotate_data_np(u: np.ndarray, rotation: int) -> np.ndarray:
  """Numpy mirror of ``Bicop::rotate_data`` (class.ipp). Counter-clockwise."""
  out = u.copy()
  if rotation == 0:
    return out
  if rotation == 90:
    out = out[:, [1, 0]]
    out[:, 1] = 1.0 - out[:, 1]
    return out
  if rotation == 180:
    return 1.0 - out
  if rotation == 270:
    out = out[:, [1, 0]]
    out[:, 0] = 1.0 - out[:, 0]
    return out
  raise ValueError(rotation)


@pytest.mark.parametrize("rotation", [0, 90, 180, 270])
def test_rotations_apply_canonical_formulas(rotation: int) -> None:
  """TorchBicop must mirror ``Bicop``'s rotation handling end-to-end.

  TLL is rotationless on the C++ side, so we can't ask ``pv.Bicop`` to fit a
  rotated TLL directly. Instead we fit a rotation-0 TLL, construct a
  ``TorchBicop`` with the same grid but an arbitrary rotation, and compare
  against the rotation formulas from ``class.ipp`` applied to the underlying
  rotation-0 evaluations.
  """
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[7, 8, 9])
  cop_tll = _fit_tll(u_fit)
  assert cop_tll.rotation == 0

  bc_zero = TorchBicop.from_bicop(cop_tll)
  bc_rot = TorchBicop(
    grid_points=bc_zero.interp_grid.grid_points,
    values=bc_zero.interp_grid.values,
    rotation=rotation,
    norm_times=0,
  )

  u = _eval_grid(400, seed=30 + rotation)
  u_t = torch.from_numpy(u)
  u_rot = _rotate_data_np(u, rotation)

  pdf_exp = cop_tll.pdf(u_rot)  # |Jacobian| = 1 for any rotation
  np.testing.assert_allclose(
    bc_rot.pdf(u_t).numpy(), pdf_exp, atol=1e-10, rtol=1e-10
  )

  cdf0 = cop_tll.cdf(u_rot)
  if rotation == 0:
    cdf_exp = cdf0
  elif rotation == 90:
    cdf_exp = u[:, 1] - cdf0
  elif rotation == 180:
    cdf_exp = cdf0 - 1.0 + u[:, 0] + u[:, 1]
  else:
    cdf_exp = u[:, 0] - cdf0
  np.testing.assert_allclose(
    bc_rot.cdf(u_t).numpy(), cdf_exp, atol=1e-10, rtol=1e-10
  )

  if rotation == 0:
    h1_exp = cop_tll.hfunc1(u_rot)
    h2_exp = cop_tll.hfunc2(u_rot)
  elif rotation == 90:
    h1_exp = cop_tll.hfunc2(u_rot)
    h2_exp = 1.0 - cop_tll.hfunc1(u_rot)
  elif rotation == 180:
    h1_exp = 1.0 - cop_tll.hfunc1(u_rot)
    h2_exp = 1.0 - cop_tll.hfunc2(u_rot)
  else:
    h1_exp = 1.0 - cop_tll.hfunc2(u_rot)
    h2_exp = cop_tll.hfunc1(u_rot)
  np.testing.assert_allclose(
    bc_rot.hfunc1(u_t).numpy(), h1_exp, atol=1e-10, rtol=1e-10
  )
  np.testing.assert_allclose(
    bc_rot.hfunc2(u_t).numpy(), h2_exp, atol=1e-10, rtol=1e-10
  )


@pytest.mark.parametrize("rotation", [0, 90, 180, 270])
def test_hinv_roundtrip_with_rotation(rotation: int) -> None:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.7]]))
  u_fit = cop.simulate(2000, seeds=[4, 5, 6])
  cop_tll = _fit_tll(u_fit)

  bc_zero = TorchBicop.from_bicop(cop_tll)
  bc = TorchBicop(
    grid_points=bc_zero.interp_grid.grid_points,
    values=bc_zero.interp_grid.values,
    rotation=rotation,
    norm_times=0,
  )

  u = _eval_grid(400, seed=70 + rotation)
  u_t = torch.from_numpy(u)

  u2 = bc.hinv1(u_t).unsqueeze(-1)
  back1 = bc.hfunc1(torch.cat([u_t[:, 0:1], u2], dim=-1)).numpy()
  np.testing.assert_allclose(back1, u[:, 1], atol=1e-5, rtol=1e-5)

  u1 = bc.hinv2(u_t).unsqueeze(-1)
  back2 = bc.hfunc2(torch.cat([u1, u_t[:, 1:2]], dim=-1)).numpy()
  np.testing.assert_allclose(back2, u[:, 0], atol=1e-5, rtol=1e-5)


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
