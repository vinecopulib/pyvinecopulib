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


def test_inverse_rosenblatt_matches_pvvinecop() -> None:
  u_fit = _simulate(d=5, n=2000, seed=3)
  cop_tll = _fit_tll_vine(u_fit)

  bc = TorchVinecop.from_vinecop(cop_tll)
  # Inverse-rosenblatt takes independent uniforms; sample fresh.
  w = _eval_grid(400, d=5, seed=13)

  out_torch = bc.inverse_rosenblatt(torch.from_numpy(w)).numpy()
  out_cpp = cop_tll.inverse_rosenblatt(w)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-5, rtol=1e-5)


def test_inverse_rosenblatt_roundtrip() -> None:
  u_fit = _simulate(d=5, n=2000, seed=4)
  cop_tll = _fit_tll_vine(u_fit)

  bc = TorchVinecop.from_vinecop(cop_tll)
  # Stay away from u ≈ 0 / 1: the conditional CDFs of the cascade become
  # extremely flat near the unit-cube corners and even the C++ round-trip
  # saturates there (we saw a 0.024 worst-case at u=0.02 on the same data).
  rng = np.random.default_rng(14)
  u_eval = rng.uniform(0.05, 0.95, size=(400, 5))
  u_t = torch.from_numpy(u_eval)
  back = bc.inverse_rosenblatt(bc.rosenblatt(u_t)).numpy()
  np.testing.assert_allclose(back, u_eval, atol=1e-5, rtol=1e-5)


def test_independent_vine_short_circuits() -> None:
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
    bc.pdf(u_t).numpy(), np.ones(200), atol=1e-12, rtol=0
  )
  np.testing.assert_allclose(
    bc.rosenblatt(u_t).numpy(), u, atol=1e-10, rtol=1e-10
  )
  np.testing.assert_allclose(
    bc.inverse_rosenblatt(u_t).numpy(), u, atol=1e-10, rtol=1e-10
  )


def test_truncated_vine_pdf() -> None:
  u_fit = _simulate(d=6, n=2000, seed=5)
  cop_tll = _fit_tll_vine(u_fit, trunc_lvl=2)
  assert cop_tll.structure.trunc_lvl == 2

  bc = TorchVinecop.from_vinecop(cop_tll)
  assert bc.trunc_lvl == 2

  u_eval = _eval_grid(300, d=6, seed=15)
  np.testing.assert_allclose(
    bc.pdf(torch.from_numpy(u_eval)).numpy(),
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
