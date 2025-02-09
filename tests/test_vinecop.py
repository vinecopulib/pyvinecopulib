import os
import shutil

import numpy as np

import pyvinecopulib as pv

from .helpers import compare_rvinestructure, compare_vinecop, random_data


def test_vinecop():
  d = 5
  n = 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  controls = pv.FitControlsVinecop(family_set=[pv.BicopFamily.gaussian])
  assert controls.family_set == [pv.BicopFamily.gaussian]
  cop = pv.Vinecop.from_data(u, controls=controls)

  # Test get_pair_copula method
  pair_copula = cop.get_pair_copula(0, 0)
  assert isinstance(pair_copula, pv.Bicop)

  # Test get_family method
  family = cop.get_family(0, 0)
  assert family == pv.BicopFamily.gaussian

  # Test get_rotation method
  rotation = cop.get_rotation(0, 0)
  assert rotation == 0

  # Test get_parameters method
  parameters = cop.get_parameters(0, 0)
  assert isinstance(parameters, np.ndarray)

  # Test get_tau method
  tau = cop.get_tau(0, 0)
  assert isinstance(tau, float)

  for method in ["pdf", "cdf"]:
    values = getattr(cop, method)(u)
    assert isinstance(values, np.ndarray)
    assert values.shape == (n,)

  # Test passing a single row of data (#169 & #170 fix)
  u1 = u[0, :].reshape(1, d)
  for method in ["pdf", "cdf"]:
    values = getattr(cop, method)(u1)
    assert isinstance(values, np.ndarray)
    assert values.shape == (1,)

  # Test simulate method
  simulated_data = cop.simulate(n)
  assert simulated_data.shape == (n, d)

  # Test loglik method
  loglik_value = cop.loglik(u)
  assert isinstance(loglik_value, float)

  # Test AIC method
  aic_value = cop.aic(u)
  assert isinstance(aic_value, float)

  # Test BIC method
  bic_value = cop.bic(u)
  assert isinstance(bic_value, float)

  # Test MBICV method
  mbicv_value = cop.mbicv(u)
  assert isinstance(mbicv_value, float)

  # Test truncate method
  cop.truncate(2)
  assert cop.trunc_lvl == 2

  # Test order and structure
  assert isinstance(cop.order, list)
  assert isinstance(cop.structure, pv.RVineStructure)

  # Test to_json and from_json
  new_cop = pv.Vinecop.from_json(cop.to_json())
  compare_vinecop(cop, new_cop)
  test_folder = "test_dump"
  os.makedirs(test_folder, exist_ok=True)
  filename = test_folder + "/test_vinecop.json"
  cop.to_file(filename)
  new_cop = pv.Vinecop.from_file(filename)
  compare_vinecop(cop, new_cop)

  # Clean up
  shutil.rmtree(test_folder)


def test_rvinestructure():
  # Test RVineStructure class
  rvine = pv.RVineStructure(5)
  assert isinstance(rvine, pv.RVineStructure)
  assert rvine.dim == 5
  assert rvine.matrix.shape == (5, 5)
  assert rvine.trunc_lvl == 4
  assert rvine.order == list(range(1, 6))

  # Should be the same as the previous one
  dvine = pv.DVineStructure(rvine.order)
  compare_rvinestructure(dvine, rvine, True)

  # Test to_json and from_json
  new_rvine = pv.RVineStructure.from_json(rvine.to_json())
  compare_rvinestructure(rvine, new_rvine)
  test_folder = "test_dump"
  os.makedirs(test_folder, exist_ok=True)
  filename = test_folder + "/test_rvine.json"
  rvine.to_file(filename)
  new_rvine = pv.RVineStructure.from_file(filename)
  compare_rvinestructure(rvine, new_rvine)

  # Clean up
  shutil.rmtree(test_folder)
