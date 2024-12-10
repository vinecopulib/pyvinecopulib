import os
import shutil

import numpy as np

import pyvinecopulib as pv


def random_data(d=5, n=1000):
  # Simulate some data
  np.random.seed(1234)  # seed for the random generator
  mean = np.random.normal(size=d)  # mean vector
  cov = np.random.normal(size=(d, d))  # covariance matrix
  cov = np.dot(cov.transpose(), cov)  # make it non-negative definite
  x = np.random.multivariate_normal(mean, cov, n)
  return x

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
  test_folder = "test_dump"
  os.makedirs(test_folder, exist_ok=True)
  filename = test_folder + "/test_vinecop.json"
  cop.to_file(filename)
  new_cop = pv.Vinecop.from_file(filename)
  assert new_cop.dim == cop.dim
  assert new_cop.trunc_lvl == cop.trunc_lvl
  assert new_cop.order == cop.order

  # Clean up
  shutil.rmtree(test_folder)
