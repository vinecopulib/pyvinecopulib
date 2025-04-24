import numpy as np

import pyvinecopulib as pv

from .helpers import compare_vinecop, random_data


def test_vinecop(test_dump_folder: str) -> None:
  d = 5
  n = 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  controls = pv.FitControlsVinecop(family_set=[pv.gaussian])
  assert controls.family_set == [pv.gaussian]
  cop = pv.Vinecop.from_data(u, controls=controls)

  # Test get_pair_copula method
  for t in range(1, d):
    for e in range(d - t - 1):
      pair_copula = cop.get_pair_copula(t, e)
      assert isinstance(pair_copula, pv.Bicop)

      # Test get_family method
      family = cop.get_family(0, 0)
      assert family == pv.gaussian

      # Test get_rotation method
      rotation = cop.get_rotation(0, 0)
      assert rotation == 0

      # Test get_parameters method
      parameters = cop.get_parameters(0, 0)
      assert isinstance(parameters, np.ndarray)
      assert parameters.shape == (1, 1)
      assert -1 < parameters[0, 0] < 1

      # Test get_tau method
      tau = cop.get_tau(0, 0)
      assert isinstance(tau, float)

  for method in ["pdf", "cdf"]:
    values = getattr(cop, method)(u)
    assert isinstance(values, np.ndarray)
    assert values.shape == (n,)
    assert values.dtype == np.float64

  for method in ["rosenblatt", "inverse_rosenblatt"]:
    values = getattr(cop, method)(u)
    assert isinstance(values, np.ndarray)
    assert values.shape == (n, d)
    assert values.dtype == np.float64

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
  assert set(cop.order) == set(range(1, d + 1))
  assert isinstance(cop.structure, pv.RVineStructure)
  matrix = cop.matrix
  assert isinstance(matrix, np.ndarray)
  assert matrix.shape == (d, d)
  assert matrix.dtype == np.uint64
  assert np.all(np.logical_and(matrix >= 0, matrix <= d))

  # Test to_json and from_json
  new_cop = pv.Vinecop.from_json(cop.to_json())
  compare_vinecop(cop, new_cop)
  filename = test_dump_folder + "/test_vinecop.json"
  cop.to_file(filename)
  new_cop = pv.Vinecop.from_file(filename)
  compare_vinecop(cop, new_cop)
