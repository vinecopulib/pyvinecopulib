import os
import shutil

import numpy as np
import pytest
import pyvinecopulib as pv


def test_bicop():
  bicop = pv.Bicop()

  # Test default initialization
  assert bicop.family == pv.BicopFamily.indep
  assert bicop.rotation == 0
  assert bicop.parameters.shape == (0, 0)
  assert bicop.var_types == ["c", "c"]

  # Test initialization with arguments
  data = [[0.1, 0.2], [0.3, 0.4]]
  controls = pv.FitControlsBicop()
  bicop = pv.Bicop(data, controls)

  assert bicop.family == pv.BicopFamily.indep
  assert bicop.rotation == 0
  assert bicop.parameters.shape == (0, 0)
  assert bicop.var_types == ["c", "c"]

  # Test to_json method
  test_folder = "test_dump"
  os.makedirs(test_folder, exist_ok=True)
  filename = test_folder + "/test_bicop.json"
  bicop.to_json(filename)
  assert os.path.exists(filename)
  new_bicop = pv.Bicop(filename)
  assert bicop.family == new_bicop.family
  assert bicop.rotation == new_bicop.rotation
  assert bicop.parameters.shape == new_bicop.parameters.shape
  assert bicop.var_types == new_bicop.var_types

  # Test properties
  bicop = pv.Bicop(family=pv.BicopFamily.gumbel, rotation=90)
  bicop.rotation = 0
  assert bicop.rotation == 0
  with pytest.raises(RuntimeError):
    bicop.rotation = 45

  bicop.parameters = [3]
  assert bicop.parameters.shape == (1, 1)

  bicop.var_types = ["d", "d"]
  assert bicop.var_types == ["d", "d"]

  # Test read-only properties
  assert isinstance(bicop.tau, float)
  assert bicop.npars == 1
  with pytest.raises(AttributeError):
    bicop.npars = 2

  # Test loglik method
  bicop.var_types = ["c", "c"]
  u = [[0.1, 0.2], [0.3, 0.4]]
  loglik = bicop.loglik(u)
  assert isinstance(loglik, float)

  # Test aic method
  aic = bicop.aic(u)
  assert isinstance(aic, float)

  # Test bic method
  bic = bicop.bic(u)
  assert isinstance(bic, float)

  # Test mbic method
  psi0 = 0.9
  mbic = bicop.mbic(u, psi0)
  assert isinstance(mbic, float)

  # Test __repr__ method
  assert isinstance(repr(bicop), str)

  # Test str method
  assert isinstance(bicop.str(), str)

  # Test parameters_to_tau method
  parameters = [[0.5, 0.6], [0.7, 0.8]]
  tau = bicop.parameters_to_tau(parameters)
  assert isinstance(tau, float)

  # Test tau_to_parameters method
  tau = 0.5
  parameters = bicop.tau_to_parameters(tau)
  assert isinstance(parameters, np.ndarray)

  # Test parameters_lower_bounds method
  lower_bounds = bicop.parameters_lower_bounds()
  assert isinstance(lower_bounds, np.ndarray)

  # Test parameters_upper_bounds method
  upper_bounds = bicop.parameters_upper_bounds()
  assert isinstance(upper_bounds, np.ndarray)

  for method in ["pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"]:
    values = getattr(bicop, method)(u)
    assert isinstance(values, np.ndarray)
    assert values.shape == (2,)

  # Test simulate method
  n = 100
  qrng = False
  seeds = []
  samples = bicop.simulate(n, qrng, seeds)
  assert samples.shape == (n, 2)

  # Test fit method
  data = [[0.1, 0.2], [0.3, 0.4]]
  controls = pv.FitControlsBicop()
  bicop.fit(data, controls)

  # Test select method
  data = [[0.1, 0.2], [0.3, 0.4]]
  controls = pv.FitControlsBicop()
  bicop.select(data, controls)

  # Clean up
  shutil.rmtree(test_folder)
