import pickle

import numpy as np
import pytest

import pyvinecopulib as pv

from .helpers import (
  compare_bicop,
  compare_kde1d,
  compare_properties,
  compare_rvinestructure,
  compare_vinecop,
  random_data,
)


def _custom_criterion(data: np.ndarray, weights: np.ndarray) -> float:
  # Module-level (hence picklable) custom tree criterion.
  return float(pv.utils.wdm(data[:, 0], data[:, 1], "tau"))


def test_fitcontrolsbicop() -> None:
  original_controls = pv.FitControlsBicop()
  original_controls.family_set = pv.families.itau
  original_controls.parametric_method = "itau"

  # Serialize the object
  serialized = pickle.dumps(original_controls)

  # Deserialize the object
  deserialized_controls = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  attrs = [
    "family_set",
    "parametric_method",
    "nonparametric_method",
    "nonparametric_mult",
    "nonparametric_grid_size",
    "selection_criterion",
    "weights",
    "psi0",
    "preselect_families",
    "num_threads",
  ]
  compare_properties(original_controls, deserialized_controls, attrs)


def test_fitcontrolsvinecop() -> None:
  # Create an instance of FitControlsVinecop with some configuration
  original_controls = pv.FitControlsVinecop()
  original_controls.family_set = pv.families.itau
  original_controls.parametric_method = "itau"

  # Serialize the object
  serialized = pickle.dumps(original_controls)

  # Deserialize the object
  deserialized_controls = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  attrs = [
    "family_set",
    "parametric_method",
    "nonparametric_method",
    "weights",
    "nonparametric_mult",
    "nonparametric_grid_size",
    "trunc_lvl",
    "tree_criterion",
    "tree_criterion_function",
    "threshold",
    "selection_criterion",
    "psi0",
    "preselect_families",
    "select_trunc_lvl",
    "select_threshold",
    "select_families",
    "show_trace",
    "num_threads",
    "tree_algorithm",
    "allow_rotations",
    "seeds",
  ]
  compare_properties(original_controls, deserialized_controls, attrs)


def test_fitcontrolsvinecop_custom_criterion() -> None:
  # A module-level custom criterion round-trips by reference; the getter
  # returns the same object after unpickling.
  original_controls = pv.FitControlsVinecop(tree_criterion="custom")
  original_controls.tree_criterion_function = _custom_criterion

  deserialized_controls = pickle.loads(pickle.dumps(original_controls))
  assert deserialized_controls.tree_criterion == "custom"
  assert deserialized_controls.tree_criterion_function is _custom_criterion

  # A non-picklable callable (lambda) raises the standard pickling error.
  original_controls.tree_criterion_function = lambda data, weights: 0.0
  with pytest.raises((pickle.PicklingError, AttributeError)):
    pickle.dumps(original_controls)


def test_bicop() -> None:
  original_bicop = pv.Bicop(pv.families.gaussian)
  original_bicop.parameters = np.array([[0.5]])

  # Serialize the object
  serialized = pickle.dumps(original_bicop)

  # Deserialize the object
  deserialized_bicop = pickle.loads(serialized)

  # Assert that the deserialized object's properties match the original
  compare_bicop(original_bicop, deserialized_bicop)


def test_rvinestructure() -> None:
  # Create an instance of RVineStructure with some configuration
  original_structure = pv.RVineStructure.sample(5)

  # Serialize the object
  serialized = pickle.dumps(original_structure)

  # Deserialize the object
  deserialized_structure = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  compare_rvinestructure(original_structure, deserialized_structure)


def test_kde1d() -> None:
  # Test with unfitted Kde1d object first
  original_kde = pv.core.Kde1d(
    xmin=-5.0,
    xmax=5.0,
    type="continuous",
    multiplier=1.5,
    degree=1,
    bandwidth=0.1,
    grid_size=100,
  )

  # Serialize the unfitted object
  serialized = pickle.dumps(original_kde)

  # Deserialize the object
  deserialized_kde = pickle.loads(serialized)

  # Assert that the deserialized object's properties match the original (unfitted)
  compare_kde1d(original_kde, deserialized_kde)

  # Now test with fitted model
  np.random.seed(1234)
  x = np.random.normal(0, 1, 100)
  original_kde.fit(x)

  # Serialize the fitted object
  serialized_fitted = pickle.dumps(original_kde)

  # Deserialize the fitted object
  deserialized_fitted = pickle.loads(serialized_fitted)

  # Assert that the deserialized fitted object's properties match the original
  compare_kde1d(original_kde, deserialized_fitted)


def test_vinecop() -> None:
  d = 5
  n = 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  controls = pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  assert controls.family_set == [pv.families.gaussian]
  original_cop = pv.Vinecop.from_data(u, controls=controls)

  # Serialize the object
  serialized = pickle.dumps(original_cop)

  # Deserialize the object
  deserialized_cop = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  compare_vinecop(original_cop, deserialized_cop)
