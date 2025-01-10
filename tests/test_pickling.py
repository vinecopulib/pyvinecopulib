import pickle

import numpy as np

import pyvinecopulib as pv

from .helpers import (
  compare_bicop,
  compare_properties,
  compare_rvinestructure,
  compare_vinecop,
  random_data,
)


def test_fitcontrolsbicop():
  original_controls = pv.FitControlsBicop()

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
    "selection_criterion",
    "weights",
    "psi0",
    "preselect_families",
    "num_threads",
  ]
  compare_properties(original_controls, deserialized_controls, attrs)


def test_fitcontrolsvinecop():
  # Create an instance of FitControlsVinecop with some configuration
  original_controls = pv.FitControlsVinecop()

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
    "trunc_lvl",
    "tree_criterion",
    "threshold",
    "selection_criterion",
    "psi0",
    "preselect_families",
    "select_trunc_lvl",
    "select_threshold",
    "select_families",
    "show_trace",
    "num_threads",
    "mst_algorithm",
  ]
  compare_properties(original_controls, deserialized_controls, attrs)


def test_bicop():
  original_bicop = pv.Bicop(pv.BicopFamily.gaussian)
  original_bicop.parameters = np.array([[0.5]])

  # Serialize the object
  serialized = pickle.dumps(original_bicop)

  # Deserialize the object
  deserialized_bicop = pickle.loads(serialized)

  # Assert that the deserialized object's properties match the original
  compare_bicop(original_bicop, deserialized_bicop)


def test_rvinestructure():
  # Create an instance of RVineStructure with some configuration
  original_structure = pv.RVineStructure.simulate(5)

  # Serialize the object
  serialized = pickle.dumps(original_structure)

  # Deserialize the object
  deserialized_structure = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  compare_rvinestructure(original_structure, deserialized_structure)


def test_vinecop():
  d = 5
  n = 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  controls = pv.FitControlsVinecop(family_set=[pv.BicopFamily.gaussian])
  assert controls.family_set == [pv.BicopFamily.gaussian]
  original_cop = pv.Vinecop.from_data(u, controls=controls)

  # Serialize the object
  serialized = pickle.dumps(original_cop)

  # Deserialize the object
  deserialized_cop = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  compare_vinecop(original_cop, deserialized_cop)
