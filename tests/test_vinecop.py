import os

import numpy as np
import pytest

import pyvinecopulib as pv

from .helpers import compare_vinecop, random_data


def _wdm_tau(data: np.ndarray, weights: np.ndarray) -> float:
  # Reproduce the built-in "tau" edge weight (signed Kendall's tau). The
  # selector applies abs() and the sqrt(freq) factor on top, so a custom
  # function returning signed tau yields the same structure as "tau".
  if weights.size:
    return float(pv.utils.wdm(data[:, 0], data[:, 1], "tau", weights=weights))
  return float(pv.utils.wdm(data[:, 0], data[:, 1], "tau"))


def _raising_criterion(data: np.ndarray, weights: np.ndarray) -> float:
  raise ValueError("custom criterion failure")


def test_vinecop(unique_json_path: str) -> None:
  d = 5
  n = 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  controls = pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  assert controls.family_set == [pv.families.gaussian]
  cop = pv.Vinecop.from_data(u, controls=controls)

  # Test get_pair_copula method
  for t in range(1, d):
    for e in range(d - t - 1):
      pair_copula = cop.get_pair_copula(t, e)
      assert isinstance(pair_copula, pv.Bicop)

      # Test get_family method
      family = cop.get_family(0, 0)
      assert family == pv.families.gaussian

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
  filename = os.fspath(unique_json_path)
  cop.to_file(filename)
  new_cop = pv.Vinecop.from_file(filename)
  compare_vinecop(cop, new_cop)


def test_custom_criterion_matches_builtin_tau() -> None:
  # A custom function returning signed Kendall's tau must reproduce the
  # structure selected by the built-in "tau" criterion exactly.
  d, n = 5, 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  cop_tau = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(tree_criterion="tau")
  )

  controls = pv.FitControlsVinecop(tree_criterion="custom")
  controls.tree_criterion_function = _wdm_tau
  cop_custom = pv.Vinecop.from_data(u, controls=controls)

  assert np.array_equal(cop_tau.matrix, cop_custom.matrix)
  assert cop_tau.order == cop_custom.order


def test_tree_criterion_custom_accepted() -> None:
  # Regression that the submodule bump exposes the "custom" criterion.
  controls = pv.FitControlsVinecop()
  controls.tree_criterion = "custom"
  assert controls.tree_criterion == "custom"


def test_tree_criterion_function_roundtrip_property() -> None:
  controls = pv.FitControlsVinecop()
  assert controls.tree_criterion_function is None
  controls.tree_criterion_function = _wdm_tau
  # The getter hands back the original Python callable.
  assert controls.tree_criterion_function is _wdm_tau
  controls.tree_criterion_function = None
  assert controls.tree_criterion_function is None


def test_custom_criterion_unset_raises() -> None:
  d, n = 5, 300
  u = pv.to_pseudo_obs(random_data(d, n))
  controls = pv.FitControlsVinecop(tree_criterion="custom")
  with pytest.raises(RuntimeError):
    pv.Vinecop.from_data(u, controls=controls)


def test_custom_criterion_ignored_when_not_custom() -> None:
  # A function set while tree_criterion != "custom" is never called.
  d, n = 5, 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  controls = pv.FitControlsVinecop(tree_criterion="tau")
  controls.tree_criterion_function = _raising_criterion  # must be ignored
  cop_custom = pv.Vinecop.from_data(u, controls=controls)

  cop_tau = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(tree_criterion="tau")
  )
  assert np.array_equal(cop_tau.matrix, cop_custom.matrix)


def test_custom_criterion_exception_propagates() -> None:
  # A Python exception raised inside the criterion surfaces from from_data.
  d, n = 5, 1000
  u = pv.to_pseudo_obs(random_data(d, n))
  controls = pv.FitControlsVinecop(tree_criterion="custom", num_threads=1)
  controls.tree_criterion_function = _raising_criterion
  with pytest.raises(ValueError):
    pv.Vinecop.from_data(u, controls=controls)


def test_custom_criterion_multithread_smoke() -> None:
  # Calls acquire the GIL, so num_threads > 1 must still give the same result.
  d, n = 5, 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  c1 = pv.FitControlsVinecop(tree_criterion="custom", num_threads=1)
  c1.tree_criterion_function = _wdm_tau
  cop1 = pv.Vinecop.from_data(u, controls=c1)

  c2 = pv.FitControlsVinecop(tree_criterion="custom", num_threads=2)
  c2.tree_criterion_function = _wdm_tau
  cop2 = pv.Vinecop.from_data(u, controls=c2)

  assert np.array_equal(cop1.matrix, cop2.matrix)


def _check_triangular(arr: object, d: int, leaf_check) -> None:
  """Validate a ``[tree][edge]`` nested list: tree ``i`` has ``d - 1 - i``
  entries, each passing ``leaf_check``. isinstance narrowing at every level
  keeps the (untyped) dict/list contents type-checkable."""
  assert isinstance(arr, list)
  for i, tree in enumerate(arr):
    assert isinstance(tree, list)
    assert len(tree) == d - 1 - i
    for leaf in tree:
      leaf_check(leaf)


def test_vinecop_pdf_full() -> None:
  d, n = 4, 200
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  )

  # keep_all=True (default): dict with the density + per-edge triangular arrays.
  full = cop.pdf_full(u)
  assert isinstance(full, dict)
  assert set(full) == {
    "pdf",
    "pdf_edges",
    "hfunc1",
    "hfunc2",
    "hfunc1_sub",
    "hfunc2_sub",
  }
  assert isinstance(full["pdf"], np.ndarray) and full["pdf"].shape == (n,)
  np.testing.assert_allclose(full["pdf"], cop.pdf(u), rtol=1e-10, atol=1e-12)

  # Triangular arrays are nested lists [tree][edge], with d - 1 edges in tree 0
  # and one fewer per subsequent tree. Per-edge densities are always length n;
  # h-functions are only filled where the cascade needs them (unneeded edges are
  # returned empty).
  def _is_len_n_vector(vec: object) -> None:
    assert isinstance(vec, np.ndarray) and vec.shape == (n,)

  def _is_edge_vector(vec: object) -> None:
    assert isinstance(vec, np.ndarray) and vec.shape in ((n,), (0,))

  _check_triangular(full["pdf_edges"], d, _is_len_n_vector)
  for key in ("hfunc1", "hfunc2", "hfunc1_sub", "hfunc2_sub"):
    _check_triangular(full[key], d, _is_edge_vector)

  # keep_all=False: only the density.
  simple = cop.pdf_full(u, keep_all=False)
  assert set(simple) == {"pdf"}
  np.testing.assert_allclose(simple["pdf"], cop.pdf(u), rtol=1e-10, atol=1e-12)

  # num_threads must not change the result.
  np.testing.assert_allclose(
    cop.pdf_full(u, num_threads=2)["pdf"], full["pdf"], rtol=1e-12, atol=1e-14
  )


def test_vinecop_scores_and_hessian() -> None:
  d, n = 4, 200
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  )

  # scores: one row per observation, one column per parameter.
  scores = cop.scores(u)
  assert isinstance(scores, np.ndarray) and scores.ndim == 2
  assert scores.shape[0] == n
  p = scores.shape[1]
  assert p == round(cop.npars)

  # scores_cov / hessian_avg: (p, p) matrices; scores_cov is symmetric.
  cov = cop.scores_cov(u)
  assert cov.shape == (p, p)
  np.testing.assert_allclose(cov, cov.T, rtol=1e-10, atol=1e-12)
  assert cop.hessian_avg(u).shape == (p, p)

  # step_wise=False and threading run and are consistent.
  assert cop.scores(u, step_wise=False).shape == (n, p)
  np.testing.assert_allclose(
    cop.scores(u, num_threads=2), scores, rtol=1e-10, atol=1e-12
  )

  # Per-observation hessian: nested lists [tree][edge], each a list of matrices.
  def _is_matrix_list(leaf: object) -> None:
    assert isinstance(leaf, list)
    for mat in leaf:
      assert isinstance(mat, np.ndarray) and mat.ndim == 2

  _check_triangular(cop.hessian(u), d, _is_matrix_list)


def test_struct_array_accessors_and_factory() -> None:
  d, n = 5, 200
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  )
  structure = cop.structure

  for natural_order in (False, True):
    # get_struct_array returns a [tree][edge] nested list matching the
    # per-entry accessor.
    arr = structure.get_struct_array(natural_order=natural_order)
    assert isinstance(arr, list) and len(arr) == d - 1
    for t, tree in enumerate(arr):
      assert isinstance(tree, list)
      assert len(tree) == d - 1 - t
      for e, entry in enumerate(tree):
        assert entry == structure.struct_array(
          t, e, natural_order=natural_order
        )

    # The vine's convenience accessor matches its structure's.
    assert cop.get_struct_array(natural_order=natural_order) == arr

    # from_struct_array(order, struct_array) rebuilds the same structure.
    rebuilt = pv.RVineStructure.from_struct_array(
      structure.order, arr, natural_order=natural_order
    )
    assert np.array_equal(rebuilt.matrix, structure.matrix)

  # Rows that do not form a triangular array are rejected.
  with pytest.raises(RuntimeError):
    pv.RVineStructure.from_struct_array(structure.order, [[1, 2], [3, 4]])
