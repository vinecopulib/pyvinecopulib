import copy
import itertools
import os
import pickle

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

  # Test sample method
  simulated_data = cop.sample(n)
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

  # CBOR round-trip: a ``.cbor`` filename selects the binary format
  # (vinecopulib#684).
  cbor_filename = filename.removesuffix(".json") + ".cbor"
  cop.to_file(cbor_filename)
  compare_vinecop(cop, pv.Vinecop.from_file(cbor_filename))
  with open(cbor_filename, "rb") as f:
    assert f.read(1) != b"{"


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


def test_custom_criterion_rejected_when_not_custom() -> None:
  # A function set while tree_criterion != "custom" would be silently unused,
  # so it is rejected up front rather than ignored.
  d, n = 5, 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  controls = pv.FitControlsVinecop(tree_criterion="tau")
  controls.tree_criterion_function = _raising_criterion
  with pytest.raises(RuntimeError, match="tree_criterion"):
    pv.Vinecop.from_data(u, controls=controls)


def test_custom_criterion_exception_propagates() -> None:
  # A Python exception raised inside the criterion surfaces from from_data.
  d, n = 5, 1000
  u = pv.to_pseudo_obs(random_data(d, n))
  controls = pv.FitControlsVinecop(tree_criterion="custom", num_threads=1)
  controls.tree_criterion_function = _raising_criterion
  with pytest.raises(ValueError):
    pv.Vinecop.from_data(u, controls=controls)


def test_custom_criterion_multithread_smoke() -> None:
  # A custom criterion is evaluated on the thread that entered select(), so
  # num_threads > 1 must still give the same result.
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
  for key in ("hfunc1", "hfunc2"):
    _check_triangular(full[key], d, _is_edge_vector)

  # The left-limit h-functions are only computed for models with discrete
  # variables; for an all-continuous model they come back empty
  # (vinecopulib#692).
  assert full["hfunc1_sub"] == []
  assert full["hfunc2_sub"] == []

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

  # scores_cov / hessian: (p, p) matrices; scores_cov is symmetric.
  cov = cop.scores_cov(u)
  assert cov.shape == (p, p)
  np.testing.assert_allclose(cov, cov.T, rtol=1e-10, atol=1e-12)
  assert cop.hessian(u).shape == (p, p)

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

  _check_triangular(cop.hessian_full(u), d, _is_matrix_list)


def test_vinecop_gradient() -> None:
  d, n = 4, 200
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  )
  p = round(cop.npars)

  # The gradient is the observation-average of the scores, for both step_wise
  # modes, and threading must not change the result.
  for step_wise in (True, False):
    grad = cop.gradient(u, step_wise=step_wise)
    assert isinstance(grad, np.ndarray) and grad.shape == (p,)
    np.testing.assert_allclose(
      grad,
      cop.scores(u, step_wise=step_wise).mean(axis=0),
      rtol=1e-10,
      atol=1e-12,
    )
  np.testing.assert_allclose(
    cop.gradient(u, num_threads=2), cop.gradient(u), rtol=1e-10, atol=1e-12
  )


_SCORES_FULL_CACHE_KEYS = (
  "pdf_edges",
  "logpdf_deriv_pars",
  "hfunc1_deriv_pars",
  "hfunc2_deriv_pars",
  "logpdf_deriv_u1",
  "logpdf_deriv_u2",
  "hfunc1_deriv_u1",
  "hfunc2_deriv_u2",
)


def test_vinecop_scores_full() -> None:
  d, n = 4, 200
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  )
  p = round(cop.npars)

  # keep_all=True (default): dict with the scores + per-edge caches.
  for step_wise in (True, False):
    res = cop.scores_full(u, step_wise=step_wise)
    assert isinstance(res, dict)
    assert set(res) == {"scores", *_SCORES_FULL_CACHE_KEYS}
    assert res["scores"].shape == (n, p)
    np.testing.assert_allclose(
      res["scores"],
      cop.scores(u, step_wise=step_wise),
      rtol=1e-10,
      atol=1e-12,
    )

  # The step-wise path only fills the per-parameter log-density derivatives
  # (the step-wise scores themselves); the full gradient additionally fills
  # the caches consumed by the chain rule through the vine.
  def _is_vector_list(leaf: object) -> None:
    assert isinstance(leaf, list)
    for vec in leaf:
      assert isinstance(vec, np.ndarray) and vec.shape == (n,)

  res_sw = cop.scores_full(u, step_wise=True)
  _check_triangular(res_sw["logpdf_deriv_pars"], d, _is_vector_list)
  assert res_sw["pdf_edges"] == []

  res_full = cop.scores_full(u, step_wise=False)
  _check_triangular(res_full["logpdf_deriv_pars"], d, _is_vector_list)
  assert res_full["pdf_edges"] != []

  # keep_all=False: only the scores.
  simple = cop.scores_full(u, keep_all=False)
  assert set(simple) == {"scores"}

  # Threading must not change the result.
  np.testing.assert_allclose(
    cop.scores_full(u, num_threads=2)["scores"],
    res_sw["scores"],
    rtol=1e-12,
    atol=1e-14,
  )


def test_vinecop_scores_discrete_finite_differences() -> None:
  # Models with discrete variables fall back to finite differences: scores /
  # gradient / hessian still work, and the scores_full caches stay empty.
  rng = np.random.RandomState(7)
  n = 200
  x = rng.binomial(4, 0.5, size=n)
  pmf = np.array([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0
  cdf = np.cumsum(pmf)
  u1 = cdf[x]
  u1_sub = np.where(x > 0, cdf[x - 1], 0.0)
  u2 = rng.uniform(size=n)
  u = np.column_stack([u1, u2, u1_sub, u2])

  bicop = pv.Bicop(
    family=pv.families.gaussian,
    parameters=np.array([[0.5]]),
    var_types=["d", "c"],
  )
  cop = pv.Vinecop.from_structure(
    matrix=np.array([[1, 1], [2, 0]]),
    pair_copulas=[[bicop]],
    var_types=["d", "c"],
  )

  scores = cop.scores(u)
  assert scores.shape == (n, 1)
  assert np.isfinite(scores).all()
  assert cop.hessian(u).shape == (1, 1)

  # Step-wise: the per-edge scores cache is filled (with the discrete edge's
  # finite-difference values); the full-gradient-only caches stay empty.
  res = cop.scores_full(u)
  np.testing.assert_allclose(res["scores"], scores, rtol=1e-12)
  assert res["logpdf_deriv_pars"] != []
  for key in _SCORES_FULL_CACHE_KEYS:
    if key != "logpdf_deriv_pars":
      assert res[key] == []

  # Full gradient: whole-vine finite differences, all caches empty.
  res_full = cop.scores_full(u, step_wise=False)
  assert np.isfinite(res_full["scores"]).all()
  for key in _SCORES_FULL_CACHE_KEYS:
    assert res_full[key] == []


def test_vinecop_scores_rejects_nonparametric() -> None:
  rng = np.random.RandomState(11)
  u = pv.to_pseudo_obs(rng.normal(size=(300, 2)))
  cop = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.tll])
  )
  for call in (
    cop.scores,
    cop.gradient,
    cop.scores_full,
    cop.hessian,
    cop.hessian_full,
  ):
    with pytest.raises(RuntimeError, match="parametric"):
      call(u)


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


def test_vinecop_sample_conditional() -> None:
  # Conditional sampling given fixed values of a subset of variables
  # (vinecopulib#696). The conditioning variables are the last k of the vine
  # order; their natural columns in the output reproduce u_cond exactly.
  d, n = 5, 700
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  )
  order = [int(v) for v in cop.structure.order]

  for k in (1, 2):
    u_cond = np.random.RandomState(k).uniform(0.05, 0.95, size=(50, k))
    sim = cop.sample_conditional(u_cond, seeds=[1, 2, 3])
    assert sim.shape == (50, d)
    cond_cols = [v - 1 for v in order[-k:]]
    np.testing.assert_allclose(
      sim[:, cond_cols], u_cond, rtol=1e-12, atol=1e-12
    )
    # Deterministic given seeds; all columns are valid uniforms.
    sim2 = cop.sample_conditional(u_cond, seeds=[1, 2, 3])
    np.testing.assert_allclose(sim, sim2, rtol=1e-12, atol=1e-14)
    assert np.all((sim >= 0.0) & (sim <= 1.0))


def test_vinecop_sample_conditional_discrete() -> None:
  # Discrete conditioning variables are supported: u_cond carries an extra
  # left-limit F(x^-) column per discrete conditioner, and the conditioning
  # column of the output lands within the atom [F(x^-), F(x)] (vinecopulib#696).
  rng = np.random.RandomState(11)
  n = 800
  x = rng.binomial(4, 0.5, size=n)
  pmf = np.array([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0
  cdf = np.cumsum(pmf)
  u1, u1_sub = cdf[x], np.where(x > 0, cdf[x - 1], 0.0)
  cont = pv.to_pseudo_obs(rng.normal(size=(n, 2)))
  # var_types = [d, c, c]; data is the main columns then the discrete left-limit.
  data = np.column_stack([u1, cont, u1_sub])

  # Force the discrete variable (1) to the order tail so we can condition on it.
  fc = pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  fc.conditioning_set = [1]
  cop = pv.Vinecop.from_data(data, var_types=["d", "c", "c"], controls=fc)
  assert [int(v) for v in cop.structure.order][-1] == 1

  # Condition var 1 on the atom x = 2: u_cond = [F(2), F(1)] (value, left-limit).
  m = 20
  u_cond = np.tile([[cdf[2], cdf[1]]], (m, 1))
  sim = cop.sample_conditional(u_cond, seeds=[1, 2, 3])
  assert sim.shape == (m, 3)
  # var 1 is natural column 0; the draw lands within the conditioned atom.
  col1 = sim[:, 0]
  assert np.all((col1 >= cdf[1] - 1e-9) & (col1 <= cdf[2] + 1e-9))


def test_vinecop_reorient() -> None:
  # reorient relabels to an equivalent vine whose order tail equals the given
  # set, without refitting (value-preserving: pdf / loglik invariant) (#697).
  d, n = 5, 700
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(u)
  ll0, pdf0 = cop.loglik(u), cop.pdf(u)
  order = [int(v) for v in cop.structure.order]

  # Exercise a GENUINE reorientation: search for an admissible tail that is not
  # already the current suffix (a suffix hits the no-op early return and would
  # test nothing). Inadmissible sets raise, exercising the rejection path.
  did_reorient = did_reject = False
  for k in (1, 2):
    for cand in itertools.combinations(range(1, d + 1), k):
      if set(cand) == set(order[-k:]):  # same tail set -> (near) no-op, skip
        continue
      vc = copy.deepcopy(cop)
      try:
        vc.reorient(list(cand))
      except RuntimeError:
        did_reject = True  # set is not an admissible sampling-order tail
        continue
      new_order = [int(v) for v in vc.structure.order]
      assert set(new_order[-k:]) == set(cand)  # requested set placed at tail
      assert new_order != order  # order genuinely changed (not a no-op)
      # Value-preserving: the reoriented vine has the same density / loglik.
      np.testing.assert_allclose(vc.loglik(u), ll0, rtol=1e-9, atol=1e-10)
      np.testing.assert_allclose(vc.pdf(u), pdf0, rtol=1e-8, atol=1e-10)
      did_reorient = True
  assert did_reorient, "expected at least one admissible non-suffix tail"
  assert did_reject, "expected at least one inadmissible set to be rejected"

  # Out-of-range variables are rejected by input validation.
  with pytest.raises((RuntimeError, ValueError)):
    copy.deepcopy(cop).reorient([d + 1])


def test_vinecop_get_trees_and_roundtrip() -> None:
  # Vinecop.get_trees(): readable [tree][edge] dicts carrying the fitted
  # pair-copulas; RVineStructure.get_trees() + from_trees() round-trip
  # matrix-exactly (vinecopulib#698).
  d, n = 5, 700
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(
    u,
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.gaussian, pv.families.clayton]
    ),
  )

  trees = cop.get_trees()
  assert isinstance(trees, list) and len(trees) == d - 1
  for t, tree in enumerate(trees):
    assert len(tree) == d - 1 - t
    for edge in tree:
      assert set(edge) == {"conditioned", "conditioning", "pair_copula"}
      a, b = edge["conditioned"]
      assert isinstance(a, int) and isinstance(b, int)
      assert isinstance(edge["conditioning"], list)
      assert isinstance(edge["pair_copula"], pv.Bicop)

  # get_trees() carries the FITTED pair copulas (not independence placeholders):
  # the family and Kendall's-tau multisets match the vine's own accessors.
  ft_fams = sorted(str(e["pair_copula"].family) for tr in trees for e in tr)
  vine_fams = sorted(str(f) for tr in cop.families for f in tr)
  assert ft_fams == vine_fams
  ft_taus = sorted(round(e["pair_copula"].tau, 10) for tr in trees for e in tr)
  vine_taus = sorted(round(t, 10) for tr in cop.taus for t in tr)
  assert ft_taus == vine_taus

  # Structure round-trip via from_trees (matrix-exactness is covered thoroughly
  # in test_rvinestructure_get_trees_faithful_roundtrip).
  s = cop.structure
  assert pv.RVineStructure.from_trees(s.dim, s.get_trees()) == s
  # __eq__ / __ne__ are genuine equality: a truncated copy compares unequal.
  st = copy.deepcopy(s)
  st.truncate(1)
  assert st.trunc_lvl == 1
  assert s == s and s != st


def test_vinecop_per_observation_parameters() -> None:
  # Per-observation-parameter overloads (vinecopulib#699): passing the fitted
  # parameter vector broadcast over rows reproduces the fixed-parameter result.
  d, n = 4, 600
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(
    u,
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.student], preselect_families=False
    ),
  )
  # Flatten fitted parameters in scores()' (tree, edge, parameter) column order.
  flat = np.concatenate(
    [np.asarray(edge).ravel() for tree in cop.parameters for edge in tree]
  )
  assert flat.shape[0] == cop.scores(u).shape[1]
  pars = np.tile(flat, (n, 1))

  np.testing.assert_allclose(
    cop.pdf(u, parameters=pars), cop.pdf(u), rtol=1e-8, atol=1e-10
  )
  np.testing.assert_allclose(
    cop.loglik(u, parameters=pars), cop.loglik(u), rtol=1e-8, atol=1e-10
  )
  for method in ("scores", "gradient", "hessian", "scores_cov"):
    fn = getattr(cop, method)
    np.testing.assert_allclose(
      fn(u, parameters=pars), fn(u), rtol=1e-8, atol=1e-10
    )
  np.testing.assert_allclose(
    cop.pdf_full(u, parameters=pars)["pdf"], cop.pdf(u), rtol=1e-8, atol=1e-10
  )
  np.testing.assert_allclose(
    cop.scores_full(u, parameters=pars)["scores"],
    cop.scores(u),
    rtol=1e-8,
    atol=1e-10,
  )
  # hessian_full (per-observation Hessians) also honors the parameters overload.
  hf_p, hf_0 = cop.hessian_full(u, parameters=pars), cop.hessian_full(u)
  np.testing.assert_allclose(hf_p[0][0], hf_0[0][0], rtol=1e-8, atol=1e-10)

  # Existing positional signatures still work (parameters appended last).
  np.testing.assert_allclose(cop.pdf(u, 2), cop.pdf(u), rtol=1e-12)
  np.testing.assert_allclose(
    cop.scores(u, False, 2), cop.scores(u, step_wise=False), rtol=1e-10
  )


def test_fit_controls_vinecop_conditioning_set() -> None:
  # conditioning_set property + pickle round-trip, and conditioning-aware
  # selection places the set at the tail of the vine order (vinecopulib#697).
  fc = pv.FitControlsVinecop()
  assert fc.conditioning_set == []
  fc.conditioning_set = [2, 3]
  assert fc.conditioning_set == [2, 3]
  fc2 = pickle.loads(pickle.dumps(fc))
  assert fc2.conditioning_set == [2, 3]

  d, n = 5, 700
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(u, controls=fc)
  order = [int(v) for v in cop.structure.order]
  assert set(order[-2:]) == {2, 3}


def test_inverse_rosenblatt_is_thread_safe() -> None:
  # vinecopulib#712 removed the var_types round-trip that mutated the model
  # during evaluation, so concurrent calls on one object are now safe.
  import concurrent.futures

  d, n = 4, 200
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(u)
  expected = cop.inverse_rosenblatt(u)

  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    futures = [pool.submit(cop.inverse_rosenblatt, u) for _ in range(8)]
    results = [f.result() for f in futures]

  for got in results:
    np.testing.assert_array_equal(got, expected)


def test_from_data_without_structure_or_matrix() -> None:
  # vinecopulib made the data constructor's `matrix` argument required, since
  # `Vinecop(data)` was ambiguous. `from_data` passes one explicitly, so
  # omitting both structure and matrix must still select a structure.
  d, n = 4, 300
  u = pv.to_pseudo_obs(random_data(d, n))
  cop = pv.Vinecop.from_data(u)
  assert cop.dim == d
  assert set(int(v) for v in cop.structure.order) == set(range(1, d + 1))


def test_vinecop_get_trees_feeds_rvinestructure_from_trees() -> None:
  # ``Vinecop.get_trees()`` labels each edge with a mapping rather than the
  # triple ``RVineStructure.get_trees()`` uses; both spell the same
  # decomposition, so both must reassemble the same structure.
  rng = np.random.default_rng(3)
  d = 5
  u = pv.to_pseudo_obs(
    rng.standard_normal((400, d)) @ rng.standard_normal((d, d))
  )
  vine = pv.Vinecop.from_data(u)

  assert (
    pv.RVineStructure.from_trees(vine.dim, vine.get_trees()) == vine.structure
  )


def test_vinecop_pair_copulas_is_settable() -> None:
  rng = np.random.default_rng(4)
  d = 4
  u = pv.to_pseudo_obs(
    rng.standard_normal((400, d)) @ rng.standard_normal((d, d))
  )
  vine = pv.Vinecop.from_data(u)
  fitted = vine.pair_copulas

  vine.pair_copulas = [
    [pv.Bicop.from_family(pv.families.indep) for _ in tree] for tree in fitted
  ]
  np.testing.assert_allclose(vine.pdf(u), 1.0)

  vine.pair_copulas = fitted
  assert not np.allclose(vine.pdf(u), 1.0)

  with pytest.raises(RuntimeError):
    vine.pair_copulas = [[pv.Bicop.from_family(pv.families.indep)] * 9]


def test_fit_controls_vinecop_from_bicop_controls() -> None:
  bicop_controls = pv.FitControlsBicop(
    family_set=[pv.families.gaussian], num_threads=2, psi0=0.5
  )
  controls = pv.FitControlsVinecop.from_bicop_controls(
    bicop_controls, trunc_lvl=2, threshold=0.1
  )

  assert controls.family_set == [pv.families.gaussian]
  assert controls.num_threads == 2
  assert controls.psi0 == 0.5
  assert controls.trunc_lvl == 2
  assert controls.threshold == 0.1
  assert controls.bicop_controls.family_set == [pv.families.gaussian]

  controls.bicop_controls = pv.FitControlsBicop(
    family_set=[pv.families.clayton]
  )
  assert controls.family_set == [pv.families.clayton]


def _cop_and_data(d: int = 5, n: int = 400):
  u = pv.to_pseudo_obs(random_data(d, n))
  return pv.Vinecop.from_data(u), u


def test_rosenblatt_conditioning_set_order_tail_is_identity() -> None:
  # Conditioning on the set the vine order already ends with is the identity:
  # upstream evaluates it on the original representation.
  cop, u = _cop_and_data()
  tail = [int(v) for v in cop.structure.order][-2:]

  np.testing.assert_array_equal(
    cop.rosenblatt(u, conditioning_set=tail), cop.rosenblatt(u)
  )
  w = cop.rosenblatt(u)
  np.testing.assert_array_equal(
    cop.inverse_rosenblatt(w, conditioning_set=tail),
    cop.inverse_rosenblatt(w),
  )


def test_rosenblatt_conditioning_set_round_trip() -> None:
  # The inverse transform undoes the forward one under the same set.
  cop, u = _cop_and_data()
  tail = [int(v) for v in cop.structure.order][-2:]
  w = cop.rosenblatt(u, conditioning_set=tail)
  back = cop.inverse_rosenblatt(w, conditioning_set=tail)
  np.testing.assert_allclose(back, u, rtol=1e-8, atol=1e-8)


def test_rosenblatt_conditioning_set_matches_reorient() -> None:
  # The argument applies the same relabeling as `reorient`, without mutating
  # the model. Not every set is an admissible sampling-order tail, so fit a
  # vine that admits this one.
  cs = [2, 3]
  controls = pv.FitControlsVinecop()
  controls.conditioning_set = cs
  u = pv.to_pseudo_obs(random_data(5, 400))
  cop = pv.Vinecop.from_data(u, controls=controls)

  order_before = [int(v) for v in cop.structure.order]
  with_arg = cop.rosenblatt(u, conditioning_set=cs)

  # The keyword form evaluates through a view and leaves the model alone;
  # `reorient` performs the same relabeling in place.
  assert [int(v) for v in cop.structure.order] == order_before
  cop.reorient(cs)
  np.testing.assert_allclose(
    with_arg, cop.rosenblatt(u), rtol=1e-10, atol=1e-10
  )


def test_rosenblatt_second_positional_is_still_num_threads() -> None:
  # C++ puts `conditioning_set` in position 2; Python keeps `num_threads`
  # there, so `rosenblatt(u, 4)` must not silently become a conditioning set.
  cop, u = _cop_and_data(d=4, n=200)
  np.testing.assert_allclose(
    cop.rosenblatt(u, 4), cop.rosenblatt(u), rtol=1e-12
  )
  np.testing.assert_allclose(
    cop.inverse_rosenblatt(u, 2), cop.inverse_rosenblatt(u), rtol=1e-12
  )


@pytest.mark.parametrize(
  "bad", [[], [1, 1], [0], [99], [1, 2, 3, 4, 5]], ids=str
)
def test_rosenblatt_rejects_inadmissible_conditioning_set(bad) -> None:
  cop, u = _cop_and_data()
  with pytest.raises(RuntimeError):
    cop.rosenblatt(u, conditioning_set=bad)


def test_rosenblatt_conditioning_set_rejects_truncated_vine() -> None:
  # The reorientation needs every tree, so a truncated vine cannot serve an
  # arbitrary conditioning set even when the set itself is admissible.
  u = pv.to_pseudo_obs(random_data(5, 400))
  cop = pv.Vinecop.from_data(u, controls=pv.FitControlsVinecop(trunc_lvl=2))

  with pytest.raises(RuntimeError):
    cop.rosenblatt(u, conditioning_set=[1, 2])
  with pytest.raises(RuntimeError):
    cop.inverse_rosenblatt(u, conditioning_set=[1, 2])


def _cop_conditioned_on(cs: list[int], d: int = 5, n: int = 400):
  """A vine whose order tail is ``cs``, so ``cs`` is an admissible set."""
  controls = pv.FitControlsVinecop()
  controls.conditioning_set = cs
  u = pv.to_pseudo_obs(random_data(d, n))
  return pv.Vinecop.from_data(u, controls=controls), u


def test_sample_conditional_explicit_set_matches_implicit() -> None:
  # Given in the tail's own order, the explicit form reproduces the implicit
  # one exactly (vinecopulib#729).
  cs = [2, 3]
  cop, u = _cop_conditioned_on(cs)
  tail = [int(v) for v in cop.structure.order][-2:]
  u_cond = u[:20, [t - 1 for t in tail]]

  implicit = cop.sample_conditional(u_cond, seeds=[1, 2])
  explicit = cop.sample_conditional(u_cond, seeds=[1, 2], conditioning_set=tail)
  np.testing.assert_array_equal(implicit, explicit)


def test_sample_conditional_explicit_set_column_mapping() -> None:
  # The two forms map u_cond's columns differently: implicitly by the order
  # tail, explicitly by the given set. Reversing the set must permute the
  # columns it consumes, not be ignored.
  cs = [2, 3]
  cop, u = _cop_conditioned_on(cs)
  tail = [int(v) for v in cop.structure.order][-2:]
  u_cond = u[:20, [t - 1 for t in tail]]

  forward = cop.sample_conditional(u_cond, seeds=[3], conditioning_set=tail)
  reversed_ = cop.sample_conditional(
    u_cond[:, ::-1], seeds=[3], conditioning_set=tail[::-1]
  )
  np.testing.assert_allclose(forward, reversed_, rtol=1e-10, atol=1e-10)


def test_sample_conditional_does_not_mutate_the_model() -> None:
  cs = [2, 3]
  cop, u = _cop_conditioned_on(cs)
  tail = [int(v) for v in cop.structure.order][-2:]
  order_before = [int(v) for v in cop.structure.order]
  cop.sample_conditional(u[:10, [t - 1 for t in tail]], conditioning_set=tail)
  assert [int(v) for v in cop.structure.order] == order_before


def test_sample_conditional_rejects_inadmissible_set() -> None:
  cop, u = _cop_and_data()
  with pytest.raises(RuntimeError):
    cop.sample_conditional(u[:10, :2], conditioning_set=[])


def test_sample_conditional_discrete_set_takes_a_left_limit_column() -> None:
  # A discrete conditioning variable is described by two columns, its cdf and
  # its left limit, so ``u_cond`` is wider than the conditioning set.
  rng = np.random.default_rng(11)
  d, n, m = 4, 600, 32

  continuous = rng.standard_normal((n, d)) @ rng.standard_normal((d, d))
  counts = np.floor(4.0 * pv.to_pseudo_obs(continuous[:, [d - 1]]))
  data = np.column_stack([continuous[:, : d - 1], counts])

  # Expanded layout: d columns, then a left-limit column per discrete variable.
  ecdf = pv.to_pseudo_obs(data)
  left = pv.to_pseudo_obs(np.column_stack([counts - 1.0]))
  u = np.column_stack([ecdf, left])

  var_types = ["c"] * (d - 1) + ["d"]
  controls = pv.FitControlsVinecop()
  controls.conditioning_set = [d]
  cop = pv.Vinecop.from_data(u, var_types=var_types, controls=controls)
  assert int(cop.structure.order[-1]) == d

  u_cond = np.column_stack([u[:m, d - 1], u[:m, d]])
  drawn = cop.sample_conditional(u_cond, seeds=[1, 2, 3])

  assert drawn.shape == (m, d)
  # The conditioning column is reproduced, not resampled.
  np.testing.assert_allclose(drawn[:, -1], u_cond[:, 0], rtol=1e-12)

  # Dropping the left-limit column leaves the model under-specified.
  with pytest.raises(RuntimeError):
    cop.sample_conditional(u_cond[:, :1], seeds=[1, 2, 3])
