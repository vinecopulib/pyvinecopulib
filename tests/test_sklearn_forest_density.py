"""Tests for VineForestDensity."""

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("pandas")
pytest.importorskip("joblib")
pytest.importorskip("scipy")

from pyvinecopulib.sklearn import VineDensity, VineForestDensity  # noqa: E402


@pytest.fixture
def fitted_forest_density(sample_array_data):
  X, _, _ = sample_array_data
  forest = VineForestDensity(
    base_params={"batch_size": 50},
    n_vines=3,
    n_jobs=1,
    val_fraction=0.2,
    seed=42,
  )
  forest.fit(X)
  return forest, X


def test_fit_properties(fitted_forest_density):
  forest, _ = fitted_forest_density
  assert hasattr(forest, "_estimators")
  assert hasattr(forest, "_estimators_logliks")
  assert len(forest._estimators) > 0
  # n_vines + (1 if add_dissmann).
  assert len(forest._estimators) <= 4


def test_schema_inference_array(sample_array_data):
  """All base estimators should agree on the auto-inferred schema."""
  X, _, _ = sample_array_data
  forest = VineForestDensity(
    base_params={"batch_size": 100},
    n_vines=2,
    n_jobs=1,
    val_fraction=0.3,
    seed=42,
  )
  forest.fit(X)
  for estimator in forest._estimators:
    assert estimator._schema["kde1d_types"] == ["continuous", "continuous"]


def test_schema_propagation_dataframe(sample_dataframe_data):
  X_df, expected_expanded_cols = sample_dataframe_data
  forest = VineForestDensity(n_vines=2, n_jobs=1, val_fraction=0.2)
  forest.fit(X_df)

  expected_types = [
    "continuous",
    "continuous",
    "discrete",
    "discrete",
    "continuous",
  ]
  for estimator in forest._estimators:
    assert estimator._schema["kde1d_types"] == expected_types
    if hasattr(estimator, "_expanded_columns"):
      assert estimator._expanded_columns == expected_expanded_cols


def test_score_samples_shapes_and_types(fitted_forest_density):
  forest, X = fitted_forest_density
  X_test = X[:20]
  log_scores = forest.score_samples(X_test)
  assert isinstance(log_scores, np.ndarray)
  assert log_scores.shape == (20,)
  assert np.all(np.isfinite(log_scores))


def test_pdf_shapes_and_types(fitted_forest_density):
  forest, X = fitted_forest_density
  X_test = X[:20]
  densities = forest.pdf(X_test)
  assert isinstance(densities, np.ndarray)
  assert densities.shape == (20,)
  assert np.all(densities > 0)
  assert np.all(np.isfinite(densities))


def test_pdf_log_consistency(fitted_forest_density):
  forest, X = fitted_forest_density
  X_test = X[:20]
  log_scores = forest.score_samples(X_test)
  densities = forest.pdf(X_test)
  np.testing.assert_allclose(densities, np.exp(log_scores), rtol=1e-10)


def test_score_method(fitted_forest_density):
  forest, X = fitted_forest_density
  X_test = X[:20]
  score = forest.score(X_test)
  log_scores = forest.score_samples(X_test)
  assert isinstance(score, (float, np.floating))
  np.testing.assert_allclose(score, np.mean(log_scores))


@pytest.mark.parametrize("n_samples", [1, 5, 10])
def test_sample_method(fitted_forest_density, n_samples):
  forest, _ = fitted_forest_density
  samples = forest.sample(n_samples)
  assert isinstance(samples, np.ndarray)
  assert samples.shape == (n_samples, 2)
  assert np.all(np.isfinite(samples))


def test_sample_from_estimator(fitted_forest_density):
  forest, _ = fitted_forest_density
  n_estimators = len(forest._estimators)
  for i in range(n_estimators):
    samples = forest.sample_from_estimator(i, 5)
    assert isinstance(samples, np.ndarray)
    assert samples.shape == (5, 2)
    assert np.all(np.isfinite(samples))

  with pytest.raises(ValueError):
    forest.sample_from_estimator(n_estimators, 1)


def test_ensemble_vs_single_estimator(sample_array_data):
  """Forest predictions should correlate with a single-vine baseline."""
  X, _, _ = sample_array_data
  X_test = X[:50]

  single = VineDensity(batch_size=50)
  single.fit(X)
  single_scores = single.score_samples(X_test)

  forest = VineForestDensity(
    base_params={"batch_size": 50},
    n_vines=4,
    n_jobs=1,
    val_fraction=0.3,
    seed=123,
  )
  forest.fit(X)
  forest_scores = forest.score_samples(X_test)

  assert np.all(np.isfinite(single_scores))
  assert np.all(np.isfinite(forest_scores))
  correlation = np.corrcoef(single_scores, forest_scores)[0, 1]
  assert correlation > 0.3


def test_val_fraction_zero(sample_array_data):
  X, _, _ = sample_array_data
  forest = VineForestDensity(n_vines=3, val_fraction=0.0, n_jobs=1)
  forest.fit(X)
  assert len(forest._estimators) > 0


def test_best_only_option(sample_array_data):
  X, _, _ = sample_array_data

  forest_all = VineForestDensity(
    n_vines=5, best_only=False, n_jobs=1, val_fraction=0.2
  )
  forest_all.fit(X)

  forest_best = VineForestDensity(
    n_vines=5, best_only=True, n_jobs=1, val_fraction=0.2, seed=42
  )
  forest_best.fit(X)

  assert len(forest_best._estimators) <= len(forest_all._estimators)
  assert len(forest_best._estimators) == 1


@pytest.mark.parametrize("vines_sampling", ["uniform", "local"])
def test_reproducibility(sample_array_data, vines_sampling):
  X, _, _ = sample_array_data
  X_test = X[:20]

  def create_and_fit():
    forest = VineForestDensity(
      n_vines=3,
      vines_sampling=vines_sampling,
      seed=42,
      n_jobs=1,
      val_fraction=0.2,
    )
    forest.fit(X)
    return forest.score_samples(X_test)

  scores1 = create_and_fit()
  scores2 = create_and_fit()
  np.testing.assert_allclose(scores1, scores2, rtol=1e-12)


def test_uniform_vs_local(sample_array_data):
  X, _, _ = sample_array_data
  X_test = X[:20]

  def create_and_fit(vines_sampling):
    forest = VineForestDensity(
      n_vines=3,
      vines_sampling=vines_sampling,
      seed=42,
      n_jobs=1,
      val_fraction=0.2,
    )
    forest.fit(X)
    return forest.score_samples(X_test)

  scores_uniform = create_and_fit("uniform")
  scores_local = create_and_fit("local")
  assert not np.allclose(scores_uniform, scores_local, rtol=1e-10)


def test_parallel_vs_sequential(sample_array_data):
  X, _, _ = sample_array_data
  X_test = X[:20]

  forest_seq = VineForestDensity(n_vines=3, seed=42, n_jobs=1, val_fraction=0.2)
  forest_seq.fit(X)
  scores_seq = forest_seq.score_samples(X_test)

  forest_par = VineForestDensity(n_vines=3, seed=42, n_jobs=2, val_fraction=0.2)
  forest_par.fit(X)
  scores_par = forest_par.score_samples(X_test)

  np.testing.assert_allclose(scores_seq, scores_par, rtol=1e-10)
