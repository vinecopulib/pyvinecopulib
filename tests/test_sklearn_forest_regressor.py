"""Tests for VineForestRegressor."""

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("pandas")
pytest.importorskip("joblib")
pytest.importorskip("scipy")

from pyvinecopulib.sklearn import VineForestRegressor, VineRegressor  # noqa: E402


@pytest.fixture
def regression_setup(regression_data):
  X, y, true_coef, _ = regression_data
  X_train, X_test = X[:200], X[200:]
  y_train, y_test = y[:200], y[200:]

  def true_mean(x):
    return x @ true_coef

  return X_train, X_test, y_train, y_test, true_mean, true_coef


@pytest.fixture
def fitted_forest_regressor(regression_setup):
  X_train, X_test, y_train, y_test, true_mean, _ = regression_setup
  forest = VineForestRegressor(
    base_params={"mean": True, "batch_size": 50},
    n_vines=3,
    n_jobs=1,
    val_fraction=0.2,
    seed=42,
  )
  forest.fit(X_train, y_train)
  return forest, X_train, X_test, y_train, y_test, true_mean


def test_fit_properties(fitted_forest_regressor):
  forest, _, _, _, _, _ = fitted_forest_regressor
  assert hasattr(forest, "_estimators")
  assert hasattr(forest, "_estimators_logliks")
  assert len(forest._estimators) > 0
  assert len(forest._estimators) <= 4


def test_schema_inference_array(regression_setup):
  X_train, _, y_train, _, _, _ = regression_setup
  forest = VineForestRegressor(
    base_params={"mean": True, "batch_size": 100},
    n_vines=2,
    n_jobs=1,
    val_fraction=0.3,
    seed=42,
  )
  forest.fit(X_train, y_train)
  for estimator in forest._estimators:
    assert estimator._schema["kde1d_types"] == ["continuous", "continuous"]


def test_schema_propagation_dataframe(sample_dataframe_data, random_state):
  X_df, expected_expanded_cols = sample_dataframe_data
  y = random_state.randn(len(X_df))

  forest = VineForestRegressor(
    base_params={"mean": True}, n_vines=2, n_jobs=1, val_fraction=0.2
  )
  forest.fit(X_df, y)

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


def test_predict_shapes_and_types(fitted_forest_regressor):
  forest, _, X_test, _, _, _ = fitted_forest_regressor
  predictions = forest.predict(X_test)
  assert isinstance(predictions, np.ndarray)
  assert predictions.shape == (len(X_test),)
  assert np.all(np.isfinite(predictions))


def test_predict_with_quantiles(regression_setup):
  X_train, X_test, y_train, _, _, _ = regression_setup
  forest = VineForestRegressor(
    base_params={"mean": True, "quantiles": [0.1, 0.9]},
    n_vines=3,
    n_jobs=1,
    val_fraction=0.2,
  )
  forest.fit(X_train, y_train)
  predictions = forest.predict(X_test)

  assert predictions.shape == (len(X_test), 3)
  q10_pred = predictions[:, 1]
  q90_pred = predictions[:, 2]
  assert np.mean(q10_pred <= q90_pred) > 0.8


def test_prediction_accuracy(fitted_forest_regressor):
  forest, _, X_test, _, _, true_mean = fitted_forest_regressor
  pred_mean = forest.predict(X_test)
  true_mean_test = true_mean(X_test)

  correlation = np.corrcoef(pred_mean, true_mean_test)[0, 1]
  assert correlation > 0.8

  rmse = np.sqrt(np.mean((pred_mean - true_mean_test) ** 2))
  assert rmse < 0.5


def test_ensemble_vs_single_estimator(regression_setup):
  X_train, X_test, y_train, _, _, _ = regression_setup

  single = VineRegressor(mean=True, batch_size=50)
  single.fit(X_train, y_train)
  single_pred = single.predict(X_test)

  forest = VineForestRegressor(
    base_params={"mean": True, "batch_size": 50},
    n_vines=3,
    n_jobs=1,
    val_fraction=0.2,
    seed=42,
  )
  forest.fit(X_train, y_train)
  forest_pred = forest.predict(X_test)

  # Forest averages multiple structures; outputs should differ.
  assert not np.allclose(single_pred, forest_pred, rtol=0.1)
  correlation = np.corrcoef(single_pred, forest_pred)[0, 1]
  assert correlation > 0.5


def test_val_fraction_zero(regression_setup):
  X_train, _, y_train, _, _, _ = regression_setup
  forest = VineForestRegressor(
    base_params={"mean": True}, n_vines=3, val_fraction=0.0, n_jobs=1
  )
  forest.fit(X_train, y_train)
  assert len(forest._estimators) > 0


def test_best_only_option(regression_setup):
  X_train, _, y_train, _, _, _ = regression_setup

  forest_all = VineForestRegressor(
    base_params={"mean": True},
    n_vines=5,
    best_only=False,
    n_jobs=1,
    val_fraction=0.2,
  )
  forest_all.fit(X_train, y_train)

  forest_best = VineForestRegressor(
    base_params={"mean": True},
    n_vines=5,
    best_only=True,
    n_jobs=1,
    val_fraction=0.2,
    seed=42,
  )
  forest_best.fit(X_train, y_train)

  assert len(forest_best._estimators) <= len(forest_all._estimators)
  assert len(forest_best._estimators) == 1


@pytest.mark.parametrize("vines_sampling", ["uniform", "local"])
def test_reproducibility(regression_setup, vines_sampling):
  X_train, X_test, y_train, _, _, _ = regression_setup

  def create_and_fit():
    forest = VineForestRegressor(
      base_params={"mean": True},
      n_vines=3,
      vines_sampling=vines_sampling,
      seed=42,
      n_jobs=1,
      val_fraction=0.2,
    )
    forest.fit(X_train, y_train)
    return forest.predict(X_test)

  pred1 = create_and_fit()
  pred2 = create_and_fit()
  np.testing.assert_allclose(pred1, pred2, rtol=1e-12)


def test_uniform_vs_local(regression_setup):
  X_train, X_test, y_train, _, _, _ = regression_setup

  def create_and_fit(vines_sampling):
    forest = VineForestRegressor(
      base_params={"mean": True},
      n_vines=3,
      vines_sampling=vines_sampling,
      seed=42,
      n_jobs=1,
      val_fraction=0.2,
    )
    forest.fit(X_train, y_train)
    return forest.predict(X_test)

  pred_uniform = create_and_fit("uniform")
  pred_local = create_and_fit("local")
  assert not np.allclose(pred_uniform, pred_local, rtol=1e-10)


def test_normalize_weights_propagated(regression_setup):
  """The forest sets _normalize_weights=False on every base estimator."""
  X_train, _, y_train, _, _, _ = regression_setup
  forest = VineForestRegressor(
    base_params={"mean": True}, n_vines=2, n_jobs=1, val_fraction=0.2
  )
  forest.fit(X_train, y_train)
  for estimator in forest._estimators:
    assert estimator._normalize_weights is False


def test_wrong_dimensions(fitted_forest_regressor):
  forest, _, X_test, _, _, _ = fitted_forest_regressor
  X_wrong = X_test[:, :1]
  with pytest.raises(ValueError):
    forest.predict(X_wrong)


def test_loglik_estimator_method(fitted_forest_regressor):
  forest, _, X_test, _, y_test, _ = fitted_forest_regressor
  estimator = forest._estimators[0]
  loglik = forest._loglik_estimator(estimator, X_test, y_test)
  assert isinstance(loglik, np.ndarray)
  assert loglik.shape == (len(X_test),)
  assert np.all(np.isfinite(loglik))


def test_loglik_estimator_requires_y(fitted_forest_regressor):
  forest, _, X_test, _, _, _ = fitted_forest_regressor
  estimator = forest._estimators[0]
  with pytest.raises(ValueError, match="requires y"):
    forest._loglik_estimator(estimator, X_test, None)


def test_parallel_vs_sequential(regression_setup):
  X_train, X_test, y_train, _, _, _ = regression_setup

  forest_seq = VineForestRegressor(
    base_params={"mean": True}, n_vines=3, seed=42, n_jobs=1, val_fraction=0.2
  )
  forest_seq.fit(X_train, y_train)
  pred_seq = forest_seq.predict(X_test)

  forest_par = VineForestRegressor(
    base_params={"mean": True}, n_vines=3, seed=42, n_jobs=2, val_fraction=0.2
  )
  forest_par.fit(X_train, y_train)
  pred_par = forest_par.predict(X_test)

  np.testing.assert_allclose(pred_seq, pred_par, rtol=1e-10)
