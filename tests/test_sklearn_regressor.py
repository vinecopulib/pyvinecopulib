"""Tests for VineRegressor class."""

import numpy as np
import pytest

pytest.importorskip("sklearn")

from pyvinecopulib.sklearn import VineRegressor  # noqa: E402


@pytest.fixture
def regression_setup(regression_data):
  """Setup regression data and split."""
  X, y, true_coef, noise_std = regression_data
  X_train, X_test = X[:200], X[200:]
  y_train, y_test = y[:200], y[200:]

  def true_mean(x):
    return x @ true_coef

  return X_train, X_test, y_train, y_test, true_mean, true_coef


@pytest.fixture
def fitted_regressor(regression_setup):
  """Fitted VineRegressor for testing."""
  X_train, X_test, y_train, y_test, true_mean, true_coef = regression_setup
  regressor = VineRegressor(
    mean=True
  )  # Only mean, no quantiles for simpler testing
  regressor.fit(X_train, y_train)
  return regressor, X_train, X_test, y_train, y_test, true_mean


def test_fit_properties(fitted_regressor):
  """Test properties after fitting."""
  regressor, _, _, _, _, _ = fitted_regressor
  assert hasattr(regressor, "_vine")
  assert hasattr(regressor, "_x_kde1d")
  assert hasattr(regressor, "_y_kde1d")
  assert hasattr(regressor, "_y_train")
  assert regressor.n_features_in_ == 2


def test_predict_mean_only(regression_setup):
  """Test prediction with mean only."""
  X_train, X_test, y_train, y_test, _, _ = regression_setup
  regressor = VineRegressor(mean=True)
  regressor.fit(X_train, y_train)
  pred_mean = regressor.predict(X_test)

  assert isinstance(pred_mean, np.ndarray)
  assert pred_mean.shape == (len(X_test),)
  assert np.all(np.isfinite(pred_mean))


def test_predict_quantiles_only(regression_setup):
  """Test prediction with quantiles only."""
  X_train, X_test, y_train, y_test, _, _ = regression_setup
  regressor = VineRegressor(mean=False, quantiles=[0.1, 0.5, 0.9])
  regressor.fit(X_train, y_train)
  pred_quant = regressor.predict(X_test)

  assert isinstance(pred_quant, np.ndarray)
  assert pred_quant.shape == (len(X_test), 3)

  # Check quantile ordering
  assert np.all(pred_quant[:, 0] <= pred_quant[:, 1])  # Q10 <= Q50
  assert np.all(pred_quant[:, 1] <= pred_quant[:, 2])  # Q50 <= Q90


def test_predict_mean_and_quantiles(regression_setup):
  """Test prediction with both mean and quantiles."""
  X_train, X_test, y_train, y_test, _, _ = regression_setup
  regressor = VineRegressor(mean=True, quantiles=[0.1, 0.9])
  regressor.fit(X_train, y_train)
  pred_both = regressor.predict(X_test)

  assert pred_both.shape == (len(X_test), 3)  # mean + 2 quantiles


def test_prediction_accuracy(fitted_regressor):
  """Test prediction accuracy against true mean."""
  regressor, _, X_test, _, y_test, true_mean = fitted_regressor

  pred_mean = regressor.predict(X_test)
  true_mean_test = true_mean(X_test)

  # Check correlation between predicted and true means
  correlation = np.corrcoef(pred_mean, true_mean_test)[0, 1]
  assert correlation > 0.8, "Predicted means should correlate with true means"

  # Check RMSE is reasonable
  rmse = np.sqrt(np.mean((pred_mean - true_mean_test) ** 2))
  assert rmse < 0.5, "RMSE should be reasonable"


def test_quantile_predictions(regression_setup):
  """Test quantile prediction properties."""
  X_train, X_test, y_train, y_test, _, _ = regression_setup

  regressor_quant = VineRegressor(mean=False, quantiles=[0.1, 0.5, 0.9])
  regressor_quant.fit(X_train, y_train)
  pred_quant = regressor_quant.predict(X_test)

  # Median should be close to mean prediction
  regressor_mean = VineRegressor(mean=True)
  regressor_mean.fit(X_train, y_train)
  pred_mean = regressor_mean.predict(X_test)

  median_pred = pred_quant[:, 1]  # 0.5 quantile
  correlation = np.corrcoef(median_pred, pred_mean)[0, 1]
  assert correlation > 0.9, "Median should be close to mean"


def test_wrong_dimensions(fitted_regressor):
  """Test error handling for wrong dimensions."""
  regressor, _, X_test, _, y_test, _ = fitted_regressor

  # Wrong number of features
  X_wrong = X_test[:, :1]  # Only 1 feature instead of 2
  with pytest.raises(ValueError):
    regressor.predict(X_wrong)


def test_invalid_configuration(regression_setup):
  """Test invalid regressor configurations surface at fit time."""
  X_train, _, y_train, _, _, _ = regression_setup
  # Neither mean nor quantiles enabled
  est = VineRegressor(mean=False, quantiles=None)
  with pytest.raises(ValueError):
    est.fit(X_train, y_train)


@pytest.mark.parametrize(
  "kwargs, exc",
  [
    ({"mean": "yes"}, TypeError),
    ({"use_grid": 1}, TypeError),
    ({"batch_size": 0}, ValueError),
    ({"batch_size": -3}, ValueError),
    ({"batch_size": 1.5}, TypeError),
    ({"quantiles": [0.0, 0.5]}, ValueError),
    ({"quantiles": [0.5, 1.0]}, ValueError),
    ({"quantiles": []}, ValueError),
  ],
)
def test_constructor_validation(kwargs, exc, regression_setup):
  """Bad parameters surface at ``fit`` time, per the sklearn dev guide."""
  X_train, _, y_train, _, _, _ = regression_setup
  est = VineRegressor(**kwargs)
  with pytest.raises((exc, ValueError, TypeError)):
    est.fit(X_train, y_train)


@pytest.mark.parametrize("use_grid", [True, False])
def test_parameter_variations(regression_setup, use_grid):
  """Test different parameter combinations."""
  X_train, X_test, y_train, _, _, _ = regression_setup

  regressor = VineRegressor(mean=True, use_grid=use_grid)
  regressor.fit(X_train, y_train)
  pred = regressor.predict(X_test)

  assert isinstance(pred, np.ndarray)
  assert pred.shape == (len(X_test),)
  assert np.all(np.isfinite(pred))


def test_normalize_weights_parameter(regression_setup):
  """``normalize_weights`` parameter toggles row-wise weight normalization."""
  X_train, X_test, y_train, _, _, _ = regression_setup

  reg_default = VineRegressor(mean=True).fit(X_train, y_train)
  pred_default = reg_default.predict(X_test)

  reg_raw = VineRegressor(mean=True, normalize_weights=False).fit(
    X_train, y_train
  )
  pred_raw = reg_raw.predict(X_test)

  # With normalized weights, rows sum to 1 and the prediction is a
  # convex combination of training responses; without normalization it
  # picks up the absolute scale of the copula density, so outputs
  # generally differ.
  assert not np.allclose(pred_default, pred_raw)


def test_copula_marginal_density_single_covariate(regression_setup):
  """``_copula_marginal_density`` recovers :math:`c_X \\equiv 1` in 2-d.

  Integrating a bivariate copula density over one of its arguments is
  exactly one, so a fit with a single covariate pins the Simpson
  quadrature, the batching loop, the even-``n_grid`` fix-up, and the
  ``log=True`` branch in one shot. This method has no in-tree caller,
  so without this test it would be uncovered.
  """
  X_train, X_test, y_train, _, _, _ = regression_setup
  reg = VineRegressor(mean=True, batch_size=7).fit(X_train[:, :1], y_train)

  c_x = reg._copula_marginal_density(X_test[:20, :1], n_grid=200)
  assert c_x.shape == (20,)
  assert np.all(c_x > 0)
  assert np.all(np.isfinite(c_x))
  # The identity is exact; the deviation is quadrature plus TLL
  # boundary error, which is a few 1e-3 in the bulk and up to ~0.15 for
  # a point near the edge of the support.
  assert abs(np.median(c_x) - 1.0) < 1e-2
  assert np.allclose(c_x, 1.0, atol=0.25)

  log_c = reg._copula_marginal_density(X_test[:20, :1], n_grid=200, log=True)
  assert np.allclose(log_c, np.log(c_x))


def test_vine_regressor_is_a_regressor_to_sklearn() -> None:
  # `RegressorMixin` has to precede `VineBase` in the MRO, or
  # `BaseEstimator.__sklearn_tags__` -- which does not call `super()` -- wins
  # and the estimator type is never set. Meta-estimators check this tag before
  # they will accept an estimator at all.
  from sklearn.base import is_regressor
  from sklearn.ensemble import StackingRegressor

  assert is_regressor(VineRegressor())

  rng = np.random.default_rng(0)
  X = rng.normal(size=(80, 3))
  y = X @ np.array([1.0, 2.0, -1.0]) + 0.1 * rng.normal(size=80)
  stacked = StackingRegressor(estimators=[("v", VineRegressor())], cv=2)
  assert stacked.fit(X, y).predict(X[:5]).shape == (5,)


def test_vine_regressor_predict_keeps_the_sample_axis() -> None:
  # Only the output axis collapses for a single output. A one-row X used to
  # come back as a 0-d scalar, so any caller predicting row by row got a
  # different rank than the same call on two rows.
  rng = np.random.default_rng(1)
  X = rng.normal(size=(80, 3))
  y = X @ np.array([1.0, 2.0, -1.0]) + 0.1 * rng.normal(size=80)

  mean_only = VineRegressor().fit(X, y)
  assert mean_only.predict(X[:1]).shape == (1,)
  assert mean_only.predict(X[:4]).shape == (4,)

  multi = VineRegressor(quantiles=[0.25, 0.75]).fit(X, y)
  assert multi.predict(X[:1]).shape == (1, 3)
  assert multi.predict(X[:4]).shape == (4, 3)
