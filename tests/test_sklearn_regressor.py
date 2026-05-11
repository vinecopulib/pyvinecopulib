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


def test_pdf_shapes_and_types(fitted_regressor):
  """Test pdf method shapes and types."""
  regressor, _, X_test, _, y_test, _ = fitted_regressor

  pdf_vals = regressor.pdf(X_test, y_test)
  log_pdf_vals = regressor.pdf(X_test, y_test, log=True)

  # Check shapes and types
  assert isinstance(pdf_vals, np.ndarray)
  assert isinstance(log_pdf_vals, np.ndarray)
  assert pdf_vals.shape == (len(X_test),)
  assert log_pdf_vals.shape == (len(X_test),)

  # Check positivity and consistency
  assert np.all(pdf_vals > 0)
  np.testing.assert_allclose(pdf_vals, np.exp(log_pdf_vals), rtol=1e-10)


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

  # Wrong y dimension for pdf
  y_wrong = np.random.randn(len(X_test), 2)  # Should be 1D
  with pytest.raises(ValueError):
    regressor.pdf(X_test, y_wrong)


def test_invalid_configuration():
  """Test invalid regressor configurations."""
  # Neither mean nor quantiles enabled
  with pytest.raises(ValueError):
    VineRegressor(mean=False, quantiles=None)


@pytest.mark.parametrize("use_grid", [True, False])
@pytest.mark.parametrize("normalize_weights", [True, False])
def test_parameter_variations(regression_setup, use_grid, normalize_weights):
  """Test different parameter combinations."""
  X_train, X_test, y_train, y_test, _, _ = regression_setup

  regressor = VineRegressor(
    mean=True, use_grid=use_grid, normalize_weights=normalize_weights
  )
  regressor.fit(X_train, y_train)
  pred = regressor.predict(X_test)

  assert isinstance(pred, np.ndarray)
  assert pred.shape == (len(X_test),)
  assert np.all(np.isfinite(pred))


def test_copula_only_pdf(fitted_regressor):
  """Test that regressor PDF uses copula-only mode."""
  regressor, _, X_test, _, y_test, _ = fitted_regressor

  # Get PDF from regressor (should be copula-only)
  copula_pdf = regressor.pdf(X_test, y_test)

  # Get copula-only from base method
  copula_direct = regressor._pdf_samples(
    X_test, y=y_test, log=False, copula_only=True
  )

  # Should be the same
  np.testing.assert_allclose(copula_pdf, copula_direct, rtol=1e-12)
