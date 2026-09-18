"""Tests for VineRegressor class."""

from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("pandas")

import pandas as pd
from sklearn.exceptions import NotFittedError

from pyvinecopulib.sklearn import VineRegressor

# The two shared fixtures return wide tuples; naming them keeps every
# signature below readable.
_Setup = tuple[
  np.ndarray,
  np.ndarray,
  np.ndarray,
  np.ndarray,
  Callable[[np.ndarray], np.ndarray],
  np.ndarray,
]
_Fitted = tuple[
  VineRegressor,
  np.ndarray,
  np.ndarray,
  np.ndarray,
  np.ndarray,
  Callable[[np.ndarray], np.ndarray],
]
_RegressionData = tuple[np.ndarray, np.ndarray, np.ndarray, float]


@pytest.fixture
def regression_setup(
  regression_data: _RegressionData,
) -> _Setup:
  """Setup regression data and split."""
  X, y, true_coef, _noise_std = regression_data
  X_train, X_test = X[:200], X[200:]
  y_train, y_test = y[:200], y[200:]

  def true_mean(x: np.ndarray) -> np.ndarray:
    return x @ true_coef

  return X_train, X_test, y_train, y_test, true_mean, true_coef


@pytest.fixture
def fitted_regressor(
  regression_setup: _Setup,
) -> _Fitted:
  """Fitted VineRegressor for testing."""
  X_train, X_test, y_train, y_test, true_mean, _true_coef = regression_setup
  regressor = VineRegressor(
    mean=True
  )  # Only mean, no quantiles for simpler testing
  regressor.fit(X_train, y_train)
  return regressor, X_train, X_test, y_train, y_test, true_mean


def test_fit_properties(
  fitted_regressor: _Fitted,
) -> None:
  """Test properties after fitting."""
  regressor, _, _, _, _, _ = fitted_regressor
  assert hasattr(regressor, "_vine")
  assert regressor.n_features_in_ == 2
  # The joint distribution is over (Y, X), the response first.
  assert regressor.distribution_.dim == 3


def test_predict_mean_only(
  regression_setup: _Setup,
) -> None:
  """Test prediction with mean only."""
  X_train, X_test, y_train, _y_test, _, _ = regression_setup
  regressor = VineRegressor(mean=True)
  regressor.fit(X_train, y_train)
  pred_mean = regressor.predict(X_test)

  assert isinstance(pred_mean, np.ndarray)
  assert pred_mean.shape == (len(X_test),)
  assert np.all(np.isfinite(pred_mean))


def test_predict_quantiles_only(
  regression_setup: _Setup,
) -> None:
  """Test prediction with quantiles only."""
  X_train, X_test, y_train, _y_test, _, _ = regression_setup
  regressor = VineRegressor(mean=False, quantiles=[0.1, 0.5, 0.9])
  regressor.fit(X_train, y_train)
  pred_quant = regressor.predict(X_test)

  assert isinstance(pred_quant, np.ndarray)
  assert pred_quant.shape == (len(X_test), 3)

  # Check quantile ordering
  assert np.all(pred_quant[:, 0] <= pred_quant[:, 1])  # Q10 <= Q50
  assert np.all(pred_quant[:, 1] <= pred_quant[:, 2])  # Q50 <= Q90


def test_quantile_regressor_score_is_rejected(
  regression_setup: _Setup,
) -> None:
  """R² is defined for the mean prediction, not a quantile matrix."""
  X_train, X_test, y_train, y_test, _, _ = regression_setup
  est = VineRegressor(mean=False, quantiles=[0.25, 0.75]).fit(X_train, y_train)
  with pytest.raises(ValueError, match="mean predictions"):
    est.score(X_test, y_test)


def test_predict_mean_and_quantiles(
  regression_setup: _Setup,
) -> None:
  """Test prediction with both mean and quantiles."""
  X_train, X_test, y_train, _y_test, _, _ = regression_setup
  regressor = VineRegressor(mean=True, quantiles=[0.1, 0.9])
  regressor.fit(X_train, y_train)
  pred_both = regressor.predict(X_test)

  assert pred_both.shape == (len(X_test), 3)  # mean + 2 quantiles


def test_prediction_accuracy(
  fitted_regressor: _Fitted,
) -> None:
  """Test prediction accuracy against true mean."""
  regressor, _, X_test, _, _y_test, true_mean = fitted_regressor

  pred_mean = regressor.predict(X_test)
  true_mean_test = true_mean(X_test)

  # Check correlation between predicted and true means
  correlation = np.corrcoef(pred_mean, true_mean_test)[0, 1]
  assert correlation > 0.8, "Predicted means should correlate with true means"

  # Check RMSE is reasonable
  rmse = np.sqrt(np.mean((pred_mean - true_mean_test) ** 2))
  assert rmse < 0.5, "RMSE should be reasonable"


def test_quantile_predictions(
  regression_setup: _Setup,
) -> None:
  """Test quantile prediction properties."""
  X_train, X_test, y_train, _y_test, _, _ = regression_setup

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


def test_wrong_dimensions(
  fitted_regressor: _Fitted,
) -> None:
  """Test error handling for wrong dimensions."""
  regressor, _, X_test, _, _y_test, _ = fitted_regressor

  # Wrong number of features
  X_wrong = X_test[:, :1]  # Only 1 feature instead of 2
  with pytest.raises(ValueError):
    regressor.predict(X_wrong)


def test_invalid_configuration(
  regression_setup: _Setup,
) -> None:
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
    ({"n_nodes": 1}, ValueError),
    ({"n_nodes": 12.5}, TypeError),
    ({"batch_size": 0}, ValueError),
    ({"batch_size": -3}, ValueError),
    ({"batch_size": 1.5}, TypeError),
    ({"quantiles": [0.0, 0.5]}, ValueError),
    ({"quantiles": [0.5, 1.0]}, ValueError),
    ({"quantiles": []}, ValueError),
  ],
)
def test_constructor_validation(
  kwargs: dict[str, Any],
  exc: type[Exception],
  regression_setup: _Setup,
) -> None:
  """Bad parameters surface at ``fit`` time, per the sklearn dev guide."""
  X_train, _, y_train, _, _, _ = regression_setup
  est = VineRegressor(**kwargs)
  with pytest.raises((exc, ValueError, TypeError)):
    est.fit(X_train, y_train)


@pytest.mark.parametrize("use_grid", [True, False])
def test_parameter_variations(
  regression_setup: _Setup,
  use_grid: bool,
) -> None:
  """Test different parameter combinations."""
  X_train, X_test, y_train, _, _, _ = regression_setup

  regressor = VineRegressor(mean=True, use_grid=use_grid)
  regressor.fit(X_train, y_train)
  pred = regressor.predict(X_test)

  assert isinstance(pred, np.ndarray)
  assert pred.shape == (len(X_test),)
  assert np.all(np.isfinite(pred))


def test_normalize_weights_parameter(
  regression_setup: _Setup,
) -> None:
  """``normalize_weights`` parameter toggles row-wise weight normalization."""
  X_train, X_test, y_train, _, _, _ = regression_setup

  reg_default = VineRegressor(mean=True).fit(X_train, y_train)
  pred_default = reg_default.predict(X_test)

  reg_raw = VineRegressor(mean=True, normalize_weights=False).fit(
    X_train, y_train
  )
  pred_raw = reg_raw.predict(X_test)

  # The conditional mean is `sum(w y) / sum(w)`, so it does not depend on the
  # weights' scale: the flag changes what `conditional_weights` returns -- what
  # a caller combining several vines normalizes across -- not the prediction.
  # This test used to assert the opposite, pinning a mean that was scaled by
  # the weight total.
  np.testing.assert_allclose(pred_default, pred_raw, rtol=1e-10, atol=1e-10)
  w_default = reg_default.conditional_weights(X_test[:3])
  w_raw = reg_raw.conditional_weights(X_test[:3])
  np.testing.assert_allclose(w_default.sum(axis=1), 1.0)
  assert not np.allclose(w_default, w_raw)


def testcopula_marginal_density_single_covariate(
  regression_setup: _Setup,
) -> None:
  r"""``copula_marginal_density`` recovers :math:`c_X \\equiv 1` in 2-d.

  Integrating a bivariate copula density over one of its arguments is
  exactly one, so a fit with a single covariate pins the Simpson
  quadrature, the batching loop, the even-``n_grid`` fix-up, and the
  ``log=True`` branch in one shot. This method has no in-tree caller,
  so without this test it would be uncovered.
  """
  X_train, X_test, y_train, _, _, _ = regression_setup
  reg = VineRegressor(mean=True, batch_size=7).fit(X_train[:, :1], y_train)

  c_x = reg.copula_marginal_density(X_test[:20, :1], n_grid=200)
  assert c_x.shape == (20,)
  assert np.all(c_x > 0)
  assert np.all(np.isfinite(c_x))
  # The identity is exact; the deviation is quadrature plus TLL
  # boundary error, which is a few 1e-3 in the bulk and up to ~0.15 for
  # a point near the edge of the support.
  assert abs(np.median(c_x) - 1.0) < 1e-2
  assert np.allclose(c_x, 1.0, atol=0.25)

  log_c = reg.copula_marginal_density(X_test[:20, :1], n_grid=200, log=True)
  assert np.allclose(log_c, np.log(c_x))


def test_regressor_dataframe_inputs_are_validated(
  regression_setup: _Setup,
) -> None:
  """The regressor accepts DataFrames at fit and marginal-density time."""
  X_train, X_test, y_train, _, _, _ = regression_setup
  train = pd.DataFrame({"x": X_train[:, 0]})
  query = pd.DataFrame({"x": X_test[:20, 0]})
  reg = VineRegressor(mean=True, batch_size=7).fit(train, y_train.tolist())

  density = reg.copula_marginal_density(query, n_grid=200)
  assert density.shape == (len(query),)
  assert np.all(np.isfinite(density))

  with pytest.raises(ValueError, match="Column names/order do not match"):
    reg.copula_marginal_density(query.rename(columns={"x": "other"}))


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


@pytest.mark.parametrize("use_grid", [True, False])
def test_normalize_weights_does_not_move_the_conditional_mean(
  use_grid: bool,
  regression_data: _RegressionData,
) -> None:
  """The mean is `sum(w y) / sum(w)`, whatever the weights' scale.

  It used to be a plain dot product, so `normalize_weights=False` returned a
  mean scaled by the weight total -- around 100x with the default grid -- while
  the quantile columns of the same matrix stayed correct. The flag exists so a
  caller combining several vines can normalize once, across all of them.
  """
  X, y, _, _ = regression_data
  preds = [
    VineRegressor(mean=True, use_grid=use_grid, normalize_weights=nw)
    .fit(X, y)
    .predict(X[:40])
    for nw in (True, False)
  ]
  np.testing.assert_allclose(preds[0], preds[1], rtol=1e-10, atol=1e-10)
  assert np.std(preds[0]) < 5.0 * np.std(y)


def test_conditional_weights_are_the_estimator_before_it_is_summarized(
  regression_data: _RegressionData,
) -> None:
  """`predict` is `conditional_weights` against `y_nodes_`, and says so.

  The weights are what a caller combining several fitted vines works with, so
  they are reachable and the response nodes they are attached to are a fitted
  attribute. Batching is internal: a query longer than `batch_size` comes back
  as one matrix, in the order it was asked in.
  """
  X, y, _, _ = regression_data
  est = VineRegressor(mean=True, batch_size=3).fit(X, y)
  X_test = X[:7]

  w = est.conditional_weights(X_test)
  assert w.shape == (7, est.y_nodes_.shape[0])
  np.testing.assert_allclose(
    w @ est.y_nodes_ / w.sum(axis=1),
    est.predict(X_test),
    rtol=1e-12,
    atol=0.0,
  )
  # Row 0 of a seven-row query is row 0 of a one-row one: a batch is written
  # at the offset it was read at, not at a batch counter.
  np.testing.assert_allclose(
    est.conditional_weights(X_test[:1])[0], w[0], rtol=0.0, atol=0.0
  )


def test_conditional_weights_needs_a_fit(
  regression_data: _RegressionData,
) -> None:
  """Every post-fit method checks first, and this is one."""
  X, _, _, _ = regression_data
  with pytest.raises(NotFittedError):
    VineRegressor().conditional_weights(X[:2])
