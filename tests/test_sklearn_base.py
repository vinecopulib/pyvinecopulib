"""
Tests for VineBase functionality.
"""

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("pandas")

import pandas as pd  # noqa: E402

from sklearn.exceptions import DataConversionWarning  # noqa: E402

from pyvinecopulib.sklearn import VineDensity, VineRegressor  # noqa: E402
from pyvinecopulib.sklearn._base import expand_factors  # noqa: E402


def test_expand_factors_numeric_unchanged() -> None:
  """Test that numeric columns remain unchanged."""
  df = pd.DataFrame({"x1": [1.0, 2.0, 3.0], "x2": [4, 5, 6]})
  result: pd.DataFrame = expand_factors(df)

  pd.testing.assert_frame_equal(result, df)


def test_expand_factors_ordered_categorical_unchanged() -> None:
  """Test that ordered categoricals remain unchanged."""
  df = pd.DataFrame(
    {"cat": pd.Categorical(["low", "med", "high"], ordered=True)}
  )
  result: pd.DataFrame = expand_factors(df)
  pd.testing.assert_frame_equal(result, df)


def test_expand_factors_unordered_categorical_expanded() -> None:
  """Test that unordered categoricals are expanded into dummies."""
  df = pd.DataFrame({"color": pd.Categorical(["red", "blue", "red", "green"])})
  result = expand_factors(df)

  # Should have 2 columns (drop first level 'blue')
  assert result.shape == (4, 2)
  assert list(result.columns) == ["color_green", "color_red"]

  # Check values
  assert result["color_red"].tolist() == [1, 0, 1, 0]
  assert result["color_green"].tolist() == [0, 0, 0, 1]

  # Check that dummies are ordered categorical
  assert isinstance(result["color_red"].dtype, pd.CategoricalDtype)
  assert result["color_red"].dtype.ordered


def test_expand_factors_mixed_columns() -> None:
  """Test DataFrame with mixed column types."""
  df = pd.DataFrame(
    {
      "numeric": [1.0, 2.0, 3.0],
      "ordered": pd.Categorical(["a", "b", "c"], ordered=True),
      "unordered": pd.Categorical(["x", "y", "x"]),
    }
  )
  result = expand_factors(df)

  # Should have 3 columns: numeric, ordered, unordered_y
  assert result.shape == (3, 3)
  assert list(result.columns) == ["numeric", "ordered", "unordered_y"]

  # Check numeric unchanged
  assert result["numeric"].tolist() == [1.0, 2.0, 3.0]

  # Check dummy
  assert result["unordered_y"].tolist() == [0, 1, 0]


def test_vinebase_array_input_validation(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """Test array input validation and processing."""
  X, _, _ = sample_array_data
  density: VineDensity = VineDensity()

  # Test successful array processing
  X_processed = density._validate_input(X, reset=True)
  assert isinstance(X_processed, np.ndarray)
  assert X_processed.shape == X.shape
  assert density.n_features_in_ == 2
  assert density.schema_ is not None
  assert len(density.schema_["kde1d_types"]) == 2
  assert all(t == "continuous" for t in density.schema_["kde1d_types"])


def test_dataframe_prediction_does_not_mutate_caller_categories() -> None:
  """Schema validation recategorizes a private copy, never the caller's frame."""
  train = pd.DataFrame(
    {"x": [0.0, 1.0, 2.0, 3.0], "c": pd.Categorical(["a", "b", "a", "b"])}
  )
  est = VineDensity().fit(train)
  query = pd.DataFrame({"x": [0.5, 1.5], "c": pd.Categorical(["b", "a"])})
  before = query.copy(deep=True)
  est.score_samples(query)
  pd.testing.assert_frame_equal(query, before)
  assert est.n_features_in_ == 2
  assert est.n_model_features_ == 2


def test_vinebase_dataframe_expansion(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """Test DataFrame expansion and schema creation."""
  X_df, expected_expanded_cols = sample_dataframe_data
  density: VineDensity = VineDensity()

  # Test DataFrame processing via the canonical sklearn-style helper
  X_processed = density._validate_input(X_df, reset=True)
  assert isinstance(X_processed, np.ndarray)
  assert X_processed.shape[0] == len(X_df)
  assert X_processed.shape[1] == len(expected_expanded_cols)

  # Check that original columns are stored under the canonical sklearn name
  assert list(density.feature_names_in_) == list(X_df.columns)
  assert density._expanded_columns == expected_expanded_cols

  # Check schema creation
  expected_types = [
    "continuous",
    "continuous",
    "discrete",
    "discrete",
    "continuous",
  ]
  assert density.schema_ is not None
  assert density.schema_["kde1d_types"] == expected_types


def test_vinebase_dataframe_prediction_validation(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """Test DataFrame validation during prediction."""
  X_df, expected_expanded_cols = sample_dataframe_data
  density = VineDensity()
  density.fit(X_df)

  # Test successful prediction with matching DataFrame
  X_test = X_df.iloc[:10].copy()
  X_processed = density._validate_input(X_test, reset=False)
  assert isinstance(X_processed, np.ndarray)
  assert X_processed.shape == (10, len(expected_expanded_cols))

  # Test error with wrong columns
  X_wrong = X_test.drop("cont1", axis=1)
  with pytest.raises(ValueError, match="Column names/order do not match"):
    density._validate_input(X_wrong, reset=False)


def test_vinebase_schema_attribute(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """Pre-set ``schema_`` overrides the auto-inferred schema."""
  X, _, _ = sample_array_data
  schema = {"kde1d_types": ["continuous", "discrete"]}
  density = VineDensity()
  density.schema_ = schema

  density._validate_input(X, reset=True)
  assert density.schema_ is not None
  assert density.schema_["kde1d_types"] == ["continuous", "discrete"]

  # Schema length mismatch raises.
  density_wrong = VineDensity()
  density_wrong.schema_ = {"kde1d_types": ["continuous"]}  # Too short
  with pytest.raises(ValueError):
    density_wrong._validate_input(X, reset=True)

  # So does a pre-set bounds list of the wrong length.
  density_bounds = VineDensity()
  density_bounds.schema_ = {"bounds": [(0.0, 1.0)]}  # Too short
  with pytest.raises(ValueError, match="bounds"):
    density_bounds._validate_input(X, reset=True)


def test_vinebase_marginal_fitting(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """Test marginal distribution fitting."""
  X, _, _ = sample_array_data
  density = VineDensity()
  X_processed = density._validate_input(X, reset=True)
  assert isinstance(X_processed, np.ndarray)
  density._fit_marginals(X_processed)

  # Check that marginals are fitted
  assert len(density._x_margins) == 2
  assert all(margin.is_fitted for margin in density._x_margins)


def test_vinebase_pseudoobservations(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """Test pseudo-observation transformation."""
  X, _, _ = sample_array_data
  density = VineDensity()
  X_processed = density._validate_input(X, reset=True)
  assert isinstance(X_processed, np.ndarray)
  density._fit_marginals(X_processed)

  U = density._to_u_scale(X_processed)

  # Check shape and range
  assert U.shape == X_processed.shape
  assert np.all(U >= 0)
  assert np.all(U <= 1)

  # Check that values are approximately uniform
  for j in range(U.shape[1]):
    # Kolmogorov-Smirnov test would be ideal, but just check basic distribution
    assert 0.1 < np.mean(U[:, j]) < 0.9


def test_vinebase_pdf_samples_unified_method(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
  regression_data: tuple[np.ndarray, np.ndarray, np.ndarray, float],
) -> None:
  """Test the unified _pdf_samples method."""
  X, _, _ = sample_array_data
  X_reg, y_reg, _, _ = regression_data

  # Test density estimation case
  density = VineDensity()
  density.fit(X)

  X_test = X[:20]

  # Test different modes
  log_full = density._pdf_samples(X_test, log=True, copula_only=False)
  log_copula = density._pdf_samples(X_test, log=True, copula_only=True)
  pdf_full = density._pdf_samples(X_test, log=False, copula_only=False)
  pdf_copula = density._pdf_samples(X_test, log=False, copula_only=True)

  # Check shapes and types
  assert isinstance(log_full, np.ndarray)
  assert log_full.shape == (20,)
  assert np.allclose(pdf_full, np.exp(log_full))
  assert np.allclose(pdf_copula, np.exp(log_copula))

  # Test joint density case (using VineRegressor-like setup)
  from pyvinecopulib.sklearn import VineRegressor

  regressor = VineRegressor()
  regressor.fit(X_reg[:200], y_reg[:200])

  X_test_reg = X_reg[200:220]
  y_test_reg = y_reg[200:220]

  # Test joint density
  log_joint = regressor._pdf_samples(
    X_test_reg, y=y_test_reg, log=True, copula_only=False
  )
  copula_joint = regressor._pdf_samples(
    X_test_reg, y=y_test_reg, log=True, copula_only=True
  )

  assert isinstance(log_joint, np.ndarray)
  assert log_joint.shape == (20,)

  assert isinstance(copula_joint, np.ndarray)
  assert copula_joint.shape == (20,)


def test_non_finite_input_raises_rather_than_crashing() -> None:
  # The marginal estimator is C++ and reads NaN as a segmentation fault, so
  # the guard has to be on this side of the boundary.
  rng = np.random.default_rng(0)
  X = rng.normal(size=(60, 3))
  y = X @ np.array([1.0, 2.0, -1.0])

  bad_X = X.copy()
  bad_X[0, 0] = np.nan
  with pytest.raises(ValueError, match="NaN"):
    VineDensity().fit(bad_X)

  bad_y = y.copy()
  bad_y[3] = np.inf
  with pytest.raises(ValueError, match="infinity|inf"):
    VineRegressor().fit(X, bad_y)


def test_column_vector_y_is_raveled_with_a_warning() -> None:
  rng = np.random.default_rng(1)
  X = rng.normal(size=(60, 3))
  y = (X @ np.array([1.0, 2.0, -1.0])).reshape(-1, 1)

  with pytest.warns(DataConversionWarning):
    model = VineRegressor().fit(X, y)
  assert model.predict(X[:4]).shape == (4,)
