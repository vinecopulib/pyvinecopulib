"""
Tests for VineBase functionality.
"""

from typing import Any

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("pandas")

import pandas as pd
from sklearn.exceptions import DataConversionWarning

from pyvinecopulib.sklearn import VineDensity, VineRegressor
from pyvinecopulib.sklearn._base import expand_factors


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
  assert len(density.schema_["var_types"]) == 2
  assert all(t == "c" for t in density.schema_["var_types"])


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
  expected_types = ["c", "c", "d", "d", "c"]
  assert density.schema_ is not None
  assert density.schema_["var_types"] == expected_types


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
  schema = {"var_types": ["c", "d"]}
  density = VineDensity()
  density.schema_ = schema

  density._validate_input(X, reset=True)
  assert density.schema_ is not None
  assert density.schema_["var_types"] == ["c", "d"]

  # Schema length mismatch raises.
  density_wrong = VineDensity()
  density_wrong.schema_ = {"var_types": ["c"]}  # Too short
  with pytest.raises(ValueError):
    density_wrong._validate_input(X, reset=True)

  # So does a pre-set bounds list of the wrong length.
  density_bounds = VineDensity()
  density_bounds.schema_ = {"supports": [(0.0, 1.0)]}  # Too short
  with pytest.raises(ValueError, match="supports"):
    density_bounds._validate_input(X, reset=True)


def test_vinebase_marginal_fitting(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """Test marginal distribution fitting."""
  X, _, _ = sample_array_data
  density = VineDensity()
  X_processed = density._validate_input(X, reset=True)
  assert isinstance(X_processed, np.ndarray)
  density._resolve_runtime_state()
  density._fit_distribution(X_processed)

  # Check that marginals are fitted
  assert len(density._x_margins) == 2
  # `is_fitted` is an optional capability on `MarginLike`, so it is read the
  # way the library reads it.
  assert all(
    getattr(margin, "is_fitted", False) for margin in density._x_margins
  )


def test_vinebase_pseudoobservations(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """Test pseudo-observation transformation."""
  X, _, _ = sample_array_data
  density = VineDensity()
  X_processed = density._validate_input(X, reset=True)
  assert isinstance(X_processed, np.ndarray)
  density._resolve_runtime_state()
  density._fit_distribution(X_processed)

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
  with pytest.raises(ValueError, match=r"infinity|inf"):
    VineRegressor().fit(X, bad_y)


def test_column_vector_y_is_raveled_with_a_warning() -> None:
  rng = np.random.default_rng(1)
  X = rng.normal(size=(60, 3))
  y = (X @ np.array([1.0, 2.0, -1.0])).reshape(-1, 1)

  with pytest.warns(DataConversionWarning):
    model = VineRegressor().fit(X, y)
  assert model.predict(X[:4]).shape == (4,)


def test_fit_resets_the_schema_a_previous_fit_derived(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """A second `fit` must not validate against the first fit's schema.

  `_validate_input(reset=True)` read a leftover `schema_`, so refitting the
  same estimator on an array raised instead of refitting -- even at the same
  width -- while a clone of it fitted fine.
  """
  X_df, _ = sample_dataframe_data
  est = VineDensity().fit(X_df)
  rng = np.random.RandomState(0)
  est.fit(rng.normal(size=(80, 3)))
  assert est.n_features_in_ == 3
  assert est.n_model_features_ == 3
  assert not hasattr(est, "feature_names_in_")


def test_an_array_is_read_as_the_modeled_layout() -> None:
  """After an expanding fit, an array is the wide layout and only that.

  `sample` emits the modeled width, so the estimator's own output has to be a
  legal input to its own density. The public width has no reading as an array:
  which columns are the levels of which factor is what the frame carries.
  """
  rng = np.random.RandomState(0)
  train = pd.DataFrame(
    {
      "x": rng.normal(size=60),
      "c": pd.Categorical(rng.choice(["a", "b", "c"], size=60)),
    }
  )
  est = VineDensity().fit(train)
  assert (est.n_features_in_, est.n_model_features_) == (2, 3)

  assert est.score_samples(est.sample(5)).shape == (5,)
  with pytest.raises(ValueError, match="expecting 3 features"):
    est.score_samples(rng.normal(size=(5, 2)))


def test_a_caller_preset_schema_survives_fit() -> None:
  """The pre-settable `schema_` hook still overrides the array default."""
  est = VineDensity()
  est.schema_ = {
    "var_types": ["d", "c"],
    "supports": [None] * 2,
  }
  rng = np.random.RandomState(0)
  est.fit(
    np.column_stack([rng.poisson(3, 80).astype(float), rng.normal(size=80)])
  )
  assert est.schema_["var_types"] == ["d", "c"]


def test_a_preset_schema_and_a_dataframe_is_refused_not_discarded() -> None:
  """A frame states its own types, so the two declarations are a conflict.

  The frame's won silently, which drops exactly what a caller pre-sets
  `schema_` to say: a `"zi"` column, or a bound, neither of which a pandas
  dtype can express. The column was then modeled as continuous.
  """
  rng = np.random.RandomState(0)
  frame = pd.DataFrame(
    {
      "zeros": np.where(rng.random(80) < 0.4, 0.0, rng.gamma(2.0, size=80)),
      "plain": rng.normal(size=80),
    }
  )
  est = VineDensity()
  est.schema_ = {"var_types": ["zi", "c"], "supports": [None] * 2}
  with pytest.raises(ValueError, match="states its own variable types"):
    est.fit(frame)
  # And the route the message names does honor it.
  clean = VineDensity()
  clean.schema_ = {"var_types": ["zi", "c"], "supports": [None] * 2}
  clean.fit(frame.to_numpy())
  assert clean.schema_["var_types"] == ["zi", "c"]


def test_dataframe_after_an_array_fit_is_not_reported_unfitted() -> None:
  """A fitted estimator must never raise `NotFittedError`."""
  rng = np.random.RandomState(0)
  est = VineDensity().fit(rng.normal(size=(80, 2)))
  frame = pd.DataFrame(rng.normal(size=(4, 2)), columns=["a", "b"])
  assert est.score_samples(frame).shape == (4,)
  # A categorical column was never modeled as one, so that is refused by name
  # rather than read as codes.
  frame["b"] = pd.Categorical(["x", "y", "x", "y"])
  with pytest.raises(ValueError, match="fitted without feature names"):
    est.score_samples(frame)


def test_an_unseen_category_is_refused_not_silently_recoded(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """`set_categories` maps an unseen level to NaN, which expands to all zeros.

  That row is indistinguishable from the reference level, so an unseen
  category used to return the reference level's density with no warning.
  """
  X_df, _ = sample_dataframe_data
  est = VineDensity().fit(X_df)
  cat_col = next(
    c for c in X_df.columns if isinstance(X_df[c].dtype, pd.CategoricalDtype)
  )
  bad = X_df.head(3).copy()
  levels = [*list(X_df[cat_col].cat.categories), "!unseen!"]
  bad[cat_col] = pd.Categorical(["!unseen!"] * 3, categories=levels)
  with pytest.raises(ValueError, match="not seen during fit"):
    est.score_samples(bad)


def test_array_like_input_is_accepted() -> None:
  """sklearn's convention is that any array-like is valid input."""
  rng = np.random.RandomState(0)
  rows = rng.normal(size=(60, 2)).tolist()
  est = VineDensity().fit(rows)
  assert est.n_features_in_ == 2
  assert est.score_samples(rows[:3]).shape == (3,)
  with pytest.raises(ValueError, match="array-like of floats"):
    VineDensity().fit([["a", "b"], ["c", "d"]])


def test_the_default_controls_are_resolvable_before_fitting() -> None:
  """Which controls an estimator would fit with depends only on `__init__`.

  It is answered from `distribution`, so it holds on an estimator nothing has
  fitted; reading the fitted `distribution_class_` made it answerable only
  *after* a fit, which is the wrong way round for a question about defaults --
  and raised `AttributeError` for anyone who asked earlier.
  """
  from sklearn.utils.validation import check_is_fitted

  import pyvinecopulib as pv

  for est in (VineDensity(), VineRegressor(quantiles=[0.5])):
    controls: Any = est._default_copula_controls()
    assert controls is not None
    assert list(controls.family_set) == [pv.families.tll]
    assert controls.trunc_lvl == 20
    # Asking must not leave the estimator looking fitted.
    with pytest.raises(Exception, match="not fitted"):
      check_is_fitted(est)


def test_the_default_controls_come_from_the_vine_class_declaration() -> None:
  """The default is built from `vinecop_class.controls_class`, not a name here.

  A lane whose controls have no family to choose gets `None` instead, and the
  vine resolves its own default -- which is what retired the `is pv.Vinecop`
  identity check this used to dispatch on.
  """
  torch = pytest.importorskip("torch")
  del torch
  import pyvinecopulib as pv
  from pyvinecopulib.core import Vinedist
  from pyvinecopulib.torch import TorchVinedist

  vinecop_class: Any = Vinedist.vinecop_class
  assert vinecop_class.controls_class is pv.FitControlsVinecop
  assert VineDensity()._default_copula_controls() is not None
  # The torch vine fits TLL grids and nothing else, so there is no family set
  # to narrow and no density default to express.
  assert (
    VineDensity(distribution=TorchVinedist)._default_copula_controls() is None
  )
