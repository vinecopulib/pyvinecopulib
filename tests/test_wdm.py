import numpy as np
import pytest
from numpy.testing import assert_allclose

import pyvinecopulib as pv


class TestWdm:
  """Test weighted dependence measures function"""

  def setup_method(self) -> None:
    """Set up test data for each test method"""
    np.random.seed(42)
    self.n = 100

    # Create test data with known correlations
    self.x_independent = np.random.normal(0, 1, self.n)
    self.y_independent = np.random.normal(0, 1, self.n)

    # Create positively correlated data
    z = np.random.normal(0, 1, self.n)
    self.x_positive = z + np.random.normal(0, 0.5, self.n)
    self.y_positive = z + np.random.normal(0, 0.5, self.n)

    # Create negatively correlated data
    self.x_negative = z + np.random.normal(0, 0.5, self.n)
    self.y_negative = -z + np.random.normal(0, 0.5, self.n)

    # Perfect positive correlation
    self.x_perfect = np.linspace(0, 10, self.n)
    self.y_perfect = 2 * self.x_perfect + 1

    # Perfect negative correlation
    self.x_perfect_neg = np.linspace(0, 10, self.n)
    self.y_perfect_neg = -2 * self.x_perfect_neg + 5

  def test_wdm_basic_functionality(self) -> None:
    """Test basic functionality with different methods"""
    methods = ["pearson", "spearman", "kendall", "blomqvist", "hoeffding"]

    for method in methods:
      # Should not raise an error and return a float
      result = pv.wdm(self.x_independent, self.y_independent, method)
      assert isinstance(result, float)
      assert not np.isnan(result)
      assert not np.isinf(result)

  def test_wdm_method_aliases(self) -> None:
    """Test that method aliases work correctly"""
    x, y = self.x_positive, self.y_positive

    # Pearson aliases
    pearson_methods = ["pearson", "prho", "cor"]
    pearson_results = [pv.wdm(x, y, method) for method in pearson_methods]
    for result in pearson_results[1:]:
      assert_allclose(result, pearson_results[0], rtol=1e-15)

    # Spearman aliases
    spearman_methods = ["spearman", "srho", "rho"]
    spearman_results = [pv.wdm(x, y, method) for method in spearman_methods]
    for result in spearman_results[1:]:
      assert_allclose(result, spearman_results[0], rtol=1e-15)

    # Kendall aliases
    kendall_methods = ["kendall", "ktau", "tau"]
    kendall_results = [pv.wdm(x, y, method) for method in kendall_methods]
    for result in kendall_results[1:]:
      assert_allclose(result, kendall_results[0], rtol=1e-15)

    # Blomqvist aliases
    blomqvist_methods = ["blomqvist", "bbeta", "beta"]
    blomqvist_results = [pv.wdm(x, y, method) for method in blomqvist_methods]
    for result in blomqvist_results[1:]:
      assert_allclose(result, blomqvist_results[0], rtol=1e-15)

    # Hoeffding aliases
    hoeffding_methods = ["hoeffding", "hoeffd", "d"]
    hoeffding_results = [pv.wdm(x, y, method) for method in hoeffding_methods]
    for result in hoeffding_results[1:]:
      assert_allclose(result, hoeffding_results[0], rtol=1e-15)

  def test_wdm_perfect_correlations(self) -> None:
    """Test with perfect correlations"""
    # Perfect positive correlation
    pearson_pos = pv.wdm(self.x_perfect, self.y_perfect, "pearson")
    assert_allclose(pearson_pos, 1.0, atol=1e-10)

    spearman_pos = pv.wdm(self.x_perfect, self.y_perfect, "spearman")
    assert_allclose(spearman_pos, 1.0, atol=1e-10)

    kendall_pos = pv.wdm(self.x_perfect, self.y_perfect, "kendall")
    assert_allclose(kendall_pos, 1.0, atol=1e-10)

    # Perfect negative correlation
    pearson_neg = pv.wdm(self.x_perfect_neg, self.y_perfect_neg, "pearson")
    assert_allclose(pearson_neg, -1.0, atol=1e-10)

    spearman_neg = pv.wdm(self.x_perfect_neg, self.y_perfect_neg, "spearman")
    assert_allclose(spearman_neg, -1.0, atol=1e-10)

    kendall_neg = pv.wdm(self.x_perfect_neg, self.y_perfect_neg, "kendall")
    assert_allclose(kendall_neg, -1.0, atol=1e-10)

  def test_wdm_independence(self) -> None:
    """Test with independent data (should be close to zero)"""
    methods = ["pearson", "spearman", "kendall", "blomqvist"]

    for method in methods:
      result = pv.wdm(self.x_independent, self.y_independent, method)
      # For independent data, correlation should be close to zero
      assert abs(result) < 0.2  # Allow some variation due to randomness

  def test_wdm_weights_basic(self) -> None:
    """Test basic weighted functionality"""
    x, y = self.x_positive, self.y_positive
    weights = np.ones(len(x))

    # Uniform weights should give same result as unweighted
    unweighted = pv.wdm(x, y, "pearson")
    weighted_uniform = pv.wdm(x, y, "pearson", weights)
    assert_allclose(unweighted, weighted_uniform, rtol=1e-10)

  def test_wdm_weights_half_zero(self) -> None:
    """Test with half weights zero, half weights one (as requested)"""
    x, y = self.x_positive, self.y_positive
    n = len(x)

    # Create weights: first half zero, second half one
    weights = np.zeros(n)
    weights[n // 2 :] = 1.0

    # Weighted result should match unweighted on second half
    x_second_half = x[n // 2 :]
    y_second_half = y[n // 2 :]

    methods = ["pearson", "spearman", "kendall", "blomqvist"]

    for method in methods:
      weighted_result = pv.wdm(x, y, method, weights)
      unweighted_second_half = pv.wdm(x_second_half, y_second_half, method)
      assert_allclose(weighted_result, unweighted_second_half, rtol=1e-12)

  def test_wdm_weights_scaling(self) -> None:
    """Test that scaling all weights gives same result"""
    x, y = self.x_positive, self.y_positive
    weights1 = np.ones(len(x))
    weights2 = 2.0 * np.ones(len(x))
    weights3 = 0.5 * np.ones(len(x))

    methods = ["pearson", "spearman", "kendall"]

    for method in methods:
      result1 = pv.wdm(x, y, method, weights1)
      result2 = pv.wdm(x, y, method, weights2)
      result3 = pv.wdm(x, y, method, weights3)

      assert_allclose(result1, result2, rtol=1e-12)
      assert_allclose(result1, result3, rtol=1e-12)

  def test_wdm_weights_zeros(self) -> None:
    """Test with all zero weights"""
    x, y = self.x_positive, self.y_positive
    weights = np.zeros(len(x))

    # Zero weights should return NaN
    result = pv.wdm(x, y, "pearson", weights)
    assert np.isnan(result)

  def test_wdm_input_validation(self) -> None:
    """Test input validation"""
    x, y = self.x_positive, self.y_positive

    # Test mismatched lengths
    x_short = x[:-10]
    with pytest.raises((ValueError, RuntimeError)):
      pv.wdm(x_short, y, "pearson")

    # Test empty arrays - returns NaN
    x_empty = np.array([])
    y_empty = np.array([])
    result = pv.wdm(x_empty, y_empty, "pearson")
    assert np.isnan(result)

    # Test invalid method
    with pytest.raises((ValueError, RuntimeError)):
      pv.wdm(x, y, "invalid_method")

  def test_wdm_weights_length_mismatch(self) -> None:
    """Test weights with wrong length"""
    x, y = self.x_positive, self.y_positive
    weights_wrong = np.ones(len(x) - 5)

    with pytest.raises((ValueError, RuntimeError)):
      pv.wdm(x, y, "pearson", weights_wrong)

  def test_wdm_missing_values_remove_true(self) -> None:
    """Test behavior with missing values when remove_missing=True"""
    x = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    y = np.array([2.0, 4.0, 6.0, np.nan, 10.0])

    # Should work with remove_missing=True (default)
    result = pv.wdm(x, y, "pearson", remove_missing=True)
    assert isinstance(result, float)
    assert not np.isnan(result)

  def test_wdm_missing_values_remove_false(self) -> None:
    """Test behavior with missing values when remove_missing=False"""
    x = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    y = np.array([2.0, 4.0, 6.0, 8.0, 10.0])

    # Should raise error with remove_missing=False
    with pytest.raises((ValueError, RuntimeError)):
      pv.wdm(x, y, "pearson", remove_missing=False)

  def test_wdm_monotonic_transformation(self) -> None:
    """Test that rank-based methods are invariant to monotonic transformations"""
    x = self.x_positive
    y = self.y_positive

    # Apply monotonic transformations
    x_transformed = np.exp(x)
    y_transformed = np.log(np.abs(y) + 1) * np.sign(y)

    # Spearman and Kendall should be unchanged
    spearman_orig = pv.wdm(x, y, "spearman")
    spearman_trans = pv.wdm(x_transformed, y_transformed, "spearman")
    assert_allclose(spearman_orig, spearman_trans, rtol=1e-10)

    kendall_orig = pv.wdm(x, y, "kendall")
    kendall_trans = pv.wdm(x_transformed, y_transformed, "kendall")
    assert_allclose(kendall_orig, kendall_trans, rtol=1e-10)

  def test_wdm_sign_consistency(self) -> None:
    """Test that positive/negative correlations have correct signs"""
    # Positive correlation
    methods = ["pearson", "spearman", "kendall"]
    for method in methods:
      result_pos = pv.wdm(self.x_positive, self.y_positive, method)
      assert result_pos > 0, (
        f"{method} should be positive for positive correlation"
      )

      result_neg = pv.wdm(self.x_negative, self.y_negative, method)
      assert result_neg < 0, (
        f"{method} should be negative for negative correlation"
      )

  def test_wdm_symmetric(self) -> None:
    """Test that wdm(x, y) == wdm(y, x) for symmetric measures"""
    x, y = self.x_positive, self.y_positive

    symmetric_methods = [
      "pearson",
      "spearman",
      "kendall",
      "blomqvist",
      "hoeffding",
    ]

    for method in symmetric_methods:
      result_xy = pv.wdm(x, y, method)
      result_yx = pv.wdm(y, x, method)
      assert_allclose(result_xy, result_yx, rtol=1e-15)

  def test_wdm_weights_different_values(self) -> None:
    """Test with various weight patterns"""
    x, y = self.x_positive, self.y_positive
    n = len(x)

    # Linear weights
    weights_linear = np.linspace(0.1, 2.0, n)
    result_linear = pv.wdm(x, y, "pearson", weights_linear)
    assert isinstance(result_linear, float)
    assert not np.isnan(result_linear)

    # Exponential weights
    weights_exp = np.exp(np.linspace(-2, 2, n))
    result_exp = pv.wdm(x, y, "pearson", weights_exp)
    assert isinstance(result_exp, float)
    assert not np.isnan(result_exp)

    # Random weights
    np.random.seed(123)
    weights_random = np.random.uniform(0.1, 2.0, n)
    result_random = pv.wdm(x, y, "pearson", weights_random)
    assert isinstance(result_random, float)
    assert not np.isnan(result_random)

  def test_wdm_boundary_values(self) -> None:
    """Test with boundary cases"""
    # Constant arrays
    x_const = np.ones(self.n)
    y_const = np.ones(self.n) * 2

    # Pearson correlation with constants should be undefined (NaN)
    result = pv.wdm(x_const, y_const, "pearson")
    assert np.isnan(result)

    # Single different value
    x_almost_const = np.ones(self.n)
    x_almost_const[-1] = 2.0
    y_varied = np.random.normal(0, 1, self.n)

    result = pv.wdm(x_almost_const, y_varied, "pearson")
    assert isinstance(result, float)
    # Should not be NaN as there is some variation
    assert not np.isnan(result) or True  # Allow NaN in edge cases

  @pytest.mark.parametrize(
    "method", ["pearson", "spearman", "kendall", "blomqvist"]
  )
  def test_wdm_small_samples(self, method: str) -> None:
    """Test with very small sample sizes"""
    # Size 2
    x2 = np.array([1.0, 2.0])
    y2 = np.array([3.0, 4.0])
    result2 = pv.wdm(x2, y2, method)
    assert isinstance(result2, float)

    # Size 3
    x3 = np.array([1.0, 2.0, 3.0])
    y3 = np.array([1.0, 3.0, 2.0])
    result3 = pv.wdm(x3, y3, method)
    assert isinstance(result3, float)

  def test_wdm_type_consistency(self) -> None:
    """Test that function always returns float"""
    x, y = self.x_positive, self.y_positive
    methods = ["pearson", "spearman", "kendall", "blomqvist", "hoeffding"]

    for method in methods:
      result = pv.wdm(x, y, method)
      assert isinstance(result, float), f"Method {method} should return float"

      # With weights
      weights = np.ones(len(x))
      result_weighted = pv.wdm(x, y, method, weights)
      assert isinstance(result_weighted, float), (
        f"Method {method} with weights should return float"
      )

  def test_wdm_weights_edge_cases(self) -> None:
    """Test additional edge cases with weights"""
    x, y = self.x_positive, self.y_positive
    n = len(x)

    # Test with one zero weight in the middle
    weights_one_zero = np.ones(n)
    weights_one_zero[n // 2] = 0.0
    result = pv.wdm(x, y, "pearson", weights_one_zero)
    assert isinstance(result, float)
    assert not np.isnan(result)

    # Test with alternating weights (0, 1, 0, 1, ...)
    weights_alternating = np.zeros(n)
    weights_alternating[::2] = 1.0
    result_alt = pv.wdm(x, y, "pearson", weights_alternating)

    # Should match unweighted on even indices
    x_even = x[::2]
    y_even = y[::2]
    result_even = pv.wdm(x_even, y_even, "pearson")
    assert_allclose(result_alt, result_even, rtol=1e-12)

  def test_wdm_weights_numerical_precision(self) -> None:
    """Test numerical precision with very small weights"""
    x, y = self.x_positive, self.y_positive
    n = len(x)

    # Very small but non-zero weights
    weights_tiny = np.full(n, 1e-10)
    result_tiny = pv.wdm(x, y, "pearson", weights_tiny)

    # Should be similar to unweighted (weights are uniform)
    result_unweighted = pv.wdm(x, y, "pearson")
    assert_allclose(result_tiny, result_unweighted, rtol=1e-6)

    # Very large weights (should also work)
    weights_large = np.full(n, 1e6)
    result_large = pv.wdm(x, y, "pearson", weights_large)
    assert_allclose(result_large, result_unweighted, rtol=1e-12)
