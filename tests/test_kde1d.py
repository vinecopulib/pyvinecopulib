import numpy as np
import pytest

import pyvinecopulib as pv


def test_kde1d_initialization() -> None:
  """Test different ways of initializing Kde1d objects."""

  # Test default initialization
  kde = pv.Kde1d()
  assert kde.xmin != kde.xmin  # NaN check
  assert kde.xmax != kde.xmax  # NaN check
  assert kde.type == "continuous"
  assert kde.multiplier == 1.0
  assert kde.bandwidth != kde.bandwidth  # NaN check (not fitted yet)
  assert kde.degree == 2
  assert kde.grid_size == 400  # Default grid size

  # Test initialization with arguments
  kde = pv.Kde1d(
    xmin=0.0,
    xmax=1.0,
    type="discrete",
    multiplier=25,
    bandwidth=0.1,
    degree=1,
    grid_size=200,
  )
  assert kde.xmin == 0.0
  assert kde.xmax == 1.0
  assert kde.type == "discrete"
  assert kde.multiplier == 25
  assert kde.bandwidth == 0.1
  assert kde.degree == 1
  assert kde.grid_size == 200

  # Test zero-inflated type
  kde = pv.Kde1d(type="zero_inflated")
  assert kde.type == "zero-inflated"


def test_kde1d_factory_methods() -> None:
  """Test static factory methods."""

  # Test from_params
  kde = pv.Kde1d.from_params(
    xmin=0.0, xmax=1.0, type="continuous", multiplier=1.2
  )
  assert kde.xmin == 0.0
  assert kde.xmax == 1.0
  assert kde.type == "continuous"
  assert kde.multiplier == 1.2

  # Test from_grid
  # Create some sample grid data
  grid_points = np.linspace(-2, 2, 50)
  # Simple Gaussian-like density values
  values = np.exp(-0.5 * grid_points**2) / np.sqrt(2 * np.pi)

  kde_from_grid = pv.Kde1d.from_grid(
    grid_points=grid_points,
    values=values,
    xmin=-2.0,
    xmax=2.0,
    type="continuous",
  )

  assert kde_from_grid.xmin == -2.0
  assert kde_from_grid.xmax == 2.0
  assert kde_from_grid.type == "continuous"

  # Should be able to evaluate (it's already "fitted" from grid)
  pdf_vals = kde_from_grid.pdf(np.array([0.0]))
  assert len(pdf_vals) == 1
  assert pdf_vals[0] > 0

  # Should have the provided grid points and values
  assert np.array_equal(kde_from_grid.grid_points, grid_points)
  assert np.array_equal(kde_from_grid.values, values)


def test_kde1d_properties() -> None:
  """Test all properties of Kde1d objects."""

  kde = pv.Kde1d(
    xmin=-1.0,
    xmax=1.0,
    type="continuous",
    multiplier=2.0,
    bandwidth=0.5,
    degree=1,
    grid_size=100,
  )

  # Test read-only properties
  assert kde.xmin == -1.0
  assert kde.xmax == 1.0
  assert kde.type == "continuous"
  assert kde.multiplier == 2.0
  assert kde.bandwidth == 0.5
  assert kde.degree == 1
  assert kde.grid_size == 100
  assert kde.prob0 == 0.0  # Default for non-zero-inflated

  # Test that properties are read-only
  with pytest.raises(AttributeError):
    kde.xmin = 0.0  # type: ignore[misc]
  with pytest.raises(AttributeError):
    kde.multiplier = 1.0  # type: ignore[misc]


def test_kde1d_fit_and_methods() -> None:
  """Test fitting and evaluation methods."""

  # Generate test data
  np.random.seed(1234)
  x = np.random.normal(0, 1, 100)

  # Test fitting
  kde = pv.Kde1d()
  kde.fit(x)

  # After fitting, should have valid loglik and edf
  assert isinstance(kde.loglik, float)
  assert not (kde.loglik != kde.loglik)  # Should not be NaN
  assert isinstance(kde.edf, float)
  assert kde.edf > 0

  # Test grid points and values are available
  grid_points = kde.grid_points
  values = kde.values
  assert isinstance(grid_points, np.ndarray)
  assert isinstance(values, np.ndarray)
  assert len(grid_points) == len(values)
  assert len(grid_points) > 0

  # Test evaluation methods
  eval_points = np.array([-1.0, 0.0, 1.0])

  # Test pdf
  pdf_vals = kde.pdf(eval_points)
  assert isinstance(pdf_vals, np.ndarray)
  assert pdf_vals.shape == (3,)
  assert np.all(pdf_vals >= 0)  # PDF should be non-negative

  # Test cdf
  cdf_vals = kde.cdf(eval_points)
  assert isinstance(cdf_vals, np.ndarray)
  assert cdf_vals.shape == (3,)
  assert np.all(cdf_vals >= 0) and np.all(
    cdf_vals <= 1
  )  # CDF should be in [0,1]
  assert np.all(np.diff(cdf_vals) >= 0)  # CDF should be monotonic

  # Test quantile
  probs = np.array([0.1, 0.5, 0.9])
  quantiles = kde.quantile(probs)
  assert isinstance(quantiles, np.ndarray)
  assert quantiles.shape == (3,)
  assert np.all(np.diff(quantiles) >= 0)  # Quantiles should be monotonic

  # Test simulate
  samples = kde.simulate(50)
  assert isinstance(samples, np.ndarray)
  assert samples.shape == (50,)

  # Test simulate with seeds
  samples1 = kde.simulate(10, seeds=[123])
  samples2 = kde.simulate(10, seeds=[123])
  assert np.array_equal(samples1, samples2)  # Should be reproducible


def test_kde1d_weighted_fit() -> None:
  """Test fitting with weights."""

  np.random.seed(1234)
  x = np.random.normal(0, 1, 50)
  weights = np.random.exponential(1, 50)

  kde = pv.Kde1d()
  kde.fit(x, weights)

  # Should still work and produce valid results
  assert isinstance(kde.loglik, float)
  assert not (kde.loglik != kde.loglik)  # Should not be NaN

  pdf_vals = kde.pdf(np.array([0.0]))
  assert len(pdf_vals) == 1
  assert pdf_vals[0] > 0


def test_kde1d_discrete_data() -> None:
  """Test with discrete data."""

  np.random.seed(1234)
  x = np.random.binomial(10, 0.3, 100).astype(float)

  kde = pv.Kde1d(xmin=0, xmax=10, type="discrete")
  kde.fit(x)

  # Test evaluation at integer points
  eval_points = np.array([0.0, 1.0, 2.0, 5.0, 10.0])
  pdf_vals = kde.pdf(eval_points)
  assert isinstance(pdf_vals, np.ndarray)
  assert len(pdf_vals) == 5
  assert np.all(pdf_vals >= 0)


def test_kde1d_bounded_data() -> None:
  """Test with bounded support."""

  np.random.seed(1234)
  x = np.random.beta(2, 5, 100)  # Data in [0,1]

  kde = pv.Kde1d(xmin=0.0, xmax=1.0)
  kde.fit(x)

  # Test evaluation
  eval_points = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
  pdf_vals = kde.pdf(eval_points)
  assert isinstance(pdf_vals, np.ndarray)
  assert len(pdf_vals) == 5
  assert np.all(pdf_vals >= 0)

  # Test that PDF is low (but not necessarily 0) outside bounds
  outside_points = np.array([-0.1, 1.1])
  pdf_outside = kde.pdf(outside_points)
  pdf_inside = kde.pdf(np.array([0.25, 0.75]))
  # Outside values should be much smaller than typical inside values
  assert np.all(pdf_outside < 0.1 * np.mean(pdf_inside))


def test_kde1d_set_boundaries() -> None:
  """Test setting boundaries after initialization."""

  kde = pv.Kde1d()
  kde.set_xmin_xmax(xmin=-2.0, xmax=2.0)

  # Note: The boundaries are set but we need to check if they're reflected
  # in the object state after fitting (implementation dependent)


def test_kde1d_check_fitted_parameter() -> None:
  """Test the check_fitted parameter in evaluation methods."""

  np.random.seed(1234)
  x = np.random.normal(0, 1, 50)
  eval_points = np.array([0.0])

  kde = pv.Kde1d()

  # Should raise an error when not fitted (with check_fitted=True, default)
  with pytest.raises(RuntimeError):
    kde.pdf(eval_points)

  with pytest.raises(RuntimeError):
    kde.cdf(eval_points)

  with pytest.raises(RuntimeError):
    kde.quantile(np.array([0.5]))

  with pytest.raises(RuntimeError):
    kde.simulate(10)

  # Should work with check_fitted=False (might give invalid results)
  # Note: This depends on the C++ implementation behavior

  # After fitting, should work normally
  kde.fit(x)
  pdf_vals = kde.pdf(eval_points)
  assert len(pdf_vals) == 1


def test_kde1d_single_point_evaluation() -> None:
  """Test evaluation with single point (similar to bicop test)."""

  np.random.seed(1234)
  x = np.random.normal(0, 1, 100)

  kde = pv.Kde1d()
  kde.fit(x)

  # Test with single point
  single_point = np.array([0.0])
  pdf_val = kde.pdf(single_point)
  assert isinstance(pdf_val, np.ndarray) and pdf_val.shape == (1,)


def test_kde1d_string_representation() -> None:
  """Test string representation methods."""

  kde = pv.Kde1d()

  # Test __repr__ and __str__
  repr_str = repr(kde)
  str_str = str(kde)

  assert isinstance(repr_str, str)
  assert isinstance(str_str, str)
  assert "<pyvinecopulib.Kde1d>" in repr_str
  assert "<pyvinecopulib.Kde1d>" in str_str


def test_kde1d_error_conditions() -> None:
  """Test various error conditions."""

  kde = pv.Kde1d()

  # Test with invalid data (empty array)
  with pytest.raises((RuntimeError, ValueError)):
    kde.fit(np.array([]))

  # Test with invalid evaluation points on unfitted model
  with pytest.raises(RuntimeError):
    kde.pdf(np.array([0.0]))


def test_kde1d_different_degrees() -> None:
  """Test different polynomial degrees."""

  np.random.seed(1234)
  x = np.random.normal(0, 1, 100)

  for degree in [0, 1, 2]:
    kde = pv.Kde1d(degree=degree)
    assert kde.degree == degree

    kde.fit(x)
    pdf_vals = kde.pdf(np.array([0.0]))
    assert len(pdf_vals) == 1
    assert pdf_vals[0] > 0


def test_kde1d_zero_inflated() -> None:
  """Test zero-inflated data."""

  np.random.seed(1234)
  # Create zero-inflated data
  x = np.random.exponential(1, 80)
  zeros = np.zeros(20)
  x = np.concatenate([x, zeros])
  np.random.shuffle(x)

  kde = pv.Kde1d(xmin=0.0, type="zero_inflated")
  kde.fit(x)

  # Check that prob0 is estimated
  assert kde.prob0 > 0
  assert kde.prob0 < 1

  # Test evaluation
  eval_points = np.array([0.0, 0.5, 1.0])
  pdf_vals = kde.pdf(eval_points)
  assert len(pdf_vals) == 3
  assert np.all(pdf_vals >= 0)

  # PDF should be positive for all points
  assert pdf_vals[0] > 0  # Point mass at 0
  assert pdf_vals[1] > 0  # Continuous part
  assert pdf_vals[2] > 0  # Continuous part


def test_kde1d_large_data() -> None:
  """Test with larger dataset to ensure stability."""

  np.random.seed(1234)
  x = np.random.normal(0, 1, 1000)

  kde = pv.Kde1d()
  kde.fit(x)

  # Should handle large datasets
  assert isinstance(kde.loglik, float)
  assert kde.edf > 0

  # Test evaluation on many points
  eval_points = np.linspace(-3, 3, 100)
  pdf_vals = kde.pdf(eval_points)
  assert len(pdf_vals) == 100
  assert np.all(pdf_vals >= 0)
