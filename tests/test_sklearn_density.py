"""Tests for VineDensity class."""

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("scipy")

from scipy.stats import multivariate_normal  # noqa: E402

from pyvinecopulib.sklearn import VineDensity  # noqa: E402


@pytest.fixture
def multivariate_normal_dist():
  """Factory for creating multivariate normal distributions."""

  def _create_dist(mean, cov):
    return multivariate_normal(mean, cov)

  return _create_dist


@pytest.fixture
def fitted_density(sample_array_data):
  """Fitted VineDensity for testing."""
  X, _, _ = sample_array_data
  density = VineDensity(batch_size=50)
  density.fit(X)
  return density, X


def test_fit_properties(fitted_density):
  """Test properties after fitting."""
  density, _ = fitted_density
  assert hasattr(density, "_vine")
  assert hasattr(density, "_x_kde1d")
  assert density.n_features_in_ == 2
  assert len(density._x_kde1d) == 2


def test_score_samples_shapes_and_types(fitted_density):
  """Test score_samples method shapes and types."""
  density, X = fitted_density
  X_test = X[:20]
  log_scores = density.score_samples(X_test)

  # Check type and shape
  assert isinstance(log_scores, np.ndarray)
  assert log_scores.shape == (20,)
  assert np.all(np.isfinite(log_scores))


def test_pdf_shapes_and_types(fitted_density):
  """Test pdf method shapes and types."""
  density, X = fitted_density
  X_test = X[:20]
  densities = density.pdf(X_test)

  # Check type and shape
  assert isinstance(densities, np.ndarray)
  assert densities.shape == (20,)
  assert np.all(densities > 0)
  assert np.all(np.isfinite(densities))


def test_pdf_log_consistency(fitted_density):
  """Test consistency between pdf and score_samples."""
  density, X = fitted_density
  X_test = X[:20]
  log_scores = density.score_samples(X_test)
  densities = density.pdf(X_test)

  # Check that pdf = exp(log_scores)
  np.testing.assert_allclose(densities, np.exp(log_scores), rtol=1e-10)


def test_cdf_shapes_and_range(fitted_density):
  """CDF returns values in [0, 1] with the right shape."""
  density, X = fitted_density
  X_test = X[:20]
  cdf_vals = density.cdf(X_test, seeds=[1, 2, 3])

  assert isinstance(cdf_vals, np.ndarray)
  assert cdf_vals.shape == (20,)
  assert np.all(cdf_vals >= 0.0)
  assert np.all(cdf_vals <= 1.0)
  assert np.all(np.isfinite(cdf_vals))


def test_cdf_monotone_along_axis(fitted_density):
  """CDF is approximately monotone along a single-coordinate sweep."""
  density, X = fitted_density
  # Build a path: sweep x_1 from its 5th to 95th percentile, fix x_2 at the median.
  x1_lo, x1_hi = np.quantile(X[:, 0], [0.05, 0.95])
  x2_med = np.median(X[:, 1])
  X_sweep = np.column_stack(
    [np.linspace(x1_lo, x1_hi, 30), np.full(30, x2_med)]
  )
  cdf_sweep = density.cdf(X_sweep, N=20000, seeds=[42, 43, 44])
  # Should be roughly non-decreasing along the sweep; allow some MC noise.
  diffs = np.diff(cdf_sweep)
  assert (diffs > -0.05).all(), (
    f"CDF should be approx. non-decreasing, got {diffs}"
  )
  assert cdf_sweep[-1] > cdf_sweep[0]


def test_score_method(fitted_density):
  """Test score method."""
  density, X = fitted_density
  X_test = X[:20]
  score = density.score(X_test)
  log_scores = density.score_samples(X_test)

  # Check that score is mean of score_samples
  assert isinstance(score, (float, np.floating))
  np.testing.assert_allclose(score, np.mean(log_scores))


def test_density_accuracy(sample_array_data, multivariate_normal_dist):
  """Test density estimation accuracy against true distribution."""
  X, mean, cov = sample_array_data
  true_dist = multivariate_normal_dist(mean, cov)

  density = VineDensity()
  density.fit(X)

  X_test = X[:50]

  # Get estimated log-densities
  estimated_log_densities = density.score_samples(X_test)

  # Get true log-densities
  true_log_densities = true_dist.logpdf(X_test)

  # Check correlation (should be high for good estimation)
  correlation = np.corrcoef(estimated_log_densities, true_log_densities)[0, 1]
  assert correlation > 0.7, (
    "Density estimation should correlate with true density"
  )

  # Check that mean difference is reasonable
  mean_diff = np.mean(estimated_log_densities - true_log_densities)
  assert abs(mean_diff) < 2.0, "Mean density difference should be reasonable"


@pytest.mark.parametrize("n_samples", [1, 5, 10])
def test_sample_method(fitted_density, n_samples):
  """Test sample generation."""
  density, _ = fitted_density
  samples = density.sample(n_samples)

  # Check type and shape
  assert isinstance(samples, np.ndarray)
  assert samples.shape == (n_samples, 2)
  assert np.all(np.isfinite(samples))


def test_sample_distribution(fitted_density):
  """Test that samples follow approximately the right distribution."""
  density, X = fitted_density
  n_samples = 1000
  samples = density.sample(n_samples)

  # Check that sample mean is close to fitted data mean
  sample_mean = np.mean(samples, axis=0)
  data_mean = np.mean(X, axis=0)
  np.testing.assert_allclose(sample_mean, data_mean, atol=0.3)

  # Check that sample std is reasonable
  sample_std = np.std(samples, axis=0)
  data_std = np.std(X, axis=0)
  np.testing.assert_allclose(sample_std, data_std, rtol=0.5)


def test_batch_processing(sample_array_data):
  """Test that different batch sizes give same results."""
  X, _, _ = sample_array_data
  X_test = X[:50]

  # Test with different batch sizes
  results = []
  for batch_size in [1, 10, 25, 50]:
    density = VineDensity(batch_size=batch_size)
    density.fit(X)
    log_scores = density.score_samples(X_test)
    results.append(log_scores)

  # All results should be very close
  for i in range(1, len(results)):
    np.testing.assert_allclose(results[0], results[i], rtol=1e-12)
