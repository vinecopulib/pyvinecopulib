import uuid
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")


# --- Fixtures for the pyvinecopulib.sklearn estimator tests ---


@pytest.fixture
def random_state() -> np.random.RandomState:
  """Fixed random state for reproducibility across estimator tests."""
  return np.random.RandomState(42)


@pytest.fixture
def sample_array_data(
  random_state: np.random.RandomState,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  """2-D multivariate normal sample for density/regression tests."""
  n_samples = 300
  mean = np.array([1.0, -0.5])
  cov = np.array([[2.0, 0.8], [0.8, 1.5]])
  X = random_state.multivariate_normal(mean, cov, n_samples)
  return X, mean, cov


@pytest.fixture
def sample_dataframe_data(
  random_state: np.random.RandomState,
):
  """Mixed-dtype DataFrame for factor-expansion tests."""
  import pandas as pd

  n_samples = 200
  X_array = random_state.multivariate_normal(
    [0, 0], [[1, 0.5], [0.5, 1]], n_samples
  )
  X_df = pd.DataFrame(
    {
      "cont1": X_array[:, 0],
      "cont2": X_array[:, 1],
      "cat1": pd.Categorical(random_state.choice(["A", "B", "C"], n_samples)),
      "discrete": random_state.randint(0, 5, n_samples),
    }
  )
  expected_expanded_cols = ["cont1", "cont2", "cat1_B", "cat1_C", "discrete"]
  return X_df, expected_expanded_cols


@pytest.fixture
def regression_data(
  random_state: np.random.RandomState,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
  """Linear regression dataset for VineRegressor tests."""
  n_samples = 300
  X = random_state.multivariate_normal([0, 0], [[1, 0.3], [0.3, 1]], n_samples)
  true_coef = np.array([1.5, -0.8])
  noise_std = 0.2
  y = X @ true_coef + noise_std * random_state.randn(n_samples)
  return X, y, true_coef, noise_std


# # Per-test unique directory, xdist-friendly (each worker gets its own base dir)
# @pytest.fixture(scope="function")
# def test_dump_folder(
#   tmp_path_factory: pytest.TempPathFactory, request: pytest.FixtureRequest
# ) -> str:
#   worker = os.getenv("PYTEST_XDIST_WORKER", "worker0")
#   base = tmp_path_factory.mktemp(f"dump-{worker}")
#   # also isolate by test node to avoid clashes inside same worker
#   d = base / request.node.name.replace(os.sep, "_")
#   d.mkdir(parents=True, exist_ok=True)
#   return str(os.fspath(d))


# If you prefer a unique file path directly:
@pytest.fixture
def unique_json_path(tmp_path: Path, request: pytest.FixtureRequest) -> Path:
  return tmp_path / f"{request.node.name}-{uuid.uuid4().hex}.json"
