import math
import statistics
import uuid
from pathlib import Path
from typing import Any, Optional

import matplotlib
import numpy as np
import pytest
from array_api_compat import array_namespace

from pyvinecopulib.core import BicopBase

matplotlib.use("Agg")


# --- Toy conditional pair copula for the non-simplified / conditional tests ---

# The standard-normal CDF and quantile on the NumPy path come from the standard
# library rather than `scipy.special`: `scipy` is an optional extra, and the
# built wheels are tested in an environment that installs neither it nor torch,
# so a pair copula hosted by a `core` test has to work without both.
_erfc = np.vectorize(math.erfc, otypes=[float])
_norm_inv_cdf = np.vectorize(statistics.NormalDist().inv_cdf, otypes=[float])


def _std_normal_cdf(z: Any) -> Any:
  """Standard normal CDF, dispatched by array backend (torch / numpy)."""
  if type(z).__module__.split(".", 1)[0] == "torch":
    import torch

    return torch.special.ndtr(z)
  return 0.5 * _erfc(-np.asarray(z) / math.sqrt(2.0))


def _std_normal_ppf(p: Any) -> Any:
  """Standard normal quantile, dispatched by array backend (torch / numpy)."""
  if type(p).__module__.split(".", 1)[0] == "torch":
    import torch

    return torch.special.ndtri(p)
  return _norm_inv_cdf(np.asarray(p))


class GaussianBicop(BicopBase[Any]):
  """Toy conditional Gaussian pair copula (correlation depends on ``x``).

  The correlation is a Fisher-style link of a *position-weighted mean* of the
  conditioning matrix, ``rho = rho_max * tanh(scale * mean_j (j + 1) * x[:, j])``,
  so the copula is genuinely non-simplified when hosted with a
  :class:`~pyvinecopulib.core.NonSimplifiedContext`, and — because the weights
  differ per column — its output *depends on the column order of* ``x`` (used to
  pin the C1 order). With ``x=None`` the correlation is ``base_rho``. Capping at
  ``rho_max < 1`` keeps ``rho`` away from ``±1`` (where the Gaussian copula
  degenerates and the cascade's ``[1e-10, 1-1e-10]`` clamp would break the
  round-trip); normalizing by the column count keeps it bounded across the
  varying ``x_e`` widths of a vine. Array-backend-agnostic (numpy / torch); has
  closed-form ``hfunc`` / ``hinv`` so the vine round-trip is exact.
  """

  supports_batched: bool = False

  def __init__(
    self,
    *,
    scale: float = 1.0,
    base_rho: float = 0.0,
    rho_max: float = 0.8,
  ) -> None:
    self._scale = float(scale)
    self._base_rho = float(base_rho)
    self._rho_max = float(rho_max)

  def _rho(self, u: Any, x: Optional[Any]) -> Any:
    """Per-row correlation from the (position-weighted) conditioning ``x``."""
    xp = array_namespace(u)
    n = u.shape[0]
    if x is None:
      return xp.full((n,), self._base_rho, dtype=u.dtype, device=u.device)
    xa: Any = x
    k = xa.shape[1]
    # Distinct per-column weights -> the link is sensitive to the x column
    # order (pins the C1 contract), unlike a plain sum. Normalize by k and cap
    # at rho_max so rho stays well away from +-1 across varying x_e widths.
    weights = xp.arange(1, k + 1, dtype=u.dtype, device=u.device)
    z = self._scale * xp.sum(xa * weights, axis=-1) / k
    return self._rho_max * xp.tanh(z)

  def pdf(self, u: Any, x: Optional[Any] = None) -> Any:
    xp = array_namespace(u)
    uc = xp.clip(u, 1e-10, 1.0 - 1e-10)
    z1, z2 = _std_normal_ppf(uc[:, 0]), _std_normal_ppf(uc[:, 1])
    rho = self._rho(u, x)
    one_minus = 1.0 - rho * rho
    quad = 2.0 * rho * z1 * z2 - rho * rho * (z1 * z1 + z2 * z2)
    return xp.exp(quad / (2.0 * one_minus)) / xp.sqrt(one_minus)

  def hfunc1(self, u: Any, x: Optional[Any] = None) -> Any:
    # P(U2 <= u2 | U1 = u1) = Phi((z2 - rho z1) / sqrt(1 - rho^2)).
    xp = array_namespace(u)
    uc = xp.clip(u, 1e-10, 1.0 - 1e-10)
    z1, z2 = _std_normal_ppf(uc[:, 0]), _std_normal_ppf(uc[:, 1])
    rho = self._rho(u, x)
    return _std_normal_cdf((z2 - rho * z1) / xp.sqrt(1.0 - rho * rho))

  def hfunc2(self, u: Any, x: Optional[Any] = None) -> Any:
    # P(U1 <= u1 | U2 = u2) = Phi((z1 - rho z2) / sqrt(1 - rho^2)).
    xp = array_namespace(u)
    uc = xp.clip(u, 1e-10, 1.0 - 1e-10)
    z1, z2 = _std_normal_ppf(uc[:, 0]), _std_normal_ppf(uc[:, 1])
    rho = self._rho(u, x)
    return _std_normal_cdf((z1 - rho * z2) / xp.sqrt(1.0 - rho * rho))

  def hinv1(self, u: Any, x: Optional[Any] = None) -> Any:
    # Invert hfunc1 w.r.t. u2: u = [u1, p] -> z2 = rho z1 + sqrt(1-rho^2) Phi^-1(p).
    xp = array_namespace(u)
    uc = xp.clip(u, 1e-10, 1.0 - 1e-10)
    z1, zp = _std_normal_ppf(uc[:, 0]), _std_normal_ppf(uc[:, 1])
    rho = self._rho(u, x)
    return _std_normal_cdf(rho * z1 + xp.sqrt(1.0 - rho * rho) * zp)

  def hinv2(self, u: Any, x: Optional[Any] = None) -> Any:
    # Invert hfunc2 w.r.t. u1: u = [p, u2] -> z1 = rho z2 + sqrt(1-rho^2) Phi^-1(p).
    xp = array_namespace(u)
    uc = xp.clip(u, 1e-10, 1.0 - 1e-10)
    zp, z2 = _std_normal_ppf(uc[:, 0]), _std_normal_ppf(uc[:, 1])
    rho = self._rho(u, x)
    return _std_normal_cdf(rho * z2 + xp.sqrt(1.0 - rho * rho) * zp)


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
