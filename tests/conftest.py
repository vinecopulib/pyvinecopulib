import math
import statistics
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, cast

import matplotlib
import numpy as np
import pytest
from array_api_compat import array_namespace

import pyvinecopulib as pv
from pyvinecopulib.core import BicopBase, BicopLike, VinecopBase

if TYPE_CHECKING:
  import pandas as pd

matplotlib.use("Agg")


class HostedVinecop(VinecopBase[Any]):
  """A ``VinecopBase`` over a nested list of already-fitted pair copulas.

  The minimal concrete vine the extension-layer tests evaluate through. Its
  ``_sample_uniform`` is the *compiled* draw, so a vine hosting a ``Vinecop``'s
  own pair copulas reproduces that vine's sampling bit for bit and the two can
  be compared directly.

  Parameters
  ----------
  pairs : list of list of BicopLike
      Pair copulas indexed ``[tree][edge]``.
  structure : RVineStructure
      The structure to evaluate along.
  var_types : list of str, optional
      Per-variable types; ``None`` means all continuous.
  context : ConditioningContext, optional
      Per-edge conditioning policy; ``None`` is the simplified default.
  """

  def __init__(
    self,
    pairs: Any,
    structure: Any,
    var_types: Optional[list[str]] = None,
    context: Optional[Any] = None,
  ) -> None:
    self._pairs = pairs
    self._bind_vine(structure, context, var_types=var_types)

  def set_pair_copulas(self, pair_copulas: Any) -> None:
    self._pairs = pair_copulas

  def get_pair_copula(self, tree: int, edge: int) -> BicopLike[Any]:
    return self._pairs[tree][edge]

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> Any:
    return pv.utils.sample_uniform(n, self.d, qrng, list(seeds))


def host_vinecop(
  cop: Any, var_types: Optional[list[str]] = None
) -> HostedVinecop:
  """Host a compiled ``Vinecop``'s pair copulas in a :class:`HostedVinecop`.

  Parameters
  ----------
  cop : Vinecop
      The fitted vine whose pair copulas and structure to reuse.
  var_types : list of str, optional
      Per-variable types; ``None`` reads them off ``cop``.

  Returns
  -------
  HostedVinecop
      The same model, evaluated by the array-agnostic cascades.
  """
  pairs = [
    [cop.get_pair_copula(t, e) for e in range(cop.dim - 1 - t)]
    for t in range(cop.dim - 1)
  ]
  types = list(cop.var_types) if var_types is None else var_types
  if all(v == "c" for v in types):
    types = None
  return HostedVinecop(pairs, cop.structure, types)


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


def position_weighted_mean(x: Any, ref: Any) -> Any:
  """Mean of ``x``'s columns, weighted by 1-based column position.

  The one link both conditional pair-copula doubles push a correlation
  through, and it is shared because two spellings of it would let the two
  legs of a conditional-vine comparison disagree about what ``x`` means.
  Distinct per-column weights make anything built on it sensitive to the
  *column order* of ``x``, which is what pins the C1 conditioning contract; a
  plain sum would not. Dividing by the column count keeps it bounded across a
  vine's varying ``x_e`` widths.

  ``ref`` supplies the dtype and device the weights are built on, since ``x``
  may be an integer array while the model evaluates in floating point.
  """
  xp = array_namespace(ref)
  xa: Any = x
  k = xa.shape[1]
  weights = xp.arange(1, k + 1, dtype=ref.dtype, device=ref.device)
  return xp.sum(xa * weights, axis=-1) / k


class GaussianBicop(BicopBase[Any]):
  """Toy conditional Gaussian pair copula (correlation depends on ``x``).

  The correlation is a Fisher-style link of a *position-weighted mean* of the
  conditioning matrix, ``rho = rho_max * tanh(scale * mean_j (j + 1) * x[:, j])``,
  so the copula is actually non-simplified when hosted with a
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
    # Capped at rho_max so rho stays well away from +-1 across varying
    # x_e widths.
    z = self._scale * position_weighted_mean(x, u)
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


class MinimalBicop(BicopBase[Any]):
  """The smallest valid pair copula: independence, and nothing declared twice.

  Named for what it exercises rather than for what it models -- the library's
  own ``IndependenceBicop`` is the class to reach for outside the suite.
  Implements only the abstract surface (``pdf`` / ``hfunc1`` / ``hfunc2``), so
  ``hinv1`` / ``hinv2`` / ``cdf`` / ``flip`` come from :class:`BicopBase` --
  the two inverses numerically, the latter two as the raising stubs -- and are
  what the tests hosting it exercise. Array-backend-agnostic apart from
  ``_sample_uniform``, the one hook with no array-agnostic default.
  """

  def pdf(self, u: Any, *, x: Optional[Any] = None) -> Any:
    xp = array_namespace(u)
    return xp.ones((u.shape[0],), dtype=u.dtype, device=u.device)

  def hfunc1(self, u: Any, *, x: Optional[Any] = None) -> Any:
    return u[:, 1]

  def hfunc2(self, u: Any, *, x: Optional[Any] = None) -> Any:
    return u[:, 0]

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> Any:
    rng = np.random.default_rng(seeds[0] if seeds else 0)
    return rng.uniform(size=(n, 2))


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
) -> tuple["pd.DataFrame", list[str]]:
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


@pytest.fixture
def unique_json_path(tmp_path: Path, request: pytest.FixtureRequest) -> Path:
  return tmp_path / f"{request.node.name}-{uuid.uuid4().hex}.json"


def _cuda_available() -> bool:
  """Whether a CUDA device is usable, without requiring torch to be installed."""
  try:
    import torch
  except Exception:  # torch is an optional extra
    return False
  try:
    return bool(torch.cuda.is_available())
  except Exception:  # a half-installed driver must not break collection
    return False


_HAS_CUDA = _cuda_available()


# Parameterized over both devices. The marks sit on the *param*, not on the
# fixture body, so `-m cuda` selects and `-m "not cuda"` deselects: a mark
# added from inside the fixture would land after collection and do neither.
# A test that must pin one device overrides with
# `@pytest.mark.parametrize("device", ["cpu"])`.
@pytest.fixture(
  params=[
    "cpu",
    pytest.param(
      "cuda",
      marks=[
        pytest.mark.cuda,
        pytest.mark.skipif(not _HAS_CUDA, reason="no CUDA device available"),
      ],
    ),
  ]
)
def device(request: pytest.FixtureRequest) -> str:
  """Torch device a test runs on."""
  return cast("str", request.param)


@pytest.fixture
def count_sample() -> np.ndarray:
  """400 Poisson(4) counts, including zeros."""
  return np.random.default_rng(1).poisson(4.0, size=400).astype(float)
