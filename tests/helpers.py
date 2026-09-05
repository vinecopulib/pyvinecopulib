import contextlib
from typing import Any, Optional

import numpy as np
from numpy.typing import NDArray

from pyvinecopulib.core import MarginBase


def random_data(d: int = 5, n: int = 1000) -> NDArray[np.float64]:
  # Simulate some data
  np.random.seed(1234)  # seed for the random generator
  mean = np.random.normal(size=d)  # mean vector
  cov = np.random.normal(size=(d, d))  # covariance matrix
  cov = np.dot(cov.transpose(), cov)  # make it non-negative definite
  x = np.random.multivariate_normal(mean, cov, n)
  return x


def compare_properties(
  obj1: Any, obj2: Any, attrs: list[str], subclass: bool = False
) -> None:
  if subclass:
    assert issubclass(type(obj1), type(obj2)), (
      "Objects must be of the same type"
    )
  else:
    assert type(obj1) is type(obj2), "Objects must be of the same type"
  for attr in attrs:
    val1 = getattr(obj1, attr)
    val2 = getattr(obj2, attr)
    if isinstance(val1, np.ndarray):
      assert isinstance(val2, np.ndarray) and np.array_equal(val1, val2), (
        f"Mismatch in {attr}: {val1} != {val2}"
      )
    else:
      assert val1 == val2, f"Mismatch in {attr}: {val1} != {val2}"


def compare_bicop(cop1: Any, cop2: Any) -> None:
  attrs = ["family", "rotation", "parameters", "var_types"]
  compare_properties(cop1, cop2, attrs)


def compare_vinecop(cop1: Any, cop2: Any) -> None:
  attrs = ["dim", "trunc_lvl", "order", "matrix"]
  compare_properties(cop1, cop2, attrs)

  d = cop1.dim
  trunc_lvl = cop1.trunc_lvl
  for t in range(trunc_lvl):
    for e in range(d - t - 1):
      bicop1 = cop1.get_pair_copula(t, e)
      bicop2 = cop2.get_pair_copula(t, e)
      compare_bicop(bicop1, bicop2)


def compare_kde1d(kde1: Any, kde2: Any) -> None:
  is_fitted1 = kde1.is_fitted
  is_fitted2 = kde2.is_fitted

  # Both should have same fitted status
  assert is_fitted1 == is_fitted2, (
    f"Fitted status mismatch: {is_fitted1} != {is_fitted2}"
  )

  # Always compare basic properties
  attrs = ["xmin", "xmax", "type"]

  if is_fitted1:
    # A fitted estimator is more than its evaluation grid: its controls and
    # attained diagnostics govern inspection and any later refit.
    attrs += [
      "prob0",
      "grid_points",
      "values",
      "multiplier",
      "bandwidth",
      "degree",
      "grid_size",
      "boundary_repair",
      "edf",
      "n_parameters",
    ]
  else:
    # For unfitted models: compare fitting parameters. Every one of them has
    # to be here -- a key the state carries but this list skips is a key that
    # can stop round-tripping without a test noticing.
    attrs += [
      "multiplier",
      "bandwidth",
      "degree",
      "grid_size",
      "boundary_repair",
    ]

  compare_properties(kde1, kde2, attrs)
  if is_fitted1:
    assert kde1.loglik() == kde2.loglik()


def compare_rvinestructure(
  struct1: Any, struct2: Any, subclass: bool = False
) -> None:
  attrs = ["dim", "order", "trunc_lvl", "matrix"]
  compare_properties(struct1, struct2, attrs, subclass)


def assert_on_device(
  module: Any,
  device: str,
  *outputs: Any,
  extra: Any = (),
) -> None:
  """Assert every registered tensor and every output sits on ``device``.

  Walks ``named_parameters()`` and ``named_buffers()`` recursively, then
  every tensor reachable in ``outputs``. ``extra`` covers state that is
  deliberately unregistered -- notably ``TorchVinecop._batched``, installed
  via ``object.__setattr__`` so it stays out of ``state_dict()`` and is
  therefore invisible to ``named_buffers()``.

  Compares ``device.type``, so ``"cuda"`` matches ``"cuda:0"``.
  """
  import torch

  want = torch.device(device).type
  for name, t in list(module.named_parameters()) + list(module.named_buffers()):
    assert t.device.type == want, (
      f"{type(module).__name__}.{name} on {t.device}, expected {want}"
    )

  def walk(obj: Any, path: str) -> None:
    if isinstance(obj, torch.Tensor):
      assert obj.device.type == want, f"{path} on {obj.device}, expected {want}"
    elif isinstance(obj, (list, tuple)):
      for i, v in enumerate(obj):
        walk(v, f"{path}[{i}]")
    elif isinstance(obj, dict):
      for k, v in obj.items():
        walk(v, f"{path}[{k!r}]")
    elif isinstance(obj, torch.nn.Module):
      assert_on_device(obj, device)

  for i, out in enumerate(outputs):
    walk(out, f"output[{i}]")
  for i, ex in enumerate(
    extra if isinstance(extra, (list, tuple)) else [extra]
  ):
    if ex is not None:
      walk(ex, f"extra[{i}]")


class TransferCounts:
  """Device<->host movements tallied by :func:`count_transfers`."""

  def __init__(self) -> None:
    self.d2h = 0
    self.h2d = 0
    self.scalars = 0
    self.ops: list[str] = []

  def assert_no_d2h(self, what: str) -> None:
    assert self.d2h == 0 and self.scalars == 0, (
      f"{what} moved data off the device: {self.d2h} copies, "
      f"{self.scalars} scalar reads via {sorted(set(self.ops))}"
    )

  def __repr__(self) -> str:
    return (
      f"TransferCounts(d2h={self.d2h}, h2d={self.h2d}, scalars={self.scalars})"
    )


@contextlib.contextmanager
def count_transfers(device: str) -> Any:
  """Count device<->host movement inside the block.

  Tallies ``aten::_local_scalar_dense`` -- what ``.item()``, ``float(t)``,
  ``int(t)`` and ``bool(t)`` all lower to, so unlike patching ``Tensor.item``
  this sees the last three -- and cross-device ``copy_`` / ``_to_copy``.

  A no-op on CPU, where there is no device to leave.
  """
  import torch
  from torch.utils._python_dispatch import TorchDispatchMode

  counts = TransferCounts()
  if torch.device(device).type != "cuda":
    yield counts
    return

  class _Counter(TorchDispatchMode):
    def __torch_dispatch__(
      self,
      func: Any,
      types: Any,
      args: "tuple[Any, ...]" = (),
      kwargs: "dict[str, Any] | None" = None,
    ) -> Any:
      kwargs = kwargs or {}
      name = str(func)
      if "_local_scalar_dense" in name:
        counts.scalars += 1
        counts.ops.append(name)
      out = func(*args, **kwargs)
      if "copy" in name:
        src = args[1] if len(args) > 1 else None
        dst = args[0] if args else None
        st = getattr(src, "device", None)
        dt = getattr(dst, "device", None)
        if st is not None and dt is not None and st.type != dt.type:
          if st.type == "cuda":
            counts.d2h += 1
          else:
            counts.h2d += 1
          counts.ops.append(name)
      return out

  with _Counter():
    yield counts


class AtomicMargin(MarginBase[NDArray[np.float64]]):
  """A margin whose mass is not a density, with an atomic inverse CDF.

  Every shipped margin has a density, and every one that declares atoms
  declares `var_type="d"` as well, so two guards have no shipped exercise: the
  refusal of a density-less margin, and the quadrature's degeneration to a
  weighted sum over atoms. This double supplies both -- `cdf` / `icdf` are the
  empirical ones and `pdf` declines.
  """

  def __init__(self) -> None:
    self._sorted: Optional[NDArray[np.float64]] = None

  @property
  def is_fitted(self) -> bool:
    return self._sorted is not None

  @property
  def support(self) -> tuple[float, float]:
    assert self._sorted is not None
    return (float(self._sorted[0]), float(self._sorted[-1]))

  def fit(
    self,
    y: NDArray[np.float64],
    /,
    controls: object = None,
    *,
    x: Optional[NDArray[np.float64]] = None,
    weights: Optional[NDArray[np.float64]] = None,
  ) -> "AtomicMargin":
    self._sorted = np.sort(np.asarray(y, dtype=float).ravel())
    return self

  def pdf(
    self, y: NDArray[np.float64], /, *, x: Optional[NDArray[np.float64]] = None
  ) -> NDArray[np.float64]:
    raise NotImplementedError(
      "an atomic distribution has no density with respect to Lebesgue measure"
    )

  def cdf(
    self, y: NDArray[np.float64], /, *, x: Optional[NDArray[np.float64]] = None
  ) -> NDArray[np.float64]:
    assert self._sorted is not None
    ranks = np.searchsorted(self._sorted, np.asarray(y, dtype=float), "right")
    return np.asarray(ranks / (self._sorted.size + 1.0), dtype=float)

  def icdf(
    self, p: NDArray[np.float64], /, *, x: Optional[NDArray[np.float64]] = None
  ) -> NDArray[np.float64]:
    assert self._sorted is not None
    n = self._sorted.size
    idx = np.clip(
      np.ceil(np.asarray(p, dtype=float) * (n + 1.0)) - 1.0, 0, n - 1
    )
    return self._sorted[idx.astype(int)]


def widen(value: object) -> Any:
  """Widen a value to ``Any``.

  ``nn.Module.__getattr__`` hands back a ``Tensor | Parameter | Module`` union
  that no static checker can narrow, and ``Vinedist.margins`` is typed by the
  margin protocol rather than by the concrete class. Both are read through this.
  """
  return value


def run_without(package: str, body: str) -> None:
  """Assert ``body`` succeeds in a fresh interpreter that cannot import
  ``package``.

  Whether a module avoids importing an optional dependency can only be checked
  by making the import fail and seeing whether anything notices, and that needs
  a separate interpreter. ``body`` signals its own verdict with ``sys.exit``: 0
  for the behavior under test, anything else for a miss.

  Parameters
  ----------
  package : str
      Top-level package name to block.
  body : str
      Script to run once the block is installed. ``sys`` is already imported.

  Returns
  -------
  None
  """
  import subprocess
  import sys as _sys

  # `level` decides whether the name is the third-party package at all: a
  # relative import resolves inside the importing package, so
  # `from .scipy import ...` arrives here as name="scipy", level=1 and must
  # not be blocked. Ignoring `level` would make a submodule unimportable
  # whenever it shares a name with the package being blocked.
  preamble = (
    "import sys, builtins\n"
    "real = builtins.__import__\n"
    "def blocked(name, globals=None, locals=None, fromlist=(), level=0):\n"
    f"  if level == 0 and name.split('.')[0] == {package!r}:\n"
    f"    raise ImportError('no {package}')\n"
    "  return real(name, globals, locals, fromlist, level)\n"
    "builtins.__import__ = blocked\n"
  )
  result = subprocess.run(  # noqa: S603
    [_sys.executable, "-c", preamble + body], capture_output=True, text=True
  )
  assert result.returncode == 0, result.stdout + result.stderr
