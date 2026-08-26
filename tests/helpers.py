import contextlib
from typing import Any

import numpy as np
from numpy.typing import NDArray


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


def compare_kde1d(kde1: Any, kde2: Any, from_grid: bool = False) -> None:
  # Check if models are fitted by testing grid_points
  is_fitted1 = len(kde1.grid_points) > 0
  is_fitted2 = len(kde2.grid_points) > 0

  # Both should have same fitted status
  assert is_fitted1 == is_fitted2, (
    f"Fitted status mismatch: {is_fitted1} != {is_fitted2}"
  )

  # Always compare basic properties
  attrs = ["xmin", "xmax", "type"]

  if is_fitted1:
    # For fitted models: compare grid data and prob0
    attrs += ["prob0", "grid_points", "values"]
  else:
    # For unfitted models: compare fitting parameters
    attrs += ["multiplier", "bandwidth", "degree"]

  compare_properties(kde1, kde2, attrs)


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
