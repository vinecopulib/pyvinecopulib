"""Clamp bounds that keep copula arguments strictly inside the unit square.

The cascades clamp every h-function and distribution-function value away
from ``0`` and ``1``: a ``0`` or a ``1`` reaching a downstream normal
quantile is an infinity, and an argument outside ``[0, 1]`` extrapolates
off the interpolation grid.

The bound has to be representable in the working precision to do that.
``1 - 1e-10`` is not: in ``float32`` it rounds to exactly ``1.0``, so the
clamp admits the value it exists to exclude. :func:`trim_bounds` therefore
derives the upper bound from the dtype, and returns the historical
``float64`` pair unchanged so that precision's results are unmoved.
"""

from typing import Any, Tuple

_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10


def trim_bounds(xp: Any, dtype: Any) -> Tuple[float, float]:
  """Clamp bounds for ``dtype``, strictly inside ``(0, 1)``.

  Parameters
  ----------
  xp : module
      Array namespace exposing ``finfo`` (NumPy or PyTorch).
  dtype : dtype
      Floating dtype the values are held in.

  Returns
  -------
  tuple of float
      ``(lo, hi)`` with ``0 < lo < hi < 1`` in ``dtype``. For ``float64``
      this is ``(1e-10, 1 - 1e-10)``.
  """
  eps = float(xp.finfo(dtype).eps)
  return _TRIM_LO, min(_TRIM_HI, 1.0 - eps)


def trim(xp: Any, a: Any) -> Any:
  """Clamp ``a`` into the open unit interval at its own precision.

  Parameters
  ----------
  xp : module
      Array namespace of ``a`` (NumPy or PyTorch).
  a : array
      Values to clamp.

  Returns
  -------
  array
      ``a`` clamped to :func:`trim_bounds` for its dtype.
  """
  lo, hi = trim_bounds(xp, a.dtype)
  return xp.clip(a, lo, hi)
