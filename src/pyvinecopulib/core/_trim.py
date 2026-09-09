"""Clamp bounds that keep copula arguments strictly inside the unit square.

The cascades clamp every h-function and distribution-function value away
from ``0`` and ``1``: a ``0`` or a ``1`` reaching a downstream normal
quantile is an infinity, and an argument outside ``[0, 1]`` extrapolates
off the interpolation grid.

The bounds have to be representable in the working precision to do that.
``1 - 1e-10`` rounds to exactly ``1.0`` in ``float32``, while ``1e-10``
rounds to zero in ``float16``. ``trim_bounds`` therefore derives safe
bounds from the dtype, while returning the historical ``float64`` pair
unchanged so that precision's results are unmoved.
"""

from types import ModuleType
from typing import Any, Optional, Tuple, cast

from array_api_compat import array_namespace

from .protocols import ArrayT

_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10


def trim_bounds(xp: ModuleType, dtype: object) -> Tuple[float, float]:
  """Clamp bounds for ``dtype``, strictly inside ``(0, 1)``.

  Parameters
  ----------
  xp : module
      Array namespace exposing ``finfo`` (NumPy or PyTorch).
  dtype : object
      Floating dtype the values are held in, as the namespace spells one.

  Returns
  -------
  tuple of float
      ``(lo, hi)`` with ``0 < lo < hi < 1`` in ``dtype``. For ``float64``
      this is ``(1e-10, 1 - 1e-10)``.
  """
  eps = float(xp.finfo(dtype).eps)
  return max(_TRIM_LO, eps * eps), min(_TRIM_HI, 1.0 - eps)


def trim(a: ArrayT, xp: Optional[ModuleType] = None) -> ArrayT:
  """Clamp ``a`` into the open unit interval at its own precision.

  For a value **this library produced** -- an h-function, a distribution
  function, an interpolated integral -- where landing on ``0`` or ``1`` is
  arithmetic rounding and clamping is the only sane answer. A value a *caller*
  supplied is a different question: one that has been through a probability
  integral transform cannot legitimately be ``0`` or ``1``, so clamping it
  converts a real defect upstream into a plausible number and hides it. Refuse
  such a value rather than passing it through here.

  Parameters
  ----------
  a : array
      Values to clamp.
  xp : module, or None, optional
      Array namespace of ``a``, resolved from ``a`` when omitted. Pass it only
      where the caller already holds it: the cascades do, and the lookup is
      per-call.

  Returns
  -------
  array
      ``a`` clamped to ``trim_bounds`` for its dtype.
  """
  ns = array_namespace(a) if xp is None else xp
  lo, hi = trim_bounds(ns, cast("Any", a).dtype)
  return cast("ArrayT", ns.clip(a, lo, hi))
