"""Vectorized monotone root-finder (the backend-agnostic ``hinv`` / ``icdf`` default).

Array-backend-agnostic (numpy / torch) via the Array API: resolves the
namespace from the target array and uses only standard elementwise ops. Grad is
disabled by the *caller's* evaluation context (torch ``no_grad`` /
``nullcontext``), not here, so this stays a pure function. Arrays are typed
``Any`` per the ``pyvinecopulib.core`` typing policy (the Array API namespace is
untyped).
"""

from __future__ import annotations

import math
from typing import Any, Callable

from array_api_compat import array_namespace

__all__ = ["solve_increasing"]


def _is_finite_bound(xp: Any, bound: Any, broadcast: Any) -> bool:
  """Whether ``bound`` is finite everywhere.

  Scalars answer without touching the array namespace, which keeps the
  h-function inverses (always the unit interval) free of a device
  synchronization on every call.

  Parameters
  ----------
  xp : module
      Array namespace resolved from the target values.
  bound : float or array
      The bracket endpoint as supplied by the caller.
  broadcast : array
      The bound already broadcast against the target values.

  Returns
  -------
  bool
      ``True`` if no entry is infinite or NaN.
  """
  try:
    return math.isfinite(bound)
  except TypeError:
    return bool(xp.all(xp.isfinite(broadcast)))


def _bracket(
  xp: Any,
  f: Callable[[Any], Any],
  p: Any,
  a: Any,
  b: Any,
  finite_lo: bool,
  finite_hi: bool,
  max_expand: int,
) -> tuple[Any, Any]:
  """Widen an infinite bracket outward until it contains the solution.

  Each step doubles the interval on the open side, so the reachable range
  grows as ``2 ** max_expand`` and the default cap covers any float64
  magnitude. Endpoints the caller pinned to a finite value are never moved.

  Parameters
  ----------
  xp : module
      Array namespace.
  f : callable
      The monotone increasing function being inverted.
  p : array
      Target values.
  a, b : array
      Current bracket, already broadcast against ``p``.
  finite_lo, finite_hi : bool
      Whether the caller's ``lo`` / ``hi`` were finite (and so must be kept).
  max_expand : int
      Maximum number of doublings.

  Returns
  -------
  tuple of array
      The widened ``(a, b)``.
  """
  one = xp.ones_like(p)
  # Seed the open side one unit past the closed one, so the first doubling has
  # a finite width to work with.
  if not finite_lo:
    a = xp.where(xp.isfinite(a), a, (b - one) if finite_hi else -one)
  if not finite_hi:
    b = xp.where(xp.isfinite(b), b, (a + one) if finite_lo else one)

  for _ in range(max_expand):
    width = b - a
    grow_lo = (f(a) > p) if not finite_lo else None
    grow_hi = (f(b) < p) if not finite_hi else None
    if grow_lo is not None:
      a = xp.where(grow_lo, a - width, a)
    if grow_hi is not None:
      b = xp.where(grow_hi, b + width, b)
    done = True
    if grow_lo is not None:
      done = done and not bool(xp.any(grow_lo))
    if grow_hi is not None:
      done = done and not bool(xp.any(grow_hi))
    if done:
      break
  return a, b


def solve_increasing(
  f: Callable[[Any], Any],
  p: Any,
  *,
  lo: Any = 0.0,
  hi: Any = 1.0,
  n_iter: int = 50,
  max_expand: int = 64,
) -> Any:
  """Solve ``f(x) = p`` for ``x`` in ``[lo, hi]``, elementwise.

  ``f`` must be monotone increasing per element. Bisection with ``n_iter``
  steps gives a bracket width ``(hi - lo) * 2 ** -n_iter`` — with the default 50
  that is far below the ``[1e-10, 1 - 1e-10]`` clamp the h-functions impose, so
  the result is exact to that floor. A superlinear (ITP) upgrade can drop in
  behind this signature later without touching callers.

  An infinite ``lo`` or ``hi`` is widened outward to a finite bracket first, so
  a distribution on the whole real line can be inverted from its ``cdf`` alone.
  Finite brackets skip that search entirely and cost exactly what they always
  did.

  Parameters
  ----------
  f : callable
      Maps an array of candidates ``x`` to ``f(x)``, monotone increasing.
  p : array, shape (n,)
      Target values.
  lo, hi : float or array, optional
      Search bracket (default the unit interval). Arrays broadcast against
      ``p``, so a per-element support is allowed; infinite entries are widened.
  n_iter : int
      Number of bisection steps.
  max_expand : int
      Maximum bracket doublings when ``lo`` or ``hi`` is infinite.

  Returns
  -------
  array, shape (n,)
      Solutions clamped to ``[lo, hi]``.

  Raises
  ------
  ValueError
      If a bracket endpoint is infinite and ``max_expand`` is ``0``, since
      bisecting an infinite interval yields ``nan``.
  """
  xp = array_namespace(p)
  zeros = xp.zeros_like(p)
  a = zeros + lo
  b = zeros + hi
  clip_lo, clip_hi = a, b

  finite_lo = _is_finite_bound(xp, lo, a)
  finite_hi = _is_finite_bound(xp, hi, b)
  if not (finite_lo and finite_hi):
    if max_expand <= 0:
      raise ValueError(
        "an infinite bracket needs max_expand > 0; bisecting "
        f"[{lo}, {hi}] directly would return nan"
      )
    a, b = _bracket(xp, f, p, a, b, finite_lo, finite_hi, max_expand)

  for _ in range(n_iter):
    mid = 0.5 * (a + b)
    lower = f(mid) < p
    a = xp.where(lower, mid, a)
    b = xp.where(lower, b, mid)
  return xp.clip(0.5 * (a + b), clip_lo, clip_hi)
