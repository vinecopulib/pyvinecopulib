"""Vectorized monotone root-finder (the backend-agnostic ``hinv`` default).

Array-backend-agnostic (numpy / torch) via the Array API: resolves the
namespace from the target array and uses only standard elementwise ops. Grad is
disabled by the *caller's* evaluation context (torch ``no_grad`` /
``nullcontext``), not here, so this stays a pure function. Arrays are typed
``Any`` per the ``pyvinecopulib.core`` typing policy (the Array API namespace is
untyped).
"""

from __future__ import annotations

from typing import Any, Callable

from array_api_compat import array_namespace

__all__ = ["solve_increasing"]


def solve_increasing(
  f: Callable[[Any], Any],
  p: Any,
  *,
  lo: float = 0.0,
  hi: float = 1.0,
  n_iter: int = 50,
) -> Any:
  """Solve ``f(x) = p`` for ``x`` in ``[lo, hi]``, elementwise.

  ``f`` must be monotone increasing per element. Bisection with ``n_iter``
  steps gives a bracket width ``(hi - lo) * 2 ** -n_iter`` — with the default 50
  that is far below the ``[1e-10, 1 - 1e-10]`` clamp the h-functions impose, so
  the result is exact to that floor. A superlinear (ITP) upgrade can drop in
  behind this signature later without touching callers.

  Parameters
  ----------
  f : callable
      Maps an array of candidates ``x`` to ``f(x)``, monotone increasing.
  p : array, shape (n,)
      Target values.
  lo, hi : float
      Search bracket (default the unit interval).
  n_iter : int
      Number of bisection steps.

  Returns
  -------
  array, shape (n,)
      Solutions clamped to ``[0, 1]``.
  """
  xp = array_namespace(p)
  a = xp.full_like(p, lo)
  b = xp.full_like(p, hi)
  for _ in range(n_iter):
    mid = 0.5 * (a + b)
    lower = f(mid) < p
    a = xp.where(lower, mid, a)
    b = xp.where(lower, b, mid)
  return xp.clip(0.5 * (a + b), 0.0, 1.0)
