"""Vectorized ITP root-finder for the torch bicop submodule.

Used by :meth:`TorchBicop.hinv1` / :meth:`hinv2` (and the corresponding
cache-build paths in :class:`InterpolationGrid2D`) when
``cache_integrals=False``. Grid construction lives on
:class:`InterpolationGrid2D` itself.
"""

from __future__ import annotations

from typing import Callable

import torch


@torch.no_grad()
def solve_itp(
  fun: Callable[[torch.Tensor], torch.Tensor],
  x_a: torch.Tensor,
  x_b: torch.Tensor,
  n_iter: int = 30,
  k1: float = 0.2,
  k2: float = 2.0,
  eps: float = 1e-12,
) -> torch.Tensor:
  """Vectorized ITP root-finder for ``fun(x) = 0`` on ``[x_a, x_b]``.

  Fixed-iter port of the Interpolate-Truncate-Project method of Oliveira &
  Takahashi (2020). Each iter:

  1. **Interpolation** — false-position candidate
     ``x_f = (f_b * a - f_a * b) / (f_b - f_a)``.
  2. **Truncation** — shift ``x_f`` toward the midpoint by ``k1 *
     |b-a|^k2`` (unless ``x_f`` is already closer to the midpoint than
     that).
  3. **Projection** — clamp the truncated point to lie within
     ``r = ½|b-a| - eps`` of the midpoint, guaranteeing the new bracket
     stays strictly inside the old one.
  4. Evaluate ``f`` at the projected point; update the bracket based on
     its sign.

  Cost per iter: **one** ``fun`` evaluation (same as bisection). ITP
  converges super-linearly without needing a derivative — at fixed
  ``n_iter=30`` it consistently beats bisection-at-35-iters on smooth
  AND piecewise-linear (cached-h / trapezoidal) functions, which is
  what we actually need.

  Sign convention: ``f`` must be monotone increasing on ``[x_a, x_b]``
  (negative below the root, positive above). The bracket
  ``f(x_a) ≤ 0 ≤ f(x_b)`` is enforced via a one-shot swap if the
  initial endpoints are inverted.

  Args:
    fun: Tensor → tensor residual map.
    x_a, x_b: Initial bracket; same shape as the desired output.
    n_iter: Iteration count.
    k1, k2: Standard ITP truncation parameters (defaults match the
      Oliveira-Takahashi reference).
    eps: Projection radius slack — keeps the new ``x`` strictly inside
      the bracket so subsequent ``fun`` calls don't see the boundaries.

  Returns:
    Tensor of approximate roots, same shape as ``x_a``.
  """
  x_a = x_a.clone()
  x_b = x_b.clone()
  fa = fun(x_a)
  fb = fun(x_b)
  # Defensive orientation: ensure fa <= 0 <= fb.
  swap = fa > 0
  x_a, x_b = torch.where(swap, x_b, x_a), torch.where(swap, x_a, x_b)
  fa, fb = torch.where(swap, fb, fa), torch.where(swap, fa, fb)

  for _ in range(n_iter):
    width = x_b - x_a
    mid = 0.5 * (x_a + x_b)
    # Step 1: interpolation (false position).
    denom = (fb - fa).clamp_min(1e-300)
    x_f = (fb * x_a - fa * x_b) / denom
    # Step 2: truncation toward midpoint.
    diff = mid - x_f
    sigma = torch.sign(diff)
    delta = k1 * width.abs().pow(k2)
    x_t = torch.where(delta <= diff.abs(), x_f + sigma * delta, mid)
    # Step 3: projection to keep x_itp within r of the midpoint.
    r = (0.5 * width.abs() - eps).clamp_min(0.0)
    x_itp = torch.where((x_t - mid).abs() <= r, x_t, mid - sigma * r)
    # Step 4: evaluate and update bracket.
    fx = fun(x_itp)
    below = fx < 0.0
    x_a = torch.where(below, x_itp, x_a)
    x_b = torch.where(below, x_b, x_itp)
    fa = torch.where(below, fx, fa)
    fb = torch.where(below, fb, fx)

  return 0.5 * (x_a + x_b)
