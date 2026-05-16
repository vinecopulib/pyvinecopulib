"""Utilities for the torch bicop submodule.

Two vectorized root-finders for h-function inversion (``hinv1`` / ``hinv2``):

* ``solve_itp`` — the default. Fixed-iter Interpolate-Truncate-Project
  (Oliveira & Takahashi 2020). Combines false-position-style
  interpolation with bracket-preserving truncation and projection;
  super-linear convergence without derivatives. No host-side early-exit
  → data-independent, CUDA-Graph-friendly.

* ``solve_bisection`` — fixed 35-iter bisection, mirrors C++
  ``tools_eigen::invert_f``. Kept as a private fallback / regression
  baseline; ITP supersedes it at every call site we have.
"""

from __future__ import annotations

from typing import Callable

import torch


@torch.no_grad()
def solve_bisection(
  fun: Callable[[torch.Tensor], torch.Tensor],
  x_a: torch.Tensor,
  x_b: torch.Tensor,
  n_iter: int = 35,
) -> torch.Tensor:
  """Vectorized bisection root-finder for ``fun(x) = 0`` on ``[x_a, x_b]``.

  Mirrors ``tools_eigen::invert_f`` in the C++ pyvinecopulib library: fixed
  ``n_iter`` iterations, no host-side early-exit, all ops elementwise.

  Args:
    fun: Function from a tensor of candidate roots to a tensor of residuals.
      Must be monotone increasing on ``[x_a, x_b]`` (the sign convention is
      "residual negative below the root, positive above").
    x_a: Lower bracket; tensor of the same shape as the desired output.
    x_b: Upper bracket; same shape as ``x_a``.
    n_iter: Number of bisection iterations. The default ``35`` matches the
      C++ side and yields ``0.5**35 ≈ 6e-11`` accuracy.

  Returns:
    Tensor of approximate roots, same shape as ``x_a``.
  """
  x_a = x_a.clone()
  x_b = x_b.clone()
  for _ in range(n_iter):
    x_mid = 0.5 * (x_a + x_b)
    f_mid = fun(x_mid)
    below = f_mid < 0.0
    x_a = torch.where(below, x_mid, x_a)
    x_b = torch.where(below, x_b, x_mid)
  return x_a


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


def make_normal_grid(
  m: int, dtype: torch.dtype = torch.float64
) -> torch.Tensor:
  """Return the default vinecopulib normal grid: ``Phi(linspace(-3.25, 3.25, m))``.

  Mirrors ``KernelBicop::make_normal_grid`` in the C++ library. Endpoints are
  forced to exactly 0 and 1 to avoid extrapolation, matching the
  ``InterpolationGrid`` constructor.
  """
  z = torch.linspace(-3.25, 3.25, m, dtype=dtype)
  # Standard normal CDF via the error function.
  grid = 0.5 * (1.0 + torch.erf(z / (2.0**0.5)))
  grid[0] = 0.0
  grid[-1] = 1.0
  return grid


def make_linear_grid(
  m: int, dtype: torch.dtype = torch.float64
) -> torch.Tensor:
  """Return a uniform ``[0, 1]`` grid of size ``m``."""
  return torch.linspace(0.0, 1.0, m, dtype=dtype)
