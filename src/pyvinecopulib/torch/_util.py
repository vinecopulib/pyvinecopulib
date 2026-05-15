"""Utilities for the torch bicop submodule.

``solve_bisection`` is a fixed-iteration vectorized bisection root-finder
used to invert h-functions (``hinv1`` / ``hinv2``). It mirrors C++
``tools_eigen::invert_f``: 35 iterations, no host-side early-exit. The
fixed iteration count keeps the entire op sequence elementwise and
data-independent, which is required for ``torch.cuda.CUDAGraph`` capture.
At 35 iters the bracket-width accuracy is ``0.5**35 ≈ 6e-11``.
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
