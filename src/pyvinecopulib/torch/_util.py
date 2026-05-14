"""Utilities for the torch bicop submodule.

``solve_ITP`` is a vectorized root-finder used to invert h-functions
(``hinv1`` / ``hinv2``). It is the Interpolate-Truncate-Project method of
Oliveira & Takahashi (2020), with guaranteed average performance strictly
better than bisection. The implementation is adapted from torchvinecopulib.
"""

from __future__ import annotations

from typing import Callable

import torch

_EPS: float = 1e-12


@torch.no_grad()
def solve_ITP(
  fun: Callable[[torch.Tensor], torch.Tensor],
  x_a: torch.Tensor,
  x_b: torch.Tensor,
  epsilon: float = _EPS,
  num_iter_max: int = 31,
  k_1: float = 0.2,
) -> torch.Tensor:
  """Vectorized ITP root-finder for `fun(x) = 0` on `[x_a, x_b]`.

  Args:
    fun: Function from a tensor of candidate roots to a tensor of residuals.
    x_a: Lower bracket; tensor of the same shape as the desired output.
    x_b: Upper bracket; same shape as ``x_a``.
    epsilon: Convergence tolerance on the bracket width.
    num_iter_max: Maximum number of iterations.
    k_1: ITP truncation scaling factor.

  Returns:
    Tensor of approximate roots, same shape as ``x_a``.
  """
  y_a, y_b = fun(x_a), fun(x_b)
  x_a = torch.where(
    condition=y_b.abs() < epsilon,
    input=x_b - epsilon * num_iter_max,
    other=x_a,
  )
  x_b = torch.where(
    condition=y_a.abs() < epsilon,
    input=x_a + epsilon * num_iter_max,
    other=x_b,
  )
  y_a, y_b, x_wid = fun(x_a), fun(x_b), x_b - x_a
  eps_2 = torch.as_tensor(epsilon * 2.0, device=x_a.device, dtype=x_a.dtype)
  eps_scale = epsilon * 2.0 ** (
    (x_wid / epsilon).max().clamp_min(1.0).log2().ceil().clamp_min(1.0).int()
  )
  x_half = torch.empty_like(x_wid)
  rho = torch.empty_like(x_wid)
  sigma = torch.empty_like(x_wid)
  delta = torch.empty_like(x_wid)
  for _ in range(num_iter_max):
    if (x_wid < eps_2).all():
      break
    x_half.copy_(0.5 * (x_a + x_b))
    rho.copy_(eps_scale - 0.5 * x_wid)
    x_f = (y_b * x_a - y_a * x_b) / (y_b - y_a)
    sigma.copy_(x_half - x_f)
    delta.copy_(k_1 * x_wid.square())
    x_t = torch.where(
      condition=delta <= sigma.abs(),
      input=x_f + torch.copysign(delta, sigma),
      other=x_half,
    )
    x_itp = torch.where(
      condition=rho >= (x_t - x_half).abs(),
      input=x_t,
      other=x_half - torch.copysign(rho, sigma),
    )
    y_itp = fun(x_itp)
    idx = y_itp > 0.0
    x_b[idx], y_b[idx] = x_itp[idx], y_itp[idx]
    idx = ~idx
    x_a[idx], y_a[idx] = x_itp[idx], y_itp[idx]
    x_wid = x_b - x_a
    eps_scale *= 0.5
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
