"""PyTorch port of :class:`vinecopulib::tools_interpolation::InterpolationGrid`.

Stores a density on a tensor-product grid in [0, 1]^2 and provides bilinear
interpolation plus trapezoidal integration along one or both axes, with the
same partial-cell handling as the C++ implementation.
"""

from __future__ import annotations


import torch
from torch import Tensor

_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10
_STRIP_FLOOR: float = 1e-4


class InterpolationGrid2D(torch.nn.Module):
  """Bilinear interpolation grid for a bivariate density on `[0, 1]^2`.

  Mirrors the C++ ``vinecopulib::tools_interpolation::InterpolationGrid``:
  the same non-uniform ``grid_points`` are used along both axes, ``values``
  stores the density at the tensor-product grid, and integration uses the
  trapezoidal rule with linear interpolation inside the partial cell that
  contains the upper limit.
  """

  values: Tensor
  grid_points: Tensor

  def __init__(
    self,
    grid_points: Tensor,
    values: Tensor,
    norm_times: int = 3,
  ) -> None:
    super().__init__()
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
      raise ValueError("values must be a square 2D tensor")
    if grid_points.ndim != 1 or grid_points.shape[0] != values.shape[0]:
      raise ValueError(
        "grid_points must be 1D and match the side length of values"
      )

    grid_points = grid_points.clone()
    # Force boundary points to exactly 0 / 1 so we never extrapolate.
    grid_points[0] = 0.0
    grid_points[-1] = 1.0

    self.register_buffer("grid_points", grid_points.contiguous())
    self.register_buffer("values", values.clone().contiguous())
    self.normalize_margins(norm_times)

  # --------------------------------------------------------------------- #
  # Margin normalization (port of InterpolationGrid::normalize_margins).  #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def normalize_margins(self, times: int) -> None:
    """Renormalize ``values`` so both margins integrate to 1.

    Same algorithm as the C++ ``normalize_margins``: alternating row / column
    trapezoidal-integral divides, repeated ``times`` rounds.
    """
    if times <= 0:
      return
    dgrid = self.grid_points[1:] - self.grid_points[:-1]  # (m-1,)
    for _ in range(times):
      row_int = 0.5 * ((self.values[:, :-1] + self.values[:, 1:]) * dgrid).sum(
        dim=-1
      )
      self.values.div_(row_int.clamp_min(1e-20).unsqueeze(-1))
      col_int = 0.5 * (
        (self.values[:-1, :] + self.values[1:, :]) * dgrid.unsqueeze(-1)
      ).sum(dim=0)
      self.values.div_(col_int.clamp_min(1e-20))

  @torch.no_grad()
  def flip(self) -> None:
    """Transpose ``values`` in place (mirror of ``KernelBicop::flip``)."""
    self.values = self.values.t().contiguous()

  # --------------------------------------------------------------------- #
  # Helpers                                                                #
  # --------------------------------------------------------------------- #

  def _cell_index(self, u: Tensor) -> Tensor:
    """Cell index for each value of ``u`` (shape preserved), clamped to ``[0, m-2]``."""
    m = self.grid_points.shape[0]
    return (
      torch.searchsorted(self.grid_points, u.contiguous(), right=False) - 1
    ).clamp(0, m - 2)

  def _int_on_grid(
    self,
    upr: Tensor,
    vals: Tensor,
  ) -> Tensor:
    """Vectorized trapezoidal integral of `(grid_points, vals)` from 0 to ``upr``.

    Args:
      upr: shape ``(*B,)`` of upper limits, each broadcast-compatible with
        ``vals[..., 0]``.
      vals: shape ``(*B, m)`` of piecewise-linear values.

    Returns:
      Tensor of shape ``(*B,)`` with the trapezoidal integrals.
    """
    grid = self.grid_points
    m = grid.shape[0]
    dgrid = grid[1:] - grid[:-1]  # (m-1,)

    trap = 0.5 * (vals[..., :-1] + vals[..., 1:]) * dgrid  # (..., m-1)
    zero = torch.zeros_like(trap[..., :1])
    cumulative = torch.cat([zero, trap.cumsum(dim=-1)], dim=-1)  # (..., m)

    upr_clamped = upr.clamp(0.0, 1.0)
    cell = (
      torch.searchsorted(grid, upr_clamped.contiguous(), right=False) - 1
    ).clamp(0, m - 2)  # (...,)

    cell_exp = cell.unsqueeze(-1)
    v_k = torch.gather(vals, dim=-1, index=cell_exp).squeeze(-1)
    v_k1 = torch.gather(vals, dim=-1, index=cell_exp + 1).squeeze(-1)
    w_k = torch.gather(cumulative, dim=-1, index=cell_exp).squeeze(-1)

    g_k = grid[cell]
    g_k1 = grid[cell + 1]
    dx_cell = g_k1 - g_k
    dx = upr_clamped - g_k
    frac = dx / dx_cell
    partial = (2.0 * v_k + (v_k1 - v_k) * frac) * dx * 0.5

    # Where upr <= grid[0] the integral is exactly 0 (dx == 0 already, so
    # ``partial`` is 0 and ``w_k`` is 0 — the formula is self-correcting).
    return w_k + partial

  # --------------------------------------------------------------------- #
  # Public eval API                                                        #
  # --------------------------------------------------------------------- #

  def interpolate(self, u: Tensor) -> Tensor:
    """Bilinear interpolation of ``values`` at ``u``.

    Args:
      u: shape ``(n, 2)``, each row a query point in ``[0, 1]^2``.

    Returns:
      Tensor of shape ``(n,)`` with the interpolated densities.
    """
    if u.ndim != 2 or u.shape[1] != 2:
      raise ValueError(f"u must have shape (n, 2); got {tuple(u.shape)}")
    u = u.clamp(0.0, 1.0)
    grid = self.grid_points
    i = self._cell_index(u[:, 0])  # (n,)
    j = self._cell_index(u[:, 1])  # (n,)

    z11 = self.values[i, j]
    z12 = self.values[i, j + 1]
    z21 = self.values[i + 1, j]
    z22 = self.values[i + 1, j + 1]

    x1 = grid[i]
    x2 = grid[i + 1]
    y1 = grid[j]
    y2 = grid[j + 1]
    x = u[:, 0]
    y = u[:, 1]

    x2x = x2 - x
    y2y = y2 - y
    xx1 = x - x1
    yy1 = y - y1
    denom = (x2 - x1) * (y2 - y1)
    return (
      z11 * x2x * y2y + z21 * xx1 * y2y + z12 * x2x * yy1 + z22 * xx1 * yy1
    ) / denom

  def integrate_1d(self, u: Tensor, cond_var: int) -> Tensor:
    """Conditional integral along one axis.

    ``cond_var=1`` returns ``H1(u1, u2) = int_0^{u2} c(u1, s) ds / int_0^1 c(u1, s) ds``;
    ``cond_var=2`` returns the symmetric quantity. Output is clamped to
    ``[1e-10, 1-1e-10]``.
    """
    if cond_var not in (1, 2):
      raise ValueError(f"cond_var must be 1 or 2; got {cond_var}")
    u = u.clamp(0.0, 1.0)
    grid = self.grid_points

    if cond_var == 1:
      u_fixed = u[:, 0]
      u_free = u[:, 1]
    else:
      u_fixed = u[:, 1]
      u_free = u[:, 0]

    cell = self._cell_index(u_fixed)
    g_lo = grid[cell]
    g_hi = grid[cell + 1]
    t = ((u_fixed - g_lo) / (g_hi - g_lo)).unsqueeze(-1)  # (n, 1)

    if cond_var == 1:
      v_lo = self.values[cell, :]
      v_hi = self.values[cell + 1, :]
    else:
      v_lo = self.values.index_select(1, cell).t()
      v_hi = self.values.index_select(1, cell + 1).t()

    strip = ((1.0 - t) * v_lo + t * v_hi).clamp_min(_STRIP_FLOOR)  # (n, m)

    numer = self._int_on_grid(u_free, strip)
    denom = self._int_on_grid(torch.ones_like(u_free), strip)
    return (numer / denom).clamp(_TRIM_LO, _TRIM_HI)

  def integrate_2d(self, u: Tensor) -> Tensor:
    """Bivariate CDF: ``int_0^{u1} int_0^{u2} c(s, t) dt ds``.

    Implementation: trapezoidally integrate each row ``values[k, :]`` up to
    ``u2`` to get an ``(n, m)`` strip, then trapezoidally integrate that
    strip along ``grid_points`` up to ``u1``. Clamped to ``[1e-10, 1-1e-10]``.
    """
    u = u.clamp(0.0, 1.0)
    n = u.shape[0]
    m = self.grid_points.shape[0]

    u1 = u[:, 0]
    u2 = u[:, 1]

    # Inner pass: for each query i and grid row k, integrate values[k, :]
    # up to u2[i]. Broadcasting: upr is (n, m), vals is (1, m, m) expanded
    # to (n, m, m). Output (n, m).
    upr_inner = u2.unsqueeze(-1).expand(n, m)
    vals_inner = self.values.unsqueeze(0).expand(n, m, m)
    strip = self._int_on_grid(upr_inner, vals_inner)  # (n, m)

    # Outer pass: integrate strip[i, :] up to u1[i].
    return self._int_on_grid(u1, strip).clamp(_TRIM_LO, _TRIM_HI)

  # --------------------------------------------------------------------- #
  # Cached evaluation grids (optional, see TorchBicop(cache_integrals=…)). #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def build_caches(self) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    """Precompute ``cdf`` / ``hfunc1`` / ``hfunc2`` / ``hinv1`` / ``hinv2``
    at every grid node.

    Returns five ``(m, m)`` tensors; ``TorchBicop`` stores them as buffers
    and bilinearly interpolates them when ``cache_integrals=True``. The
    inverse-h-function caches are built by inverting the just-computed
    h-function caches via a single batched bisection over the full
    ``m^2`` grid (each call is then one lookup instead of 35 iterations
    of bisection at eval time).
    """
    from ._util import solve_itp

    m = self.grid_points.shape[0]
    gi, gj = torch.meshgrid(self.grid_points, self.grid_points, indexing="ij")
    grid_pairs = torch.stack(
      [gi.reshape(-1), gj.reshape(-1)], dim=-1
    )  # (m^2, 2)
    cdf_vals = self.integrate_2d(grid_pairs).reshape(m, m)
    h1_vals = self.integrate_1d(grid_pairs, cond_var=1).reshape(m, m)
    h2_vals = self.integrate_1d(grid_pairs, cond_var=2).reshape(m, m)

    # hinv1: for each grid point (g_i, g_j), solve hfunc1((g_i, u2)) = g_j
    # for u2. Use the just-built h1_vals cache so the bisection's inner
    # function is a fast bilinear interp.
    a = grid_pairs[:, 0:1]
    p = grid_pairs[:, 1:2]

    def _fun_hinv1(u2: Tensor) -> Tensor:
      u_e = torch.cat([a, u2], dim=-1)
      return (
        self.interp_at(h1_vals, u_e).clamp(_TRIM_LO, _TRIM_HI).unsqueeze(-1) - p
      )

    hinv1_vals = (
      solve_itp(_fun_hinv1, torch.zeros_like(p), torch.ones_like(p))
      .squeeze(-1)
      .clamp(0.0, 1.0)
      .reshape(m, m)
    )

    # hinv2: symmetric — solve hfunc2((u1, g_j)) = g_i for u1.
    b = grid_pairs[:, 1:2]
    p = grid_pairs[:, 0:1]

    def _fun_hinv2(u1: Tensor) -> Tensor:
      u_e = torch.cat([u1, b], dim=-1)
      return (
        self.interp_at(h2_vals, u_e).clamp(_TRIM_LO, _TRIM_HI).unsqueeze(-1) - p
      )

    hinv2_vals = (
      solve_itp(_fun_hinv2, torch.zeros_like(p), torch.ones_like(p))
      .squeeze(-1)
      .clamp(0.0, 1.0)
      .reshape(m, m)
    )

    return cdf_vals, h1_vals, h2_vals, hinv1_vals, hinv2_vals

  @torch.no_grad()
  def interp_at(self, cache: Tensor, u: Tensor) -> Tensor:
    """Bilinearly interpolate a precomputed ``(m, m)`` cache at points ``u``.

    Convenience for the cached-integrals path; uses the same bilinear formula
    as :meth:`interpolate` but on an arbitrary user-provided grid.
    """
    if cache.shape != self.values.shape:
      raise ValueError("cache must have the same shape as values")
    if u.ndim != 2 or u.shape[1] != 2:
      raise ValueError(f"u must have shape (n, 2); got {tuple(u.shape)}")
    u = u.clamp(0.0, 1.0)
    grid = self.grid_points
    i = self._cell_index(u[:, 0])
    j = self._cell_index(u[:, 1])

    z11 = cache[i, j]
    z12 = cache[i, j + 1]
    z21 = cache[i + 1, j]
    z22 = cache[i + 1, j + 1]
    x1 = grid[i]
    x2 = grid[i + 1]
    y1 = grid[j]
    y2 = grid[j + 1]
    x = u[:, 0]
    y = u[:, 1]
    x2x = x2 - x
    y2y = y2 - y
    xx1 = x - x1
    yy1 = y - y1
    denom = (x2 - x1) * (y2 - y1)
    return (
      z11 * x2x * y2y + z21 * xx1 * y2y + z12 * x2x * yy1 + z22 * xx1 * yy1
    ) / denom
