"""PyTorch port of :class:`vinecopulib::tools_interpolation::InterpolationGrid`.

Stores a density on a tensor-product grid in [0, 1]^2 and provides bilinear
interpolation plus trapezoidal integration along one or both axes, with the
same partial-cell handling as the C++ implementation.
"""

from __future__ import annotations

import math

import torch
from torch import Tensor

# The trapezoidal-integration / bilinear-interpolation kernels live in
# ``_batched`` (they are shape-polymorphic over leading batch dims). The
# scalar methods below are thin ``N=1`` wrappers so there is a single
# source of truth for the numerics shared with the batched vine cascade.
from ._batched import (
  _batched_cell_index,
  int_on_grid_batched,
  integrate_1d_batched,
  integrate_2d_batched,
  interp_at_batched,
  interpolate_batched,
)

_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10

_SQRT_2 = math.sqrt(2.0)
# u_lo = Phi(-3.25); the lower end of the "effective" support that the C++
# normal grid evaluates the TLL KDE over. Reused by ``make_kde_eval_points``
# below so the linear grid sees the same z-range.
_NORMAL_GRID_U_LO: float = 0.5 * (1.0 + math.erf(-3.25 / _SQRT_2))
_NORMAL_GRID_Z_LIMIT: float = 3.25

GRID_TYPES = ("normal", "linear")


def _trap_weights(grid_points: Tensor) -> Tensor:
  """Trapezoid weights of ``grid_points``, summing to 1 on ``[0, 1]``.

  ``_trap_weights(g) @ v`` is the exact integral of the piecewise-linear
  function through ``(g, v)``, because the telescoping sum collapses to
  ``g[-1] - g[0]``.

  Parameters
  ----------
  grid_points : Tensor, shape (m,), dtype float
      A strictly increasing grid whose endpoints are 0 and 1.

  Returns
  -------
  Tensor, shape (m,), dtype float
      The weights.
  """
  m = grid_points.shape[0]
  if m < 2:
    return torch.zeros(0, dtype=grid_points.dtype, device=grid_points.device)
  w = torch.empty_like(grid_points)
  w[0] = (grid_points[1] - grid_points[0]) / 2.0
  w[1:-1] = (grid_points[2:] - grid_points[:-2]) / 2.0
  w[-1] = (grid_points[-1] - grid_points[-2]) / 2.0
  return w


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
  trap_weights: Tensor

  def __init__(
    self,
    grid_points: Tensor,
    values: Tensor,
    norm_maxiter: int = 25,
    is_linear: bool = False,
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
    # Trapezoid weights of `grid_points`, so `trap_weights @ v` integrates the
    # piecewise-linear function through `(grid_points, v)` over [0, 1]. The grid
    # is immutable after construction, so this is built once; the normalization
    # and the h-function denominator both read it instead of walking the cells.
    self.register_buffer("trap_weights", _trap_weights(grid_points))
    # When ``is_linear`` is True the grid is assumed to be ``linspace(0, 1, m)``
    # (with the endpoint clamp above leaving it unchanged), so cell-finding is
    # O(1) — ``floor(u * (m - 1))`` — instead of an O(log m) ``searchsorted``.
    self._is_linear = bool(is_linear)
    self.normalize_margins(norm_maxiter)

  # --------------------------------------------------------------------- #
  # Grid construction (factories shared by all callers)                    #
  # --------------------------------------------------------------------- #

  @staticmethod
  def make_grid_points(
    grid_type: str, m: int, dtype: torch.dtype = torch.float64
  ) -> Tensor:
    r"""Builds the storage grid for a kernel-style bicop on ``[0, 1]^2``.

    Mirrors ``KernelBicop::make_normal_grid`` in the C++ library.
    Endpoints are forced to exactly ``0`` / ``1`` (same as the
    `InterpolationGrid2D` constructor) so callers never need to
    special-case the boundary.

    The ``"normal"`` grid uses :math:`u_i = \Phi(z_i)` with
    :math:`z_i` uniformly spaced on :math:`[-3.25, 3.25]` and
    :math:`\Phi` the standard-normal CDF — the natural domain of
    the TLL density. The ``"linear"`` grid uses
    :math:`u_i = i / (m - 1)`; it trades boundary distortion for
    O(1) cell-finding.

    Parameters
    ----------
    grid_type : {"normal", "linear"}
        ``"normal"`` (Phi-spaced) or ``"linear"`` (uniform on
        ``[0, 1]``).
    m : int
        Number of grid points per axis.
    dtype : torch.dtype, default=torch.float64
        Floating-point dtype.

    Returns
    -------
    Tensor, shape (m,), dtype float
        Grid points on ``[0, 1]``.
    """
    if grid_type not in GRID_TYPES:
      raise ValueError(
        f"grid_type must be one of {GRID_TYPES}; got {grid_type!r}"
      )
    if grid_type == "linear":
      return torch.linspace(0.0, 1.0, m, dtype=dtype)
    z = torch.linspace(
      -_NORMAL_GRID_Z_LIMIT, _NORMAL_GRID_Z_LIMIT, m, dtype=dtype
    )
    grid = 0.5 * (1.0 + torch.erf(z / _SQRT_2))
    grid[0] = 0.0
    grid[-1] = 1.0
    return grid

  @staticmethod
  def make_kde_eval_points(
    grid_type: str, m: int, dtype: torch.dtype = torch.float64
  ) -> Tensor:
    """U-space points used by the TLL KDE evaluator.

    For ``"normal"``: identical to :meth:`make_grid_points` — the KDE
    evaluates at the same un-forced ``Phi(linspace(-3.25, 3.25))`` grid
    that C++ uses, and the storage-grid endpoints get force-clamped to
    0 / 1 inside the :class:`InterpolationGrid2D` constructor (the density
    values are NOT recomputed, matching C++ exactly).

    For ``"linear"``: ``linspace(0, 1, m)`` clamped to
    ``[Phi(-3.25), 1 - Phi(-3.25)]`` so ``qnorm`` stays finite at the
    endpoints — the same trick as the normal grid, only the storage
    coordinate system is uniform on ``[0, 1]``.
    """
    if grid_type not in GRID_TYPES:
      raise ValueError(
        f"grid_type must be one of {GRID_TYPES}; got {grid_type!r}"
      )
    if grid_type == "normal":
      z = torch.linspace(
        -_NORMAL_GRID_Z_LIMIT, _NORMAL_GRID_Z_LIMIT, m, dtype=dtype
      )
      return 0.5 * (1.0 + torch.erf(z / _SQRT_2))
    return torch.linspace(0.0, 1.0, m, dtype=dtype).clamp(
      _NORMAL_GRID_U_LO, 1.0 - _NORMAL_GRID_U_LO
    )

  # --------------------------------------------------------------------- #
  # Margin normalization (port of InterpolationGrid::normalize_margins).  #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def normalize_margins(self, max_iter: int) -> None:
    """Renormalize ``values`` so both margins integrate to 1.

    Port of the C++ ``normalize_margins``. Each pass is the **elementwise
    geometric mean of the two ways of rescaling the grid** — rows then columns,
    and columns then rows. Both are rank-one rescalings of the same values, so
    the mean factorizes into a single rank-one scaling and is applied as one
    fused multiply. Averaging them is what leaves the two margins equally close
    to uniform and makes the pass commute with transposition exactly, so a grid
    and its flipped counterpart normalize to flipped counterparts whether or not
    the iteration has converged.

    Three details are load-bearing rather than incidental, and match the
    reference: the residual is measured on the margins *before* the scaling, so
    an already-normalized grid costs one margin computation and no scaling; the
    scaling is one fused multiply rather than two successive rescalings, which
    would round the two orders differently and lose the equivariance; and the
    transpose is materialized rather than left as a view, so both margins are
    the same reduction and transposing the grid swaps them bit for bit.

    Args:
      max_iter: maximum number of rescaling passes; ``0`` leaves the values
        untouched. Rescaling also stops as soon as both margins integrate to 1
        within ``1e-10``.
    """
    m = self.grid_points.shape[0]
    if max_iter < 1 or m < 2:
      return
    tol, min_mass = 1e-10, 1e-20
    w = self.trap_weights
    for _ in range(max_iter):
      vt = self.values.t().contiguous()
      r = (self.values @ w).clamp_min(min_mass)
      c = (vt @ w).clamp_min(min_mass)
      err = torch.maximum((r - 1.0).abs().max(), (c - 1.0).abs().max())
      if bool(err < tol):
        break
      r2 = (self.values @ (w / c)).clamp_min(min_mass)
      c2 = (vt @ (w / r)).clamp_min(min_mass)
      sr = (r * r2).sqrt().reciprocal()
      sc = (c * c2).sqrt().reciprocal()
      self.values.mul_(sr.unsqueeze(-1)).mul_(sc)

  @torch.no_grad()
  def flip(self) -> None:
    """Transpose ``values`` in place (mirror of ``KernelBicop::flip``)."""
    self.values = self.values.t().contiguous()

  # --------------------------------------------------------------------- #
  # Helpers                                                                #
  # --------------------------------------------------------------------- #

  def _cell_index(self, u: Tensor) -> Tensor:
    """Cell index for each value of ``u`` (shape preserved), clamped to ``[0, m-2]``."""
    return _batched_cell_index(self.grid_points, u, self._is_linear)

  def _int_on_grid(self, upr: Tensor, vals: Tensor) -> Tensor:
    """Vectorized trapezoidal integral of `(grid_points, vals)` from 0 to ``upr``.

    Thin wrapper over :func:`._batched.int_on_grid_batched` (the single
    source of truth). ``upr`` has shape ``(*B,)`` and ``vals`` shape
    ``(*B, m)``; returns shape ``(*B,)``.
    """
    return int_on_grid_batched(self.grid_points, upr, vals, self._is_linear)

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
    return interpolate_batched(
      self.grid_points,
      self.values.unsqueeze(0),
      u.unsqueeze(0),
      self._is_linear,
    ).squeeze(0)

  def integrate_1d(self, u: Tensor, cond_var: int) -> Tensor:
    """Conditional integral along one axis.

    ``cond_var=1`` returns ``H1(u1, u2) = int_0^{u2} c(u1, s) ds / int_0^1 c(u1, s) ds``;
    ``cond_var=2`` returns the symmetric quantity. Output is clamped to
    ``[1e-10, 1-1e-10]``. Thin ``N=1`` wrapper over
    :func:`._batched.integrate_1d_batched`.
    """
    return integrate_1d_batched(
      self.grid_points,
      self.values.unsqueeze(0),
      u.unsqueeze(0),
      cond_var,
      self._is_linear,
    ).squeeze(0)

  def integrate_2d(self, u: Tensor) -> Tensor:
    """Bivariate CDF: ``int_0^{u1} int_0^{u2} c(s, t) dt ds``.

    Trapezoidally integrate each grid row up to ``u2`` to get an
    ``(n, m)`` strip, then integrate that strip up to ``u1``,
    renormalizing by the full-strip outer integral so C(1, u2) = u2
    holds exactly (post-vinecopulib#667 C++ behavior). Clamped to
    ``[1e-10, 1-1e-10]``. Thin ``N=1`` wrapper over
    :func:`._batched.integrate_2d_batched`.
    """
    return integrate_2d_batched(
      self.grid_points,
      self.values.unsqueeze(0),
      u.unsqueeze(0),
      self._is_linear,
    ).squeeze(0)

  def inverse_integrate_1d(self, u: Tensor, cond_var: int) -> Tensor:
    """Closed-form inverse of :meth:`integrate_1d` in its free argument.

    Port of the C++ ``InterpolationGrid::inverse_integrate_1d``
    (vinecopulib#691). The conditional density along the free axis is the
    knot vector linearly interpolated at the conditioning value (floored
    at ``1e-4``, exactly like the forward pass), so the conditional cdf is
    piecewise quadratic and inverts cell by cell: cumulative trapezoidal
    masses, ``searchsorted`` for the bracketing cell, then a numerically
    stable quadratic root clamped to the cell.

    ``cond_var=1`` solves ``H1(u1, x) = p`` for ``x`` given
    ``u = [u1, p]``; ``cond_var=2`` solves ``H2(x, u2) = p`` given
    ``u = [p, u2]``. Rows with NaN inputs return NaN, mirroring the C++
    ``binaryExpr_or_nan`` wrapper.

    Args:
      u: shape ``(n, 2)``; see above for the column convention.
      cond_var: 1 or 2, the conditioning variable.

    Returns:
      Tensor of shape ``(n,)`` with the conditional quantiles in
      ``[0, 1]``.
    """
    if u.ndim != 2 or u.shape[1] != 2:
      raise ValueError(f"u must have shape (n, 2); got {tuple(u.shape)}")
    if cond_var == 1:
      cond, p = u[:, 0], u[:, 1]
    else:
      p, cond = u[:, 0], u[:, 1]
    nan_mask = torch.isnan(cond) | torch.isnan(p)
    cond = cond.nan_to_num(0.5).clamp(0.0, 1.0)
    p = p.nan_to_num(0.5).clamp(_TRIM_LO, _TRIM_HI)

    g = self.grid_points
    m = g.shape[0]

    # Knot vector: the density along the free axis, linearly interpolated
    # at the conditioning value (rows of ``values`` for cond_var=1, columns
    # for cond_var=2), guarded only against rounding, as the forward
    # strip is.
    i = self._cell_index(cond)  # (n,)
    x1, x2 = g[i], g[i + 1]
    w = ((cond - x1) / (x2 - x1)).unsqueeze(-1)  # (n, 1)
    vals = self.values if cond_var == 1 else self.values.t()
    knots = ((1.0 - w) * vals[i, :] + w * vals[i + 1, :]).clamp_min(0.0)

    # Cumulative trapezoidal masses of the (unnormalized) conditional cdf.
    dg = (g[1:] - g[:-1]).unsqueeze(0)  # (1, m - 1)
    masses = 0.5 * (knots[:, :-1] + knots[:, 1:]) * dg  # (n, m - 1)
    incl = masses.cumsum(dim=-1)
    target = (p * incl[:, -1]).unsqueeze(-1)  # (n, 1)

    # Bracketing cell: first k with cumulative mass >= target (the trimmed
    # target is strictly below the total, so k <= m - 2 always holds).
    k = torch.searchsorted(incl.contiguous(), target).clamp(0, m - 2)  # (n, 1)

    # Solve target = cum + v_k s + (v_k1 - v_k) / (2 dg_k) s^2 for
    # s in [0, dg_k]. Without the old 1e-4 knot floor `b = v_k` can be
    # exactly zero, so the stable root needs its own branch: a cell
    # carrying no mass is one the cdf is flat across, where every point
    # is a quantile and the left endpoint is the smallest.
    v_k = knots.gather(-1, k).squeeze(-1)
    v_k1 = knots.gather(-1, k + 1).squeeze(-1)
    cum = torch.where(
      (k > 0).squeeze(-1),
      incl.gather(-1, (k - 1).clamp_min(0)).squeeze(-1),
      torch.zeros_like(v_k),
    )
    k = k.squeeze(-1)
    dg_k = g[k + 1] - g[k]
    a = (v_k1 - v_k) / (2.0 * dg_k)
    b = v_k
    c = cum - target.squeeze(-1)
    disc = (b * b - 4.0 * a * c).clamp_min(0.0)
    denom = b + disc.sqrt()
    safe_b = torch.where(b == 0.0, torch.ones_like(b), b)
    safe_d = torch.where(denom == 0.0, torch.ones_like(denom), denom)
    s = torch.where(
      denom <= 0.0,
      torch.zeros_like(denom),
      torch.where(a.abs() < 1e-300, -c / safe_b, 2.0 * (-c) / safe_d),
    )
    s = torch.minimum(s.clamp_min(0.0), dg_k)
    out = g[k] + s
    return torch.where(nan_mask, torch.full_like(out, torch.nan), out)

  # --------------------------------------------------------------------- #
  # Cached evaluation grids (optional, see TorchBicop(cache_integrals=…)). #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def build_caches(self) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    """Precompute ``cdf`` / ``hfunc1`` / ``hfunc2`` / ``hinv1`` / ``hinv2``
    at every grid node.

    Returns five ``(m, m)`` tensors; ``TorchBicop`` stores them as buffers
    and bilinearly interpolates them when ``cache_integrals=True``. The
    inverse-h-function caches evaluate the exact closed-form
    :meth:`inverse_integrate_1d` at the grid nodes (each call is then one
    lookup instead of a root-finding pass at eval time).
    """
    m = self.grid_points.shape[0]
    gi, gj = torch.meshgrid(self.grid_points, self.grid_points, indexing="ij")
    grid_pairs = torch.stack(
      [gi.reshape(-1), gj.reshape(-1)], dim=-1
    )  # (m^2, 2)
    cdf_vals = self.integrate_2d(grid_pairs).reshape(m, m)
    h1_vals = self.integrate_1d(grid_pairs, cond_var=1).reshape(m, m)
    h2_vals = self.integrate_1d(grid_pairs, cond_var=2).reshape(m, m)

    # hinv1: for each grid point (g_i, g_j), solve hfunc1((g_i, u2)) = g_j
    # for u2; hinv2 is the symmetric problem. Both are exact closed-form
    # conditional quantiles of the interpolated density.
    hinv1_vals = (
      self.inverse_integrate_1d(grid_pairs, cond_var=1)
      .clamp(0.0, 1.0)
      .reshape(m, m)
    )
    hinv2_vals = (
      self.inverse_integrate_1d(grid_pairs, cond_var=2)
      .clamp(0.0, 1.0)
      .reshape(m, m)
    )

    return cdf_vals, h1_vals, h2_vals, hinv1_vals, hinv2_vals

  def interp_at(self, cache: Tensor, u: Tensor) -> Tensor:
    """Bilinearly interpolate a precomputed ``(m, m)`` cache at points ``u``.

    Convenience for the cached-integrals path; uses the same bilinear formula
    as :meth:`interpolate` but on an arbitrary user-provided grid.
    """
    if cache.shape != self.values.shape:
      raise ValueError("cache must have the same shape as values")
    if u.ndim != 2 or u.shape[1] != 2:
      raise ValueError(f"u must have shape (n, 2); got {tuple(u.shape)}")
    return interp_at_batched(
      self.grid_points, cache.unsqueeze(0), u.unsqueeze(0), self._is_linear
    ).squeeze(0)
