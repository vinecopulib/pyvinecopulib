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
    is_linear: bool = False,
    norm_tol: float | None = None,
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
    # When ``is_linear`` is True the grid is assumed to be ``linspace(0, 1, m)``
    # (with the endpoint clamp above leaving it unchanged), so cell-finding is
    # O(1) — ``floor(u * (m - 1))`` — instead of an O(log m) ``searchsorted``.
    self._is_linear = bool(is_linear)
    self.normalize_margins(norm_times, tol=norm_tol)

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
  def normalize_margins(self, times: int, tol: float | None = None) -> None:
    """Renormalize ``values`` so both margins integrate to 1.

    Same algorithm as the C++ ``normalize_margins``: alternating row /
    column trapezoidal-integral divides, repeated up to ``times`` rounds.

    Args:
      times: maximum number of normalization rounds. ``times=3`` matches
        the C++ TLL pipeline byte-for-byte.
      tol: optional convergence tolerance. When provided, iteration stops
        as soon as the largest absolute deviation of the row+col
        integrals from 1 drops below ``tol``. Default ``None`` preserves
        the C++ "fixed-budget" semantics — useful for parity-sensitive
        callers. Setting e.g. ``tol=1e-9`` together with a generous
        ``times=50`` reproduces IPFP-to-convergence.
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
      if tol is not None:
        err = max(
          (row_int - 1.0).abs().max().item(),
          (col_int - 1.0).abs().max().item(),
        )
        if err < tol:
          break

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
    renormalising by the full-strip outer integral so C(1, u2) = u2
    holds exactly (post-vinecopulib#667 C++ behaviour). Clamped to
    ``[1e-10, 1-1e-10]``. Thin ``N=1`` wrapper over
    :func:`._batched.integrate_2d_batched`.
    """
    return integrate_2d_batched(
      self.grid_points,
      self.values.unsqueeze(0),
      u.unsqueeze(0),
      self._is_linear,
    ).squeeze(0)

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
    return interp_at_batched(
      self.grid_points, cache.unsqueeze(0), u.unsqueeze(0), self._is_linear
    ).squeeze(0)
