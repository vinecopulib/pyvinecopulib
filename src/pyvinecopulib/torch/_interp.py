"""PyTorch port of :class:`vinecopulib::tools_interpolation::InterpolationGrid`.

Stores a density on a tensor-product grid in [0, 1]^2 and provides bilinear
interpolation plus trapezoidal integration along one or both axes, with the
same partial-cell handling as the C++ implementation.
"""

from __future__ import annotations

import math
from typing import Optional

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
  interpolate_batched,
  _MIN_MASS,
  _hfunc_from_cells,
  _locate,
  _trim,
)


#: Above this grid side the passes below stay on the device: the launch
#: overhead they save no longer outweighs moving the grid to the host.
_HOST_NORMALIZE_MAX_M: int = 256

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
  _dgrid: Tensor

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
    # A density is nonnegative, and the exact prefix tables rely on it: every
    # integral they serve is a sum of nonnegative terms, so no cancellation can
    # amplify a rounding error. A negative node would break that silently.
    if bool((values < 0.0).any()):
      raise ValueError("values must be nonnegative; it is a density grid")

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
    self.register_buffer("_dgrid", grid_points[1:] - grid_points[:-1])
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
    grid_type: str,
    m: int,
    dtype: torch.dtype = torch.float64,
    device: Optional[torch.device] = None,
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
      return torch.linspace(0.0, 1.0, m, dtype=dtype, device=device)
    z = torch.linspace(
      -_NORMAL_GRID_Z_LIMIT,
      _NORMAL_GRID_Z_LIMIT,
      m,
      dtype=dtype,
      device=device,
    )
    grid = 0.5 * (1.0 + torch.erf(z / _SQRT_2))
    grid[0] = 0.0
    grid[-1] = 1.0
    return grid

  @staticmethod
  def make_kde_eval_points(
    grid_type: str,
    m: int,
    dtype: torch.dtype = torch.float64,
    device: Optional[torch.device] = None,
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
        -_NORMAL_GRID_Z_LIMIT,
        _NORMAL_GRID_Z_LIMIT,
        m,
        dtype=dtype,
        device=device,
      )
      return 0.5 * (1.0 + torch.erf(z / _SQRT_2))
    return torch.linspace(0.0, 1.0, m, dtype=dtype, device=device).clamp(
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
    # A storage grid is `m x m` with `m` in the tens, and this runs up to
    # `max_iter` passes of a dozen reductions over it. On an accelerator that
    # is a few hundred kernel launches to move a few thousand numbers, which
    # costs far more than the arithmetic; below `_HOST_NORMALIZE_MAX_M` the
    # round trip to the host is cheaper than the launches it saves.
    values, w = self.values, self.trap_weights
    on_host = values.device.type != "cpu" and m <= _HOST_NORMALIZE_MAX_M
    if on_host:
      values, w = values.cpu(), w.cpu()

    tol, min_mass = 1e-10, 1e-20
    for _ in range(max_iter):
      vt = values.t().contiguous()
      r = (values @ w).clamp_min(min_mass)
      c = (vt @ w).clamp_min(min_mass)
      err = torch.maximum((r - 1.0).abs().max(), (c - 1.0).abs().max())
      if bool(err < tol):
        break
      r2 = (values @ (w / c)).clamp_min(min_mass)
      c2 = (vt @ (w / r)).clamp_min(min_mass)
      sr = (r * r2).sqrt().reciprocal()
      sc = (c * c2).sqrt().reciprocal()
      values.mul_(sr.unsqueeze(-1)).mul_(sc)

    if on_host:
      self.values.copy_(values.to(self.values.device))

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
    p = _trim(p.nan_to_num(0.5))

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

  def build_caches(self) -> tuple[Tensor, Tensor, Tensor]:
    """Precompute the three prefix-integral tables the exact cache runs on.

    Returns ``(sy, sx, p)``, each ``(m, m)``:

    * ``sy[i, j] = int_0^{g_j} chat(g_i, t) dt`` -- cumulative along the second
      argument, one row per grid line of the first;
    * ``sx[i, j] = int_0^{g_i} chat(s, g_j) ds`` -- the transpose situation;
    * ``p[i, j]`` -- the double integral over ``[0, g_i] x [0, g_j]``.

    Every one is a cumulative trapezoid, which is **exact** here: ``chat`` is
    bilinear, so along a fixed grid line it is piecewise linear, and its
    integral against ``s`` is again piecewise linear across cells. That is what
    lets :meth:`cdf_cached` and :meth:`hfunc_cached` reconstruct the integral at
    an arbitrary point in O(1) *without approximating it*.

    Three cumulative sums over ``(m, m)``. Built inside the graph, so an exact
    gradient with respect to ``values`` survives; they are registered as
    buffers, so
    ``TorchBicop._tables`` rebuilds them when ``values`` starts tracking grad
    afterwards.

    Returns
    -------
    tuple of Tensor
        The three ``(m, m)`` tables.
    """
    inc = 0.5 * (self.values[:, :-1] + self.values[:, 1:]) * self._dgrid
    sy = torch.cat([torch.zeros_like(inc[:, :1]), inc.cumsum(dim=1)], dim=1)
    incx = (
      0.5
      * (self.values[:-1, :] + self.values[1:, :])
      * self._dgrid.unsqueeze(-1)
    )
    sx = torch.cat([torch.zeros_like(incx[:1, :]), incx.cumsum(dim=0)], dim=0)
    incp = 0.5 * (sy[:-1, :] + sy[1:, :]) * self._dgrid.unsqueeze(-1)
    p = torch.cat([torch.zeros_like(incp[:1, :]), incp.cumsum(dim=0)], dim=0)
    return sy, sx, p

  def cdf_cached(self, u: Tensor, sy: Tensor, sx: Tensor, p: Tensor) -> Tensor:
    """The exact distribution function at ``u``, in O(1) per point.

    The integral over ``[0, u1] x [0, u2]`` splits into the whole cells below
    and left of the query, the strip partial in the first argument, the strip
    partial in the second, and the corner cell partial in both. Each of the four
    is a closed form in the tables, because the ``u2``-partial of a grid line is
    a fixed linear combination of two columns of ``values`` -- so integrating it
    over the first argument reads ``sx`` rather than needing its own table.

    The result carries the same ``* u2 / total`` renormalization and
    ``[1e-10, 1-1e-10]`` clamp as :meth:`integrate_2d`, so it is that function's
    value and not a different definition of it.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Query points.
    sy, sx, p : Tensor, shape (m, m), dtype float
        The tables from :meth:`build_caches`.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Distribution values.
    """
    g = self.grid_points
    m = g.shape[0]
    uu = u.clamp(0.0, 1.0)
    u1, u2 = uu[:, 0], uu[:, 1]
    ic = self._cell_index(u1)
    jc = self._cell_index(u2)
    dx1 = u1 - g[ic]
    f1 = dx1 / (g[ic + 1] - g[ic])
    dx2 = u2 - g[jc]
    f2 = dx2 / (g[jc + 1] - g[jc])
    # the u2-partial of a grid line, as alpha * column jc + beta * column jc+1
    al = dx2 / 2.0 * (2.0 - f2)
    be = dx2 / 2.0 * f2

    def s_partial(v0: Tensor, v1: Tensor) -> Tensor:
      """Partial cell of a piecewise-linear function of the first argument."""
      return (2.0 * v0 + (v1 - v0) * f1) * dx1 / 2.0

    out = p[ic, jc]
    out = out + s_partial(sy[ic, jc], sy[ic + 1, jc])
    out = out + al * sx[ic, jc] + be * sx[ic, jc + 1]
    out = out + al * s_partial(self.values[ic, jc], self.values[ic + 1, jc])
    out = out + be * s_partial(
      self.values[ic, jc + 1], self.values[ic + 1, jc + 1]
    )
    # the same expression at u1 = 1, where both first-argument partials vanish
    last = torch.full_like(jc, m - 1)
    total = p[last, jc] + al * sx[last, jc] + be * sx[last, jc + 1]
    return _trim(out * u2 / total.clamp_min(_MIN_MASS))

  def _interval_weights(self, lo: Tensor, hi: Tensor) -> Tensor:
    """Nonnegative quadrature weights for ``int_lo^hi`` on the grid.

    Returns ``w`` with ``w @ v == int_lo^hi f(t) dt`` for the piecewise-linear
    ``f`` through ``(grid_points, v)``, exactly. Every entry is nonnegative,
    because the hat functions are, which is the property :meth:`rect_mass`
    needs: a difference of two cumulative integrals would be exact too, but
    would carry the larger one's rounding error into a result of order
    ``hi - lo``.

    Parameters
    ----------
    lo, hi : Tensor, shape (n,), dtype float
        Interval endpoints, clamped to ``[0, 1]``; ``hi < lo`` gives zero.

    Returns
    -------
    Tensor, shape (n, m), dtype float
        The weights.
    """
    g = self.grid_points
    m = g.shape[0]
    a = lo.clamp(0.0, 1.0)
    b = hi.clamp(0.0, 1.0).clamp_min(a)
    ka, kb = self._cell_index(a), self._cell_index(b)

    def cell_pair(k: Tensor, s0: Tensor, d: Tensor) -> tuple[Tensor, Tensor]:
      """Weights on nodes ``k`` / ``k + 1`` for the sub-cell ``[s0, s0 + d]``.

      Parameterized by the *width* rather than by the upper end, so a narrow
      interval never forms it as a difference of two numbers of order one --
      which would put the ``1 / w`` amplification back into the weights.
      """
      h = g[k + 1] - g[k]
      q = 0.5 * d * (2.0 * s0 + d)
      return h * (d - q), h * q

    # Whole cells strictly between the two partial ones: the trapezoid rule is
    # exact for a piecewise-linear integrand, so each contributes h / 2 to both
    # of its nodes.
    idx = torch.arange(m - 1, device=g.device)
    inner = (idx[None, :] > ka[:, None]) & (idx[None, :] < kb[:, None])
    half = torch.where(inner, 0.5 * self._dgrid.expand(a.shape[0], m - 1), 0.0)
    zc = torch.zeros_like(half[:, :1])
    w = torch.cat([zc, half], dim=1) + torch.cat([half, zc], dim=1)

    zero = torch.zeros_like(a)
    same = ka == kb
    ha, hb = g[ka + 1] - g[ka], g[kb + 1] - g[kb]
    s0a = (a - g[ka]) / ha
    a0, a1 = cell_pair(ka, s0a, torch.where(same, (b - a) / ha, 1.0 - s0a))
    b0, b1 = cell_pair(kb, zero, (b - g[kb]) / hb)
    b0 = torch.where(same, zero, b0)
    b1 = torch.where(same, zero, b1)
    for col, val in ((ka, a0), (ka + 1, a1), (kb, b0), (kb + 1, b1)):
      w = w.scatter_add(1, col.unsqueeze(1), val.unsqueeze(1))
    return w

  def _raw_mass(self, a1: Tensor, b1: Tensor, a2: Tensor, b2: Tensor) -> Tensor:
    """The mass of the interpolant over ``[a1, b1] x [a2, b2]``, exactly.

    Nonnegative weights against a nonnegative grid, so every term of the sum is
    nonnegative and the result carries no cancellation whatever the rectangle's
    width. This is the density's own mass; :meth:`rect_mass` wraps it in the
    renormalization the distribution function applies.
    """
    wx = self._interval_weights(a1, b1)
    wy = self._interval_weights(a2, b2)
    return (wx * (wy @ self.values.t())).sum(dim=1)

  def rect_mass(self, a1: Tensor, b1: Tensor, a2: Tensor, b2: Tensor) -> Tensor:
    """The exact probability of ``(a1, b1] x (a2, b2]``, without the cancellation.

    The value the four-corner difference of :meth:`cdf_cached` defines, arranged
    so that almost none of it cancels. Differencing those four values turns an
    absolute error ``eps`` into ``~4 eps / (w1 w2)`` in the rectangle's widths;
    this route amplifies by ``1 / w2`` alone, one power instead of two. On a
    ``1.2e-4``-wide rectangle that is ``2.9e-12`` against the difference's
    ``8.7e-9``.

    ``cdf`` renormalizes each line of the grid by its own total, so the
    probability is not simply the mass. Writing ``lam(y) = y / M(1, y)`` for that
    factor, ``M`` for the mass and ``R`` for the rectangle's own mass, the
    four-corner difference is ``lam(b2) * R + (lam(b2) - lam(a2)) * S``, where
    ``S`` is the mass of ``(a1, b1] x (0, a2]``. ``R`` and ``S`` are sums of
    nonnegative terms and do not cancel at all; only the ``lam`` difference does,
    and it multiplies a term of order ``w1`` rather than one of order one, which
    is where the second power goes.

    Parameters
    ----------
    a1, b1, a2, b2 : Tensor, shape (n,), dtype float
        Rectangle bounds per query. An empty or inverted interval gives zero.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Rectangle probabilities.
    """
    zero = torch.zeros_like(a1)
    one = torch.ones_like(a1)
    lam_b = b2 / self._raw_mass(zero, one, zero, b2).clamp_min(_MIN_MASS)
    lam_a = a2 / self._raw_mass(zero, one, zero, a2).clamp_min(_MIN_MASS)
    return lam_b * self._raw_mass(a1, b1, a2, b2) + (
      lam_b - lam_a
    ) * self._raw_mass(a1, b1, zero, a2)

  def hfunc_cached(self, u: Tensor, cond_var: int, sy: Tensor) -> Tensor:
    """The exact conditional distribution function at ``u``, in O(1) per point.

    ``chat`` at a fixed first argument is the linear interpolation of the two
    bracketing grid lines, so both the partial and the total integral along the
    free argument are that same interpolation of the corresponding entries of
    ``sy`` -- no quadrature at evaluation time.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Query points, ``[u_cond, u_free]`` for ``cond_var=1`` and the other way
        round for ``cond_var=2``, matching :meth:`integrate_1d`.
    cond_var : int
        1 or 2, the argument held fixed.
    sy : Tensor, shape (m, m), dtype float
        The first table from :meth:`build_caches` for ``cond_var=1``; its
        transpose-situation twin is built from the transposed values, so
        ``TorchBicop`` passes the one that matches.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Conditional distribution values.
    """
    uu = u.clamp(0.0, 1.0)
    cond = uu[:, 0] if cond_var == 1 else uu[:, 1]
    free = uu[:, 1] if cond_var == 1 else uu[:, 0]
    ic, w, _ = _locate(self.grid_points, cond, self._is_linear)
    jc, frac, dx = _locate(self.grid_points, free, self._is_linear)
    vals = self.values if cond_var == 1 else self.values.t()
    # A batch of one, through the kernel the stacked path uses. The two are
    # pinned to agree to the last bit, and the surest way to keep them there
    # is for there to be only one of them.
    return _hfunc_from_cells(
      vals.unsqueeze(0).contiguous(),
      sy.unsqueeze(0).contiguous(),
      ic.unsqueeze(0),
      w.unsqueeze(0),
      jc.unsqueeze(0),
      frac.unsqueeze(0),
      dx.unsqueeze(0),
    ).squeeze(0)
