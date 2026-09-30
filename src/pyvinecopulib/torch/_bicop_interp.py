"""PyTorch port of :class:`vinecopulib::tools_interpolation::InterpolationGrid`.

Stores a density on a tensor-product grid in [0, 1]^2 and provides bilinear
interpolation plus trapezoidal integration along one or both axes, with the
same partial-cell handling as the C++ implementation.
"""

from __future__ import annotations

import math

import torch
from torch import Tensor
from torch.utils.weak import WeakTensorKeyDictionary

from ..core._trim import trim

# The trapezoidal-integration / bilinear-interpolation kernels live in
# ``_batched`` (they are shape-polymorphic over leading batch dims). The
# scalar methods below are thin ``N=1`` wrappers so there is a single
# source of truth for the numerics shared with the batched vine cascade.
from ._placement import TENSOR_NS
from ._vinecop_batched import (
  _MIN_MASS,
  _batched_cell_index,
  _hfunc_from_cells,
  _locate,
  cond_interval_mass_batched,
  int_on_grid_batched,
  integrate_1d_batched,
  integrate_2d_batched,
  interpolate_batched,
  inverse_integrate_1d_batched,
  mass_tables,
  rect_mass_batched,
  trap_weights,
)

#: Each grid's mass tables, keyed weakly on its ``values`` tensor together with
#: that tensor's version counter; see ``InterpolationGrid2D._mass``.
_MASS_CACHE: WeakTensorKeyDictionary = WeakTensorKeyDictionary()

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
    # Two points are the minimum a cell needs, and the cell widths divide every
    # interpolation and every prefix table: a single point leaves no cell and
    # produced NaN, while a repeated or decreasing point divides by zero or by a
    # negative width and produced a plausible-looking wrong number.
    if grid_points.shape[0] < 2:
      raise ValueError(
        "grid_points must contain at least two points, "
        f"got {grid_points.shape[0]}"
      )
    if not bool((grid_points[1:] > grid_points[:-1]).all()):
      raise ValueError("grid_points must be strictly increasing")
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
    self.register_buffer("trap_weights", trap_weights(grid_points))
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
    device: torch.types.Device = None,
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
    device: torch.types.Device = None,
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

    Three details are required rather than incidental, and match the
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
    ``cond_var=2`` returns the symmetric quantity. Output is clamped
    strictly inside ``[0, 1]``. Thin ``N=1`` wrapper over
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
    """Bivariate CDF: the grid's mass ``M`` over ``[0, u1] x [0, u2]``, rescaled.

    ``M(u1, u2) f(u1) g(u2) M(1, 1)``, with ``f(x) = x / M(x, 1)`` and
    ``g(y) = y / M(1, y)``, so that both margins are exactly uniform and the
    function is symmetric in its arguments, as ``InterpolationGrid`` defines
    it. Clamped strictly inside ``[0, 1]``. Thin ``N=1`` wrapper over
    :func:`._batched.integrate_2d_batched`.
    """
    return integrate_2d_batched(
      self.grid_points,
      self.values.unsqueeze(0),
      u.unsqueeze(0),
      self._is_linear,
    ).squeeze(0)

  def inverse_integrate_1d(
    self, u: Tensor, cond_var: int, cum: Tensor | None = None
  ) -> Tensor:
    """Closed-form inverse of :meth:`integrate_1d` in its free argument.

    Port of the C++ ``InterpolationGrid::inverse_integrate_1d``
    (vinecopulib#691). The conditional density along the free axis is the
    knot vector linearly interpolated at the conditioning value, so the
    conditional cdf is piecewise quadratic and inverts cell by cell.

    ``cond_var=1`` solves ``H1(u1, x) = p`` for ``x`` given
    ``u = [u1, p]``; ``cond_var=2`` solves ``H2(x, u2) = p`` given
    ``u = [p, u2]``. Rows with NaN inputs return NaN, mirroring the C++
    ``binaryExpr_or_nan`` wrapper.

    Args:
      u: shape ``(n, 2)``; see above for the column convention.
      cond_var: 1 or 2, the conditioning variable.
      cum: the prefix-integral table matching ``cond_var``, shape
        ``(m, m)``, or ``None`` to quadrature the conditional cumulative
        from the knots instead.

    Returns:
      Tensor of shape ``(n,)`` with the conditional quantiles in
      ``[0, 1]``.
    """
    if u.ndim != 2 or u.shape[1] != 2:
      raise ValueError(f"u must have shape (n, 2); got {tuple(u.shape)}")
    # A batch of one, through the kernel the stacked path uses -- the same
    # reason `interpolate` and `hfunc_cached` delegate.
    return inverse_integrate_1d_batched(
      self.grid_points,
      self.values.unsqueeze(0),
      u.unsqueeze(0),
      cond_var,
      self._is_linear,
      None if cum is None else cum.unsqueeze(0),
    ).squeeze(0)

  # --------------------------------------------------------------------- #
  # Cached evaluation grids (optional, see TorchTllBicop(cache_integrals=…)). #
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
    ``TorchTllBicop._tables`` rebuilds them when ``values`` starts tracking grad
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

    The result carries the same rescaling of both margins and the same
    domain clamp as :meth:`integrate_2d`, so it is that function's value and
    not a different definition of it.

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
    # the same expression at u1 = 1, where both first-argument partials vanish,
    # and at u2 = 1, where both second-argument ones do
    last = torch.full_like(jc, m - 1)
    m_y = p[last, jc] + al * sx[last, jc] + be * sx[last, jc + 1]
    m_x = p[ic, last] + s_partial(sy[ic, last], sy[ic + 1, last])
    fx = u1 / m_x.clamp_min(_MIN_MASS)
    gy = u2 / m_y.clamp_min(_MIN_MASS)
    inside = (u1 > 0) & (u2 > 0)
    cdf = torch.where(inside, p[-1, -1] * ((fx * gy) * out), 0.0)
    return trim(cdf, TENSOR_NS)

  def _mass(self) -> Tensor:
    """This grid's mass tables, rebuilt only when the grid has changed.

    Held outside the module, keyed weakly on the grid tensor and checked
    against its version counter, so an in-place edit, a replacement or a device
    move each rebuild them, and they reach neither a pickle nor
    ``state_dict``. Rebuilt in the graph, uncached, whenever a gradient in the
    grid is being taken, as ``TorchTllBicop._tables`` is.

    Returns
    -------
    Tensor, shape (1, m * m, 7), dtype float
        From :func:`mass_tables`.
    """
    values = self.values
    if torch.is_grad_enabled() and values.requires_grad:
      return mass_tables(self.grid_points, values.unsqueeze(0))
    hit: tuple[int, Tensor] | None = _MASS_CACHE.get(values)
    if hit is None or hit[0] != values._version:
      hit = (
        values._version,
        mass_tables(self.grid_points, values.unsqueeze(0)),
      )
      _MASS_CACHE[values] = hit
    return hit[1]

  def cond_interval_mass(
    self, u_cond: Tensor, lo: Tensor, hi: Tensor, cond_var: int
  ) -> Tensor:
    """Probability that the free argument falls in ``(lo, hi]``, given the other.

    A batch of one through :func:`cond_interval_mass_batched`, which the
    stacked vine cascade calls too, so the two routes cannot drift.

    Parameters
    ----------
    u_cond : Tensor, shape (n,), dtype float
        The argument held fixed.
    lo, hi : Tensor, shape (n,), dtype float
        Bounds in the free argument, in either order.
    cond_var : int
        1 or 2, the argument held fixed, as for :meth:`integrate_1d`.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Conditional probabilities.
    """
    return cond_interval_mass_batched(
      self.grid_points,
      self._mass(),
      u_cond.unsqueeze(0),
      lo.unsqueeze(0),
      hi.unsqueeze(0),
      cond_var,
      self._is_linear,
    ).squeeze(0)

  def rect_mass(self, a1: Tensor, b1: Tensor, a2: Tensor, b2: Tensor) -> Tensor:
    """The exact probability of ``(a1, b1] x (a2, b2]``, without the cancellation.

    A batch of one through :func:`rect_mass_batched`, where the arrangement
    that avoids the four-corner cancellation is written down. On a
    ``1.2e-4``-wide rectangle it errs by ``2.8e-15`` against exact truth,
    where the difference errs by ``1.8e-8``.

    Parameters
    ----------
    a1, b1, a2, b2 : Tensor, shape (n,), dtype float
        Rectangle bounds per query, in either order. An empty rectangle gives
        zero.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Rectangle probabilities.
    """
    return rect_mass_batched(
      self.grid_points,
      self._mass(),
      a1.unsqueeze(0),
      b1.unsqueeze(0),
      a2.unsqueeze(0),
      b2.unsqueeze(0),
      self._is_linear,
    ).squeeze(0)

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
        ``TorchTllBicop`` passes the one that matches.

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
