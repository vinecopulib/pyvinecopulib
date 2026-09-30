"""Batched pair-copula primitives for ``TorchVinecop``.

Stacks the per-pair ``InterpolationGrid2D`` state at one tree level into
``(N, m, m)`` tensors so the whole level fires one batched bilinear interp
(or one batched trapezoidal integration) instead of N separate calls.

Exposes :func:`interpolate_batched`, :func:`int_on_grid_batched`,
:func:`integrate_1d_batched`, :func:`integrate_2d_batched`,
:func:`inverse_integrate_1d_batched`, :func:`rect_mass_batched` and
:func:`cond_interval_mass_batched` — the ``(N, m, m)`` analogs of the
unbatched operations in :mod:`._interp`, the last two being what a discrete
slot reads an atom's probability off — plus :class:`BatchedTreeLevel`,
:class:`BatchedWave` and :class:`BatchedVine`, which stage one tree level
(resp. one wave of the inverse cascade, resp. an entire vine) of stacked
grids and wire-up tensors.

Intentionally side-by-side with :mod:`._interp` rather than rewriting it:
the legacy / lazy paths stay untouched so any regression is bisectable
to this file.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import torch
from torch import Tensor

from ..core._trim import trim
from ..core.bicop_base import DELTA_MIN
from ..core.extend import NotBatchable
from ..pyvinecopulib_ext import RVineStructure
from ._placement import TENSOR_NS

if TYPE_CHECKING:
  from .vinecop import TorchVinecop

#: Guard on a conditional total mass, so a zero-mass grid line cannot 0/0.
#: Floor for a *mass* the cascades divide by -- a renormalizing integral or a
#: rectangle's probability. Not `trim`'s own lower bound: that is a **domain**
#: bound on copula arguments (1e-10), and using it here clamped a denominator
#: ten orders of magnitude early. Both sites reachable only where the numerator
#: vanishes too, so nothing moved; the constant was simply the wrong one.
_MIN_MASS: float = 1e-20

#: Peak a mixed level's evaluation aims to stay under, the budget the TLL
#: fit's grid blocking uses too.
_DISCRETE_MEM_BUDGET_BYTES: int = 256 * 1024 * 1024

#: Values a mixed level's evaluation holds live per (pair, row) at its peak:
#: measured 8.2 KB at float64 with the density, nearly all of it the stencil
#: gathers of the stacked probability calls, and rounded up.
_DISCRETE_VALUES_PER_QUERY: int = 1100


# --------------------------------------------------------------------------- #
# Batched bilinear interpolation                                               #
# --------------------------------------------------------------------------- #


def _batched_cell_index(
  grid_points: Tensor, u: Tensor, is_linear: bool = False
) -> Tensor:
  """Per-element cell index, clamped to ``[0, m-2]``. Preserves ``u``'s shape.

  When ``is_linear`` is True the grid is assumed to be ``linspace(0, 1, m)``
  and the index is computed as ``floor(u * (m - 1))`` — O(1) vs the O(log m)
  ``searchsorted`` of the default path.
  """
  m = grid_points.shape[0]
  if is_linear:
    return (u * (m - 1)).long().clamp(0, m - 2)
  return (
    torch.searchsorted(grid_points, u.contiguous(), right=False) - 1
  ).clamp(0, m - 2)


def _locate(
  grid_points: Tensor, x: Tensor, is_linear: bool
) -> tuple[Tensor, Tensor, Tensor]:
  """Locate ``x`` on the grid: ``(cell, weight, offset)``.

  ``weight`` is the position within the cell and ``offset`` the distance from
  its lower edge. Both the density lookup and the two h-functions need exactly
  this triple for each argument, so it is computed once per argument and
  shared rather than three times over.

  ``x`` is already clamped to the grid and ``cell`` to ``[0, m-2]``, so the
  ratio needs no clamp of its own.
  """
  cell = _batched_cell_index(grid_points, x, is_linear)
  lo = grid_points[cell]
  off = x - lo
  return cell, off / (grid_points[cell + 1] - lo), off


def _bilinear(
  values: Tensor, i: Tensor, j: Tensor, wx: Tensor, wy: Tensor
) -> Tensor:
  """Bilinear value of ``values`` at cell ``(i, j)``, offsets ``(wx, wy)``.

  Three ``lerp`` calls rather than the four-term weighted sum: identical
  arithmetic to within rounding, a third of the kernel launches, and this
  cascade is bound by launch count rather than by flops.
  """
  n_batch = values.shape[0]
  n_pts = i.shape[-1]
  rows = torch.arange(n_batch, device=values.device).unsqueeze(-1)
  rows = rows.expand(n_batch, n_pts)
  lo = torch.lerp(values[rows, i, j], values[rows, i + 1, j], wx)
  hi = torch.lerp(values[rows, i, j + 1], values[rows, i + 1, j + 1], wx)
  return torch.lerp(lo, hi, wy)


def interpolate_batched(
  grid_points: Tensor, values: Tensor, u: Tensor, is_linear: bool = False
) -> Tensor:
  """Batched bilinear interpolation.

  Args:
    grid_points: shape ``(m,)``, shared across all pairs.
    values: shape ``(N, m, m)``, one grid per pair.
    u: shape ``(N, n, 2)``, queries per pair (in the unrotated frame —
      the rotation must already be applied to ``values``).

  Returns:
    Tensor of shape ``(N, n)``.
  """
  if values.ndim != 3:
    raise ValueError(f"values must be 3D (N, m, m); got {tuple(values.shape)}")
  if u.ndim != 3 or u.shape[-1] != 2:
    raise ValueError(f"u must be (N, n, 2); got {tuple(u.shape)}")
  if u.shape[0] != values.shape[0]:
    raise ValueError(
      f"u.shape[0]={u.shape[0]} != values.shape[0]={values.shape[0]}"
    )
  # The *closed* interval, not `trim`'s open one: this clamp keeps a grid
  # lookup inside the grid, which is a different question from the domain step
  # `trim` applies to a copula argument. Every `.clamp(0.0, 1.0)` in this
  # module is one of these, and the kernels apply `trim` on the way out.
  u = u.clamp(0.0, 1.0)
  _N, _n, _ = u.shape

  i, wx, _ = _locate(grid_points, u[..., 0], is_linear)
  j, wy, _ = _locate(grid_points, u[..., 1], is_linear)
  return _bilinear(values, i, j, wx, wy)


def int_on_grid_batched(
  grid_points: Tensor, upr: Tensor, vals: Tensor, is_linear: bool = False
) -> Tensor:
  """Vectorized trapezoidal integral of ``(grid_points, vals)`` from 0 to ``upr``.

  The function is shape-polymorphic in the same way the per-pair
  :meth:`InterpolationGrid2D._int_on_grid` is: ``upr`` of shape ``(*B,)``
  and ``vals`` of shape ``(*B, m)`` produce an output of shape ``(*B,)``.
  ``*B`` can carry leading batch dimensions (``(N, n)`` typically).
  """
  m = grid_points.shape[0]
  dgrid = grid_points[1:] - grid_points[:-1]  # (m-1,)

  trap = 0.5 * (vals[..., :-1] + vals[..., 1:]) * dgrid
  zero = torch.zeros_like(trap[..., :1])
  cumulative = torch.cat([zero, trap.cumsum(dim=-1)], dim=-1)

  upr_clamped = upr.clamp(0.0, 1.0)
  if is_linear:
    cell = (upr_clamped * (m - 1)).long().clamp(0, m - 2)
  else:
    cell = (
      torch.searchsorted(grid_points, upr_clamped.contiguous(), right=False) - 1
    ).clamp(0, m - 2)

  cell_exp = cell.unsqueeze(-1)
  v_k = torch.gather(vals, dim=-1, index=cell_exp).squeeze(-1)
  v_k1 = torch.gather(vals, dim=-1, index=cell_exp + 1).squeeze(-1)
  w_k = torch.gather(cumulative, dim=-1, index=cell_exp).squeeze(-1)

  g_k = grid_points[cell]
  g_k1 = grid_points[cell + 1]
  dx_cell = g_k1 - g_k
  dx = upr_clamped - g_k
  frac = dx / dx_cell
  partial = (2.0 * v_k + (v_k1 - v_k) * frac) * dx * 0.5
  return w_k + partial


def integrate_1d_batched(
  grid_points: Tensor,
  values: Tensor,
  u: Tensor,
  cond_var: int,
  is_linear: bool = False,
) -> Tensor:
  """Batched conditional 1-D integral.

  Args:
    grid_points: shape ``(m,)``.
    values: shape ``(N, m, m)``, precomputed pdf grids.
    u: shape ``(N, n, 2)``, queries (unrotated frame; the precomputation absorbs the
      rotation).
    cond_var: scalar in ``{1, 2}`` — kept consistent across all pairs in the
      batch because the precomputation puts every pair in the same "natural" frame.
      ``cond_var=1`` returns the h-function conditioning on ``u[..., 0]``;
      ``cond_var=2`` conditions on ``u[..., 1]``.

  Returns:
    Tensor of shape ``(N, n)``, strictly inside ``[0, 1]``.
  """
  if cond_var not in (1, 2):
    raise ValueError(f"cond_var must be 1 or 2; got {cond_var}")
  u = u.clamp(0.0, 1.0)
  _N, _n, _ = u.shape
  m = grid_points.shape[0]

  if cond_var == 1:
    u_fixed = u[..., 0]  # (N, n)
    u_free = u[..., 1]
    fixed_axis = 1  # rows of (m, m) — i.e. dim 1 of (N, m, m)
  else:
    u_fixed = u[..., 1]
    u_free = u[..., 0]
    fixed_axis = 2  # columns

  cell = _batched_cell_index(grid_points, u_fixed, is_linear)  # (N, n)
  g_lo = grid_points[cell]
  g_hi = grid_points[cell + 1]
  t = ((u_fixed - g_lo) / (g_hi - g_lo)).unsqueeze(-1)  # (N, n, 1)

  # Gather two strips of shape (N, n, m): for each (k, l) take
  #   values[k, cell[k, l], :]      (cond_var=1)
  #   values[k, :, cell[k, l]]      (cond_var=2)
  # gather along ``fixed_axis`` after expanding the cell index along the
  # other (free) axis to size m.
  strip = _cond_strip(values, cell, t, fixed_axis, m)

  number = int_on_grid_batched(grid_points, u_free, strip, is_linear)  # (N, n)
  denom = int_on_grid_batched(
    grid_points, torch.ones_like(u_free), strip, is_linear
  )  # (N, n)
  # Without the floor a grid line can carry no mass at all, so the
  # division needs its own guard.
  return trim(number / denom.clamp_min(_MIN_MASS), TENSOR_NS)


def _cond_strip(
  values: Tensor,
  cell: Tensor,
  t: Tensor,
  fixed_axis: int,
  m: int,
  floor: bool = True,
) -> Tensor:
  """The conditional density along the free axis, as ``(N, n, m)``.

  The grid line at the conditioning value, i.e. the blend of the two
  bracketing lines of ``values``. A bilinear interpolation of a nonnegative
  grid is nonnegative, so the guard only absorbs rounding: it used to floor at
  ``1e-4``, which made the h-functions not the conditional cdf of the density
  the same object reported -- by up to 7.5e-5 on a strongly dependent fit,
  where two thirds of the grid can sit below that floor.
  """
  n_batch, n = cell.shape
  if fixed_axis == 1:
    idx_lo = cell.unsqueeze(-1).expand(n_batch, n, m)
    idx_hi = (cell + 1).unsqueeze(-1).expand(n_batch, n, m)
    v_lo = values.gather(dim=1, index=idx_lo)
    v_hi = values.gather(dim=1, index=idx_hi)
  else:
    idx_lo = cell.unsqueeze(1).expand(n_batch, m, n)
    idx_hi = (cell + 1).unsqueeze(1).expand(n_batch, m, n)
    v_lo = values.gather(dim=2, index=idx_lo).transpose(1, 2)
    v_hi = values.gather(dim=2, index=idx_hi).transpose(1, 2)
  strip = torch.lerp(v_lo, v_hi, t)
  return strip.clamp_min(0.0) if floor else strip


def inverse_integrate_1d_batched(
  grid_points: Tensor,
  values: Tensor,
  u: Tensor,
  cond_var: int,
  is_linear: bool = False,
  cum: Tensor | None = None,
) -> Tensor:
  """Batched closed-form inverse of :func:`integrate_1d_batched`.

  The stacked twin of :meth:`InterpolationGrid2D.inverse_integrate_1d`. The
  conditional density along the free axis is the knot vector interpolated at
  the conditioning value, so the conditional cdf is piecewise quadratic and
  inverts cell by cell: cumulative trapezoidal masses, ``searchsorted`` for
  the bracketing cell, then a numerically stable quadratic root clamped to it.

  Unlike the forward direction there is no O(1) lookup to be had -- locating
  the bracketing cell needs the conditional cumulative along the whole free
  axis -- so the ``(N, n, m)`` strip is unavoidable and is what bounds this
  function's memory. The *cumulative* need not be quadratured, though: given
  the prefix table, blending its two bracketing lines is the same quantity as
  integrating the blended knots, and agrees to 2e-16 while costing a gather
  instead of a trapezoid and a scan.

  Parameters
  ----------
  grid_points : Tensor, shape (m,), dtype float
      The shared grid.
  values : Tensor, shape (N, m, m), dtype float
      Density grids, one per pair.
  u : Tensor, shape (N, n, 2), dtype float
      ``[u_cond, p]`` for ``cond_var=1``, ``[p, u_cond]`` for 2.
  cond_var : int
      1 or 2, the conditioning argument.
  is_linear : bool, default=False
      Whether the grid is uniform, enabling O(1) cell lookup.
  cum : Tensor, shape (N, m, m), or None, optional
      The prefix-integral table matching ``cond_var`` (``sy`` for 1, ``sx``
      for 2). Supplied, the conditional cumulative is read off it; otherwise
      it is quadratured from the knots.

  Returns
  -------
  Tensor, shape (N, n), dtype float
      Conditional quantiles in ``[0, 1]``; NaN wherever an input was NaN.
  """
  if values.ndim != 3:
    raise ValueError(f"values must be 3D (N, m, m); got {tuple(values.shape)}")
  if u.ndim != 3 or u.shape[-1] != 2:
    raise ValueError(f"u must be (N, n, 2); got {tuple(u.shape)}")
  m = grid_points.shape[0]
  cond, p = (u[..., 0], u[..., 1]) if cond_var == 1 else (u[..., 1], u[..., 0])
  nan_mask = torch.isnan(cond) | torch.isnan(p)
  cond = cond.nan_to_num(0.5).clamp(0.0, 1.0)
  p = trim(p.nan_to_num(0.5), TENSOR_NS)

  fixed_axis = 1 if cond_var == 1 else 2
  cell = _batched_cell_index(grid_points, cond, is_linear)
  g_lo = grid_points[cell]
  t = ((cond - g_lo) / (grid_points[cell + 1] - g_lo)).unsqueeze(-1)
  knots = _cond_strip(values, cell, t, fixed_axis, m)

  if cum is not None:
    incl = _cond_strip(cum, cell, t, fixed_axis, m, floor=False)[..., 1:]
  else:
    dg = grid_points[1:] - grid_points[:-1]
    incl = (0.5 * (knots[..., :-1] + knots[..., 1:]) * dg).cumsum(dim=-1)
  target = (p * incl[..., -1]).unsqueeze(-1)
  # The trimmed target is strictly below the total, so k <= m - 2 holds.
  k = torch.searchsorted(incl.contiguous(), target).clamp(0, m - 2)

  v_k = knots.gather(-1, k).squeeze(-1)
  v_k1 = knots.gather(-1, k + 1).squeeze(-1)
  below = torch.where(
    (k > 0).squeeze(-1),
    incl.gather(-1, (k - 1).clamp_min(0)).squeeze(-1),
    torch.zeros_like(v_k),
  )
  k = k.squeeze(-1)
  dg_k = grid_points[k + 1] - grid_points[k]

  # Solve target = below + v_k s + (v_k1 - v_k) / (2 dg_k) s^2 for s in
  # [0, dg_k]. `b = v_k` can be exactly zero, so the stable root needs its own
  # branch: a cell carrying no mass is one the cdf is flat across, where every
  # point is a quantile and the left endpoint is the smallest.
  a = (v_k1 - v_k) / (2.0 * dg_k)
  b = v_k
  c = below - target.squeeze(-1)
  denom = b + (b * b - 4.0 * a * c).clamp_min(0.0).sqrt()
  # Exact: these two guard a division, so what matters is whether the
  # denominator is the value that cannot be divided by, not whether it is
  # near it. A tolerance here would substitute 1.0 for a small real divisor.
  safe_b = torch.where(b == 0.0, torch.ones_like(b), b)  # noqa: RUF069
  safe_d = torch.where(
    denom == 0.0,  # noqa: RUF069
    torch.ones_like(denom),
    denom,
  )
  s = torch.where(
    denom <= 0.0,
    torch.zeros_like(denom),
    torch.where(a.abs() < 1e-300, -c / safe_b, 2.0 * (-c) / safe_d),
  )
  out = grid_points[k] + torch.minimum(s.clamp_min(0.0), dg_k)
  return torch.where(nan_mask, torch.full_like(out, torch.nan), out)


def integrate_2d_batched(
  grid_points: Tensor, values: Tensor, u: Tensor, is_linear: bool = False
) -> Tensor:
  """Batched bivariate CDF (trapezoidal-trapezoidal).

  Same shape contract as :func:`integrate_1d_batched`: ``values: (N, m, m)``,
  ``u: (N, n, 2)``, returns ``(N, n)`` clamped strictly inside ``[0, 1]``.
  The result is renormalized by the full-strip outer integral so
  C(1, u2) = u2 holds exactly — matches the post-vinecopulib#667 C++
  behavior and stays in parity with the unbatched
  :meth:`InterpolationGrid2D.integrate_2d`.
  """
  u = u.clamp(0.0, 1.0)
  N, n, _ = u.shape
  m = grid_points.shape[0]

  u1 = u[..., 0]  # (N, n)
  u2 = u[..., 1]

  # Inner pass: for each (k, l) and each row r, integrate values[k, r, :] up
  # to u2[k, l]. Build upr_inner of shape (N, n, m) (broadcast u2 across the
  # row axis) and vals_inner of shape (N, n, m, m) (broadcast values across
  # the query axis).
  upr_inner = u2.unsqueeze(-1).expand(N, n, m)
  vals_inner = values.unsqueeze(1).expand(N, n, m, m)
  strip = int_on_grid_batched(
    grid_points, upr_inner, vals_inner, is_linear
  )  # (N, n, m)

  # Outer pass: integrate strip[k, l, :] up to u1[k, l] and renormalize
  # by the full-first-axis integral so C(1, u2) = u2 holds exactly.
  # Guard the degenerate `tmpint1 = 0` case (e.g. cache-building at
  # raw grid endpoints) — true CDF is then 0.
  tmpint = int_on_grid_batched(grid_points, u1, strip, is_linear)
  tmpint1 = int_on_grid_batched(
    grid_points, torch.ones_like(u1), strip, is_linear
  )
  out = torch.where(
    tmpint1 > 0,
    tmpint * u2 / tmpint1.clamp_min(_MIN_MASS),
    torch.zeros_like(tmpint),
  )
  return trim(out, TENSOR_NS)


# --------------------------------------------------------------------------- #
# Batched probabilities of rectangles and conditional intervals                #
# --------------------------------------------------------------------------- #


def trap_weights(grid_points: Tensor) -> Tensor:
  """Trapezoid weights of ``grid_points``, summing to 1 on ``[0, 1]``.

  ``trap_weights(g) @ v`` is the exact integral of the piecewise-linear
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


#: Channels of a mass table (see :func:`mass_tables`): the density grid, then
#: its three prefix-integral tables, each as a value and that value's rounding
#: error.
_V, _SY, _SY_LO, _SX, _SX_LO, _P, _P_LO = range(7)


def _two_sum(a: Tensor, b: Tensor) -> tuple[Tensor, Tensor]:
  """``a + b`` as a rounded sum and its exact rounding error (Knuth)."""
  s = a + b
  bb = s - a
  return s, (a - (s - bb)) + (b - bb)


def _prefix(x: Tensor, dim: int) -> tuple[Tensor, Tensor]:
  """Compensated prefix sums of the nonnegative ``x`` along ``dim``.

  Returns ``hi`` and ``lo`` with a leading zero, so that ``hi + lo`` is the
  running sum to within ``eps**2`` of its size rather than ``eps``. That is
  what lets a *difference* of two entries be accurate relative to itself: the
  range sum of a nonnegative sequence, read off a plain prefix table, carries
  the rounding error of the whole prefix and cancels in a low-mass range.

  The error is recovered without a sequential loop: ``s = hi[k-1] + x[k]`` and
  ``hi[k]`` approximate the same nonnegative sum, so ``s - hi[k]`` is exact by
  Sterbenz's lemma, and the per-step residuals telescope to ``T[k] - hi[k]``.
  """
  hi = x.cumsum(dim)
  n = x.shape[dim]
  zero = torch.zeros_like(hi.narrow(dim, 0, 1))
  prev = torch.cat([zero, hi.narrow(dim, 0, n - 1)], dim)
  s, e = _two_sum(prev, x)
  lo = ((s - hi) + e).cumsum(dim)
  return torch.cat([zero, hi], dim), torch.cat([zero, lo], dim)


def mass_tables(grid_points: Tensor, values: Tensor) -> Tensor:
  """What :func:`rect_mass_batched` and :func:`cond_interval_mass_batched` read.

  Per pair, the density grid and its three prefix-integral tables -- ``sy``
  along the second argument, ``sx`` along the first, and ``p`` over both, the
  tables of ``InterpolationGrid2D.build_caches`` -- summed cell by cell and
  compensated (:func:`_prefix`). A probability is then a few reads of a 4x4
  stencil around its interval's ends instead of a quadrature over every grid
  node, which is what made it ``O(m)`` per query.

  Parameters
  ----------
  grid_points : Tensor, shape (m,), dtype float
      The shared grid.
  values : Tensor, shape (N, m, m), dtype float
      Density grids, one per pair; nonnegative.

  Returns
  -------
  Tensor, shape (N, m * m, 7), dtype float
      The channels ``_V`` .. ``_P_LO``, flattened over the grid row-major.
  """
  n_batch, m, _ = values.shape
  dg = grid_points[1:] - grid_points[:-1]
  # Cell by cell: each line's trapezoids along either argument, and the
  # bilinear cell masses, all nonnegative.
  cy = 0.5 * (values[..., :-1] + values[..., 1:]) * dg
  cx = 0.5 * (values[:, :-1, :] + values[:, 1:, :]) * dg[:, None]
  cell = 0.5 * (cy[:, :-1, :] + cy[:, 1:, :]) * dg[:, None]
  sy, sy_lo = _prefix(cy, 2)
  sx, sx_lo = _prefix(cx, 1)
  rows, rows_lo = _prefix(cell, 2)
  p, p_lo = _prefix(rows, 1)
  zero = torch.zeros_like(rows_lo[:, :1, :])
  p_lo = p_lo + torch.cat([zero, rows_lo.cumsum(1)], 1)
  out = torch.stack([values, sy, sy_lo, sx, sx_lo, p, p_lo], dim=-1)
  return out.reshape(n_batch, m * m, 7)


def _pieces(
  grid_points: Tensor, lo: Tensor, hi: Tensor, is_linear: bool
) -> tuple[Tensor, Tensor, Tensor]:
  """``[lo, hi]`` as two partial cells and the whole cells between them.

  The partial cells are nonnegative weights on the four nodes
  ``(ka, ka + 1, kb, kb + 1)`` bracketing the two ends, the second pair zero
  when both ends share a cell. The whole cells run from node ``ka + 1`` to
  node ``kb``, stencil positions 1 and 2, and exist where ``kb > ka``.

  Parameters
  ----------
  grid_points : Tensor, shape (m,), dtype float
      The shared grid.
  lo, hi : Tensor, shape (*B,), dtype float
      Interval ends, clamped to ``[0, 1]``; ``hi < lo`` is empty.
  is_linear : bool
      Whether the grid is uniform.

  Returns
  -------
  nodes : Tensor, shape (*B, 4), dtype long
      The four stencil nodes.
  weights : Tensor, shape (*B, 4), dtype float
      Their partial-cell weights.
  inner : Tensor, shape (*B,), dtype bool
      Whether there are whole cells between the two partial ones.
  """
  g = grid_points
  a = lo.clamp(0.0, 1.0)
  b = hi.clamp(0.0, 1.0).clamp_min(a)
  ka = _batched_cell_index(g, a, is_linear)
  kb = _batched_cell_index(g, b, is_linear)

  def cell_pair(h: Tensor, s0: Tensor, d: Tensor) -> tuple[Tensor, Tensor]:
    """Weights on a cell's two nodes for its sub-interval ``[s0, s0 + d]``.

    Parameterized by the *width* rather than by the upper end, so a narrow
    interval never forms it as a difference of two numbers of order one --
    which would put the ``1 / w`` amplification back into the weights.
    """
    q = 0.5 * d * (2.0 * s0 + d)
    return h * (d - q), h * q

  same = ka == kb
  ha, hb = g[ka + 1] - g[ka], g[kb + 1] - g[kb]
  s0a = (a - g[ka]) / ha
  a0, a1 = cell_pair(ha, s0a, torch.where(same, (b - a) / ha, 1.0 - s0a))
  b0, b1 = cell_pair(hb, torch.zeros_like(a), (b - g[kb]) / hb)
  b0 = torch.where(same, 0.0, b0)
  b1 = torch.where(same, 0.0, b1)
  nodes = torch.stack([ka, ka + 1, kb, kb + 1], dim=-1)
  return nodes, torch.stack([a0, a1, b0, b1], dim=-1), kb > ka


def _stencil(tables: Tensor, rows: Tensor, cols: Tensor, m: int) -> Tensor:
  """Every channel at every ``(row, col)`` pair of two node sets.

  Parameters
  ----------
  tables : Tensor, shape (N, m * m, 7), dtype float
      From :func:`mass_tables`.
  rows, cols : Tensor, shape (N, n, R) and (N, n, C), dtype long
      Grid nodes along the first and the second argument.
  m : int
      Grid size.

  Returns
  -------
  Tensor, shape (7, N, n, R, C), dtype float
      The gathered channels, channel first.

  Notes
  -----
  The gather reads each node's seven channels together, which is why the
  tables keep them last; the result puts them first, in one copy, so every
  channel the kernels read is contiguous. Left channel-last, each read is a
  stride-7 view, which vectorizes on no cpu kernel.
  """
  n_batch, n, r = rows.shape
  c = cols.shape[-1]
  pos = (rows.unsqueeze(-1) * m + cols.unsqueeze(-2)).reshape(n_batch, -1)
  batch = torch.arange(n_batch, device=tables.device).unsqueeze(-1)
  st = tables[batch, pos].reshape(n_batch, n, r, c, 7)
  return st.movedim(-1, 0).contiguous()


def _whole(t: Tensor, dim: int) -> Tensor:
  """Every channel's difference between stencil nodes 2 and 1 along ``dim``.

  The whole cells of an interval run between those two nodes, so this is
  every table's mass over them at once. A difference of two rounded values is
  itself correctly rounded, so with each ``_lo`` channel's difference folded
  back in -- ``d[c] + d[c + 1]`` -- the result is accurate relative
  to itself: a range sum of nonnegative terms, read without cancellation.
  """
  return t.select(dim, 2) - t.select(dim, 1)


def cond_interval_mass_batched(
  grid_points: Tensor,
  tables: Tensor,
  u_cond: Tensor,
  lo: Tensor,
  hi: Tensor,
  cond_var: int,
  is_linear: bool = False,
) -> Tensor:
  """Probability that the free argument falls in ``(lo, hi]``, given the other.

  A conditional distribution is the grid line at ``u_cond`` over its own
  total, so this is a ratio of nonnegative sums and does not cancel at all.
  Not clamped into the open unit interval, unlike
  :func:`integrate_1d_batched`: an empty interval is exactly ``0`` and the
  whole line exactly ``1``, which is what makes the masses of a partition sum
  to one. Read off the tables as :func:`rect_mass_batched` reads them, and to
  the same accuracy.

  Parameters
  ----------
  grid_points : Tensor, shape (m,), dtype float
      The shared grid.
  tables : Tensor, shape (N, m * m, 7), dtype float
      From :func:`mass_tables`, one per pair.
  u_cond : Tensor, shape (N, n), dtype float
      The argument held fixed.
  lo, hi : Tensor, shape (N, n), dtype float
      Bounds in the free argument, in either order.
  cond_var : int
      1 or 2, the argument held fixed, as for :func:`integrate_1d_batched`.
  is_linear : bool, default=False
      Whether the grid is uniform.

  Returns
  -------
  Tensor, shape (N, n), dtype float
      Conditional probabilities.
  """
  if cond_var not in (1, 2):
    raise ValueError(f"cond_var must be 1 or 2; got {cond_var}")
  m = grid_points.shape[0]
  cell, t, _ = _locate(grid_points, u_cond.clamp(0.0, 1.0), is_linear)
  nodes, w, inner = _pieces(
    grid_points, torch.minimum(lo, hi), torch.maximum(lo, hi), is_linear
  )
  lines = torch.stack([cell, cell + 1], dim=-1)
  free = torch.cat([nodes, torch.full_like(nodes[..., :1], m - 1)], dim=-1)
  # The fixed argument indexes the lines; the free one runs along them, read
  # off `sy` when that is the second argument and `sx` when it is the first.
  if cond_var == 1:
    st, ch = _stencil(tables, lines, free, m), _SY
  else:
    st, ch = _stencil(tables, free, lines, m).transpose(-1, -2), _SX
  # (7, N, n, 2 lines, 5 free nodes): partial cells, then the whole ones.
  line_mass = (w.unsqueeze(-2) * st[_V, ..., :4]).sum(-1)
  d = _whole(st, -1)
  whole = d[ch] + d[ch + 1]
  line_mass = line_mass + torch.where(inner.unsqueeze(-1), whole, 0.0)
  total = st[ch, ..., 4] + st[ch + 1, ..., 4]
  mass = torch.lerp(line_mass[..., 0], line_mass[..., 1], t)
  return mass / torch.lerp(total[..., 0], total[..., 1], t).clamp_min(_MIN_MASS)


def rect_mass_batched(
  grid_points: Tensor,
  tables: Tensor,
  a1: Tensor,
  b1: Tensor,
  a2: Tensor,
  b2: Tensor,
  is_linear: bool = False,
) -> Tensor:
  """The exact probability of ``(a1, b1] x (a2, b2]``, without the cancellation.

  The value the four-corner difference of the cached distribution function
  defines, arranged so that almost none of it cancels. Differencing those four
  values turns an absolute error ``eps`` into ``~4 eps / (w1 w2)`` in the
  rectangle's widths; this route amplifies by ``1 / w2`` alone, one power
  instead of two.

  The distribution function renormalizes each line of the grid by its own
  total, so the probability is not simply the mass. Writing
  ``lam(y) = y / M(1, y)`` for that factor, ``R`` for the rectangle's own mass
  and ``S`` for the mass of ``(a1, b1] x (0, a2]``, the four-corner difference
  is ``lam(b2) R + (lam(b2) - lam(a2)) S``. ``R`` and ``S`` are sums of
  nonnegative terms and do not cancel at all; only the ``lam`` difference
  does, and expanded over the common denominator it cancels against
  ``b2 - a2`` rather than against one -- formed as ``lam(b2) - lam(a2)`` it
  would swamp the first term on a low-mass rectangle.

  Each mass is its partial cells, summed directly, plus its whole cells, read
  as compensated differences of the prefix tables (:func:`mass_tables`), so
  the cost per query does not depend on the grid size. The tables are exact to
  ``~1e-32`` absolute, so a mass is accurate relative to itself down to about
  ``1e-16`` and absolutely below that.

  Parameters
  ----------
  grid_points : Tensor, shape (m,), dtype float
      The shared grid.
  tables : Tensor, shape (N, m * m, 7), dtype float
      From :func:`mass_tables`, one per pair.
  a1, b1, a2, b2 : Tensor, shape (N, n), dtype float
      Rectangle bounds per query, in either order. An empty rectangle gives
      zero.
  is_linear : bool, default=False
      Whether the grid is uniform.

  Returns
  -------
  Tensor, shape (N, n), dtype float
      Rectangle probabilities.
  """
  # Rotating the data can leave a left limit above its own value.
  x0, x1 = torch.minimum(a1, b1), torch.maximum(a1, b1)
  y0, y1 = torch.minimum(a2, b2), torch.maximum(a2, b2)
  g = grid_points
  m = g.shape[0]
  nodes, w, inner = _pieces(
    g, torch.stack([x0, y0]), torch.stack([x1, y1]), is_linear
  )
  (px, py), (wx, wy), (in_x, in_y) = nodes, w, inner
  # The first argument's stencil and the last node, which is where `sx` and
  # `p` hold the masses over all of the first argument.
  rows = torch.cat([px, torch.full_like(px[..., :1], m - 1)], dim=-1)
  st = _stencil(tables, rows, py, m)  # (7, N, n, 5, 4)
  v = st[_V, ..., :4, :]
  # Whole cells in the second argument, per first-argument node, and in the
  # first, per second-argument node.
  dy = _whole(st, -1)  # (7, N, n, 5)
  dx = _whole(st, -2)  # (7, N, n, 4)

  # R, the rectangle's own mass: partial cells in both arguments, partial in
  # one and whole in the other, and whole in both.
  pp = (wx.unsqueeze(-1) * wy.unsqueeze(-2) * v).sum((-1, -2))
  x_part = (wx * (dy[_SY, ..., :4] + dy[_SY_LO, ..., :4])).sum(-1)
  y_whole = dx[_SX] + dx[_SX_LO]
  x_whole = (wy * y_whole).sum(-1)
  # The whole-by-whole block is a four-corner difference of `p`, so its first
  # stage keeps its rounding error exactly and the second cannot cancel
  # against it.
  s1, e1 = _two_sum(st[_P, ..., 2, 1:3], -st[_P, ..., 1, 1:3])
  both = (s1[..., 1] - s1[..., 0]) + (
    (e1[..., 1] - e1[..., 0]) + (dx[_P_LO, ..., 2] - dx[_P_LO, ..., 1])
  )
  r = (
    pp
    + torch.where(in_y, x_part, 0.0)
    + torch.where(in_x, x_whole + torch.where(in_y, both, 0.0), 0.0)
  )

  # `y0` falls in the stencil's first cell, so the mass below it is read at
  # nodes 0 and 1: the prefix at node 0 plus the partial cell up to `y0`,
  # whose two node weights are `cdf_cached`'s `al` / `be`.
  j0 = py[..., 0]
  s = y0.clamp(0.0, 1.0) - g[j0]
  f = s / (g[j0 + 1] - g[j0])
  al, be = s / 2.0 * (2.0 - f), s / 2.0 * f
  below_line = (
    st[_SY, ..., :4, 0]
    + st[_SY_LO, ..., :4, 0]
    + al.unsqueeze(-1) * v[..., 0]
    + be.unsqueeze(-1) * v[..., 1]
  )
  s_part = (wx * below_line).sum(-1)
  # The whole cells in the first argument below `y0`: their prefix mass up to
  # node 0, plus the same partial cell, taken across them.
  p_below = dx[_P, ..., 0] + dx[_P_LO, ..., 0]
  s_mass = s_part + torch.where(
    in_x, p_below + al * y_whole[..., 0] + be * y_whole[..., 1], 0.0
  )

  # The masses over the whole first argument, off its last node.
  sx_last = st[_SX, ..., 4, :] + st[_SX_LO, ..., 4, :]
  m_strip = (wy * sx_last).sum(-1) + torch.where(
    in_y, dy[_P, ..., 4] + dy[_P_LO, ..., 4], 0.0
  )
  p_last = st[_P, ..., 4, 0] + st[_P_LO, ..., 4, 0]
  m_below = (p_last + al * sx_last[..., 0] + be * sx_last[..., 1]).clamp_min(
    _MIN_MASS
  )

  # Formed additively, so that the two terms below are consistent: an
  # independent quadrature of the whole column would not cancel against the
  # strip it is supposed to contain.
  total = (m_below + m_strip).clamp_min(_MIN_MASS)
  dlam = ((y1 - y0) * m_below - y0 * m_strip) / (total * m_below)
  out = y1 * r / total + dlam * s_mass
  empty = ~((x1 > x0) & (y1 > y0))
  return torch.where(empty, torch.zeros_like(out), out)


# --------------------------------------------------------------------------- #
# Per-tree-level stacked state + per-vine container                            #
# --------------------------------------------------------------------------- #


def _hfunc_from_cells(
  values: Tensor,
  sy: Tensor,
  ic: Tensor,
  w: Tensor,
  jc: Tensor,
  frac: Tensor,
  dx: Tensor,
) -> Tensor:
  """Exact conditional distribution function from a cumulative table.

  The batched twin of :meth:`InterpolationGrid2D.hfunc_cached`, taking the
  grid locations already resolved by :func:`_locate`. At a fixed conditioning
  argument the interpolant is the linear blend of the two bracketing grid
  lines, so both the partial and the total integral along the free argument
  are that same blend of the corresponding entries of ``sy`` -- exact, and one
  gather per corner rather than a quadrature.

  Parameters
  ----------
  values : Tensor, shape (N, m, m), dtype float
      Density grids, oriented so that ``values[k, i, :]`` is the line at a
      fixed conditioning argument -- transposed by the caller for ``hfunc2``.
  sy : Tensor, shape (N, m, m), dtype float
      ``sy[k, i, j]`` is that line's integral up to ``grid_points[j]``.
  ic, w : Tensor, shape (N, n)
      Cell and within-cell weight of the conditioning argument.
  jc, frac, dx : Tensor, shape (N, n)
      Cell, within-cell weight and offset of the free argument.

  Returns
  -------
  Tensor, shape (N, n), dtype float
      Conditional distribution values, clamped strictly inside ``[0, 1]``.
  """
  n_batch, m, _ = values.shape
  flat_v = values.reshape(n_batch, m * m)
  flat_s = sy.reshape(n_batch, m * m)

  # Both bracketing grid lines, addressed off one base so the row stride is
  # added once rather than recomputed.
  base_lo = ic * m
  base_hi = base_lo + m
  at_lo, at_hi = base_lo + jc, base_hi + jc
  last = m - 1

  def line(at: Tensor, base: Tensor) -> Tensor:
    """Integral of one grid line from 0 to the free argument."""
    v0 = flat_v.gather(1, at)
    v1 = flat_v.gather(1, at + 1)
    # cum + (2 v0 + (v1 - v0) frac) dx / 2, with the inner blend as one lerp.
    return torch.addcmul(
      flat_s.gather(1, at), v0 + torch.lerp(v0, v1, frac), dx, value=0.5
    )

  num = torch.lerp(line(at_lo, base_lo), line(at_hi, base_hi), w)
  den = torch.lerp(
    flat_s.gather(1, base_lo + last), flat_s.gather(1, base_hi + last), w
  )
  return trim(num / den.clamp_min(_MIN_MASS), TENSOR_NS)


class BatchedTreeLevel(torch.nn.Module):
  """Stacked state for every pair-copula at one tree level of a vine.

  All buffers are registered so ``.to(device)`` / ``.to(dtype)`` move them
  together with the parent :class:`BatchedVine`.

  Grids (per pair):
  - ``values: (N, m, m)`` — pdf grid (rotation-less; TLL pair-copulas in
    pyvinecopulib always have rotation 0).
  - ``grids2, tables2: (2N, m, m) | None`` — the pdf grid and its
    cumulative-trapezoid prefix integrals along argument 2, each stacked over
    the grid and its transpose, present only when every source pair was
    constructed with ``cache_integrals=True``. The two h-functions are the
    same reduction with the arguments swapped, so stacking lets a level
    evaluate both in one call on ``2N`` pairs.

  Wiring (per pair, same across pdf / rosenblatt / inverse cascades):
  - ``col0_src: (N,) long`` — column to read for ``col0`` (= edge index).
  - ``col1_src: (N,) long`` — column to read for ``col1`` (= ``min_array - 1``).
  - ``col1_use_h1: (N,) bool`` — whether ``col1`` reads from ``hfunc1``
    (else ``hfunc2``).
  - ``needs_h1, needs_h2: (N,) bool`` — whether the cascade requires this
    pair's h-function output at the next tree.

  Variable types (per pair, for a vine with discrete variables):
  - ``disc1, disc2: (N,) bool`` — whether the pair's first / second argument
    is discrete, i.e. the slot's :meth:`~VinecopBase.pair_var_types`. The
    left-limit columns are read through the same wiring as the values, and a
    continuous argument's left limit is its own value.

  Indep handling:
  - ``is_indep: (N,) bool`` — true slots get short-circuit overrides
    (pdf=1, hfunc1=col1, hfunc2=col0). Their grids are sentinels.
  """

  # Class-level type hints so the buffers registered in __init__ are
  # statically typed as Tensors instead of nn.Module (cf. ``_sy`` in
  # TorchTllBicop, same pattern).
  values: Tensor
  grids2: Tensor | None
  tables2: Tensor | None
  is_indep: Tensor
  col0_src: Tensor
  col1_src: Tensor
  col1_use_h1: Tensor
  needs_h1: Tensor
  needs_h2: Tensor
  disc1: Tensor
  disc2: Tensor
  idx_d1: Tensor
  idx_d2: Tensor
  idx_dd: Tensor
  idx_dc: Tensor
  idx_cd: Tensor
  mass: Tensor | None

  def __init__(
    self,
    *,
    values: Tensor,
    sy: Tensor | None,
    sy_t: Tensor | None,
    is_indep: Tensor,
    col0_src: Tensor,
    col1_src: Tensor,
    col1_use_h1: Tensor,
    needs_h1: Tensor,
    needs_h2: Tensor,
    disc1: Tensor | None = None,
    disc2: Tensor | None = None,
    grid_points: Tensor | None = None,
    is_linear: bool = False,
  ) -> None:
    super().__init__()
    self.register_buffer("values", values)
    if sy is not None:
      assert sy_t is not None
      # `hfunc2` reads lines of the transposed grid where `hfunc1` reads lines
      # of this one, so keep both materialized and stacked: one call on 2N
      # pairs answers a whole tree level.
      self.register_buffer(
        "grids2", torch.cat([values, values.transpose(1, 2)], 0).contiguous()
      )
      self.register_buffer("tables2", torch.cat([sy, sy_t], 0))
    else:
      self.grids2 = None
      self.tables2 = None
    self.register_buffer("is_indep", is_indep)
    self.register_buffer("col0_src", col0_src)
    self.register_buffer("col1_src", col1_src)
    self.register_buffer("col1_use_h1", col1_use_h1)
    self.register_buffer("needs_h1", needs_h1)
    self.register_buffer("needs_h2", needs_h2)
    if disc1 is None:
      disc1 = torch.zeros_like(is_indep)
    if disc2 is None:
      disc2 = torch.zeros_like(is_indep)
    self.register_buffer("disc1", disc1)
    self.register_buffer("disc2", disc2)
    # The pairs each quotient is taken over, fixed by the structure. Their
    # sizes are kept as Python ints so a group no pair belongs to is skipped
    # before anything is launched for it.
    groups = {
      "idx_d1": disc1,
      "idx_d2": disc2,
      "idx_dd": disc1 & disc2,
      "idx_dc": disc1 & ~disc2,
      "idx_cd": ~disc1 & disc2,
    }
    self._group_sizes: dict[str, int] = {}
    for name, mask in groups.items():
      idx = torch.nonzero(mask).flatten()
      self.register_buffer(name, idx)
      self._group_sizes[name] = int(idx.numel())
    self._is_linear = bool(is_linear)
    # What a discrete slot reads its atoms' probabilities off. Built from the
    # grids once per level rather than per call, and only where a slot needs
    # it: a continuous level never reads a probability.
    if self.has_discrete:
      if grid_points is None:
        raise ValueError("a level with a discrete slot needs grid_points")
      self.register_buffer("mass", mass_tables(grid_points, values))
    else:
      self.mass = None

  @property
  def n_pairs(self) -> int:
    return int(self.values.shape[0])

  @property
  def has_discrete(self) -> bool:
    """Whether any pair at this level has a discrete argument."""
    return bool(self._group_sizes["idx_d1"] or self._group_sizes["idx_d2"])

  def gather_inputs(
    self,
    hfunc1_prev: Tensor,
    hfunc2_prev: Tensor,
    hfunc1_sub: Tensor | None = None,
    hfunc2_sub: Tensor | None = None,
  ) -> Tensor:
    """Build the per-pair input from the previous level's h-function columns.

    Lays the ``col0`` / ``col1`` selection out as a pair of gathers, plus a
    ``torch.where`` on the h1-vs-h2 source flag for ``col1``. Given the
    left-limit scratch as well, the result is the four-column
    ``[u1, u2, u1^-, u2^-]`` a discrete slot reads, gathered through the same
    wiring; a continuous argument's left limit is its own value, so a stale
    entry of the scratch is never read.

    Parameters
    ----------
    hfunc1_prev, hfunc2_prev : Tensor, shape (n, d), dtype float
        The h-function scratch.
    hfunc1_sub, hfunc2_sub : Tensor, shape (n, d), dtype float, or None
        The left-limit scratch, for a vine with discrete variables.

    Returns
    -------
    Tensor, shape (N, n, 2) or (N, n, 4), dtype float
        One input per pair.
    """

    def pick(h1: Tensor, h2: Tensor) -> tuple[Tensor, Tensor]:
      # (n, d) -> (N, n) per column.
      col0 = h2.index_select(dim=1, index=self.col0_src)
      col1 = torch.where(
        self.col1_use_h1[None, :],
        h1.index_select(dim=1, index=self.col1_src),
        h2.index_select(dim=1, index=self.col1_src),
      )
      return col0.t(), col1.t()

    col0, col1 = pick(hfunc1_prev, hfunc2_prev)
    if hfunc2_sub is None:
      return torch.stack([col0, col1], dim=-1)
    assert hfunc1_sub is not None
    sub0, sub1 = pick(hfunc1_sub, hfunc2_sub)
    sub0 = torch.where(self.disc1[:, None], sub0, col0)
    sub1 = torch.where(self.disc2[:, None], sub1, col1)
    return torch.stack([col0, col1, sub0, sub1], dim=-1)

  def _locate_both(self, grid_points: Tensor, u: Tensor) -> tuple[Tensor, ...]:
    """Grid location of both arguments: ``(i, wx, dx, j, wy, dy)``.

    The density and the two h-functions all need exactly these two triples --
    the h-functions with the roles swapped -- so resolving them once is what
    makes :meth:`pdf_h1_h2` cheaper than its three parts.
    """
    uu = u.clamp(0.0, 1.0)
    i, wx, dx = _locate(grid_points, uu[..., 0], self._is_linear)
    j, wy, dy = _locate(grid_points, uu[..., 1], self._is_linear)
    return i, wx, dx, j, wy, dy

  def _pdf_at(self, i: Tensor, j: Tensor, wx: Tensor, wy: Tensor) -> Tensor:
    raw = _bilinear(self.values, i, j, wx, wy).clamp_min(1e-20)
    return torch.where(self.is_indep[:, None], torch.ones_like(raw), raw)

  def _h1_h2_at(
    self,
    gp: Tensor,
    a: Tensor,
    loc_a: tuple[Tensor, ...],
    b: Tensor,
    loc_b: tuple[Tensor, ...],
  ) -> tuple[Tensor, Tensor]:
    """``hfunc1`` at the located query ``a`` and ``hfunc2`` at ``b``.

    ``hfunc1`` conditions on argument 1 and integrates argument 2; ``hfunc2``
    does the reverse against the transposed grid. That is the same kernel with
    the two triples swapped, so with the grids stacked it is one call on
    ``2N`` pairs instead of two on ``N`` -- worth doing in a cascade whose cost
    is the number of calls. The two queries coincide on a continuous edge; a
    discrete edge reads each h-function at a different point.
    """
    i_a, wx_a, _, j_a, wy_a, dy_a = loc_a
    i_b, wx_b, dx_b, j_b, wy_b, _ = loc_b
    if self.grids2 is None or self.tables2 is None:
      raw = torch.cat(
        [
          integrate_1d_batched(gp, self.values, a, 1, self._is_linear),
          integrate_1d_batched(gp, self.values, b, 2, self._is_linear),
        ],
        0,
      )
    else:
      raw = _hfunc_from_cells(
        self.grids2,
        self.tables2,
        torch.cat([i_a, j_b], 0),
        torch.cat([wx_a, wy_b], 0),
        torch.cat([j_a, i_b], 0),
        torch.cat([wy_a, wx_b], 0),
        torch.cat([dy_a, dx_b], 0),
      )
    # An independent pair returns its own free argument -- argument 2 for
    # `hfunc1`, argument 1 for `hfunc2` -- which is the same swap again.
    both = torch.where(
      self.is_indep.repeat(2)[:, None],
      trim(torch.cat([a[..., 1], b[..., 0]], 0), TENSOR_NS),
      raw,
    )
    return both[: self.n_pairs], both[self.n_pairs :]

  def pdf_h1_h2(
    self, grid_points: Tensor, u: Tensor
  ) -> tuple[Tensor, Tensor, Tensor]:
    """``(pdf, hfunc1, hfunc2)`` for one tree level, from one grid search.

    All three read the same two cells of the same grid -- the h-functions with
    the conditioning and free arguments swapped -- so the search, the cell
    weights and the offsets are resolved once and shared. This cascade is
    bound by kernel-launch count, so that sharing is most of the cost.
    """
    loc = self._locate_both(grid_points, u)
    i, wx, _, j, wy, _ = loc
    h1, h2 = self._h1_h2_at(grid_points, u, loc, u, loc)
    return self._pdf_at(i, j, wx, wy), h1, h2

  def h1_h2(self, grid_points: Tensor, u: Tensor) -> tuple[Tensor, Tensor]:
    """``(hfunc1, hfunc2)`` for one tree level; see :meth:`pdf_h1_h2`."""
    loc = self._locate_both(grid_points, u)
    return self._h1_h2_at(grid_points, u, loc, u, loc)

  # --- the mixed-discrete level ----------------------------------------- #
  # The stacked twin of `BicopBase`'s quotients: each output is what the
  # pair's `pdf` / `hfunc1` / `hfunc2` dispatchers return on the four-column
  # input, with the per-row `DELTA_MIN` split taken as a `where` rather than a
  # boolean index, and each quotient evaluated only on the pairs whose types
  # call for it.

  def _rects(
    self,
    gp: Tensor,
    jobs: list[tuple[str, tuple[Tensor, Tensor, Tensor, Tensor]]],
  ) -> list[Tensor]:
    """Rectangle probabilities for several pair groups, in one stacked call.

    Parameters
    ----------
    gp : Tensor, shape (m,), dtype float
        The shared grid.
    jobs : list of tuple
        ``(group, (a1, b1, a2, b2))``: the name of an index buffer and the
        level-wide ``(N, n)`` bounds, read at that group's pairs.

    Returns
    -------
    list of Tensor
        One ``(K, n)`` block per job, ``K`` the size of its group.
    """
    idx = torch.cat([getattr(self, group) for group, _ in jobs])
    a1, b1, a2, b2 = (
      torch.cat(
        [bs[k].index_select(0, getattr(self, group)) for group, bs in jobs]
      )
      for k in range(4)
    )
    assert self.mass is not None
    out = rect_mass_batched(
      gp, self.mass.index_select(0, idx), a1, b1, a2, b2, self._is_linear
    )
    # `TorchTllBicop.rect_prob`'s own answer for the independence copula.
    indep = (b1 - a1).abs() * (b2 - a2).abs()
    out = torch.where(self.is_indep.index_select(0, idx)[:, None], indep, out)
    return list(out.split([self._group_sizes[group] for group, _ in jobs], 0))

  def _interval(
    self,
    gp: Tensor,
    group: str,
    u_cond: Tensor,
    lo: Tensor,
    hi: Tensor,
    cond_var: int,
  ) -> Tensor:
    """Conditional interval probabilities at one group's pairs."""
    idx = getattr(self, group)
    u_cond, lo, hi = (t.index_select(0, idx) for t in (u_cond, lo, hi))
    assert self.mass is not None
    out = cond_interval_mass_batched(
      gp,
      self.mass.index_select(0, idx),
      u_cond,
      lo,
      hi,
      cond_var,
      self._is_linear,
    )
    # `TorchTllBicop.cond_interval_prob`'s answer for the independence copula.
    indep = (hi.clamp(0.0, 1.0) - lo.clamp(0.0, 1.0)).abs()
    return torch.where(self.is_indep.index_select(0, idx)[:, None], indep, out)

  def eval_discrete(
    self, grid_points: Tensor, u: Tensor, *, with_pdf: bool
  ) -> tuple[Tensor | None, Tensor, Tensor, Tensor, Tensor]:
    """``(pdf, hfunc1, hfunc2, hfunc1^-, hfunc2^-)`` on a four-column input.

    ``hfunc1^-`` is ``hfunc1`` with the second argument at its left limit and
    ``hfunc2^-`` is ``hfunc2`` with the first at its left limit -- the two
    values the vine's left-limit scratch propagates. Both are computed for
    every pair; the cascade writes each only where the slot's types make it
    one.

    Parameters
    ----------
    grid_points : Tensor, shape (m,), dtype float
        The shared grid.
    u : Tensor, shape (N, n, 4), dtype float
        ``[u1, u2, u1^-, u2^-]`` per pair, from :meth:`gather_inputs`.
    with_pdf : bool
        Whether to evaluate the density as well; ``None`` in its place
        otherwise.

    Returns
    -------
    tuple of Tensor
        Five ``(N, n)`` tensors, the first ``None`` unless ``with_pdf``.
    """
    # The pair dispatchers clamp their whole input, left limits included, and
    # a quotient the previous level produced is not clamped on its way out.
    u = trim(u, TENSOR_NS)
    if not self.has_discrete:
      u2c = u[..., :2]
      if with_pdf:
        pdf, h1, h2 = self.pdf_h1_h2(grid_points, u2c)
      else:
        pdf, (h1, h2) = None, self.h1_h2(grid_points, u2c)
      return pdf, h1, h2, h1, h2
    # Every output is row-wise, so rows are evaluated in blocks: the stencil
    # gathers hold ~1 KB per (pair, row) each, which on a wide level at a
    # large sample is more than a card holds at once.
    n_pairs, n = int(u.shape[0]), int(u.shape[1])
    per_row = n_pairs * _DISCRETE_VALUES_PER_QUERY * u.element_size()
    block = max(1, _DISCRETE_MEM_BUDGET_BYTES // per_row)
    if n <= block:
      return self._eval_discrete_rows(grid_points, u, with_pdf)
    parts = [
      self._eval_discrete_rows(grid_points, u[:, i : i + block], with_pdf)
      for i in range(0, n, block)
    ]
    pdfs = [q[0] for q in parts]
    pdf = None
    if with_pdf:
      pdf = torch.cat([q for q in pdfs if q is not None], 1)
    h1 = torch.cat([q[1] for q in parts], 1)
    h2 = torch.cat([q[2] for q in parts], 1)
    h1_sub = torch.cat([q[3] for q in parts], 1)
    h2_sub = torch.cat([q[4] for q in parts], 1)
    return pdf, h1, h2, h1_sub, h2_sub

  def _eval_discrete_rows(
    self, grid_points: Tensor, u: Tensor, with_pdf: bool
  ) -> tuple[Tensor | None, Tensor, Tensor, Tensor, Tensor]:
    """:meth:`eval_discrete` on one block of rows, already clamped."""
    gp, lin = grid_points, self._is_linear
    u1, u2, u1m, u2m = u.unbind(-1)
    # A continuous argument's left limit is its own value, so its midpoint is
    # the value itself, exactly: every continuous reading below is therefore
    # also the plain one for an argument with no atom.
    mid1, mid2 = 0.5 * (u1 + u1m), 0.5 * (u2 + u2m)
    at = {
      name: _locate(gp, x.clamp(0.0, 1.0), lin)
      for name, x in (
        ("u1", u1),
        ("u2", u2),
        ("u1m", u1m),
        ("u2m", u2m),
        ("mid1", mid1),
        ("mid2", mid2),
      )
    }

    def loc(x: str, y: str) -> tuple[Tensor, ...]:
      return (*at[x], *at[y])

    def pt(x: Tensor, y: Tensor) -> Tensor:
      return torch.stack([x, y], dim=-1)

    # The continuous readings: `hfunc1` at the atom's midpoint in the first
    # argument, `hfunc2` in the second, each at the value and at the left
    # limit of the other argument.
    h1_a, h2_a = self._h1_h2_at(
      gp, pt(mid1, u2), loc("mid1", "u2"), pt(u1, mid2), loc("u1", "mid2")
    )
    h1_b, h2_b = self._h1_h2_at(
      gp, pt(mid1, u2m), loc("mid1", "u2m"), pt(u1m, mid2), loc("u1m", "mid2")
    )
    delta1, delta2 = (u1 - u1m).abs(), (u2 - u2m).abs()
    zero = torch.zeros_like(u1)
    sizes = self._group_sizes

    jobs: list[tuple[str, tuple[Tensor, Tensor, Tensor, Tensor]]] = []
    if sizes["idx_d1"]:
      jobs.append(("idx_d1", (u1m, u1, zero, u2)))
    if sizes["idx_d2"]:
      jobs.append(("idx_d2", (zero, u1, u2m, u2)))
    if sizes["idx_dd"]:
      jobs.extend(
        [("idx_dd", (u1m, u1, zero, u2m)), ("idx_dd", (zero, u1m, u2m, u2))]
      )
      if with_pdf:
        jobs.append(("idx_dd", (u1m, u1, u2m, u2)))
    masses = iter(self._rects(gp, jobs))

    def replace(base: Tensor, group: str, value: Tensor) -> Tensor:
      return base.index_copy(0, getattr(self, group), value)

    def take(t: Tensor, group: str) -> Tensor:
      return t.index_select(0, getattr(self, group))

    def quotient(num: Tensor, delta: Tensor, fallback: Tensor) -> Tensor:
      # `BicopBase._quotient`: the quotient over a wide-enough atom, else the
      # derivative at its midpoint.
      wide = delta > DELTA_MIN
      safe = torch.where(wide, delta, torch.ones_like(delta))
      return torch.where(wide, num / safe, fallback).abs()

    h1, h2, h1_sub, h2_sub = h1_a, h2_a, h1_b, h2_b
    if sizes["idx_d1"]:
      h1 = replace(
        h1,
        "idx_d1",
        quotient(next(masses), take(delta1, "idx_d1"), take(h1_a, "idx_d1")),
      )
    if sizes["idx_d2"]:
      h2 = replace(
        h2,
        "idx_d2",
        quotient(next(masses), take(delta2, "idx_d2"), take(h2_a, "idx_d2")),
      )
    if sizes["idx_dd"]:
      d1, d2 = take(delta1, "idx_dd"), take(delta2, "idx_dd")
      h1_sub = replace(
        h1_sub, "idx_dd", quotient(next(masses), d1, take(h1_b, "idx_dd"))
      )
      h2_sub = replace(
        h2_sub, "idx_dd", quotient(next(masses), d2, take(h2_b, "idx_dd"))
      )
    if not with_pdf:
      return None, h1, h2, h1_sub, h2_sub

    loc_mid = loc("mid1", "mid2")
    pdf = self._pdf_at(loc_mid[0], loc_mid[3], loc_mid[1], loc_mid[4])
    base = pdf
    for group, cond, lo, hi, cond_var, delta in (
      ("idx_dc", u2, u1m, u1, 2, delta1),
      ("idx_cd", u1, u2m, u2, 1, delta2),
    ):
      if not sizes[group]:
        continue
      num = self._interval(gp, group, cond, lo, hi, cond_var)
      dg = take(delta, group)
      # `BicopBase._pdf_mixed`: strictly wider than `DELTA_MIN` is wide.
      wide = dg > DELTA_MIN
      safe = torch.where(wide, dg, torch.ones_like(dg))
      pdf = replace(
        pdf, group, torch.where(wide, num / safe, take(base, group)).abs()
      )
    if sizes["idx_dd"]:
      rect = next(masses)
      d1, d2 = take(delta1, "idx_dd"), take(delta2, "idx_dd")
      # `BicopBase._pdf_d_d`: the four regimes, which are disjoint.
      both = torch.where(d1 > d2, d1, d2) < DELTA_MIN
      only1 = (d1 < DELTA_MIN) & ~both
      only2 = (d2 < DELTA_MIN) & ~both
      wide = ~(both | only1 | only2)
      one = torch.ones_like(d1)
      val = torch.where(
        wide,
        rect / torch.where(wide, d1 * d2, one),
        torch.where(
          only1,
          take(h1_a - h1_b, "idx_dd") / torch.where(only1, d2, one),
          torch.where(
            only2,
            take(h2_a - h2_b, "idx_dd") / torch.where(only2, d1, one),
            take(base, "idx_dd"),
          ),
        ),
      )
      pdf = replace(pdf, "idx_dd", val.abs())
    return pdf, h1, h2, h1_sub, h2_sub


def inverse_waves(
  s: RVineStructure, d: int, trunc_lvl: int
) -> list[list[tuple[int, int]]]:
  """Group the inverse cascade's ``(var, tree)`` cells into parallel waves.

  ``_inverse_rosenblatt`` walks variables outward and, within each, trees
  inward. Cell ``(var, tree)`` reads ``hinv2[tree + 1, var]`` -- the same
  variable one tree further out -- and, at the same tree, either
  ``hinv2[tree, m - 1]`` or ``hfunc1[tree, m - 1]``, the latter written by
  cell ``(m - 1, tree - 1)``. Both predecessors are fixed by the structure,
  so the dependency graph is static and can be levelled once, here.

  The grouping is *not* the tree level -- each wave holds one cell from
  almost every tree -- and it is not the anti-diagonal either: with
  ``m - 1 == var + 1`` off the diagonal, which is the generic D-vine cell,
  ``(var + 1, tree - 1)`` lands on the same anti-diagonal as ``(var, tree)``.
  Levelling the actual graph is both correct and tighter than any fixed key.

  Parameters
  ----------
  s : RVineStructure
      The vine structure, read for ``min_array`` / ``struct_array``.
  d : int
      Vine dimension.
  trunc_lvl : int
      Truncation level.

  Returns
  -------
  list of list of tuple
      Cells per wave, in execution order. Every cell in a wave is
      independent of the others, so a wave is one stacked call.
  """
  deps: dict[tuple[int, int], set[tuple[int, int]]] = {}
  for var in range(d - 2, -1, -1):
    for tree in range(min(trunc_lvl - 1, d - var - 2), -1, -1):
      pred: set[tuple[int, int]] = set()
      if tree + 1 <= min(trunc_lvl - 1, d - var - 2):
        pred.add((var, tree + 1))
      m = int(s.min_array(tree, var))
      if m == int(s.struct_array(tree, var, natural_order=True)):
        pred.add((m - 1, tree))
      elif tree - 1 >= 0:
        pred.add((m - 1, tree - 1))
      deps[var, tree] = pred
  for cell in deps:
    deps[cell] &= deps.keys()

  # Longest-path level of each cell; the descending sweep is already a
  # topological order, so one pass settles it.
  depth: dict[tuple[int, int], int] = {}
  for cell in sorted(deps, key=lambda c: (-c[0], -c[1])):
    depth[cell] = 0 if not deps[cell] else 1 + max(depth[p] for p in deps[cell])
  n_waves = max(depth.values()) + 1 if depth else 0
  return [sorted(c for c in deps if depth[c] == k) for k in range(n_waves)]


class BatchedWave(torch.nn.Module):
  """One parallel wave of the inverse cascade: K pairs and their wiring.

  Structurally the same object as :class:`BatchedTreeLevel` -- a stack of
  per-pair grids plus index tensors -- but keyed on a wave rather than a tree,
  and wired to the transposed ``(trunc_lvl + 1, d, n)`` scratch the inverse
  walks instead of the ``(n, d)`` columns the forward cascades use.
  """

  values: Tensor
  sy: Tensor | None
  sx: Tensor | None
  is_indep: Tensor
  col0_src: Tensor
  col1_src: Tensor
  col1_use_h1: Tensor
  out_hinv2: Tensor
  h1_rows: Tensor
  out_hfunc1: Tensor

  def __init__(
    self,
    values: Tensor,
    sy: Tensor | None,
    sx: Tensor | None,
    is_indep: Tensor,
    col0_src: Tensor,
    col1_src: Tensor,
    col1_use_h1: Tensor,
    out_hinv2: Tensor,
    h1_rows: Tensor,
    out_hfunc1: Tensor,
    is_linear: bool,
  ) -> None:
    super().__init__()
    self.register_buffer("values", values)
    if sy is not None:
      assert sx is not None
      self.register_buffer("sy", sy)
      self.register_buffer("sx", sx)
    else:
      self.sy = None
      self.sx = None
    for name, t in (
      ("is_indep", is_indep),
      ("col0_src", col0_src),
      ("col1_src", col1_src),
      ("col1_use_h1", col1_use_h1),
      ("out_hinv2", out_hinv2),
      ("h1_rows", h1_rows),
      ("out_hfunc1", out_hfunc1),
    ):
      self.register_buffer(name, t)
    self._is_linear = bool(is_linear)

  def apply_to(
    self, grid_points: Tensor, hinv2: Tensor, hfunc1: Tensor
  ) -> None:
    """Invert this wave's pairs in place on the flattened scratch.

    ``hinv2`` / ``hfunc1`` are ``((trunc_lvl + 1) * d, n)`` views, so a cell's
    slot is one row and the whole wave is one ``index_select`` per input and
    one ``index_copy_`` per output.
    """
    col0 = hinv2.index_select(0, self.col0_src)
    col1 = torch.where(
      self.col1_use_h1[:, None],
      hfunc1.index_select(0, self.col1_src),
      hinv2.index_select(0, self.col1_src),
    )
    # Clamped as a pair's `hinv2` / `hfunc1` dispatchers clamp their input: a
    # quantile this wave's predecessors produced may sit exactly on 0 or 1.
    col0, col1 = trim(col0, TENSOR_NS), trim(col1, TENSOR_NS)
    u_e = torch.stack([col0, col1], dim=-1)

    raw = inverse_integrate_1d_batched(
      grid_points, self.values, u_e, 2, self._is_linear, self.sx
    )
    inv = torch.where(self.is_indep[:, None], col0, raw)
    hinv2.index_copy_(0, self.out_hinv2, inv)

    if self.h1_rows.numel() == 0:
      return
    rows = self.h1_rows
    u_after = trim(
      torch.stack(
        [inv.index_select(0, rows), col1.index_select(0, rows)], dim=-1
      ),
      TENSOR_NS,
    )
    vals = self.values.index_select(0, rows)
    uu = u_after.clamp(0.0, 1.0)
    i, wx, _ = _locate(grid_points, uu[..., 0], self._is_linear)
    j, wy, dy = _locate(grid_points, uu[..., 1], self._is_linear)
    if self.sy is not None:
      h = _hfunc_from_cells(
        vals, self.sy.index_select(0, rows), i, wx, j, wy, dy
      )
    else:
      h = integrate_1d_batched(grid_points, vals, u_after, 1, self._is_linear)
    h = torch.where(
      self.is_indep.index_select(0, rows)[:, None],
      u_after[..., 1],
      h,
    )
    hfunc1.index_copy_(0, self.out_hfunc1, h)


def _shared_grid(
  tvc: TorchVinecop, trunc_lvl: int, d: int
) -> tuple[Tensor, bool, Tensor, Tensor, Tensor]:
  """The grid every pair stacks on, plus an independence pair built on it.

  ``TorchTllBicop`` gives an independence pair a 2x2 sentinel grid and no prefix
  tables, because none of its own evaluations read either -- every method
  short-circuits on ``is_indep``. A stacked level does read them: ``torch.stack``
  needs one shape across the level, and one pair without tables drops the whole
  level to the on-the-fly path. So the precomputation substitutes an independence density
  built on the shared grid, which is a real ``InterpolationGrid2D`` rather than
  a hand-derived table, so it cannot drift from what the pairs beside it do.

  Parameters
  ----------
  tvc : TorchVinecop
      The vine to precompute from.
  trunc_lvl, d : int
      Its truncation level and dimension.

  Returns
  -------
  tuple
      ``(grid_points, is_linear, values, sy, sx)`` -- the shared grid and the
      independence substitute for it.

  Raises
  ------
  NotBatchable
      If two pairs that are not independence copulas disagree on the grid.
  """
  # Deferred: `_interp` imports the kernels in this module, so the dependency
  # only goes the other way at call time.
  from ._bicop_interp import InterpolationGrid2D

  ref = None
  for t in range(trunc_lvl):
    for e in range(d - t - 1):
      bc = tvc._pair_module(t, e)
      if bc.is_indep:
        continue
      if ref is None:
        ref = bc
        continue
      # Every stacked pair is read against `ref`'s knots, so agreeing on the
      # shape is not enough: two grids of one size and different spacing
      # interpolate to different functions, and the level would be evaluated
      # on the wrong one without anything raising.
      if bc.interp_grid.values.shape != ref.interp_grid.values.shape:
        raise NotBatchable(
          "batched path requires one shared grid: pair copulas differ in grid "
          f"size ({tuple(ref.interp_grid.values.shape)} vs "
          f"{tuple(bc.interp_grid.values.shape)})."
        )
      if bool(bc.interp_grid._is_linear) != bool(ref.interp_grid._is_linear):
        raise NotBatchable(
          "batched path requires one shared grid: pair copulas differ in "
          "grid spacing (is_linear "
          f"{bool(ref.interp_grid._is_linear)} vs "
          f"{bool(bc.interp_grid._is_linear)})."
        )
      if not torch.equal(
        bc.interp_grid.grid_points, ref.interp_grid.grid_points
      ):
        raise NotBatchable(
          "batched path requires one shared grid: pair copulas of equal size "
          "differ in their grid points."
        )
  if ref is None:
    ref = tvc._pair_module(0, 0)
  gp = ref.interp_grid.grid_points
  m = int(gp.shape[0])
  flat = InterpolationGrid2D(
    grid_points=gp,
    values=torch.ones((m, m), dtype=gp.dtype, device=gp.device),
    norm_maxiter=0,
    is_linear=ref.interp_grid._is_linear,
  )
  sy, sx, _ = flat.build_caches()
  return gp, bool(ref.interp_grid._is_linear), flat.values, sy, sx


class BatchedVine(torch.nn.Module):
  """All tree levels of a :class:`TorchVinecop`, stacked and precomputed.

  Built lazily by :meth:`TorchVinecop._ensure_batched` on first call to any
  batched cascade. The wire-up tensors are computed once by walking the
  ``pyvinecopulib.RVineStructure`` accessors (``min_array``,
  ``struct_array``, ``needed_hfunc1`` / ``needed_hfunc2``), so the hot
  loop reads tensors only.

  Holds two groupings, because the cascades have two. ``pdf`` and
  ``rosenblatt`` run tree by tree -- edges at one tree level are independent
  going forward -- so those read :attr:`levels`. The inverse's dependencies
  run across tree levels, so it reads :attr:`waves` instead: the longest-path
  levels of the ``(var, tree)`` graph, computed once by :func:`inverse_waves`,
  each holding one cell from almost every tree.
  """

  grid_points: Tensor

  def __init__(
    self,
    *,
    grid_points: Tensor,
    levels: list[BatchedTreeLevel],
    order: list[int],
    inverse_order: list[int],
    d: int,
    trunc_lvl: int,
    waves: list[BatchedWave] | None = None,
  ) -> None:
    super().__init__()
    waves = waves or []
    self.register_buffer("grid_points", grid_points)
    self.levels = torch.nn.ModuleList(levels)
    self.waves = torch.nn.ModuleList(waves)
    # Plain Python attrs — these don't move with .to() but they're scalars
    # / int lists.
    self.order = order
    self.inverse_order = inverse_order
    self.d = d
    self.trunc_lvl = trunc_lvl

  def wave(self, k: int) -> BatchedWave:
    """Typed accessor for inverse-cascade wave ``k``."""
    return cast("BatchedWave", self.waves[k])

  @property
  def n_waves(self) -> int:
    """Number of parallel waves the inverse cascade decomposes into."""
    return len(self.waves)

  def level(self, t: int) -> BatchedTreeLevel:
    """Typed accessor for tree level ``t`` (``self.levels[t]`` returns
    ``Module``, but every element is a :class:`BatchedTreeLevel`)."""
    return cast("BatchedTreeLevel", self.levels[t])

  @classmethod
  def from_torch_vinecop(cls, tvc: TorchVinecop) -> BatchedVine:
    """Build a ``BatchedVine`` from a fitted :class:`TorchVinecop`.

    Walks ``tvc.pair_copulas`` and ``tvc.structure`` once; precomputes per-pair
    grids; collects per-level wiring tensors.
    """
    s = tvc.structure
    d = int(tvc.d)
    trunc_lvl = int(tvc.trunc_lvl)
    # Pull a reference tensor to get device.
    device = tvc._ref_tensor().device

    # Every pair that models something shares one grid, since they come from
    # the same fit; an independence pair does not, and is substituted.
    grid_points, is_linear, flat_v, flat_sy, flat_sx = _shared_grid(
      tvc, trunc_lvl, d
    )

    levels: list[BatchedTreeLevel] = []
    for t in range(trunc_lvl):
      N_t = d - t - 1
      vals: list[Tensor] = []
      sy_list: list[Tensor | None] = []
      sy_t_list: list[Tensor | None] = []
      is_indep: list[bool] = []
      col0_src: list[int] = []
      col1_src: list[int] = []
      col1_use_h1: list[bool] = []
      needs_h1_list: list[bool] = []
      needs_h2_list: list[bool] = []
      disc1: list[bool] = []
      disc2: list[bool] = []
      all_have_cache = True

      for e in range(N_t):
        bc = tvc._pair_module(t, e)
        m = int(s.min_array(t, e))
        sarr = int(s.struct_array(t, e, natural_order=True))
        if bc.is_indep:
          vals.append(flat_v)
          sy_list.append(flat_sy)
          sy_t_list.append(flat_sx.t())
        elif bc._sy is None:
          all_have_cache = False
          vals.append(bc.interp_grid.values)
          sy_list.append(None)
          sy_t_list.append(None)
        else:
          # `_tables` rather than the buffers, so a grid that started tracking
          # grad builds its tables in-graph -- `_ensure_batched` rebuilds on a
          # grad-signature change, which is what makes that reachable.
          sy, sx, _ = bc._tables()
          vals.append(bc.interp_grid.values)
          sy_list.append(sy)
          sy_t_list.append(sx.t())
        is_indep.append(bool(bc.is_indep))
        col0_src.append(e)
        col1_src.append(m - 1)
        col1_use_h1.append(m != sarr)
        needs_h1_list.append(bool(s.needed_hfunc1(t, e)))
        needs_h2_list.append(bool(s.needed_hfunc2(t, e)))
        types = tvc.pair_var_types(t, e)
        disc1.append(types[0] == "d")
        disc2.append(types[1] == "d")

      values = torch.stack(vals, dim=0).to(device=device)
      sy: Tensor | None
      sy_t: Tensor | None
      if all_have_cache:
        sy = torch.stack(cast("list[Tensor]", sy_list), dim=0).to(device=device)
        sy_t = torch.stack(cast("list[Tensor]", sy_t_list), dim=0).to(
          device=device
        )
      else:
        sy = sy_t = None

      level = BatchedTreeLevel(
        values=values,
        sy=sy,
        sy_t=sy_t,
        is_indep=torch.tensor(is_indep, dtype=torch.bool, device=device),
        col0_src=torch.tensor(col0_src, dtype=torch.long, device=device),
        col1_src=torch.tensor(col1_src, dtype=torch.long, device=device),
        col1_use_h1=torch.tensor(col1_use_h1, dtype=torch.bool, device=device),
        needs_h1=torch.tensor(needs_h1_list, dtype=torch.bool, device=device),
        needs_h2=torch.tensor(needs_h2_list, dtype=torch.bool, device=device),
        disc1=torch.tensor(disc1, dtype=torch.bool, device=device),
        disc2=torch.tensor(disc2, dtype=torch.bool, device=device),
        grid_points=grid_points.to(device=device),
        is_linear=is_linear,
      )
      levels.append(level)

    # Inverse-cascade waves: the same per-pair slabs, regrouped by the
    # dependency levelling and wired to the transposed scratch.
    waves: list[BatchedWave] = []
    for cells in inverse_waves(s, d, trunc_lvl):
      w_vals, w_sy, w_sx, w_indep = [], [], [], []
      c0, c1, use_h1, out_inv = [], [], [], []
      h1_rows, out_h1 = [], []
      w_has_cache = True
      for slot, (var, tree) in enumerate(cells):
        bc = tvc._pair_module(tree, var)
        mv = int(s.min_array(tree, var))
        sarr = int(s.struct_array(tree, var, natural_order=True))
        if bc.is_indep:
          w_vals.append(flat_v)
          w_sy.append(flat_sy)
          w_sx.append(flat_sx)
        elif bc._sy is None:
          w_has_cache = False
          w_vals.append(bc.interp_grid.values)
          w_sy.append(None)
          w_sx.append(None)
        else:
          tbl = bc._tables()
          w_vals.append(bc.interp_grid.values)
          w_sy.append(tbl[0])
          w_sx.append(tbl[1])
        w_indep.append(bool(bc.is_indep))
        c0.append((tree + 1) * d + var)
        c1.append(tree * d + (mv - 1))
        use_h1.append(mv != sarr)
        out_inv.append(tree * d + var)
        if var < d - 1 and bool(s.needed_hfunc1(tree, var)):
          h1_rows.append(slot)
          out_h1.append((tree + 1) * d + var)
      waves.append(
        BatchedWave(
          values=torch.stack(w_vals, dim=0).to(device=device),
          sy=(
            torch.stack(cast("list[Tensor]", w_sy), dim=0).to(device=device)
            if w_has_cache
            else None
          ),
          sx=(
            torch.stack(cast("list[Tensor]", w_sx), dim=0).to(device=device)
            if w_has_cache
            else None
          ),
          is_indep=torch.tensor(w_indep, dtype=torch.bool, device=device),
          col0_src=torch.tensor(c0, dtype=torch.long, device=device),
          col1_src=torch.tensor(c1, dtype=torch.long, device=device),
          col1_use_h1=torch.tensor(use_h1, dtype=torch.bool, device=device),
          out_hinv2=torch.tensor(out_inv, dtype=torch.long, device=device),
          h1_rows=torch.tensor(h1_rows, dtype=torch.long, device=device),
          out_hfunc1=torch.tensor(out_h1, dtype=torch.long, device=device),
          is_linear=is_linear,
        )
      )

    return cls(
      grid_points=grid_points.to(device=device),
      levels=levels,
      waves=waves,
      order=list(tvc.order),
      inverse_order=list(tvc.inverse_order),
      d=d,
      trunc_lvl=trunc_lvl,
    )
