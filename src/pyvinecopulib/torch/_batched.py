"""Batched pair-copula primitives for ``TorchVinecop``.

Stacks the per-pair ``InterpolationGrid2D`` state at one tree level into
``(N, m, m)`` tensors so the whole level fires one batched bilinear interp
(or one batched trapezoidal integration) instead of N separate calls.

Exposes :func:`interpolate_batched`, :func:`int_on_grid_batched`,
:func:`integrate_1d_batched`, :func:`integrate_2d_batched` and
:func:`inverse_integrate_1d_batched` — the ``(N, m, m)`` analogs of the
unbatched operations in :mod:`._interp` — plus :class:`BatchedTreeLevel`,
:class:`BatchedWave` and :class:`BatchedVine`, which stage one tree level
(resp. one wave of the inverse cascade, resp. an entire vine) of stacked
grids and wire-up tensors.

Intentionally side-by-side with :mod:`._interp` rather than rewriting it:
the legacy / lazy backends stay untouched so any regression is bisectable
to this file.
"""

from __future__ import annotations

from typing import Optional, cast

import torch
from torch import Tensor

from ..core._trim import _TRIM_LO, trim_bounds
from ..core.vinecop_base import _NotBatchable

#: Guard on a conditional total mass, so a zero-mass grid line cannot 0/0.
_MIN_MASS: float = 1e-20


def _trim(t: Tensor) -> Tensor:
  """Clamp ``t`` into the open unit interval at its own precision."""
  lo, hi = trim_bounds(torch, t.dtype)
  return t.clamp(lo, hi)


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
      the rotation must already be baked into ``values``).

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
  u = u.clamp(0.0, 1.0)
  N, n, _ = u.shape

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
    values: shape ``(N, m, m)``, baked pdf grids.
    u: shape ``(N, n, 2)``, queries (unrotated frame; the bake absorbs the
      rotation).
    cond_var: scalar in ``{1, 2}`` — kept consistent across all pairs in the
      batch because the bake puts every pair in the same "natural" frame.
      ``cond_var=1`` returns the h-function conditioning on ``u[..., 0]``;
      ``cond_var=2`` conditions on ``u[..., 1]``.

  Returns:
    Tensor of shape ``(N, n)`` with values in ``[1e-10, 1-1e-10]``.
  """
  if cond_var not in (1, 2):
    raise ValueError(f"cond_var must be 1 or 2; got {cond_var}")
  u = u.clamp(0.0, 1.0)
  N, n, _ = u.shape
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
  return _trim(number / denom.clamp_min(_MIN_MASS))


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
  cum: Optional[Tensor] = None,
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
  is_linear : bool, optional
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
  p = _trim(p.nan_to_num(0.5))

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
  safe_b = torch.where(b == 0.0, torch.ones_like(b), b)
  safe_d = torch.where(denom == 0.0, torch.ones_like(denom), denom)
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
  ``u: (N, n, 2)``, returns ``(N, n)`` clamped to ``[1e-10, 1-1e-10]``.
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
    tmpint * u2 / tmpint1.clamp_min(_TRIM_LO),
    torch.zeros_like(tmpint),
  )
  return _trim(out)


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
      Conditional distribution values, clamped to ``[1e-10, 1-1e-10]``.
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
  return _trim(num / den.clamp_min(_MIN_MASS))


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

  Indep handling:
  - ``is_indep: (N,) bool`` — true slots get short-circuit overrides
    (pdf=1, hfunc1=col1, hfunc2=col0). Their grids are sentinels.
  """

  # Class-level type hints so the buffers registered in __init__ are
  # statically typed as Tensors instead of nn.Module (cf. ``_sy`` in
  # TorchBicop, same pattern).
  values: Tensor
  grids2: Tensor | None
  tables2: Tensor | None
  is_indep: Tensor
  col0_src: Tensor
  col1_src: Tensor
  col1_use_h1: Tensor
  needs_h1: Tensor
  needs_h2: Tensor

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
    self._is_linear = bool(is_linear)

  @property
  def n_pairs(self) -> int:
    return int(self.values.shape[0])

  def gather_inputs(self, hfunc1_prev: Tensor, hfunc2_prev: Tensor) -> Tensor:
    """Build the per-pair ``(N, n, 2)`` input from the previous level's
    h-function columns. Lays the ``col0`` / ``col1`` selection out as a
    pair of gathers, plus a ``torch.where`` on the h1-vs-h2 source flag
    for ``col1``.
    """
    # hfunc{1,2}_prev: (n, d). Gather columns col0_src / col1_src -> (n, N).
    col0 = hfunc2_prev.index_select(dim=1, index=self.col0_src)  # (n, N)
    col1_h2 = hfunc2_prev.index_select(dim=1, index=self.col1_src)
    col1_h1 = hfunc1_prev.index_select(dim=1, index=self.col1_src)
    col1 = torch.where(self.col1_use_h1[None, :], col1_h1, col1_h2)
    # Stack into (N, n, 2): permute (n, N) -> (N, n) then stack on last dim.
    u_e = torch.stack([col0.t(), col1.t()], dim=-1)
    return u_e

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
    self, gp: Tensor, u: Tensor, loc: tuple
  ) -> tuple[Tensor, Tensor]:
    """Both h-functions at one located query.

    ``hfunc1`` conditions on argument 1 and integrates argument 2; ``hfunc2``
    does the reverse against the transposed grid. That is the same kernel with
    the two triples swapped, so with the grids stacked it is one call on
    ``2N`` pairs instead of two on ``N`` -- worth doing in a cascade whose cost
    is the number of calls.
    """
    i, wx, dx, j, wy, dy = loc
    if self.grids2 is None or self.tables2 is None:
      raw = torch.cat(
        [
          integrate_1d_batched(gp, self.values, u, 1, self._is_linear),
          integrate_1d_batched(gp, self.values, u, 2, self._is_linear),
        ],
        0,
      )
    else:
      raw = _hfunc_from_cells(
        self.grids2,
        self.tables2,
        torch.cat([i, j], 0),
        torch.cat([wx, wy], 0),
        torch.cat([j, i], 0),
        torch.cat([wy, wx], 0),
        torch.cat([dy, dx], 0),
      )
    # An independent pair returns its own free argument -- argument 2 for
    # `hfunc1`, argument 1 for `hfunc2` -- which is the same swap again.
    both = torch.where(
      self.is_indep.repeat(2)[:, None],
      _trim(torch.cat([u[..., 1], u[..., 0]], 0)),
      raw.clamp(0.0, 1.0),
    )
    return both[: self.n_pairs], both[self.n_pairs :]

  def _h1_at(self, gp: Tensor, u: Tensor, loc: tuple) -> Tensor:
    return self._h1_h2_at(gp, u, loc)[0]

  def _h2_at(self, gp: Tensor, u: Tensor, loc: tuple) -> Tensor:
    return self._h1_h2_at(gp, u, loc)[1]

  def pdf(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair pdf at the stacked queries. Indep slots return 1."""
    i, wx, _, j, wy, _ = self._locate_both(grid_points, u)
    return self._pdf_at(i, j, wx, wy)

  def hfunc1(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair hfunc1.

    ``cache=True`` reconstructs the conditional distribution function exactly
    from the cumulative table in O(1); ``cache=False`` runs
    :func:`integrate_1d_batched` on ``values`` at every call. The two agree to
    floating point: they are the same integral, summed in a different order.
    """
    return self._h1_at(grid_points, u, self._locate_both(grid_points, u))

  def hfunc2(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair hfunc2; see :meth:`hfunc1`."""
    return self._h2_at(grid_points, u, self._locate_both(grid_points, u))

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
    h1, h2 = self._h1_h2_at(grid_points, u, loc)
    return self._pdf_at(i, j, wx, wy), h1, h2

  def h1_h2(self, grid_points: Tensor, u: Tensor) -> tuple[Tensor, Tensor]:
    """``(hfunc1, hfunc2)`` for one tree level; see :meth:`pdf_h1_h2`."""
    return self._h1_h2_at(grid_points, u, self._locate_both(grid_points, u))


def inverse_waves(s, d: int, trunc_lvl: int) -> list[list[tuple[int, int]]]:
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
      deps[(var, tree)] = pred
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
  sy: Optional[Tensor]
  sx: Optional[Tensor]
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
    sy: Optional[Tensor],
    sx: Optional[Tensor],
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
    u_e = torch.stack([col0, col1], dim=-1)

    raw = inverse_integrate_1d_batched(
      grid_points, self.values, u_e, 2, self._is_linear, self.sx
    ).clamp(0.0, 1.0)
    inv = torch.where(self.is_indep[:, None], _trim(col0), raw)
    hinv2.index_copy_(0, self.out_hinv2, inv)

    if self.h1_rows.numel() == 0:
      return
    rows = self.h1_rows
    u_after = torch.stack(
      [inv.index_select(0, rows), col1.index_select(0, rows)], dim=-1
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
      _trim(u_after[..., 1]),
      h.clamp(0.0, 1.0),
    )
    hfunc1.index_copy_(0, self.out_hfunc1, h)


def _shared_grid(tvc, trunc_lvl: int, d: int):
  """The grid every pair stacks on, plus an independence pair built on it.

  ``TorchBicop`` gives an independence pair a 2x2 sentinel grid and no prefix
  tables, because none of its own evaluations read either -- every method
  short-circuits on ``is_indep``. A stacked level does read them: ``torch.stack``
  needs one shape across the level, and one pair without tables drops the whole
  level to the on-the-fly path. So the bake substitutes an independence density
  built on the shared grid, which is a real ``InterpolationGrid2D`` rather than
  a hand-derived table, so it cannot drift from what the pairs beside it do.

  Parameters
  ----------
  tvc : TorchVinecop
      The vine being baked.
  trunc_lvl, d : int
      Its truncation level and dimension.

  Returns
  -------
  tuple
      ``(grid_points, is_linear, values, sy, sx)`` -- the shared grid and the
      independence substitute for it.

  Raises
  ------
  _NotBatchable
      If two pairs that are not independence copulas disagree on the grid.
  """
  # Deferred: `_interp` imports the kernels in this module, so the dependency
  # only goes the other way at call time.
  from ._interp import InterpolationGrid2D

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
        raise _NotBatchable(
          "batched path requires one shared grid: pair copulas differ in grid "
          f"size ({tuple(ref.interp_grid.values.shape)} vs "
          f"{tuple(bc.interp_grid.values.shape)})."
        )
      if bool(bc.interp_grid._is_linear) != bool(ref.interp_grid._is_linear):
        raise _NotBatchable(
          "batched path requires one shared grid: pair copulas differ in "
          "grid spacing (is_linear "
          f"{bool(ref.interp_grid._is_linear)} vs "
          f"{bool(bc.interp_grid._is_linear)})."
        )
      if not torch.equal(
        bc.interp_grid.grid_points, ref.interp_grid.grid_points
      ):
        raise _NotBatchable(
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
  """All tree levels of a :class:`TorchVinecop`, stacked and pre-baked.

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
    waves: "list[BatchedWave]" = [],
  ) -> None:
    super().__init__()
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
    return cast(BatchedWave, self.waves[k])

  @property
  def n_waves(self) -> int:
    """Number of parallel waves the inverse cascade decomposes into."""
    return len(self.waves)

  def level(self, t: int) -> BatchedTreeLevel:
    """Typed accessor for tree level ``t`` (``self.levels[t]`` returns
    ``Module``, but every element is a :class:`BatchedTreeLevel`)."""
    return cast(BatchedTreeLevel, self.levels[t])

  @classmethod
  def from_torch_vinecop(cls, tvc) -> "BatchedVine":
    """Build a ``BatchedVine`` from a fitted :class:`TorchVinecop`.

    Walks ``tvc.pair_copulas`` and ``tvc.structure`` once; bakes per-pair
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
          # grad bakes its tables in-graph -- `_ensure_batched` re-bakes on a
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

      values = torch.stack(vals, dim=0).to(device=device)
      sy: Tensor | None
      sy_t: Tensor | None
      if all_have_cache:
        sy = torch.stack(cast(list[Tensor], sy_list), dim=0).to(device=device)
        sy_t = torch.stack(cast(list[Tensor], sy_t_list), dim=0).to(
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
            torch.stack(cast(list[Tensor], w_sy), dim=0).to(device=device)
            if w_has_cache
            else None
          ),
          sx=(
            torch.stack(cast(list[Tensor], w_sx), dim=0).to(device=device)
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
