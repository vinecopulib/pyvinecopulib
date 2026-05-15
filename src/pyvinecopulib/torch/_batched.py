"""Batched pair-copula primitives for ``TorchVinecop``.

Stacks the per-pair ``InterpolationGrid2D`` state at one tree level into
``(N, m, m)`` tensors so the whole level fires one batched bilinear interp
(or one batched trapezoidal integration) instead of N separate calls.

This module owns:

* :func:`_bake_pair_grid` / :func:`_bake_pair_metadata` — rotation baking.
  For each rotation in ``{0, 90, 180, 270}`` we permute / flip the per-pair
  ``values`` (and ``cdf_cache`` / ``h1_cache`` / ``h2_cache`` if available)
  so eval-time calls become branchless: ``interpolate_batched(values_baked,
  u)`` directly returns the pdf at the *unrotated* ``u``. For the linear-in-u
  post-rotation arithmetic in :meth:`TorchBicop.cdf` / :meth:`hfunc1` /
  :meth:`hfunc2` we carry per-pair sign/offset tensors and a ``(4,)`` cdf
  correction so the bake fully absorbs the rotation.

* The batched primitives :func:`interpolate_batched`, :func:`interp_at_batched`,
  :func:`int_on_grid_batched`, :func:`integrate_1d_batched`,
  :func:`integrate_2d_batched`. These are the ``(N, m, m)`` analogs of the
  unbatched ones in :mod:`._interp`.

Intentionally side-by-side with :mod:`._interp` rather than rewriting it:
the legacy / lazy backends stay untouched so any regression is bisectable
to this file.
"""

from __future__ import annotations

from typing import cast

import torch
from torch import Tensor

from ._util import solve_bisection

# Mirror the trim constants from ``_interp`` so the batched path produces
# numerically identical outputs to the per-pair path.
_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10
_STRIP_FLOOR: float = 1e-4


# --------------------------------------------------------------------------- #
# Rotation baking                                                              #
# --------------------------------------------------------------------------- #


def _bake_one(grid: Tensor, rotation: int) -> Tensor:
  """Apply a rotation to one ``(m, m)`` grid so eval at unrotated ``u`` matches.

  Encodes the same coordinate transform as :meth:`TorchBicop._rotate_input`
  applied at the grid-index level. For all four rotations the symmetric
  normal grid (`g[m-1-i] == 1 - g[i]`) lets us swap index permutations for
  what would otherwise be reflections around 0.5.
  """
  if rotation == 0:
    return grid.contiguous()
  if rotation == 90:
    # (u1, u2) -> (u2, 1 - u1)
    # baked[i, j] = orig[j, m-1-i]
    return grid.transpose(-2, -1).flip(-2).contiguous()
  if rotation == 180:
    # (u1, u2) -> (1 - u1, 1 - u2)
    # baked[i, j] = orig[m-1-i, m-1-j]
    return grid.flip(-1).flip(-2).contiguous()
  if rotation == 270:
    # (u1, u2) -> (1 - u2, u1)
    # baked[i, j] = orig[m-1-j, i]
    return grid.transpose(-2, -1).flip(-1).contiguous()
  raise ValueError(f"rotation must be in {{0, 90, 180, 270}}; got {rotation}")


def bake_pair_grid(
  values: Tensor,
  rotation: int,
  *,
  h1_cache: Tensor | None = None,
  h2_cache: Tensor | None = None,
  cdf_cache: Tensor | None = None,
) -> dict[str, Tensor | None]:
  """Bake a single pair-copula's grids for a given rotation.

  Returns a dict with keys ``values_baked``, ``h1_baked``, ``h2_baked``,
  ``cdf_baked``. The cache keys are ``None`` when the source caches are
  ``None`` (i.e. ``cache_integrals=False``).

  For ``rotation in {90, 270}`` the eval-time ``hfunc1`` of the original
  pair calls ``_hfunc_raw(u_rot, 2)`` — i.e. it uses the *other* cache.
  We absorb that swap into the bake so that ``h1_baked`` is the cache the
  batched ``hfunc1_batched`` should query (and likewise for ``h2_baked``).
  """
  vb = _bake_one(values, rotation)

  # h1 source: which of {h1_cache, h2_cache} feeds the batched hfunc1.
  # Mirrors :meth:`TorchBicop.hfunc1` rotation dispatch (bicop.py:252-260).
  if h1_cache is None or h2_cache is None:
    h1_b = None
    h2_b = None
  else:
    if rotation in (0, 180):
      h1_src, h2_src = h1_cache, h2_cache
    else:  # 90 or 270
      h1_src, h2_src = h2_cache, h1_cache
    h1_b = _bake_one(h1_src, rotation)
    h2_b = _bake_one(h2_src, rotation)

  cdf_b = _bake_one(cdf_cache, rotation) if cdf_cache is not None else None

  return {
    "values_baked": vb,
    "h1_baked": h1_b,
    "h2_baked": h2_b,
    "cdf_baked": cdf_b,
  }


def bake_pair_metadata(
  rotation: int, dtype: torch.dtype = torch.float64
) -> dict[str, Tensor]:
  """Return the per-pair sign/offset/correction tensors for a rotation.

  Intended for the ``cache_integrals=True`` cascade. After
  :func:`interp_at_batched` on the baked cache, eval applies:
      hfunc1 = h1_sign * interp + h1_offset
      hfunc2 = h2_sign * interp + h2_offset
      cdf    = c_p * interp + c_u0 * u_orig[..., 0]
                            + c_u1 * u_orig[..., 1] + c_const
  Mirrors the post-rotation arithmetic in :meth:`TorchBicop.cdf` /
  :meth:`hfunc1` / :meth:`hfunc2`. Scalars are wrapped in 0-D tensors so the
  batched callers can simply stack them.

  Note: the ``cache_integrals=False`` cascade does **not** use this metadata
  on top of the bake. Trapezoidal integration over a baked grid is not
  invariant under coordinate-permuting bakes (partial-cell handling visits
  different boundary cells), so the cache=False path applies the rotation
  at eval time to the unrotated ``values`` instead.
  """
  one = torch.tensor(1.0, dtype=dtype)
  zero = torch.tensor(0.0, dtype=dtype)
  neg_one = torch.tensor(-1.0, dtype=dtype)
  # hfunc1: 0 ->  +1, 0;  90 ->  +1, 0; 180 -> -1, 1; 270 -> -1, 1
  # hfunc2: 0 ->  +1, 0;  90 -> -1, 1; 180 -> -1, 1; 270 ->  +1, 0
  if rotation == 0:
    h1_sign, h1_offset = one, zero
    h2_sign, h2_offset = one, zero
    cdf_corr = torch.tensor([1.0, 0.0, 0.0, 0.0], dtype=dtype)
  elif rotation == 90:
    h1_sign, h1_offset = one, zero
    h2_sign, h2_offset = neg_one, one
    # cdf_rot: u[:, 1] - p  ->  (c_p, c_u0, c_u1, c_const) = (-1, 0, 1, 0)
    cdf_corr = torch.tensor([-1.0, 0.0, 1.0, 0.0], dtype=dtype)
  elif rotation == 180:
    h1_sign, h1_offset = neg_one, one
    h2_sign, h2_offset = neg_one, one
    # cdf_rot: p - 1 + u[:, 0] + u[:, 1]  ->  (1, 1, 1, -1)
    cdf_corr = torch.tensor([1.0, 1.0, 1.0, -1.0], dtype=dtype)
  elif rotation == 270:
    h1_sign, h1_offset = neg_one, one
    h2_sign, h2_offset = one, zero
    # cdf_rot: u[:, 0] - p  ->  (-1, 1, 0, 0)
    cdf_corr = torch.tensor([-1.0, 1.0, 0.0, 0.0], dtype=dtype)
  else:
    raise ValueError(f"rotation must be in {{0, 90, 180, 270}}; got {rotation}")
  return {
    "h1_sign": h1_sign,
    "h1_offset": h1_offset,
    "h2_sign": h2_sign,
    "h2_offset": h2_offset,
    "cdf_correction": cdf_corr,
  }


# --------------------------------------------------------------------------- #
# Batched bilinear interpolation                                               #
# --------------------------------------------------------------------------- #


def _batched_cell_index(grid_points: Tensor, u: Tensor) -> Tensor:
  """Per-element cell index, clamped to ``[0, m-2]``. Preserves ``u``'s shape."""
  m = grid_points.shape[0]
  return (
    torch.searchsorted(grid_points, u.contiguous(), right=False) - 1
  ).clamp(0, m - 2)


def interpolate_batched(
  grid_points: Tensor, values: Tensor, u: Tensor
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

  i = _batched_cell_index(grid_points, u[..., 0])  # (N, n)
  j = _batched_cell_index(grid_points, u[..., 1])  # (N, n)

  N_idx = (
    torch.arange(N, device=values.device).unsqueeze(-1).expand(N, n)
  )  # (N, n)

  z11 = values[N_idx, i, j]
  z12 = values[N_idx, i, j + 1]
  z21 = values[N_idx, i + 1, j]
  z22 = values[N_idx, i + 1, j + 1]

  x1 = grid_points[i]
  x2 = grid_points[i + 1]
  y1 = grid_points[j]
  y2 = grid_points[j + 1]
  x = u[..., 0]
  y = u[..., 1]
  x2x = x2 - x
  y2y = y2 - y
  xx1 = x - x1
  yy1 = y - y1
  denom = (x2 - x1) * (y2 - y1)
  return (
    z11 * x2x * y2y + z21 * xx1 * y2y + z12 * x2x * yy1 + z22 * xx1 * yy1
  ) / denom


def interp_at_batched(grid_points: Tensor, cache: Tensor, u: Tensor) -> Tensor:
  """Bilinearly interpolate a stacked precomputed cache.

  Same shape contract as :func:`interpolate_batched`; the distinction is
  semantic — ``cache`` is a precomputed integral grid (e.g. h-func or cdf),
  whereas ``values`` in :func:`interpolate_batched` is the pdf grid.
  """
  return interpolate_batched(grid_points, cache, u)


# --------------------------------------------------------------------------- #
# Batched trapezoidal integration (``cache_integrals=False`` path)             #
# --------------------------------------------------------------------------- #


def int_on_grid_batched(
  grid_points: Tensor, upr: Tensor, vals: Tensor
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
  grid_points: Tensor, values: Tensor, u: Tensor, cond_var: int
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

  cell = _batched_cell_index(grid_points, u_fixed)  # (N, n)
  g_lo = grid_points[cell]
  g_hi = grid_points[cell + 1]
  t = ((u_fixed - g_lo) / (g_hi - g_lo)).unsqueeze(-1)  # (N, n, 1)

  # Gather two strips of shape (N, n, m): for each (k, l) take
  #   values[k, cell[k, l], :]      (cond_var=1)
  #   values[k, :, cell[k, l]]      (cond_var=2)
  # gather along ``fixed_axis`` after expanding the cell index along the
  # other (free) axis to size m.
  if fixed_axis == 1:
    idx_lo = cell.unsqueeze(-1).expand(N, n, m)
    idx_hi = (cell + 1).unsqueeze(-1).expand(N, n, m)
    v_lo = values.gather(dim=1, index=idx_lo)  # (N, n, m)
    v_hi = values.gather(dim=1, index=idx_hi)
  else:
    idx_lo = cell.unsqueeze(1).expand(N, m, n)
    idx_hi = (cell + 1).unsqueeze(1).expand(N, m, n)
    v_lo = values.gather(dim=2, index=idx_lo).transpose(1, 2)  # (N, n, m)
    v_hi = values.gather(dim=2, index=idx_hi).transpose(1, 2)

  strip = ((1.0 - t) * v_lo + t * v_hi).clamp_min(_STRIP_FLOOR)  # (N, n, m)

  numer = int_on_grid_batched(grid_points, u_free, strip)  # (N, n)
  denom = int_on_grid_batched(
    grid_points, torch.ones_like(u_free), strip
  )  # (N, n)
  return (numer / denom).clamp(_TRIM_LO, _TRIM_HI)


def integrate_2d_batched(
  grid_points: Tensor, values: Tensor, u: Tensor
) -> Tensor:
  """Batched bivariate CDF (trapezoidal-trapezoidal).

  Same shape contract as :func:`integrate_1d_batched`: ``values: (N, m, m)``,
  ``u: (N, n, 2)``, returns ``(N, n)`` clamped to ``[1e-10, 1-1e-10]``.
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
  strip = int_on_grid_batched(grid_points, upr_inner, vals_inner)  # (N, n, m)

  # Outer pass: integrate strip[k, l, :] up to u1[k, l].
  return int_on_grid_batched(grid_points, u1, strip).clamp(_TRIM_LO, _TRIM_HI)


# --------------------------------------------------------------------------- #
# Per-tree-level stacked state + per-vine container                            #
# --------------------------------------------------------------------------- #


class BatchedTreeLevel(torch.nn.Module):
  # Class-level type hints so the buffers registered in __init__ are
  # statically typed as Tensors instead of nn.Module (cf. ``_cdf_cache``
  # in TorchBicop, same pattern).
  values_baked: Tensor
  h1_baked: Tensor | None
  h2_baked: Tensor | None
  cdf_baked: Tensor | None
  h1_sign: Tensor
  h1_offset: Tensor
  h2_sign: Tensor
  h2_offset: Tensor
  cdf_correction: Tensor
  is_indep: Tensor
  col0_src: Tensor
  col1_src: Tensor
  col1_use_h1: Tensor
  needs_h1: Tensor
  needs_h2: Tensor

  """Stacked state for every pair-copula at one tree level of a vine.

  All buffers are registered so ``.to(device)`` / ``.to(dtype)`` move them
  together with the parent :class:`BatchedVine`.

  Grids (per pair):
  - ``values_baked: (N, m, m)`` — pdf grid post-rotation-bake. For fitted
    TLL vines (rotation == 0) this is exactly the source ``values``.
  - ``h1_baked, h2_baked, cdf_baked: (N, m, m) | None`` — present only when
    every source pair was constructed with ``cache_integrals=True``.

  Post-bake metadata (per pair, only meaningful in the ``cache=True`` path):
  - ``h1_sign, h1_offset, h2_sign, h2_offset: (N,)`` — folds the
    ``1 - x`` post-flips for h-functions at rotations 90 / 180 / 270.
  - ``cdf_correction: (N, 4)`` — coefficients ``(c_p, c_u0, c_u1, c_const)``
    of the affine post-rotation correction for the cdf.

  Wiring (per pair, same across pdf / rosenblatt / inverse cascades):
  - ``col0_src: (N,) long`` — column to read for ``col0`` (= edge index).
  - ``col1_src: (N,) long`` — column to read for ``col1`` (= ``min_array - 1``).
  - ``col1_use_h1: (N,) bool`` — whether ``col1`` reads from ``hfunc1`` (else
    ``hfunc2``).
  - ``needs_h1, needs_h2: (N,) bool`` — whether the cascade requires this
    pair's h-function output at the next tree.

  Indep handling:
  - ``is_indep: (N,) bool`` — true slots get short-circuit overrides
    (pdf=1, hfunc1=col1, hfunc2=col0). Their baked grids are sentinels.
  """

  def __init__(
    self,
    *,
    values_baked: Tensor,
    h1_baked: Tensor | None,
    h2_baked: Tensor | None,
    cdf_baked: Tensor | None,
    h1_sign: Tensor,
    h1_offset: Tensor,
    h2_sign: Tensor,
    h2_offset: Tensor,
    cdf_correction: Tensor,
    is_indep: Tensor,
    col0_src: Tensor,
    col1_src: Tensor,
    col1_use_h1: Tensor,
    needs_h1: Tensor,
    needs_h2: Tensor,
  ) -> None:
    super().__init__()
    self.register_buffer("values_baked", values_baked)
    # Optional caches go in via register_buffer with persistent=False so
    # their absence (None) doesn't clutter state_dict.
    if h1_baked is not None:
      self.register_buffer("h1_baked", h1_baked)
      self.register_buffer("h2_baked", h2_baked)
      self.register_buffer("cdf_baked", cdf_baked)
      self._has_cache = True
    else:
      self.h1_baked = None
      self.h2_baked = None
      self.cdf_baked = None
      self._has_cache = False
    self.register_buffer("h1_sign", h1_sign)
    self.register_buffer("h1_offset", h1_offset)
    self.register_buffer("h2_sign", h2_sign)
    self.register_buffer("h2_offset", h2_offset)
    self.register_buffer("cdf_correction", cdf_correction)
    self.register_buffer("is_indep", is_indep)
    self.register_buffer("col0_src", col0_src)
    self.register_buffer("col1_src", col1_src)
    self.register_buffer("col1_use_h1", col1_use_h1)
    self.register_buffer("needs_h1", needs_h1)
    self.register_buffer("needs_h2", needs_h2)

  @property
  def has_cache(self) -> bool:
    return self._has_cache

  @property
  def n_pairs(self) -> int:
    return int(self.values_baked.shape[0])

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

  def pdf(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair pdf at the stacked queries. Indep slots return 1."""
    raw = interpolate_batched(grid_points, self.values_baked, u).clamp_min(
      1e-20
    )
    return torch.where(self.is_indep[:, None], torch.ones_like(raw), raw)

  def log_pdf(self, grid_points: Tensor, u: Tensor) -> Tensor:
    return self.pdf(grid_points, u).log()

  def hfunc1(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair hfunc1. Cache=True uses baked caches; cache=False
    integrates ``values_baked`` on-the-fly via :func:`integrate_1d_batched`.
    For rotation=0 pairs (the fitted-TLL common case) the two paths are
    numerically equivalent; for non-zero rotations the cache=False path
    incurs ~5e-5 trapezoidal noise vs :meth:`TorchBicop.hfunc1`.
    """
    if self.h1_baked is not None:
      raw = interp_at_batched(grid_points, self.h1_baked, u)
    else:
      raw = integrate_1d_batched(grid_points, self.values_baked, u, cond_var=1)
    h = (self.h1_sign[:, None] * raw + self.h1_offset[:, None]).clamp(0.0, 1.0)
    return torch.where(
      self.is_indep[:, None], u[..., 1].clamp(_TRIM_LO, _TRIM_HI), h
    )

  def hfunc2(self, grid_points: Tensor, u: Tensor) -> Tensor:
    if self.h2_baked is not None:
      raw = interp_at_batched(grid_points, self.h2_baked, u)
    else:
      raw = integrate_1d_batched(grid_points, self.values_baked, u, cond_var=2)
    h = (self.h2_sign[:, None] * raw + self.h2_offset[:, None]).clamp(0.0, 1.0)
    return torch.where(
      self.is_indep[:, None], u[..., 0].clamp(_TRIM_LO, _TRIM_HI), h
    )

  def hinv2(self, grid_points: Tensor, col1: Tensor, p: Tensor) -> Tensor:
    """Invert hfunc2 per pair: solve ``hfunc2(stack(x, col1)) = p`` for x.

    Args:
      grid_points: shared ``(m,)``.
      col1: shape ``(N, n)`` — the "fixed" second column.
      p: shape ``(N, n)`` — target.

    Returns:
      Tensor of shape ``(N, n)`` with the bisection roots.
    """
    # solve_bisection is shape-polymorphic; we flatten to (N*n, 1) to keep
    # the closure simple.
    N, n = col1.shape
    p_flat = p.reshape(-1, 1)
    # Need to remember the per-pair rotation/cache info: re-index by k via
    # u_e of shape (N, n, 2) -> rebuild inside fun. Easiest: keep col1 in
    # (N, n) and pass through the closure as a tensor.

    def fun(x_flat: Tensor) -> Tensor:
      x = x_flat.reshape(N, n)
      u_e = torch.stack([x, col1], dim=-1)
      h = self.hfunc2(grid_points, u_e)  # (N, n)
      return (h - p).reshape(-1, 1) - 0  # explicit dtype preservation

    x_a = torch.zeros_like(p_flat)
    x_b = torch.ones_like(p_flat)
    x = solve_bisection(fun, x_a=x_a, x_b=x_b)
    out = x.reshape(N, n).clamp(0.0, 1.0)
    # Indep override: hinv2(u, col1) = u (the input dimension)
    return torch.where(self.is_indep[:, None], p.clamp(_TRIM_LO, _TRIM_HI), out)


class BatchedVine(torch.nn.Module):
  grid_points: Tensor

  """All tree levels of a :class:`TorchVinecop`, stacked and pre-baked.

  Built lazily by :meth:`TorchVinecop._ensure_batched` on first call to any
  batched cascade. The wire-up tensors are computed once by walking the
  ``pyvinecopulib.RVineStructure`` accessors (``min_array``,
  ``struct_array``, ``needed_hfunc1`` / ``needed_hfunc2``), so the hot
  loop reads tensors only.

  Used by ``pdf`` and ``rosenblatt`` with ``batched=True`` (edges at a
  fixed tree level are independent in the forward cascades, so each level
  fires one batched bicop call). ``inverse_rosenblatt + batched=True``
  falls back to the per-pair legacy cascade in v1 — the inverse DAG
  spans the full ``(var, tree)`` lattice and can't be flattened to
  tree-level waves without a global topological sort.
  """

  def __init__(
    self,
    *,
    grid_points: Tensor,
    levels: list[BatchedTreeLevel],
    order: list[int],
    inverse_order: list[int],
    d: int,
    trunc_lvl: int,
  ) -> None:
    super().__init__()
    self.register_buffer("grid_points", grid_points)
    self.levels = torch.nn.ModuleList(levels)
    # Plain Python attrs — these don't move with .to() but they're scalars
    # / int lists.
    self.order = order
    self.inverse_order = inverse_order
    self.d = d
    self.trunc_lvl = trunc_lvl

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
    # Pull a reference tensor to get device / dtype.
    ref = tvc._ref_tensor()
    device, dtype = ref.device, ref.dtype

    # The grid is shared by all pairs (same `make_normal_grid(m, dtype)`).
    grid_points = tvc._pair(0, 0).interp_grid.grid_points

    levels: list[BatchedTreeLevel] = []
    for t in range(trunc_lvl):
      N_t = d - t - 1
      values_bake: list[Tensor] = []
      h1_bake: list[Tensor | None] = []
      h2_bake: list[Tensor | None] = []
      cdf_bake: list[Tensor | None] = []
      h1_sign: list[Tensor] = []
      h1_offset: list[Tensor] = []
      h2_sign: list[Tensor] = []
      h2_offset: list[Tensor] = []
      cdf_corr: list[Tensor] = []
      is_indep: list[bool] = []
      col0_src: list[int] = []
      col1_src: list[int] = []
      col1_use_h1: list[bool] = []
      needs_h1_list: list[bool] = []
      needs_h2_list: list[bool] = []
      all_have_cache = True

      for e in range(N_t):
        bc = tvc._pair(t, e)
        rotation = int(bc.rotation)
        m = int(s.min_array(t, e))
        sarr = int(s.struct_array(t, e, natural_order=True))
        # Bake the per-pair grids.
        grids = bake_pair_grid(
          values=bc.interp_grid.values,
          rotation=rotation,
          h1_cache=bc._hfunc1_cache,
          h2_cache=bc._hfunc2_cache,
          cdf_cache=bc._cdf_cache,
        )
        meta = bake_pair_metadata(rotation, dtype=dtype)

        vb = grids["values_baked"]
        assert vb is not None  # values_baked is never None
        values_bake.append(vb)
        if grids["h1_baked"] is None:
          all_have_cache = False
        h1_bake.append(grids["h1_baked"])
        h2_bake.append(grids["h2_baked"])
        cdf_bake.append(grids["cdf_baked"])
        h1_sign.append(meta["h1_sign"])
        h1_offset.append(meta["h1_offset"])
        h2_sign.append(meta["h2_sign"])
        h2_offset.append(meta["h2_offset"])
        cdf_corr.append(meta["cdf_correction"])
        is_indep.append(bool(bc.is_indep))

        col0_src.append(e)
        col1_src.append(m - 1)
        col1_use_h1.append(m != sarr)
        needs_h1_list.append(bool(s.needed_hfunc1(t, e)))
        needs_h2_list.append(bool(s.needed_hfunc2(t, e)))

      values_baked = torch.stack(values_bake, dim=0).to(device=device)
      h1_baked: Tensor | None
      h2_baked: Tensor | None
      cdf_baked: Tensor | None
      if all_have_cache:
        h1_baked = torch.stack(cast(list[Tensor], h1_bake), dim=0).to(
          device=device
        )
        h2_baked = torch.stack(cast(list[Tensor], h2_bake), dim=0).to(
          device=device
        )
        cdf_baked = torch.stack(cast(list[Tensor], cdf_bake), dim=0).to(
          device=device
        )
      else:
        h1_baked = h2_baked = cdf_baked = None

      level = BatchedTreeLevel(
        values_baked=values_baked,
        h1_baked=h1_baked,
        h2_baked=h2_baked,
        cdf_baked=cdf_baked,
        h1_sign=torch.stack(h1_sign, dim=0).to(device=device),
        h1_offset=torch.stack(h1_offset, dim=0).to(device=device),
        h2_sign=torch.stack(h2_sign, dim=0).to(device=device),
        h2_offset=torch.stack(h2_offset, dim=0).to(device=device),
        cdf_correction=torch.stack(cdf_corr, dim=0).to(device=device),
        is_indep=torch.tensor(is_indep, dtype=torch.bool, device=device),
        col0_src=torch.tensor(col0_src, dtype=torch.long, device=device),
        col1_src=torch.tensor(col1_src, dtype=torch.long, device=device),
        col1_use_h1=torch.tensor(col1_use_h1, dtype=torch.bool, device=device),
        needs_h1=torch.tensor(needs_h1_list, dtype=torch.bool, device=device),
        needs_h2=torch.tensor(needs_h2_list, dtype=torch.bool, device=device),
      )
      levels.append(level)

    return cls(
      grid_points=grid_points.to(device=device),
      levels=levels,
      order=list(tvc.order),
      inverse_order=list(tvc.inverse_order),
      d=d,
      trunc_lvl=trunc_lvl,
    )
