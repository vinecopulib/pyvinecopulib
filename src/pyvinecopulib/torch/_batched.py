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

import torch
from torch import Tensor

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
