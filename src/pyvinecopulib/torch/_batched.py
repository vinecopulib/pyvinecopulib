"""Batched pair-copula primitives for ``TorchVinecop``.

Stacks the per-pair ``InterpolationGrid2D`` state at one tree level into
``(N, m, m)`` tensors so the whole level fires one batched bilinear interp
(or one batched trapezoidal integration) instead of N separate calls.

Exposes :func:`interpolate_batched`, :func:`interp_at_batched`,
:func:`int_on_grid_batched`, :func:`integrate_1d_batched`,
:func:`integrate_2d_batched` — the ``(N, m, m)`` analogs of the unbatched
operations in :mod:`._interp` — plus :class:`BatchedTreeLevel` and
:class:`BatchedVine`, which stage one tree level (resp. an entire vine)
of stacked grids and wire-up tensors.

Intentionally side-by-side with :mod:`._interp` rather than rewriting it:
the legacy / lazy backends stay untouched so any regression is bisectable
to this file.
"""

from __future__ import annotations

from typing import cast

import torch
from torch import Tensor

# Mirror the trim constants from ``_interp`` so the batched path produces
# numerically identical outputs to the per-pair path.
_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10
_STRIP_FLOOR: float = 1e-4


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

  i = _batched_cell_index(grid_points, u[..., 0], is_linear)  # (N, n)
  j = _batched_cell_index(grid_points, u[..., 1], is_linear)  # (N, n)

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


def interp_at_batched(
  grid_points: Tensor, cache: Tensor, u: Tensor, is_linear: bool = False
) -> Tensor:
  """Bilinearly interpolate a stacked precomputed cache.

  Same shape contract as :func:`interpolate_batched`; the distinction is
  semantic — ``cache`` is a precomputed integral grid (e.g. h-func or cdf),
  whereas ``values`` in :func:`interpolate_batched` is the pdf grid.
  """
  return interpolate_batched(grid_points, cache, u, is_linear)


# --------------------------------------------------------------------------- #
# Batched trapezoidal integration (``cache_integrals=False`` path)             #
# --------------------------------------------------------------------------- #


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

  numer = int_on_grid_batched(grid_points, u_free, strip, is_linear)  # (N, n)
  denom = int_on_grid_batched(
    grid_points, torch.ones_like(u_free), strip, is_linear
  )  # (N, n)
  return (numer / denom).clamp(_TRIM_LO, _TRIM_HI)


def integrate_2d_batched(
  grid_points: Tensor, values: Tensor, u: Tensor, is_linear: bool = False
) -> Tensor:
  """Batched bivariate CDF (trapezoidal-trapezoidal).

  Same shape contract as :func:`integrate_1d_batched`: ``values: (N, m, m)``,
  ``u: (N, n, 2)``, returns ``(N, n)`` clamped to ``[1e-10, 1-1e-10]``.
  The result is renormalised by the full-strip outer integral so
  C(1, u2) = u2 holds exactly — matches the post-vinecopulib#667 C++
  behaviour and stays in parity with the unbatched
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

  # Outer pass: integrate strip[k, l, :] up to u1[k, l] and renormalise
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
  return out.clamp(_TRIM_LO, _TRIM_HI)


# --------------------------------------------------------------------------- #
# Per-tree-level stacked state + per-vine container                            #
# --------------------------------------------------------------------------- #


class BatchedTreeLevel(torch.nn.Module):
  """Stacked state for every pair-copula at one tree level of a vine.

  All buffers are registered so ``.to(device)`` / ``.to(dtype)`` move them
  together with the parent :class:`BatchedVine`.

  Grids (per pair):
  - ``values: (N, m, m)`` — pdf grid (rotation-less; TLL pair-copulas in
    pyvinecopulib always have rotation 0).
  - ``h1_cache, h2_cache, cdf_cache: (N, m, m) | None`` — present only
    when every source pair was constructed with ``cache_integrals=True``.

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
  # statically typed as Tensors instead of nn.Module (cf. ``_cdf_cache``
  # in TorchBicop, same pattern).
  values: Tensor
  h1_cache: Tensor | None
  h2_cache: Tensor | None
  cdf_cache: Tensor | None
  hinv1_cache: Tensor | None
  hinv2_cache: Tensor | None
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
    h1_cache: Tensor | None,
    h2_cache: Tensor | None,
    cdf_cache: Tensor | None,
    hinv1_cache: Tensor | None,
    hinv2_cache: Tensor | None,
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
    if h1_cache is not None:
      self.register_buffer("h1_cache", h1_cache)
      self.register_buffer("h2_cache", h2_cache)
      self.register_buffer("cdf_cache", cdf_cache)
      self.register_buffer("hinv1_cache", hinv1_cache)
      self.register_buffer("hinv2_cache", hinv2_cache)
      self._has_cache = True
    else:
      self.h1_cache = None
      self.h2_cache = None
      self.cdf_cache = None
      self.hinv1_cache = None
      self.hinv2_cache = None
      self._has_cache = False
    self.register_buffer("is_indep", is_indep)
    self.register_buffer("col0_src", col0_src)
    self.register_buffer("col1_src", col1_src)
    self.register_buffer("col1_use_h1", col1_use_h1)
    self.register_buffer("needs_h1", needs_h1)
    self.register_buffer("needs_h2", needs_h2)
    self._is_linear = bool(is_linear)

  @property
  def has_cache(self) -> bool:
    return self._has_cache

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

  def pdf(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair pdf at the stacked queries. Indep slots return 1."""
    raw = interpolate_batched(
      grid_points, self.values, u, self._is_linear
    ).clamp_min(1e-20)
    return torch.where(self.is_indep[:, None], torch.ones_like(raw), raw)

  def hfunc1(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair hfunc1. ``cache=True`` does one bilinear interp on the
    precomputed cache; ``cache=False`` runs :func:`integrate_1d_batched`
    on ``values`` at every call.
    """
    if self.h1_cache is not None:
      raw = interp_at_batched(grid_points, self.h1_cache, u, self._is_linear)
    else:
      raw = integrate_1d_batched(
        grid_points, self.values, u, cond_var=1, is_linear=self._is_linear
      )
    h = raw.clamp(0.0, 1.0)
    return torch.where(
      self.is_indep[:, None], u[..., 1].clamp(_TRIM_LO, _TRIM_HI), h
    )

  def hfunc2(self, grid_points: Tensor, u: Tensor) -> Tensor:
    if self.h2_cache is not None:
      raw = interp_at_batched(grid_points, self.h2_cache, u, self._is_linear)
    else:
      raw = integrate_1d_batched(
        grid_points, self.values, u, cond_var=2, is_linear=self._is_linear
      )
    h = raw.clamp(0.0, 1.0)
    return torch.where(
      self.is_indep[:, None], u[..., 0].clamp(_TRIM_LO, _TRIM_HI), h
    )

  def hinv1(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair hinv1. Requires ``cache_integrals=True`` — the precomputed
    ``hinv1_cache`` collapses the C++ bisection cascade to one bilinear
    interp per call. Without the cache the batched cascade can still
    invert via bisection in the caller; this method raises so the caller
    fails fast rather than silently routing through a slow path."""
    if self.hinv1_cache is None:
      raise RuntimeError(
        "BatchedTreeLevel.hinv1 requires cache_integrals=True; build the "
        "TorchBicop / TorchVinecop with cache_integrals=True to populate "
        "the hinv1 cache."
      )
    raw = interp_at_batched(
      grid_points, self.hinv1_cache, u, self._is_linear
    ).clamp(0.0, 1.0)
    return torch.where(
      self.is_indep[:, None], u[..., 1].clamp(_TRIM_LO, _TRIM_HI), raw
    )

  def hinv2(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Per-pair hinv2. See :meth:`hinv1` for the cache requirement."""
    if self.hinv2_cache is None:
      raise RuntimeError(
        "BatchedTreeLevel.hinv2 requires cache_integrals=True; build the "
        "TorchBicop / TorchVinecop with cache_integrals=True to populate "
        "the hinv2 cache."
      )
    raw = interp_at_batched(
      grid_points, self.hinv2_cache, u, self._is_linear
    ).clamp(0.0, 1.0)
    return torch.where(
      self.is_indep[:, None], u[..., 0].clamp(_TRIM_LO, _TRIM_HI), raw
    )


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
    # Pull a reference tensor to get device.
    device = tvc._ref_tensor().device

    # The grid is shared by all pairs (built once via
    # `InterpolationGrid2D.make_grid_points`).
    grid_points = tvc._pair(0, 0).interp_grid.grid_points
    # The grid type is also shared: all pair-copulas in a TorchVinecop come
    # from the same fit pipeline, so they all use the same storage grid.
    is_linear = bool(tvc._pair(0, 0).interp_grid._is_linear)

    levels: list[BatchedTreeLevel] = []
    for t in range(trunc_lvl):
      N_t = d - t - 1
      vals: list[Tensor] = []
      h1_list: list[Tensor | None] = []
      h2_list: list[Tensor | None] = []
      cdf_list: list[Tensor | None] = []
      hinv1_list: list[Tensor | None] = []
      hinv2_list: list[Tensor | None] = []
      is_indep: list[bool] = []
      col0_src: list[int] = []
      col1_src: list[int] = []
      col1_use_h1: list[bool] = []
      needs_h1_list: list[bool] = []
      needs_h2_list: list[bool] = []
      all_have_cache = True

      for e in range(N_t):
        bc = tvc._pair(t, e)
        m = int(s.min_array(t, e))
        sarr = int(s.struct_array(t, e, natural_order=True))
        vals.append(bc.interp_grid.values)
        if bc._hfunc1_cache is None:
          all_have_cache = False
        h1_list.append(bc._hfunc1_cache)
        h2_list.append(bc._hfunc2_cache)
        cdf_list.append(bc._cdf_cache)
        hinv1_list.append(bc._hinv1_cache)
        hinv2_list.append(bc._hinv2_cache)
        is_indep.append(bool(bc.is_indep))
        col0_src.append(e)
        col1_src.append(m - 1)
        col1_use_h1.append(m != sarr)
        needs_h1_list.append(bool(s.needed_hfunc1(t, e)))
        needs_h2_list.append(bool(s.needed_hfunc2(t, e)))

      values = torch.stack(vals, dim=0).to(device=device)
      h1_cache: Tensor | None
      h2_cache: Tensor | None
      cdf_cache: Tensor | None
      hinv1_cache: Tensor | None
      hinv2_cache: Tensor | None
      if all_have_cache:
        h1_cache = torch.stack(cast(list[Tensor], h1_list), dim=0).to(
          device=device
        )
        h2_cache = torch.stack(cast(list[Tensor], h2_list), dim=0).to(
          device=device
        )
        cdf_cache = torch.stack(cast(list[Tensor], cdf_list), dim=0).to(
          device=device
        )
        hinv1_cache = torch.stack(cast(list[Tensor], hinv1_list), dim=0).to(
          device=device
        )
        hinv2_cache = torch.stack(cast(list[Tensor], hinv2_list), dim=0).to(
          device=device
        )
      else:
        h1_cache = h2_cache = cdf_cache = hinv1_cache = hinv2_cache = None

      level = BatchedTreeLevel(
        values=values,
        h1_cache=h1_cache,
        h2_cache=h2_cache,
        cdf_cache=cdf_cache,
        hinv1_cache=hinv1_cache,
        hinv2_cache=hinv2_cache,
        is_indep=torch.tensor(is_indep, dtype=torch.bool, device=device),
        col0_src=torch.tensor(col0_src, dtype=torch.long, device=device),
        col1_src=torch.tensor(col1_src, dtype=torch.long, device=device),
        col1_use_h1=torch.tensor(col1_use_h1, dtype=torch.bool, device=device),
        needs_h1=torch.tensor(needs_h1_list, dtype=torch.bool, device=device),
        needs_h2=torch.tensor(needs_h2_list, dtype=torch.bool, device=device),
        is_linear=is_linear,
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
