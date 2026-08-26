"""One batched cascade over a set of fitted vines.

:class:`BatchedVineEnsemble` evaluates ``M`` vines that share a dimension
and an interpolation grid in a single stacked cascade rather than ``M``.
The shape recurs whenever a search or a resampling scheme produces many
candidate structures on one dataset -- random-search structure selection,
cross-validation over structures, bootstrap ensembles, sensitivity
analyses over a candidate pool -- and evaluating them one at a time
multiplies a launch-bound cascade's cost by ``M`` while the arithmetic per
launch stays tiny.

The stacking is a *reordering*: every pair copula is evaluated on exactly
the values it would have seen alone and no term is regrouped, so the
result agrees with looping over the vines' own batched cascades to the
last bit -- exactly, on x86-64; within one unit in the last place on
arm64 and Windows, where an elementwise kernel may take a different
vectorized path for a stacked tensor than for a single vine's.
"""

from __future__ import annotations

import math
from typing import Any, Optional, Sequence, cast

import torch
from torch import Tensor

from ..core._trim import trim
from ._batched import BatchedTreeLevel, _shared_grid
from .vinecop import TorchVinecop

#: Target size, in elements, of one chunk's per-level working set --
#: ``chunk_size * (d - 1) * n``. The cascade materializes on the order of
#: forty ``(pairs, n)``-shaped temporaries at its peak, so throughput rises
#: with the chunk until that working set stops fitting in cache and then
#: falls again. Measured on one consumer GPU; :attr:`chunk_size` is
#: overridable precisely because the knee belongs to the card.
_CHUNK_TARGET_ELEMENTS: int = 150_000


class _LevelWrites(torch.nn.Module):
  """Which rows of the h-function state one tree level writes.

  The per-vine cascades gate their two writes with a boolean mask and a
  ``where`` against the current value, which writes every column of the
  level and puts the old value back wherever the mask is false. Selecting
  the rows once at bake time and scattering only those writes the same
  values to the same rows, and costs one kernel per level less. The
  indices are unique by construction -- row ``i * d + e`` over a rectangle
  of vines and edges -- which is what makes the scatter well-defined.

  Parameters
  ----------
  needs_h1, needs_h2 : Tensor, shape (pairs,), dtype bool
      Whether a later tree reads each pair's first / second h-function.
  out_rows : Tensor, shape (pairs,), dtype long
      The state row each pair writes, which is also the row it reads its
      first argument from.
  """

  h1_src: Tensor
  h1_dst: Tensor
  h2_src: Tensor
  h2_dst: Tensor

  def __init__(
    self, needs_h1: Tensor, needs_h2: Tensor, out_rows: Tensor
  ) -> None:
    super().__init__()
    h1_src = needs_h1.nonzero().squeeze(-1)
    h2_src = needs_h2.nonzero().squeeze(-1)
    self.register_buffer("h1_src", h1_src)
    self.register_buffer("h1_dst", out_rows.index_select(0, h1_src))
    self.register_buffer("h2_src", h2_src)
    self.register_buffer("h2_dst", out_rows.index_select(0, h2_src))


class _EnsembleChunk(torch.nn.Module):
  """The stacked levels of one chunk of vines, plus its cascades.

  One chunk is one cascade at one set of shapes, which is what keeps a
  compiled ensemble to a single code object no matter how many vines it
  holds.

  Parameters
  ----------
  levels : list of BatchedTreeLevel
      One stacked level per tree, each spanning every vine in the chunk.
  writes : list of _LevelWrites
      The rows each level scatters into, one entry per level.
  seed : Tensor, shape (n_vines * d,), dtype long
      Which column of the input seeds each row of the state, in each
      vine's own sampling order.
  out_idx : Tensor, shape (n_vines * d,), dtype long
      Which row of the state holds each vine's output for each variable.
  n_vines : int
      How many vines this chunk covers.
  d : int
      Their shared dimension.
  """

  seed: Tensor
  out_idx: Tensor

  def __init__(
    self,
    levels: list[BatchedTreeLevel],
    writes: list[_LevelWrites],
    seed: Tensor,
    out_idx: Tensor,
    n_vines: int,
    d: int,
  ) -> None:
    super().__init__()
    self.levels = torch.nn.ModuleList(levels)
    self.writes = torch.nn.ModuleList(writes)
    self.register_buffer("seed", seed)
    self.register_buffer("out_idx", out_idx)
    self.n_vines = n_vines
    self.d = d
    self.trunc_lvl = len(levels)

  def level(self, t: int) -> BatchedTreeLevel:
    return cast(BatchedTreeLevel, self.levels[t])

  def write(self, t: int) -> _LevelWrites:
    return cast(_LevelWrites, self.writes[t])

  def _seed_state(self, u: Tensor) -> Tensor:
    """The ``(n_vines * d, n)`` h-function state, in each vine's order."""
    return u.t().index_select(0, self.seed)

  def pdf(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Every vine's density at the same rows.

    Parameters
    ----------
    grid_points : Tensor, shape (m,)
        The shared interpolation grid.
    u : Tensor, shape (n, d)
        Prepped copula-scale observations.

    Returns
    -------
    Tensor, shape (n_vines, n)
        Row ``i`` is this chunk's ``i``-th vine's density.
    """
    n = u.shape[0]
    hfunc2 = self._seed_state(u)
    # Zeros rather than `empty`: `gather_rows` evaluates both the h1 and
    # the h2 source for `col1` and selects between them, so it reads
    # `hfunc1` at tree 0, where the per-edge cascade never touches it.
    hfunc1 = torch.zeros_like(hfunc2)
    pdf = torch.ones(
      (self.n_vines, n), dtype=hfunc2.dtype, device=hfunc2.device
    )
    for t in range(self.trunc_lvl):
      lvl, w = self.level(t), self.write(t)
      u_e = lvl.gather_rows(hfunc1, hfunc2)
      pdf_e, h1_e, h2_e = lvl.pdf_h1_h2(grid_points, u_e)
      # Vine-major slots, so each vine's block is contiguous and its
      # product runs over the same axis, in the same order, as it would
      # alone. Reducing dim 1 of `(M, N_t, n)` also leaves the same
      # iterator geometry as the per-vine `(N_t, n)` reduce over dim 0,
      # which is what makes this bit-identical rather than merely equal.
      pdf = pdf * pdf_e.reshape(self.n_vines, -1, n).prod(dim=1)
      if w.h1_src.numel():
        hfunc1.index_copy_(0, w.h1_dst, h1_e.index_select(0, w.h1_src))
      if w.h2_src.numel():
        hfunc2.index_copy_(0, w.h2_dst, h2_e.index_select(0, w.h2_src))
    return pdf

  def rosenblatt(self, grid_points: Tensor, u: Tensor) -> Tensor:
    """Every vine's Rosenblatt transform of the same rows.

    Parameters
    ----------
    grid_points : Tensor, shape (m,)
        The shared interpolation grid.
    u : Tensor, shape (n, d)
        Prepped copula-scale observations.

    Returns
    -------
    Tensor, shape (n_vines, n, d)
        Block ``i`` is this chunk's ``i``-th vine's transform.
    """
    n = u.shape[0]
    hfunc2 = self._seed_state(u)
    hfunc1 = hfunc2.clone()
    for t in range(self.trunc_lvl):
      lvl, w = self.level(t), self.write(t)
      u_e = lvl.gather_rows(hfunc1, hfunc2)
      h1_e, h2_e = lvl.h1_h2(grid_points, u_e)
      # Every live edge writes `hfunc2`, gated or not: the transform is
      # read back off those rows.
      hfunc2.index_copy_(0, lvl.col0_src, h2_e)
      if w.h1_src.numel():
        hfunc1.index_copy_(0, w.h1_dst, h1_e.index_select(0, w.h1_src))
    out = hfunc2.index_select(0, self.out_idx)
    return trim(torch, out.view(self.n_vines, self.d, n).transpose(1, 2))


def _longs(values: list[int], device: torch.device) -> Tensor:
  return torch.tensor(values, dtype=torch.long, device=device)


def _bools(values: list[bool], device: torch.device) -> Tensor:
  return torch.tensor(values, dtype=torch.bool, device=device)


def _require(condition: bool, message: str) -> None:
  if not condition:
    raise ValueError(message)


def _check_homogeneous(vines: list[TorchVinecop]) -> tuple[int, int]:
  """Refuse a set that cannot be stacked; return its ``(d, trunc_lvl)``.

  Every one of these is a ``ValueError`` rather than the ``_NotBatchable``
  the per-vine bake raises, because the two mean different things: that
  one is a signal to the dispatch layer, which answers it by running the
  per-edge cascade instead. An ensemble has no such fallback -- looping
  the vines would silently restore both the ``M``-fold launch count and
  the compiled-variant ceiling this class exists to remove -- so a set it
  cannot stack is a caller error, named by index.
  """
  _require(bool(vines), "BatchedVineEnsemble needs at least one vine")
  for k, v in enumerate(vines):
    _require(
      isinstance(v, TorchVinecop),
      f"vines[{k}] has type {type(v).__name__}, not TorchVinecop",
    )
  ref = vines[0]
  d, trunc_lvl = int(ref.d), int(ref.trunc_lvl)
  for k, v in enumerate(vines):
    _require(int(v.d) == d, f"vines[{k}] has dimension {v.d}, vines[0] has {d}")
    _require(
      int(v.trunc_lvl) == trunc_lvl,
      f"vines[{k}] is truncated at {v.trunc_lvl}, vines[0] at {trunc_lvl}; "
      "group the set by trunc_lvl and build one ensemble per group",
    )
    _require(
      not v._n_discrete,
      f"vines[{k}] has discrete variables; the stacked levels carry no "
      "distribution function, which a discrete edge's h-functions are "
      "difference quotients of",
    )
    _require(
      not v._context.assembles_conditioning,
      f"vines[{k}] is non-simplified: its pairs read their conditioning "
      "set, which the stacked cascade does not carry",
    )
    for t in range(trunc_lvl):
      for e in range(d - t - 1):
        pair = v._pair_module(t, e)
        _require(
          getattr(pair, "supports_batched", False),
          f"vines[{k}] pair ({t}, {e}) exposes no grid internals "
          "(supports_batched=False)",
        )
        grid = pair.interp_grid
        _require(
          not (grid.values.requires_grad or grid.grid_points.requires_grad),
          f"vines[{k}] pair ({t}, {e}) tracks gradients on its grid. An "
          "ensemble copies the grids once at construction, so it cannot "
          "see a later change to them; evaluate a learnable vine through "
          "its own cascade",
        )
  return d, trunc_lvl


def _ensemble_grid(vines: list[TorchVinecop], d: int, trunc_lvl: int):
  """The one grid the whole set stacks on, plus its independence pair.

  ``_shared_grid`` establishes this within a vine, where every pair came
  out of one fit and so shares one grid object; across vines it is not
  enough, because it compares only the grid's *shape*. Two vines fitted
  at the same ``grid_size`` on a normal and on a linear grid have grids of
  identical shape and different coordinates, and the ensemble evaluates
  all of them against one. So compare the points themselves, exactly: the
  grid is a deterministic function of ``(grid_size, grid_type, dtype)``,
  which makes anything short of equality a mis-stack rather than a
  rounding difference.
  """
  per_vine = [_shared_grid(v, trunc_lvl, d) for v in vines]
  # A fully independent vine has no real pair, so `_shared_grid` falls back
  # to a 2x2 sentinel; validating that against a fitted grid would refuse a
  # legitimate set. Take the reference from the first vine that has one.
  ref_k = next(
    (
      k
      for k, v in enumerate(vines)
      if any(
        not v._pair_module(t, e).is_indep
        for t in range(trunc_lvl)
        for e in range(d - t - 1)
      )
    ),
    0,
  )
  gp0, lin0 = per_vine[ref_k][0], per_vine[ref_k][1]
  for k, (gp, lin, *_) in enumerate(per_vine):
    if k == ref_k:
      continue
    _require(
      gp.dtype == gp0.dtype,
      f"vines[{k}] is {gp.dtype}, vines[{ref_k}] is {gp0.dtype}",
    )
    _require(
      gp.device == gp0.device,
      f"vines[{k}] is on {gp.device}, vines[{ref_k}] on {gp0.device}",
    )
    _require(
      lin == lin0 and gp.shape == gp0.shape and bool(torch.equal(gp, gp0)),
      f"vines[{k}] disagrees with vines[{ref_k}] on the interpolation "
      "grid; an ensemble evaluates every vine against one grid",
    )
  return per_vine[ref_k]


def _pair_slabs(bc, flat_v: Tensor, flat_sy: Tensor, flat_sx: Tensor):
  """``(values, sy, sy_t)`` for one pair, or the independence substitute.

  ``sy`` is ``None`` when the pair carries no prefix tables, which drops
  its whole level to the on-the-fly integral.
  """
  if bc.is_indep:
    return flat_v, flat_sy, flat_sx.t()
  if bc._sy is None:
    return bc.interp_grid.values, None, None
  sy, sx, _ = bc._tables()
  return bc.interp_grid.values, sy, sx.t()


def _resolve_chunk_size(chunk_size: Optional[int], n_vines: int, d: int) -> int:
  if chunk_size is not None:
    _require(chunk_size >= 1, f"chunk_size must be >= 1; got {chunk_size}")
    return min(int(chunk_size), n_vines)
  # `n` is not known until a call, so size the chunk on the pair axis alone
  # and take a mid-range row count as the reference point.
  per_vine = max(d - 1, 1) * 2000
  target = max(1, _CHUNK_TARGET_ELEMENTS // per_vine)
  return min(n_vines, 1 << int(math.floor(math.log2(target))))


class BatchedVineEnsemble(torch.nn.Module):
  """One batched cascade over a whole set of fitted vines.

  Evaluates ``M`` vines that share a dimension, a truncation level and an
  interpolation grid in one stacked cascade instead of ``M``. Each tree
  level fires a single call over all of that level's pair copulas across
  all the vines, so the kernel-launch count per evaluation stops scaling
  with ``M`` -- and the set is *one* compiled code object, where a Python
  loop over ``M`` vines is ``M`` of them and torch keeps only
  ``torch._dynamo.config.cache_size_limit`` (eight by default) before
  silently falling back to eager.

  The stacking is a reordering, not a reformulation: nothing is regrouped,
  so results agree with looping ``pdf(u, batched=True)`` /
  ``rosenblatt(u, batched=True)`` per vine to the last bit. On x86-64 that
  is exact equality. On arm64 and Windows a few elements can differ by one
  unit in the last place, because the stacked call hands an elementwise
  kernel ``M`` times as many rows and the kernel may vectorize a different
  way at a different size; the difference is bounded by rounding, not by
  the number of vines.

  A snapshot, not a view. The grids are copied into stacked buffers at
  construction, so a vine changed afterwards needs a new ensemble; a vine
  whose grid tracks gradients is refused for that reason, since the copy
  could not see a later change to it. Gradients with respect to ``u`` flow
  as usual.

  ``inverse_rosenblatt``, ``sample`` and ``cdf`` are deliberately absent:
  they are per-vine operations with per-vine randomness, and the inverse
  cascade's dependency graph does not level by tree.

  Parameters
  ----------
  vines : sequence of TorchVinecop
      Fitted vines agreeing on dimension, truncation level, dtype, device,
      grid points and grid spacing, on having no discrete variable and no
      conditioning context, and on whether their pairs carry cached
      integrals. Structures may differ freely -- that is the point.
  chunk_size : int, optional
      How many vines one stacked cascade covers. ``None`` picks a size
      from the dimension; a larger chunk saves launches and a smaller one
      keeps the per-level working set in cache, so throughput peaks in
      between. Every full chunk shares one set of shapes, which is what
      keeps a compiled ensemble to one or two variants regardless of
      ``M``.

  Raises
  ------
  ValueError
      If ``vines`` is empty, or any vine disagrees with the first on the
      properties above, is discrete or conditional, hosts a pair without
      grid internals, or tracks gradients on a grid. The message names the
      offending index.

  See Also
  --------
  TorchVinecop : The single-vine evaluator these are built from.

  Examples
  --------
  See ``examples/09_torch_backend.ipynb`` for a worked ensemble.
  """

  grid_points: Tensor

  def __init__(
    self,
    vines: Sequence[TorchVinecop],
    chunk_size: Optional[int] = None,
  ) -> None:
    super().__init__()
    vine_list = list(vines)
    d, trunc_lvl = _check_homogeneous(vine_list)
    gp, is_linear, flat_v, flat_sy, flat_sx = _ensemble_grid(
      vine_list, d, trunc_lvl
    )
    self.d = d
    self.trunc_lvl = trunc_lvl
    self.structures = tuple(v.structure for v in vine_list)
    self.register_buffer("grid_points", gp)
    self._compile_cascades = False
    self._compiled: dict[str, Any] = {}

    size = _resolve_chunk_size(chunk_size, len(vine_list), d)
    self.chunk_size = size
    self.chunks = torch.nn.ModuleList(
      [
        _build_chunk(
          vine_list[i : i + size],
          d,
          trunc_lvl,
          is_linear,
          flat_v,
          flat_sy,
          flat_sx,
        )
        for i in range(0, len(vine_list), size)
      ]
    )

  # --- introspection ----------------------------------------------------

  @property
  def n_vines(self) -> int:
    """How many vines the ensemble holds.

    Returns
    -------
    int
        The size of the set this was built from.
    """
    return len(self.structures)

  def __len__(self) -> int:
    return self.n_vines

  def extra_repr(self) -> str:
    return (
      f"n_vines={self.n_vines}, d={self.d}, trunc_lvl={self.trunc_lvl}, "
      f"chunk_size={self.chunk_size}"
    )

  # --- compilation ------------------------------------------------------

  @property
  def compile_cascades(self) -> bool:
    """Whether the stacked cascades run through :func:`torch.compile`.

    Off by default: compilation is worth tens of seconds of Inductor at
    the first call on each input shape, which pays for itself in a
    cascade called repeatedly and not in a single evaluation. Set it at
    any point -- compilation is lazy and per cascade -- and note that the
    whole ensemble is one code object, so the per-code-object cap on
    compiled variants applies once here rather than once per vine.

    Returns
    -------
    bool
        Whether compilation is currently switched on.
    """
    return self._compile_cascades

  @compile_cascades.setter
  def compile_cascades(self, value: bool) -> None:
    self._compile_cascades = bool(value)

  def _compile(self, base: Any) -> Any:
    """Compile one chunk cascade, on CUDA through CUDA graphs.

    Kept in step with ``TorchVinecop._compile_cascade``, which does the
    same thing for a single vine: ``reduce-overhead`` replays the cascade
    as one graph, and because the replay writes into a buffer the next one
    reuses, the result is copied out before the caller sees it.
    """
    if self.grid_points.device.type != "cuda":
      return torch.compile(base, dynamic=False)
    inner = torch.compile(base, dynamic=False, mode="reduce-overhead")

    def graphed(chunk: Any, gp: Tensor, u: Tensor) -> Tensor:
      torch.compiler.cudagraph_mark_step_begin()
      return cast(Tensor, inner(chunk, gp, u)).clone()

    return graphed

  def _cascade(self, name: str) -> Any:
    base = getattr(_EnsembleChunk, name)
    if not self._compile_cascades:
      return base
    fn = self._compiled.get(name)
    if fn is None:
      fn = self._compile(base)
      self._compiled[name] = fn
    return fn

  # --- evaluation -------------------------------------------------------

  def _prep(self, u: Any) -> Tensor:
    """Coerce and trim one input block, exactly as a vine would."""
    ref = self.grid_points
    u_t = torch.as_tensor(u, dtype=ref.dtype, device=ref.device)
    if u_t.ndim != 2 or int(u_t.shape[1]) != self.d:
      raise ValueError(
        f"u must have shape (n, {self.d}); got {tuple(u_t.shape)}"
      )
    return trim(torch, u_t)

  def _run(self, name: str, u: Any, cat_dim: int) -> Tensor:
    u_p = self._prep(u)
    fn = self._cascade(name)
    parts = [fn(chunk, self.grid_points, u_p) for chunk in self.chunks]
    return parts[0] if len(parts) == 1 else torch.cat(parts, dim=cat_dim)

  def pdf(self, u: Any) -> Tensor:
    """Every vine's density at the same rows.

    Parameters
    ----------
    u : Tensor, shape (n, d)
        Copula-scale observations, shared by every vine in the set. The
        compact layout only: a single vine also accepts the expanded
        ``(n, 2d)`` form and drops its left-limit block, which an ensemble
        -- being continuous by construction -- asks the caller to slice
        instead.

    Returns
    -------
    Tensor, shape (M, n)
        Row ``k`` is ``vines[k]``'s density.
    """
    return self._run("pdf", u, 0)

  def rosenblatt(self, u: Any) -> Tensor:
    """Every vine's Rosenblatt transform of the same rows.

    Parameters
    ----------
    u : Tensor, shape (n, d)
        Copula-scale observations, shared by every vine in the set.

    Returns
    -------
    Tensor, shape (M, n, d)
        Block ``k`` is ``vines[k]``'s transform.
    """
    return self._run("rosenblatt", u, 0)

  # --- torch plumbing ---------------------------------------------------

  def _apply(self, fn, *args, **kwargs):
    # Compiled code is specialized on the tensors it was traced with, so a
    # device or dtype move invalidates it. The stacked buffers themselves
    # move with the module: unlike a vine's lazily-derived bake, they are
    # the ensemble's own state and there is nothing cheaper to rebuild
    # them from.
    self._compiled = {}
    return super()._apply(fn, *args, **kwargs)

  def __getstate__(self) -> dict:
    """The picklable state: everything except the compiled-cascade cache.

    Compiled callables do not pickle, and they are a pure cache -- the
    unpickled ensemble recompiles on demand. Delegates to
    ``nn.Module.__getstate__`` rather than copying ``__dict__``, which is
    what drops ``_compiled_call_impl`` as well.

    Returns
    -------
    dict
        The instance state, with an empty compiled-cascade cache.
    """
    state = dict(super().__getstate__())
    state["_compiled"] = {}
    return state


def _build_chunk(
  vines: list[TorchVinecop],
  d: int,
  trunc_lvl: int,
  is_linear: bool,
  flat_v: Tensor,
  flat_sy: Tensor,
  flat_sx: Tensor,
) -> _EnsembleChunk:
  """Stack one chunk's pair copulas, tree level by tree level.

  Slots run vine-major -- vine ``i``'s edge ``e`` at ``i * N_t + e`` -- and
  the wiring vectors become absolute row indices into an
  ``(n_vines * d, n)`` state, ``i * d + <column>``. That ordering is not
  cosmetic: it is what leaves each vine's block of the level contiguous
  and its density product reducing over the same axis, in the same order,
  as it would on its own.
  """
  device = flat_v.device
  m_vines = len(vines)
  levels: list[BatchedTreeLevel] = []
  writes: list[_LevelWrites] = []

  for t in range(trunc_lvl):
    n_t = d - t - 1
    vals: list[Tensor] = []
    sy_list: list[Tensor | None] = []
    sy_t_list: list[Tensor | None] = []
    is_indep: list[bool] = []
    col0_src: list[int] = []
    col1_src: list[int] = []
    col1_use_h1: list[bool] = []
    needs_h1: list[bool] = []
    needs_h2: list[bool] = []
    have_cache: list[bool] = []

    for i, v in enumerate(vines):
      s = v.structure
      for e in range(n_t):
        bc = v._pair_module(t, e)
        value, sy, sy_t = _pair_slabs(bc, flat_v, flat_sy, flat_sx)
        vals.append(value)
        sy_list.append(sy)
        sy_t_list.append(sy_t)
        have_cache.append(sy is not None)
        mn = int(s.min_array(t, e))
        sarr = int(s.struct_array(t, e, natural_order=True))
        is_indep.append(bool(bc.is_indep))
        col0_src.append(i * d + e)
        col1_src.append(i * d + (mn - 1))
        col1_use_h1.append(mn != sarr)
        needs_h1.append(bool(s.needed_hfunc1(t, e)))
        needs_h2.append(bool(s.needed_hfunc2(t, e)))

    # All or nothing per level, as the per-vine bake is: the cached and
    # the on-the-fly integrals agree to rounding rather than exactly, so
    # letting one uncached vine drag the others onto the quadrature path
    # would break their parity against their own cascades.
    _require(
      all(have_cache) or not any(have_cache),
      f"the set disagrees on cached integrals at tree level {t}: "
      f"{sum(have_cache)} of {len(have_cache)} pairs carry prefix tables. "
      "Build the vines with one cache_integrals setting",
    )
    cached = all(have_cache)

    out_rows = _longs(col0_src, device)
    h1_mask = _bools(needs_h1, device)
    h2_mask = _bools(needs_h2, device)
    levels.append(
      BatchedTreeLevel(
        values=torch.stack(vals, dim=0).to(device=device),
        sy=(
          torch.stack(cast("list[Tensor]", sy_list), 0).to(device=device)
          if cached
          else None
        ),
        sy_t=(
          torch.stack(cast("list[Tensor]", sy_t_list), 0).to(device=device)
          if cached
          else None
        ),
        is_indep=_bools(is_indep, device),
        col0_src=out_rows,
        col1_src=_longs(col1_src, device),
        col1_use_h1=_bools(col1_use_h1, device),
        needs_h1=h1_mask,
        needs_h2=h2_mask,
        is_linear=is_linear,
      )
    )
    writes.append(_LevelWrites(h1_mask, h2_mask, out_rows))

  seed = _longs([int(v.order[j]) - 1 for v in vines for j in range(d)], device)
  out_idx = _longs(
    [
      i * d + int(v.inverse_order[j])
      for i, v in enumerate(vines)
      for j in range(d)
    ],
    device,
  )
  return _EnsembleChunk(levels, writes, seed, out_idx, m_vines, d)


__all__ = ["BatchedVineEnsemble"]
