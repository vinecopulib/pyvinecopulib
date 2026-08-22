"""PyTorch evaluator for a fitted R-vine copula.

Wraps a fitted ``Vinecop`` (with *Transformed Local
Likelihood* pair copulas — ``tll``) and
exposes ``pdf`` / ``cdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` /
``sample`` on top of :class:`TorchBicop` for every pair copula. The
whole evaluation chain stays in PyTorch, so the vine can move to GPU
with ``.to("cuda")`` and be composed with autograd-aware downstream
code.

The two routes to a fitted vine are :meth:`TorchVinecop.from_vinecop`
(lift a fitted ``Vinecop``) and
:meth:`TorchVinecop.from_data` (fit directly from pseudo-observations
in pure PyTorch, selecting or reusing a structure). Fit-time controls
live on :class:`FitControlsTorchVinecop`, which mirrors
``FitControlsVinecop``.

The cascade is a direct port of the
``Vinecop`` tree-by-tree h-function chain in
``Vinecop.pdf()`` / ``Vinecop.rosenblatt()`` /
``Vinecop.inverse_rosenblatt()``: dense ``(n, d)`` scratch matrices,
fixed natural-order traversal, byte-for-byte agreement with
``Vinecop`` on the same fit. Every cascade additionally accepts
``batched=True`` to fire one stacked bicop call per group of pair
copulas -- a tree level for ``pdf`` / ``rosenblatt``, a level of the
dependency graph for ``inverse_rosenblatt``.

Variables with atoms are supported throughout — declare them with
``var_types`` and pass the left-limit columns, as ``Vinecop`` takes them. The
pair copulas remain continuous interpolation grids; the mixed-discrete surface
comes from :class:`~pyvinecopulib.core.DiscretePair`.

See Also
--------
pyvinecopulib.core.Vinecop : Reference vine copula.
FitControlsTorchVinecop : Fit-time controls.
"""

from __future__ import annotations

from typing import Any, Optional, cast

import numpy as np
import torch
from torch import Tensor

from ..core import (
  BicopLike,
  ConditioningContext,
  DiscretePair,
  VinecopBase,
)
from ..core._discrete import continuous_view
from ..core.vinecop_base import _NotBatchable
from ..pyvinecopulib_ext import (
  RVineStructure,
  Vinecop,
  indep as _INDEP_FAMILY,
  tll as _TLL_FAMILY,
)
from ..utils import sample_uniform
from ._batched import BatchedVine
from .controls import FitControlsTorchVinecop
from .bicop import TorchBicop


class TorchVinecop(VinecopBase[torch.Tensor], torch.nn.Module):
  """PyTorch R-vine copula evaluator built on ``TorchBicop``.

  Mirrors the public surface of ``Vinecop`` — ``pdf`` /
  ``rosenblatt`` / ``inverse_rosenblatt`` / ``sample`` — but keeps
  the entire evaluation chain in PyTorch so the vine can move to
  GPU with ``.to(device)`` and compose with autograd-aware
  downstream code. Single batch, no threading. Build a ``TorchVinecop`` via
  ``from_data()`` (pure-torch TLL fit) or ``from_vinecop()`` (lifts a fitted
  ``Vinecop``).

  Discrete variables are declared with ``var_types``: the stored pair copulas
  stay continuous grids, and an edge that sees an atom is evaluated through
  :class:`~pyvinecopulib.core.DiscretePair`, which builds the mixed-discrete
  density and h-functions from the grid's own ``pdf`` / ``cdf`` / ``hfunc1`` /
  ``hfunc2``. The batched fast path declines on such a vine.

  The cascade is a dense ``(n, d)``-scratch port of the ``Vinecop``
  evaluator, byte-for-byte parity with ``Vinecop``. Every cascade accepts
  ``batched=True`` to fire a single stacked bicop call per group of pair
  copulas: a tree level for ``pdf`` / ``rosenblatt``, and -- since the
  inverse's dependencies run across trees -- a level of the dependency
  graph for ``inverse_rosenblatt``. The default ``batched=None`` resolves
  to ``True`` on CUDA (3–7x faster) and ``False`` on CPU.

  Parameters
  ----------
  pair_copulas : list of list of TorchBicop
      Indexed ``[tree][edge]`` and shaped like
      ``Vinecop``'s pair-copula layout: tree 0 has ``d - 1`` edges, tree 1
      has ``d - 2``, etc., up to ``trunc_lvl``. May also be passed
      as a `torch.nn.ModuleList` of `torch.nn.ModuleList`.
  structure : RVineStructure
      Vine structure whose accessors (``min_array``,
      ``struct_array``, ``needed_hfunc1`` / ``needed_hfunc2``)
      describe how to walk the trees.
  context : ConditioningContext or None, default=None
      Conditioning-context policy that assembles each pair copula's
      ``x`` per edge. ``None`` uses
      :class:`~pyvinecopulib.core.SimplifiedContext` (an unconditional /
      simplified vine).
  var_types : list of str or None, default=None
      Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete), in
      variable order; ``None`` means all continuous.
  """

  d: int
  trunc_lvl: int
  pair_copulas: torch.nn.ModuleList

  def __init__(
    self,
    pair_copulas: list[list[TorchBicop]],
    structure: RVineStructure,
    *,
    context: Optional[ConditioningContext] = None,
    var_types: Optional[list[str]] = None,
  ) -> None:
    # Initialize nn.Module explicitly: TorchVinecop also subclasses VinecopBase
    # (a Protocol-derived ABC), whose __init__ chain would otherwise shadow
    # nn.Module's under super().
    torch.nn.Module.__init__(self)
    # Install structure + context + variable types + derived order arrays
    # (VinecopBase hook; context=None resolves to SimplifiedContext).
    self._bind_vine(structure, context, var_types=var_types)

    expected_lens = [self.d - 1 - t for t in range(self.trunc_lvl)]
    if len(pair_copulas) != self.trunc_lvl:
      raise ValueError(
        f"pair_copulas has {len(pair_copulas)} trees, expected "
        f"trunc_lvl={self.trunc_lvl}"
      )
    for t, (row, expected) in enumerate(zip(pair_copulas, expected_lens)):
      if len(row) != expected:
        raise ValueError(
          f"pair_copulas tree {t} has {len(row)} edges, expected {expected}"
        )

    self.pair_copulas = torch.nn.ModuleList(
      [torch.nn.ModuleList(list(row)) for row in pair_copulas]
    )
    # `self._batched` (lazy grid-batched state) is initialized to None by
    # `_bind_vine`; `_apply` clears it so device moves re-bake it.
    self._compile_cascades = False
    self._compiled: dict[str, Any] = {}

  # --------------------------------------------------------------------- #
  # Constructor                                                            #
  # --------------------------------------------------------------------- #

  @staticmethod
  def _resolve_cache_integrals(
    cache_integrals: Optional[bool], var_types: list[str]
  ) -> bool:
    """Whether to precompute the prefix tables. ``None`` resolves to ``True``.

    ``var_types`` no longer enters the decision. It used to: a discrete edge
    reads its density from *differences* over an atom's width, and the cache
    was a bilinear interpolation of a table -- accurate enough to evaluate, not
    to difference, at 38% maximum relative error on a ``("d","d")`` density. The
    tables are exact now, so differencing them is no longer the lossy step it
    was and there is nothing left to refuse. A discrete edge still differences
    the distribution function rather than calling ``rect_mass``, which is more
    accurate but would break the torch-to-C++ cascade parity that ``Bicop``'s
    own quotients define.

    Parameters
    ----------
    cache_integrals : bool or None
        What the caller asked for; ``None`` means "whatever suits this vine".
    var_types : list of str
        The variables' types. Retained so callers need not branch on it.

    Returns
    -------
    bool
        The effective setting.
    """
    del var_types
    return True if cache_integrals is None else bool(cache_integrals)

  @classmethod
  def from_vinecop(
    cls,
    cop: Vinecop,
    cache_integrals: Optional[bool] = None,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchVinecop":
    """Lifts a fitted ``Vinecop`` into a ``TorchVinecop``.

    Each pair copula is lifted via ``TorchBicop.from_bicop()``, so the
    resulting cascade matches the ``Vinecop`` evaluator to within
    bilinear-interpolation precision on the shared grid.

    Parameters
    ----------
    cop : Vinecop
        A fitted ``Vinecop`` whose pair copulas are all
        TLL or independence families. Discrete variables are not
        supported.
    cache_integrals : bool or None, default=None
        Forwarded to ``TorchBicop.from_bicop()`` for every pair copula. ``None``
        resolves to ``True``, or to ``False`` when ``cop`` has a discrete
        variable, where an explicit ``True`` raises instead.
    device : torch.device or None, default=None
        Placement of the underlying tensors.
    dtype : torch.dtype, default=torch.float64
        Precision of the underlying tensors.

    Returns
    -------
    TorchVinecop
        A ``TorchVinecop`` mirroring ``cop``.
    """
    var_types = list(cop.var_types)
    cache_integrals = cls._resolve_cache_integrals(cache_integrals, var_types)
    if not cop.pair_copulas:
      # An empty compiled store means independence everywhere, whatever the
      # structure's depth, so route it to the fill rather than to the
      # constructor's shape check.
      return cls.from_structure(
        structure=cop.structure,
        var_types=var_types,
        device=device,
        dtype=dtype,
      )
    pair_copulas_torch: list[list[TorchBicop]] = []
    for tree_idx, row in enumerate(cop.pair_copulas):
      tree_list: list[TorchBicop] = []
      for edge_idx, b in enumerate(row):
        if b.family == _TLL_FAMILY:
          bc = TorchBicop.from_bicop(
            b,
            cache_integrals=cache_integrals,
            device=device,
            dtype=dtype,
          )
        elif b.family == _INDEP_FAMILY:
          bc = TorchBicop(device=device, dtype=dtype)
        else:
          raise ValueError(
            f"TorchVinecop only supports tll and indep pair copulas; "
            f"pair_copulas[{tree_idx}][{edge_idx}] has family={b.family!r}"
          )
        tree_list.append(bc)
      pair_copulas_torch.append(tree_list)

    return cls(
      pair_copulas=pair_copulas_torch,
      structure=cop.structure,
      var_types=var_types,
    )

  @classmethod
  def from_structure(
    cls,
    structure: Optional[RVineStructure] = None,
    matrix: Optional[np.ndarray] = None,
    pair_copulas: list[list[TorchBicop]] = [],
    var_types: list[str] = [],
    *,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchVinecop":
    """Builds a ``TorchVinecop`` from a structure (or matrix) and pair-copulas.

    Parameters
    ----------
    structure : RVineStructure or None, default=None
        Vine structure. Provide either this or ``matrix``.
    matrix : ndarray, shape (d, d), dtype int, or None, default=None
        RVine structure matrix. Provide either this or
        ``structure``.
    pair_copulas : list of list of TorchBicop, default=[]
        Indexed ``[tree][edge]`` with tree ``t`` containing
        ``d - 1 - t`` edges. An empty list fills every edge with
        the independence copula.
    var_types : list of str, default=[]
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete); empty
        means all continuous. A discrete variable makes the cascades read its
        left limit too, and wraps the pair copulas that see it in
        :class:`~pyvinecopulib.core.DiscretePair`.
    device : torch.device or None, default=None
        Placement of the independence pair copulas that fill
        missing edges.
    dtype : torch.dtype, default=torch.float64
        Precision of the underlying tensors.

    Returns
    -------
    TorchVinecop
        A ``TorchVinecop``.
    """
    if (structure is None) == (matrix is None):
      raise ValueError("Provide exactly one of `structure` or `matrix`.")
    if structure is None:
      structure = RVineStructure.from_matrix(np.asarray(matrix))

    d = int(structure.dim)
    if var_types and len(var_types) != d:
      raise ValueError(f"var_types has {len(var_types)} entries, expected {d}")

    trunc_lvl = int(structure.trunc_lvl)
    if not pair_copulas:
      # Independence vine: TorchBicop() defaults to the independence copula.
      pair_copulas = [
        [TorchBicop(device=device, dtype=dtype) for _ in range(d - 1 - t)]
        for t in range(trunc_lvl)
      ]
    return cls(
      pair_copulas=pair_copulas,
      structure=structure,
      var_types=list(var_types) or None,
    )

  @classmethod
  def from_data(
    cls,
    u,
    structure: Optional[RVineStructure] = None,
    controls: Optional[FitControlsTorchVinecop] = None,
    var_types: list[str] = [],
  ) -> "TorchVinecop":
    """Fit a pure-PyTorch TLL vine on ``u``.

    When ``structure`` is ``None`` the R-vine structure and its pair copulas
    are selected and fit natively in torch by
    :meth:`~pyvinecopulib.core.VinecopBase.select` — the array-agnostic
    analog of :class:`~pyvinecopulib.core.Vinecop`'s Dissmann / Wilson
    selection (edge weights use Kendall's tau via ``wdm``).

    When a ``structure`` is given it is fixed and the pair copulas are fit on
    it tree by tree (:meth:`~pyvinecopulib.core.VinecopBase.fit`):
    at each ``(tree, edge)`` the pair of pseudo-obs columns is collected by
    the same rule as the pdf / rosenblatt cascade, a
    :class:`~pyvinecopulib.torch.TorchBicop` is fit via
    :meth:`~pyvinecopulib.torch.TorchBicop.from_data`, and ``hfunc1`` /
    ``hfunc2`` propagate forward.

    Either way the result matches Vinecop's TLL fit to machine precision
    when ``cache_integrals=False``. Structure selection honors
    ``controls.trunc_lvl`` / ``tree_criterion`` / ``threshold`` /
    ``tree_algorithm`` / ``seeds``.

    Parameters
    ----------
    u : ndarray or Tensor, shape (n, d), dtype float
        Pseudo-observations.
    structure : RVineStructure or None, default=None
        Vine skeleton. ``None`` selects one natively in torch (continuous
        variables only).
    controls : FitControlsTorchVinecop or None, default=None
        Pair-copula fit controls bundled with the vine-level structure-selection
        and placement / cascade knobs. ``None`` defaults to TLL on a 30x30
        normal-spaced grid, float64, ``mst_prim`` with ``trunc_lvl=20``.
    var_types : list of str, default=[]
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete); empty
        means all continuous. Given, they also fix the dimension, so ``u`` may
        carry the extra left-limit columns.

    Returns
    -------
    TorchVinecop
        A fitted ``TorchVinecop``.

    Raises
    ------
    ValueError
        If ``u`` is not 2-D, or if ``structure``'s dimension disagrees with the
        one ``var_types`` or ``u`` implies.

    Notes
    -----
    With ``structure=None`` the selected vine — its variable order, R-vine
    matrix encoding, and reused pair copulas — reproduces
    :class:`~pyvinecopulib.core.Vinecop`'s selection on the same data.
    """
    if controls is None:
      controls = FitControlsTorchVinecop()

    eff_dtype = controls.dtype if controls.dtype is not None else torch.float64
    eff_device = controls.device
    u_t = torch.as_tensor(u, dtype=eff_dtype, device=eff_device)
    if u_t.ndim != 2:
      raise ValueError(f"u must be 2-D; got shape {tuple(u_t.shape)}")
    # With left-limit columns present `u` is wider than the vine, so it is
    # `var_types` that fixes the dimension.
    d = len(var_types) if var_types else int(u_t.shape[1])

    bc_controls = controls.bicop_controls
    cache_integrals = cls._resolve_cache_integrals(
      controls.cache_integrals, list(var_types)
    )

    def fit_edge(
      tree: int,
      edge: int,
      u_e: Tensor,
      x_e: Optional[Tensor],
      var_types: Any = ("c", "c"),
    ) -> BicopLike[Tensor]:
      # `var_types` here is *this edge's* two types, which the fit engines pass
      # by that keyword; the vine's own list is the enclosing argument.
      # Simplified (unconditional) TLL fit — x_e is None here.
      del tree, edge, x_e
      bc = TorchBicop.from_data(
        u_e,
        bc_controls,
        cache_integrals=cache_integrals,
        device=u_t.device,
        dtype=eff_dtype,
        var_types=list(var_types),
      )
      # A discrete edge propagates through the mixed-discrete surface, which is
      # also what the next tree's four-column input is built from.
      if "d" not in var_types:
        return bc
      return DiscretePair(bc, (var_types[0], var_types[1]))

    if structure is None:
      # Select the structure natively in torch, reusing the pairs fit during
      # selection (reoriented onto their slots via TorchBicop.flip) — exactly
      # what Vinecop's selector does, so no re-fit is needed. Kendall's tau
      # via wdm needs a host copy; detach so grad-tracking tensors are
      # accepted.
      structure, pairs = cls.select(
        u_t,
        fit_edge,
        trunc_lvl=controls.trunc_lvl,
        tree_criterion=controls.tree_criterion,
        threshold=controls.threshold,
        tree_algorithm=controls.tree_algorithm,
        seeds=list(controls.seeds),
        to_numpy=lambda t: t.detach().cpu().numpy(),
        var_types=list(var_types) or None,
        conditioning_set=list(controls.conditioning_set) or None,
      )
    else:
      if int(structure.dim) != d:
        raise ValueError(
          f"structure.dim={structure.dim} does not match the vine dimension {d}"
        )
      # Fixed structure: fit the pairs tree by tree along it
      # (SimplifiedContext -> x_e=None).
      pairs = cls.fit(
        structure, u_t, fit_edge, var_types=list(var_types) or None
      )
    # Store the continuous grids; `_get_pair_copula` re-wraps a discrete edge,
    # so the ModuleList holds only real nn.Modules.
    modules = [[continuous_view(p) for p in row] for row in pairs]
    out = cls(
      pair_copulas=cast("list[list[TorchBicop]]", modules),
      structure=structure,
      var_types=list(var_types) or None,
    )
    out.compile_cascades = controls.compile
    return out

  # --------------------------------------------------------------------- #
  # Helpers                                                                #
  # --------------------------------------------------------------------- #

  def _pair_module(self, tree: int, edge: int) -> TorchBicop:
    """The stored (always continuous) pair copula at ``(tree, edge)``.

    The two-level :class:`torch.nn.ModuleList` is typed as ``Module`` after
    indexing, so :func:`cast` is used to recover the concrete element type
    for the static checker.
    """
    return cast(
      TorchBicop, cast(torch.nn.ModuleList, self.pair_copulas[tree])[edge]
    )

  def _get_pair_copula(self, tree: int, edge: int) -> BicopLike[Tensor]:
    """The pair copula the cascades evaluate at ``(tree, edge)``.

    The stored modules are continuous grids, so an edge with a discrete
    variable is wrapped in :class:`~pyvinecopulib.core.DiscretePair`, which
    supplies the mixed-discrete surface from the grid's ``pdf`` / ``cdf`` /
    ``hfunc1`` / ``hfunc2``. Keeping the wrapper out of the ``ModuleList``
    keeps ``state_dict`` / ``.to()`` / pickling over the real parameters only.
    """
    pair = self._pair_module(tree, edge)
    types = self.pair_var_types(tree, edge)
    if "d" not in types:
      return pair
    return DiscretePair(pair, types)

  def _ref_tensor(self) -> Tensor:
    """A registered buffer we can crib dtype/device from."""
    # Every TorchBicop registers its interpolation grid; reuse the first.
    return self._pair_module(0, 0).interp_grid.values

  def _default_batched(self) -> bool:
    """Pick a sensible ``batched`` default based on the fitted device.

    Empirical finding from ``scripts/bench_torch_vinecop.py`` (PR
    after #216, d ∈ {5, 10, 20} × n ∈ {200, 1000, 2000, 10000}):

    * On CUDA, ``batched=True`` is **3–7× faster** at every (d, n)
      because per-call kernel-launch overhead dominates the non-batched
      cascade.
    * On CPU torch, ``batched=False`` wins once ``n ≥ 2000`` (the
      regime that actually matters); ``batched=True`` is only faster at
      small ``n`` where overhead dominates either way.

    So we default to ``True`` on CUDA, ``False`` on CPU. Users can
    still override explicitly via ``batched=True`` / ``False``.
    """
    return self._ref_tensor().device.type == "cuda"

  def _prep(self, u: Tensor, name: str, *, values_only: bool = False) -> Tensor:
    # Coerce foreign input onto this vine's dtype/device, then let the base
    # class own the layout: which widths are admissible depends on `var_types`,
    # not on the array namespace.
    ref = self._ref_tensor()
    u = torch.as_tensor(u, dtype=ref.dtype, device=ref.device)
    return super()._prep(u, name, values_only=values_only)

  @property
  def compile_cascades(self) -> bool:
    """Whether the batched cascades run through :func:`torch.compile`.

    Off by default. Set it at any point -- the compilation is lazy, happens
    once per cascade and input shape, and only affects how the same cascade is
    executed, so a vine can be flipped either way between calls.
    :meth:`from_data` sets it from ``controls.compile``.

    Worth it on CUDA, where the eager cascade is bound by kernel-launch count
    rather than by arithmetic: Inductor fuses each tree level's elementwise
    chain into a handful of kernels. Not worth it for a single evaluation --
    the first call at each shape pays tens of seconds of compilation -- and
    the compiled result agrees with the eager one to floating point rather
    than exactly.

    Returns
    -------
    bool
        Whether compilation is enabled.
    """
    return self._compile_cascades

  @compile_cascades.setter
  def compile_cascades(self, value: bool) -> None:
    self._compile_cascades = bool(value)

  def _cascade(self, name: str) -> Any:
    """The named batched cascade, compiled if :attr:`compile_cascades` is set."""
    base = getattr(super(), name)
    if not self._compile_cascades:
      return base
    fn = self._compiled.get(name)
    if fn is None:
      fn = self._compile_cascade(base)
      self._compiled[name] = fn
    return fn

  def _compile_cascade(self, base: Any) -> Any:
    """Compile one cascade, on CUDA through CUDA graphs.

    ``reduce-overhead`` replays the whole cascade as one graph, which is worth
    up to another 2.5x on top of plain compilation where the launch count is
    what binds -- but it writes its result into a buffer the next replay
    reuses, so what comes back is copied out before the caller sees it. That
    copy is one row-length tensor against a cascade of hundreds of kernels.

    Parameters
    ----------
    base : callable
        The uncompiled cascade.

    Returns
    -------
    callable
        The compiled cascade, with the same signature.
    """
    if self._ref_tensor().device.type != "cuda":
      return torch.compile(base, dynamic=False)
    inner = torch.compile(base, dynamic=False, mode="reduce-overhead")

    def graphed(u: Tensor) -> Tensor:
      torch.compiler.cudagraph_mark_step_begin()
      return cast(Tensor, inner(u)).clone()

    return graphed

  def _pdf_batched(self, u: Tensor) -> Tensor:
    return cast(Tensor, self._cascade("_pdf_batched")(u))

  def _rosenblatt_batched(self, u: Tensor) -> Tensor:
    return cast(Tensor, self._cascade("_rosenblatt_batched")(u))

  def _inverse_rosenblatt_batched(self, u: Tensor) -> Tensor:
    return cast(Tensor, self._cascade("_inverse_rosenblatt_batched")(u))

  def __getstate__(self) -> dict:
    """The picklable state: everything except the compiled-cascade cache.

    Compiled callables do not pickle, and they are a pure cache -- the
    unpickled vine recompiles on demand. ``nn.Module`` grew a ``__getstate__``
    in torch 2.1 and ``object`` one in Python 3.11; this project floors at
    torch 2.0 and Python 3.10, so neither can be assumed to exist.

    Returns
    -------
    dict
        The instance state, with an empty compiled-cascade cache.
    """
    inherited = getattr(super(), "__getstate__", None)
    state = dict(self.__dict__ if inherited is None else inherited())
    state["_compiled"] = {}
    return state

  def _apply(self, fn, *args, **kwargs):
    # `.to()`, `.cuda()`, `.cpu()` all route through `_apply`. The
    # BatchedVine container holds buffers — `super()._apply` would move
    # them, but we drop the whole structure so it gets re-baked from the
    # (already-moved) source pair_copulas on next use; that keeps the
    # wire-up tensors aligned with the destination dtype/device.
    self._batched = None
    # Compiled code is specialized on the tensors it was traced with; the
    # guards would recompile anyway, so drop the stale entries.
    self._compiled = {}
    return super()._apply(fn, *args, **kwargs)

  def _grad_signature(self) -> tuple[bool, ...]:
    """Which of the pairs' grids currently track gradients.

    The grids only: they are the tensors a caller is told to flip
    (``values.requires_grad_(True)`` is how a fitted density becomes a learned
    one), and walking every buffer instead costs more for no reachable gain --
    the derived integral caches are built under ``no_grad`` and are documented
    as constants with respect to the grid.

    Returns
    -------
    tuple of bool
        Two entries per pair copula, in tree-then-edge order.
    """
    out: list[bool] = []
    # The same pairs `_build_batched` bakes, in the same order.
    for tree in range(self.trunc_lvl):
      for edge in range(self.d - tree - 1):
        grid = getattr(self._get_pair_copula(tree, edge), "interp_grid", None)
        if grid is None:
          continue
        out.append(bool(grid.values.requires_grad))
        out.append(bool(grid.grid_points.requires_grad))
    return tuple(out)

  def _ensure_batched(self) -> Any:
    """The batched state, re-baked when grad tracking has changed under it.

    The bake copies each pair's grid into a stacked tensor, and
    ``requires_grad_`` mutates a flag in place, so flipping it afterwards leaves
    the copy behind. It does not raise where it is read: the batched cascade
    returns a value detached from the grid, and the error surfaces as torch's
    generic "does not require grad" from ``backward()`` — or, with something else
    in the graph, as a silently missing term. Device moves already invalidate
    through ``_apply``; this is the same invalidation for the other thing a
    caller changes after fitting.

    Returns
    -------
    object
        The memoized batched state.
    """
    signature = self._grad_signature()
    if (
      self._batched is not None
      and getattr(self, "_grad_signature_at_bake", None) != signature
    ):
      object.__setattr__(self, "_batched", None)
      # A compiled cascade was traced against the grids the stale bake holds.
      object.__setattr__(self, "_compiled", {})
    out = super()._ensure_batched()
    object.__setattr__(self, "_grad_signature_at_bake", signature)
    return out

  # --------------------------------------------------------------------- #
  # VinecopBase hooks: RNG for sample + grad control                     #
  # --------------------------------------------------------------------- #

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> Tensor:
    """Draw ``(n, d)`` base uniforms on the fitted grid's dtype/device.

    Pseudo-random via ``torch.rand`` (the first seed seeds a fresh
    ``torch.Generator``), or quasi-random via
    ``pyvinecopulib.utils.sample_uniform`` when ``qrng=True``.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    qrng : bool
        Draw a low-discrepancy (quasi-random) sequence instead.
    seeds : list of int
        RNG seeds (see above).

    Returns
    -------
    Tensor, shape (n, d), dtype float
        Base uniforms in ``[0, 1)``.
    """
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    if qrng:
      u_np = sample_uniform(n, self.d, qrng=True, seeds=list(seeds))
      return torch.as_tensor(u_np, dtype=dtype, device=device)
    gen: Optional[torch.Generator] = None
    if seeds:
      gen = torch.Generator(device=device).manual_seed(int(seeds[0]))
    return torch.rand(n, self.d, generator=gen, dtype=dtype, device=device)

  def _eval_context(self):
    """Disable autograd for ``inverse_rosenblatt`` / ``sample`` / ``cdf``.

    Returns
    -------
    torch.autograd.grad_mode.no_grad
        A ``torch.no_grad()`` context manager.
    """
    return torch.no_grad()

  # ====================================================================== #
  # Batched fast path (`batched=True`)                                       #
  # ====================================================================== #
  #
  # The batched *cascade loops* live on VinecopBase (array-agnostic). This hook
  # supplies the TLL/grid-specific state they run on: a lazily-built BatchedVine
  # (stacked / pre-baked per-tree-level grids + caches).

  def _build_batched(self) -> "BatchedVine":
    """Bake the grid-batched state from this vine's ``TorchBicop`` pairs.

    Returns
    -------
    BatchedVine
        Stacked per-tree-level grids / caches for the batched cascades.

    Raises
    ------
    _NotBatchable
        If any pair lacks the grid/cache internals the batched path needs
        (``supports_batched`` is ``False`` — e.g. a non-``TorchBicop`` nn.Module
        pair), or if any variable is discrete: the batched level has no
        distribution-function lookup, which the mixed-discrete h-functions need.
        The dispatch layer catches it and falls back to the non-batched cascade.
    """
    if self._n_discrete:
      raise _NotBatchable(
        "batched path is continuous-only: the stacked per-level grids carry no "
        "distribution function, which a discrete edge's h-functions are "
        "difference quotients of"
      )
    if not all(
      getattr(self._pair_module(t, e), "supports_batched", False)
      for t in range(self.trunc_lvl)
      for e in range(self.d - t - 1)
    ):
      raise _NotBatchable(
        "batched path requires every pair to expose grid/cache internals "
        "(supports_batched=True); this vine has a non-grid pair copula."
      )
    return BatchedVine.from_torch_vinecop(self)
