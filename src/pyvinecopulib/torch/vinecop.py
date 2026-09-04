"""PyTorch R-vine copula on ``TorchBicop`` pair copulas.

``TorchVinecop`` is the PyTorch member of the ``VinecopBase`` family: it hosts
one pair copula per edge and inherits the whole evaluator surface -- ``pdf`` /
``cdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` / ``sample``, ``loglik`` and
``plot`` -- together with the estimator surface ``fit`` / ``select`` /
``from_data``. Every step stays in PyTorch, so a vine moves to a GPU with
``.to("cuda")``, differentiates under autograd, and composes with other torch
models; on the same fit it agrees with ``Vinecop`` to floating-point tolerance.

Routes to a vine: ``TorchVinecop.from_data()`` fits one on
pseudo-observations, selecting a structure unless handed one;
``TorchVinecop.from_vinecop()`` lifts a fitted ``Vinecop``;
``TorchVinecop.from_structure()`` assembles one from a structure and pair
copulas; and the inherited ``fit`` / ``select`` re-estimate a vine in place.
Fit configuration travels as a ``FitControlsTorchVinecop``.

Two evaluation choices are not part of the model: the batched fast path, which
fires one stacked pair-copula call per group of edges and is resolved per call,
and ``compile_cascades``, which decides whether those batched cascades run
through :func:`torch.compile`.

Variables with atoms are declared with ``var_types`` and passed in the layouts
``Vinecop`` takes. The stored pair copulas remain continuous interpolation
grids; the mixed-discrete surface comes from
:class:`~pyvinecopulib.core.DiscretePair`.

See Also
--------
pyvinecopulib.core.Vinecop : Reference vine copula.
pyvinecopulib.core.VinecopBase : The array-agnostic base.
TorchBicop : The pair copulas this hosts.
FitControlsTorchVinecop : Fit-time controls.
"""

from __future__ import annotations

from typing import Any, ClassVar, Optional, Sequence, cast

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
from ..core._independence import IndependencePair
from ..core.vinecop_base import FitEdge, _NotBatchable
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
  """PyTorch R-vine copula on ``TorchBicop`` pair copulas.

  A ``VinecopBase`` whose pair copulas are density grids and whose cascades --
  ``pdf`` / ``cdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` / ``sample`` --
  stay in PyTorch, so they run on the device the vine sits on, carry
  gradients, and agree with ``Vinecop`` on the same fit to floating-point
  tolerance. Also a ``torch.nn.Module``, so ``.to(device)``, ``state_dict``
  and pickling reach the whole vine. Input is coerced onto the vine's own
  dtype and device, so a caller may hand over the array type at hand;
  ``num_threads`` is accepted for parity with ``Vinecop`` and ignored.

  The members that make it concrete, each resting on the one before it:

  - ``get_pair_copula(tree, edge)`` returns the pair copula the cascades
    evaluate at a position. The stored modules are continuous grids, so an
    edge with a discrete variable comes back wrapped in
    :class:`~pyvinecopulib.core.DiscretePair`. The wrapper itself is not
    stored, which is what keeps ``state_dict`` / ``.to()`` / pickling over
    real ``torch.nn.Module`` parameters only.
  - ``set_pair_copulas(pair_copulas)`` stores fitted pairs, which is what
    lets the inherited ``fit`` and ``select`` install what they fitted and
    hand back ``self``.
  - ``bicop_class`` is ``TorchBicop``. Naming it is what lets
    ``TorchVinecop.from_data()`` fit with no pair-fitting callback, and lets
    structure selection -- which runs through ``VinecopBase.select()``, on
    tensors rather than through a ``Vinecop`` -- refuse a pair copula it could
    not reorient before it reads any data.

  Two further choices sit alongside the fitted model. The batched fast path
  fires a single stacked pair-copula call per group of edges: a tree level for
  ``pdf`` / ``rosenblatt``, and -- the inverse's dependencies running across
  trees -- a level of the dependency graph for ``inverse_rosenblatt``. It is
  not a control: ``batched=None`` is resolved from the vine's device on every
  call, to ``True`` on CUDA and ``False`` elsewhere, and any call may name it
  explicitly. A vine with a discrete variable declines the batched path, its
  stacked per-level grids carrying no distribution function; it does not
  decline the integral cache, which reconstructs the integral exactly.
  ``compile_cascades`` is the other: whether the batched cascades run through
  :func:`torch.compile`.

  Parameters
  ----------
  pair_copulas : list of list of TorchBicop
      The pair copulas, indexed ``[tree][edge]`` and shaped as ``Vinecop``
      lays them out: tree ``t`` holds ``d - 1 - t`` edges, up to the
      structure's ``trunc_lvl``. A ``torch.nn.ModuleList`` of
      ``torch.nn.ModuleList`` is accepted too.
  structure : RVineStructure
      The vine structure to evaluate along.
  context : ConditioningContext or None, default=None
      Per-edge policy assembling each pair copula's ``x``. ``None`` uses
      :class:`~pyvinecopulib.core.SimplifiedContext` -- an unconditional,
      simplified vine.
  var_types : list of str or None, default=None
      Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete), in
      variable order; ``None`` means all continuous.

  Raises
  ------
  ValueError
      If ``pair_copulas`` does not have one row per tree, or a row does not
      have one entry per edge of that tree.
  """

  d: int
  trunc_lvl: int
  pair_copulas: torch.nn.ModuleList

  # The pair copula this vine fits, so `from_data` needs no callback and
  # selection can check `flip` before reading the data.
  bicop_class: ClassVar[Optional[type]] = TorchBicop

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
    # A vine truncated at zero has no pair copulas, so there is no grid to read
    # dtype/device from. This buffer always exists and `.to()` moves it, which
    # is what makes an independence vine evaluable at all.
    ref = (
      pair_copulas[0][0].interp_grid.values
      if self.trunc_lvl > 0
      else torch.empty(0)
    )
    self.register_buffer(
      "_device_ref", torch.empty(0, dtype=ref.dtype, device=ref.device)
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

    ``var_types`` does not enter the decision: the prefix tables reconstruct
    the integral exactly rather than approximately, so an edge that reads its
    density from differences over an atom's width can difference them safely.
    Such an edge still differences the distribution function rather than
    calling ``rect_mass``, which is more accurate but would break the cascade
    parity with ``Vinecop`` that its own difference quotients define.

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
    """Lift a fitted ``Vinecop`` into a ``TorchVinecop``.

    The result hosts one ``TorchBicop`` per pair copula, on the same grids, so
    it agrees with ``cop`` to floating-point tolerance.

    Parameters
    ----------
    cop : Vinecop
        A fitted vine whose pair copulas are all of the ``tll`` or ``indep``
        family, in any continuous, discrete or mixed variable layout.
    cache_integrals : bool or None, default=None
        Whether each pair copula precomputes the prefix tables its ``cdf`` and
        h-functions read. ``None`` resolves to ``True`` for every variable
        layout.
    device : torch.device or None, default=None
        Placement of the underlying tensors.
    dtype : torch.dtype, default=torch.float64
        Precision of the underlying tensors.

    Returns
    -------
    TorchVinecop
        A ``TorchVinecop`` mirroring ``cop``.

    Raises
    ------
    ValueError
        If a pair copula of ``cop`` is of neither family.
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
    """Build a ``TorchVinecop`` from a structure and pair copulas.

    Parameters
    ----------
    structure : RVineStructure or None, default=None
        The vine structure. Provide either this or ``matrix``.
    matrix : ndarray, shape (d, d), dtype int, or None, default=None
        R-vine structure matrix. Provide either this or ``structure``.
    pair_copulas : list of list of TorchBicop, default=[]
        The pair copulas, indexed ``[tree][edge]`` with tree ``t`` holding
        ``d - 1 - t`` edges. Empty fills every edge with the independence
        copula.
    var_types : list of str, default=[]
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete), in
        variable order; empty means all continuous. A discrete variable makes
        the cascades read its left limit too, and the pair copulas that see it
        are evaluated through :class:`~pyvinecopulib.core.DiscretePair`.
    device : torch.device or None, default=None
        Placement of the independence pair copulas that fill missing edges.
    dtype : torch.dtype, default=torch.float64
        Precision of the independence pair copulas that fill missing edges.

    Returns
    -------
    TorchVinecop
        A vine on that structure.

    Raises
    ------
    ValueError
        If neither or both of ``structure`` and ``matrix`` are given, or if
        ``var_types`` is non-empty and does not have one entry per variable.
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
      var_types=list(var_types or []) or None,
    )

  @classmethod
  def from_data(
    cls,
    u,
    /,
    # `structure` and `controls` are typed `Any` so this stays a widening of
    # `VinecopBase.from_data`, which accepts any `ControlsLike`: narrowing a
    # parameter is what a subclass may not do. The accepted objects are an
    # `RVineStructure` and a `FitControlsTorchVinecop`, as documented below.
    structure: Optional[Any] = None,
    var_types: Optional[list[str]] = None,
    controls: Optional[Any] = None,
    *,
    fit_edge: Optional[FitEdge] = None,
  ) -> "TorchVinecop":
    """Fit a vine to pseudo-observations, in PyTorch throughout.

    The factory counterpart of ``select``, for when there is no vine yet. With
    ``structure=None`` the structure is chosen from the data and comes back
    with the pair copulas selected along it, honoring the ``trunc_lvl`` /
    ``tree_criterion`` / ``threshold`` / ``tree_algorithm`` / ``seeds`` /
    ``conditioning_set`` settings on ``controls``. A supplied ``structure`` is
    taken as given and only its pair copulas are fitted, ``threshold`` still
    leaving an edge below it independent. Either way the result reproduces a
    ``Vinecop`` TLL fit on the same data: the selected structure down to its
    matrix encoding, the density to floating-point tolerance.

    Parameters
    ----------
    u : ndarray or Tensor, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]``.
    structure : RVineStructure or None, default=None
        A fixed structure. Selected from the data when ``None``.
    var_types : list of str, or None, default=None
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete);
        ``None`` means all continuous. Given, they also fix the dimension, so
        ``u`` may carry the extra left-limit columns.
    controls : FitControlsTorchVinecop or None, default=None
        Fit configuration for both halves of the fit -- the structure
        selection the vine runs and the pair-copula fits its edges run -- plus
        placement, precision, and the cascade variants. ``None`` defaults to
        TLL on a 30x30 normal-spaced grid, float64, and ``mst_prim`` with
        ``trunc_lvl=20``.
    fit_edge : callable, or None, default=None
        ``(tree, edge, u_e, x_e) -> BicopLike``, fitting one edge's pair
        copula in place of the built-in TLL fit; an edge with a discrete
        variable additionally receives its own ``var_types`` by keyword.
        ``None`` fits a ``TorchBicop`` per edge, which is the usual case.

    Returns
    -------
    TorchVinecop
        The fitted vine.

    Raises
    ------
    ValueError
        If ``u`` is not 2-D, or if ``structure``'s dimension disagrees with the
        one ``var_types`` or ``u`` implies.

    See Also
    --------
    pyvinecopulib.core.VinecopBase.fit : Refit an existing vine's pairs.
    pyvinecopulib.core.VinecopBase.select : Reselect its structure and pairs.
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

    # The vine's controls are pair controls: `FitControlsTorchVinecop`
    # derives from `FitControlsTorchBicop`, as its core counterparts do.
    bc_controls = controls
    cache_integrals = cls._resolve_cache_integrals(
      controls.cache_integrals, list(var_types or [])
    )

    def fit_edge_tll(
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

    if fit_edge is not None:
      # A caller who brings their own pair fitter overrides the built-in TLL
      # one; this is also what makes the signature a widening of
      # `VinecopBase.from_data`, whose `fit_edge` is required.
      pair_fitter: FitEdge = fit_edge
    else:
      pair_fitter = fit_edge_tll

    # `None` resolves per device, as the evaluation cascade's `batched` does:
    # the per-level fitter buys launch amortization, which cpu has none of.
    batched_fit = controls.batched_fit
    if batched_fit is None:
      batched_fit = u_t.device.type == "cuda"

    def fit_level(
      tree: int, u_level: Tensor, types: list[tuple[str, str]]
    ) -> Sequence[BicopLike[Tensor]]:
      """Fit a whole continuous tree level in one call.

      Parameters
      ----------
      tree : int
          Tree index; unused, the fit needing no structural context.
      u_level : Tensor, shape (N_t, n, 2)
          The level's edges, stacked in ascending edge order.
      types : list of tuple of str
          Each edge's variable types; unused, a level reaching here being
          continuous throughout.

      Returns
      -------
      sequence of BicopLike
          One fitted pair copula per edge, in the same order.
      """
      del tree, types  # a level reaching here is continuous and simplified
      return TorchBicop.from_data_batched(
        u_level,
        bc_controls,
        cache_integrals=cache_integrals,
        device=u_t.device,
        dtype=eff_dtype,
      )

    level_hook = fit_level if batched_fit else None

    if structure is None:
      # Select the structure natively in torch, reusing the pairs fit during
      # selection (reoriented onto their slots via TorchBicop.flip) — exactly
      # what Vinecop's selector does, so no re-fit is needed. Kendall's tau
      # via wdm needs a host copy; detach so grad-tracking tensors are
      # accepted.
      structure, pairs = cls._select_parts(
        u_t,
        pair_fitter,
        fit_level=level_hook,
        trunc_lvl=controls.trunc_lvl,
        tree_criterion=controls.tree_criterion,
        threshold=controls.threshold,
        tree_algorithm=controls.tree_algorithm,
        seeds=list(controls.seeds),
        var_types=list(var_types or []) or None,
        conditioning_set=list(controls.conditioning_set) or None,
      )
    else:
      if int(structure.dim) != d:
        raise ValueError(
          f"structure.dim={structure.dim} does not match the vine dimension {d}"
        )
      # Fixed structure: fit the pairs tree by tree along it
      # (SimplifiedContext -> x_e=None).
      pairs = cls._fit_parts(
        structure,
        u_t,
        pair_fitter,
        var_types=list(var_types or []) or None,
        fit_level=level_hook,
        tree_criterion=controls.tree_criterion,
        threshold=controls.threshold,
      )
    # Store the continuous grids; `get_pair_copula` re-wraps a discrete edge,
    # so the ModuleList holds only real nn.Modules. A thresholded edge arrives
    # from `select` as a `core.IndependencePair`, which is not one -- it becomes
    # the grid that *is* independence, whose `pdf` is exactly 1 and whose
    # h-functions are exactly the identity, so it stores, moves and pickles like
    # any other pair.
    modules = [
      [
        # A thresholded edge arrives from the engines as a
        # `core.IndependencePair`, which is not an `nn.Module`. The
        # no-argument `TorchBicop` *is* the independence copula -- a 2x2
        # sentinel that short-circuits every method on `is_indep`, exactly
        # rather than to rounding -- so it needs no grid of its own and
        # cannot disagree with its siblings about one. `u_t.device`, as
        # `fit_edge` uses: `controls.device` is `None` whenever the caller
        # let the data carry the placement.
        TorchBicop(device=u_t.device, dtype=eff_dtype)
        if isinstance(p, IndependencePair)
        else continuous_view(p)
        for p in row
      ]
      for row in pairs
    ]
    out = cls(
      pair_copulas=cast("list[list[TorchBicop]]", modules),
      structure=structure,
      var_types=list(var_types or []) or None,
    )
    out.compile_cascades = controls.compile
    return out

  # --------------------------------------------------------------------- #
  # Helpers                                                                #
  # --------------------------------------------------------------------- #

  def set_pair_copulas(self, pair_copulas: "list[list[Any]]") -> None:
    """Install fitted pairs, so ``fit`` and ``select`` can return ``self``.

    Parameters
    ----------
    pair_copulas : list of list of BicopLike
        Fitted pairs indexed ``[tree][edge]``, as the fit engines return them.
        A pair carrying discrete variables is stored as its continuous grid
        and re-wrapped on read; an edge left independent by ``threshold`` is
        stored as the independence ``TorchBicop``.

    Returns
    -------
    None
    """
    ref = cast("Tensor", self._buffers["_device_ref"])
    # Only real `nn.Module`s go in the `ModuleList`: a discrete edge is
    # re-wrapped on read by `get_pair_copula`, and a thresholded edge arrives
    # as a `core.IndependencePair`, which the no-argument `TorchBicop` -- the
    # independence copula exactly rather than to rounding -- stands in for.
    self.pair_copulas = torch.nn.ModuleList(
      [
        torch.nn.ModuleList(
          [
            TorchBicop(device=ref.device, dtype=ref.dtype)
            if isinstance(pair, IndependencePair)
            else continuous_view(pair)
            for pair in row
          ]
        )
        for row in pair_copulas
      ]
    )

  def _pair_module(self, tree: int, edge: int) -> TorchBicop:
    """The stored (always continuous) pair copula at ``(tree, edge)``."""
    return cast(
      TorchBicop, cast(torch.nn.ModuleList, self.pair_copulas[tree])[edge]
    )

  def get_pair_copula(self, tree: int, edge: int) -> BicopLike[Tensor]:
    """The pair copula the cascades evaluate at ``(tree, edge)``.

    The stored modules are continuous grids, so an edge with a discrete
    variable comes back wrapped in
    :class:`~pyvinecopulib.core.DiscretePair`, which supplies the
    mixed-discrete surface. Keeping the wrapper out of the stored
    ``torch.nn.ModuleList`` is what keeps ``state_dict`` / ``.to()`` /
    pickling over the real parameters only.

    Parameters
    ----------
    tree : int
        Tree index (``0``-based).
    edge : int
        Edge index within the tree (``0``-based).

    Returns
    -------
    BicopLike
        The pair copula hosted at that position, wrapped for a discrete edge.
    """
    pair = self._pair_module(tree, edge)
    types = self.pair_var_types(tree, edge)
    if "d" not in types:
      return pair
    return DiscretePair(pair, types)

  def get_extra_state(self) -> dict[str, Any]:
    """Return the non-tensor model identity for ``state_dict`` round-trips.

    Returns
    -------
    dict
        The structure and variable types, which no tensor carries.
    """
    return {
      "version": 1,
      "structure": self.structure.to_json(),
      "var_types": list(self.var_types),
    }

  def set_extra_state(self, state: Any) -> None:
    """Check the non-tensor model identity in a ``state_dict``.

    The pair-copula grids are the only tensors a ``state_dict`` holds, so
    nothing in it identifies the vine they hang on: a checkpoint from a
    differently structured vine carries exactly the same keys, and is refused
    here rather than loaded into a model it does not describe.

    Parameters
    ----------
    state : dict
        State returned by ``TorchVinecop.get_extra_state()``.

    Raises
    ------
    RuntimeError
        If the checkpoint's structure or variable types differ from this
        module's, or its version is not recognized.
    """
    if not isinstance(state, dict) or state.get("version") != 1:
      raise RuntimeError("unsupported TorchVinecop state-dict version")
    if state["structure"] != self.structure.to_json():
      raise RuntimeError(
        "state_dict was saved from a vine with a different structure; "
        "rebuild the module from that structure before loading it"
      )
    if list(state["var_types"]) != list(self.var_types):
      raise RuntimeError(
        "state_dict was saved from a vine with different var_types "
        f"({list(state['var_types'])} vs {list(self.var_types)})"
      )

  def _ref_tensor(self) -> Tensor:
    """A registered buffer to read dtype and device from."""
    # Every TorchBicop registers its interpolation grid, so prefer the first --
    # but a vine truncated at zero has none, and falls back to the buffer the
    # constructor registers for exactly that case.
    if self.trunc_lvl > 0:
      return self._pair_module(0, 0).interp_grid.values
    ref = self._buffers["_device_ref"]
    assert ref is not None
    return ref

  def _default_batched(self) -> bool:
    """Whether ``batched`` defaults to ``True``, read from the vine's device.

    Returns
    -------
    bool
        ``True`` on CUDA, where the batched cascade is markedly faster because
        per-call kernel-launch overhead dominates the edge-at-a-time one, and
        ``False`` elsewhere, where it wins only at sample sizes small enough
        for the overhead to dominate either way.
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

    Off by default, and settable at any point: compilation is lazy and only
    changes how a cascade is executed, so a vine can be flipped either way
    between calls. ``TorchVinecop.from_data()`` sets it from
    ``controls.compile``. A compiled cascade agrees with the eager one to
    floating point rather than exactly.

    Worth it on CUDA for a cascade evaluated repeatedly, where the eager path
    is bound by kernel-launch count rather than by arithmetic. Not worth it
    for a single evaluation: the first call at each input shape pays tens of
    seconds of compilation.

    Returns
    -------
    bool
        Whether compilation is enabled.

    Warnings
    --------
    Torch caps how many compiled variants of one code object it keeps
    (``torch._dynamo.config.cache_size_limit``, 8 by default), and each vine
    is a variant, as is each input shape. A process that compiles more than
    that falls back to eager, which shows up as this flag doing nothing; raise
    the cap if you genuinely need many compiled vines at once.
    """
    return self._compile_cascades

  @compile_cascades.setter
  def compile_cascades(self, value: bool) -> None:
    self._compile_cascades = bool(value)

  def _cascade(self, name: str) -> Any:
    """The named batched cascade, compiled if ``compile_cascades`` is set."""
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

    The graph replay writes its result into a buffer the next replay reuses,
    so what comes back is copied out before the caller sees it.

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
    unpickled vine recompiles on demand.

    Returns
    -------
    dict
        The instance state, with an empty compiled-cascade cache.
    """
    state = dict(super().__getstate__())
    state["_compiled"] = {}
    return state

  def load_state_dict(self, *args: Any, **kwargs: Any) -> Any:
    """Load parameters and buffers, dropping anything derived from them.

    Parameters
    ----------
    *args, **kwargs
        Forwarded to :meth:`torch.nn.Module.load_state_dict`.

    Returns
    -------
    torch.nn.modules.module._IncompatibleKeys
        Whatever the base implementation returns.
    """
    # The stacked bake and the compiled cascades are copies of the grids, not
    # views of them, so a load that replaces the grids leaves both answering
    # from the old density. `_apply` drops them for the same reason.
    out = super().load_state_dict(*args, **kwargs)
    self._batched = None
    self._compiled = {}
    return out

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
    one), and the derived integral caches are constants with respect to them.

    Returns
    -------
    tuple of bool
        Two entries per pair copula, in tree-then-edge order.
    """
    out: list[bool] = []
    # The same pairs `_build_batched` bakes, in the same order.
    for tree in range(self.trunc_lvl):
      for edge in range(self.d - tree - 1):
        grid = getattr(self.get_pair_copula(tree, edge), "interp_grid", None)
        if grid is None:
          continue
        out.append(bool(grid.values.requires_grad))
        out.append(bool(grid.grid_points.requires_grad))
    return tuple(out)

  def _ensure_batched(self) -> Any:
    """The batched state, re-baked when grad tracking has changed under it.

    A bake is a copy of each pair's grid, which goes stale in two ways a
    device move does not cover. ``requires_grad_`` flips a flag in place, so
    a bake made before it is left behind; and a bake made under ``no_grad``
    -- as ``sample`` / ``cdf`` / ``inverse_rosenblatt`` are evaluated -- holds
    detached copies even where the grids themselves track grad, so it is
    redone once, for the first call that needs the graph. Neither would raise
    where it is read: the batched result would simply arrive detached.

    Returns
    -------
    object
        The memoized batched state.
    """
    signature = self._grad_signature()
    wants_graph = torch.is_grad_enabled() and any(signature)
    baked = getattr(self, "_bake_signature", None)
    if self._batched is not None and (
      baked is None or baked[0] != signature or (wants_graph and not baked[1])
    ):
      object.__setattr__(self, "_batched", None)
      # A compiled cascade was traced against the grids the stale bake holds.
      object.__setattr__(self, "_compiled", {})
    fresh = self._batched is None
    out = super()._ensure_batched()
    if fresh:
      object.__setattr__(self, "_bake_signature", (signature, wants_graph))
    return out

  # --------------------------------------------------------------------- #
  # VinecopBase hooks: RNG for sample + grad control                     #
  # --------------------------------------------------------------------- #

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> Tensor:
    """Draw ``(n, d)`` base uniforms on the fitted grid's dtype/device.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    qrng : bool
        Draw a low-discrepancy (quasi-random) sequence instead.
    seeds : list of int
        RNG seeds; only the first is read for the pseudo-random draw.

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
        Stacked per-tree-level grids and caches for the batched cascades.

    Raises
    ------
    _NotBatchable
        If any variable is discrete -- the stacked per-level grids carry no
        distribution function, which a discrete edge's h-functions are
        difference quotients of -- or if any pair lacks the grid internals the
        batched path reads (``supports_batched`` is ``False``). The dispatch
        layer catches it and falls back to the non-batched cascade.
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
