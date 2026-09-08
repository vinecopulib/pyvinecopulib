"""PyTorch R-vine copula on ``TorchTllBicop`` pair copulas.

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
TorchTllBicop : The pair copulas this hosts.
FitControlsTorchVinecop : Fit-time controls.
"""

from __future__ import annotations

from typing import (
  TYPE_CHECKING,
  Any,
  Callable,
  ClassVar,
  Optional,
  Self,
  Sequence,
  Union,
  cast,
)

import numpy as np
import torch
from torch import Tensor

from ..core import (
  BicopLike,
  ConditioningContext,
  ControlsLike,
  DiscretePair,
  VinecopBase,
)
from ..core._discrete import continuous_view
from ..core.independence import IndependencePair
from ..core._validation import reject_covariates
from ..core.vinecop_base import FitEdge, FitLevel, _NotBatchable
from ..pyvinecopulib_ext import (
  RVineStructure,
  Vinecop,
  indep as _INDEP_FAMILY,
  tll as _TLL_FAMILY,
)
from ..utils import sample_uniform
from ._batched import BatchedVine
from ._placement import TensorPlacementMixin, reference_tensor
from .controls import FitControlsTorchVinecop
from .tll_bicop import TorchTllBicop

if TYPE_CHECKING:
  from torch.nn.modules.module import _IncompatibleKeys


def _placement_of(pair: torch.nn.Module) -> Tensor:
  """A tensor to read dtype and device from, for any pair copula.

  Every pair a vine can hold is an ``nn.Module``, so it carries at least one
  floating-point parameter or buffer -- which is all a placement probe needs.
  Reading a grid-specific attribute instead would make the vine hold only grid
  pairs, and a vine is a container: what it hosts is the caller's business.

  Parameters
  ----------
  pair : torch.nn.Module
      The pair copula to read the placement from.

  Returns
  -------
  Tensor
      A tensor carrying the pair's dtype and device.

  Raises
  ------
  TypeError
      If the pair registers no floating-point parameter or buffer, so its
      placement cannot be determined.
  """
  ref = reference_tensor(pair)
  if ref is None:
    raise TypeError(
      f"{type(pair).__name__} registers no floating-point parameter or "
      "buffer, so a vine cannot read its dtype and device. Register the "
      "tensors it evaluates with, as `TorchTllBicop` registers its grid."
    )
  return ref


class TorchVinecop(
  TensorPlacementMixin, VinecopBase[torch.Tensor], torch.nn.Module
):
  """PyTorch R-vine copula on ``TorchTllBicop`` pair copulas.

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
  - ``bicop_class`` is ``TorchTllBicop``. Naming it is what lets
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
  pair_copulas : list of list of torch.nn.Module
      The pair copulas, indexed ``[tree][edge]`` and shaped as ``Vinecop``
      lays them out: tree ``t`` holds ``d - 1 - t`` edges, up to the
      structure's ``trunc_lvl``. A ``torch.nn.ModuleList`` of
      ``torch.nn.ModuleList`` is accepted too.

      Each pair must satisfy ``BicopLike`` *and* be a ``torch.nn.Module`` --
      the first to be evaluated, the second so the vine owns it as a child and
      one ``.to(device)`` moves everything. Any such pair works, not only
      ``TorchTllBicop``: subclass ``BicopBase``, define ``pdf`` / ``hfunc1`` /
      ``hfunc2``, register whatever tensors it evaluates with, and a vine will
      host it -- including one whose parameters an optimizer learns.
  structure : RVineStructure
      The vine structure to evaluate along.
  context : ConditioningContext, or None, optional
      Per-edge policy assembling each pair copula's ``x``. ``None`` uses
      :class:`~pyvinecopulib.core.SimplifiedContext` -- an unconditional,
      simplified vine.
  var_types : list of str, or None, optional
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
  bicop_class: ClassVar[Optional[type]] = TorchTllBicop

  def __init__(
    self,
    pair_copulas: Sequence[Sequence[torch.nn.Module]],
    structure: RVineStructure,
    *,
    context: Optional[ConditioningContext[Tensor]] = None,
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
      _placement_of(pair_copulas[0][0])
      if self.trunc_lvl > 0
      else torch.empty(0)
    )
    self.register_buffer(
      "_device_ref", torch.empty(0, dtype=ref.dtype, device=ref.device)
    )
    # `self._batched` (lazy grid-batched state) is initialized to None by
    # `_bind_vine`; `_apply` clears it so device moves rebuild it.
    self._compile_cascades = False
    self._compiled: dict[str, Callable[[Tensor], Tensor]] = {}

  # --------------------------------------------------------------------- #
  # Constructor                                                            #
  # --------------------------------------------------------------------- #

  @staticmethod
  def _resolve_cache_integrals(cache_integrals: Optional[bool]) -> bool:
    """Whether to precompute the prefix tables. ``None`` resolves to ``True``.

    A variable's type does not enter the decision: the prefix tables reconstruct
    the integral exactly rather than approximately, so an edge that reads its
    density from differences over an atom's width can difference them safely.
    Such an edge still differences the distribution function rather than
    calling ``rect_mass``, which is more accurate but would break the cascade
    parity with ``Vinecop`` that its own difference quotients define.

    Parameters
    ----------
    cache_integrals : bool, or None, optional
        What the caller asked for; ``None`` means "whatever suits this vine".

    Returns
    -------
    bool
        The effective setting.
    """
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

    The result hosts one ``TorchTllBicop`` per pair copula, on the same grids, so
    it agrees with ``cop`` to floating-point tolerance.

    Parameters
    ----------
    cop : Vinecop
        A fitted vine whose pair copulas are all of the ``tll`` or ``indep``
        family, in any continuous, discrete or mixed variable layout.
    cache_integrals : bool, or None, optional
        Whether each pair copula precomputes the prefix tables its ``cdf`` and
        h-functions read. ``None`` resolves to ``True`` for every variable
        layout.
    device : torch.device, or None, optional
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
    cache_integrals = cls._resolve_cache_integrals(cache_integrals)
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
    pair_copulas_torch: list[list[TorchTllBicop]] = []
    for tree_idx, row in enumerate(cop.pair_copulas):
      tree_list: list[TorchTllBicop] = []
      for edge_idx, b in enumerate(row):
        if b.family == _TLL_FAMILY:
          bc = TorchTllBicop.from_bicop(
            b,
            cache_integrals=cache_integrals,
            device=device,
            dtype=dtype,
          )
        elif b.family == _INDEP_FAMILY:
          bc = TorchTllBicop(device=device, dtype=dtype)
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
    pair_copulas: list[list[TorchTllBicop]] = [],
    var_types: list[str] = [],
    *,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchVinecop":
    """Build a ``TorchVinecop`` from a structure and pair copulas.

    Parameters
    ----------
    structure : RVineStructure, or None, optional
        The vine structure. Provide either this or ``matrix``.
    matrix : ndarray, shape (d, d), dtype int, or None, optional
        R-vine structure matrix. Provide either this or ``structure``.
    pair_copulas : list of list of TorchTllBicop, default=[]
        The pair copulas, indexed ``[tree][edge]`` with tree ``t`` holding
        ``d - 1 - t`` edges. Empty fills every edge with the independence
        copula.
    var_types : list of str, default=[]
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete), in
        variable order; empty means all continuous. A discrete variable makes
        the cascades read its left limit too, and the pair copulas that see it
        are evaluated through :class:`~pyvinecopulib.core.DiscretePair`.
    device : torch.device, or None, optional
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
      # Independence vine: TorchTllBicop() defaults to the independence copula.
      pair_copulas = [
        [TorchTllBicop(device=device, dtype=dtype) for _ in range(d - 1 - t)]
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
    u: Union[np.ndarray, Tensor],
    /,
    # Declared as the base does -- narrowing a parameter is what an override
    # may not do -- while the lane's own fields are read off a local below.
    controls: Optional[ControlsLike] = None,
    *,
    structure: Optional[RVineStructure] = None,
    var_types: Optional[list[str]] = None,
    x: Optional[Tensor] = None,
    fit_edge: Optional[FitEdge] = None,
    fit_level: Optional[FitLevel] = None,
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
    controls : FitControlsTorchVinecop, or None, optional
        Fit configuration for both halves of the fit -- the structure
        selection the vine runs and the pair-copula fits its edges run -- plus
        placement, precision, and the cascade variants. ``None`` defaults to
        TLL on a 30x30 normal-spaced grid, float64, and ``mst_prim`` with
        ``trunc_lvl=20``.
    structure : RVineStructure, or None, optional
        A fixed structure. Selected from the data when ``None``.
    var_types : list of str, or None, optional
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete);
        ``None`` means all continuous. Given, they also fix the dimension, so
        ``u`` may carry the extra left-limit columns.
    x : Tensor, or None, optional
        Refused. A TLL grid is an unconditional density, so fitting one while
        ignoring covariates would return a different model than the caller
        asked for; a conditional pair copula is fitted through ``fit_edge``.
    fit_edge : callable, or None, optional
        ``(tree, edge, u_e, x_e) -> BicopLike``, fitting one edge's pair
        copula in place of the built-in TLL fit; an edge with a discrete
        variable additionally receives its own ``var_types`` by keyword.
        ``None`` fits a ``TorchTllBicop`` per edge, which is the usual case.
    fit_level : callable, or None, optional
        ``(tree, u_level, types) -> Sequence[BicopLike]``, fitting a whole tree
        level at once, in place of the built-in batched TLL fit -- the
        level-wise counterpart of ``fit_edge``. ``None`` resolves it from
        ``controls.batched_fit``.

    Returns
    -------
    TorchVinecop
        The fitted vine.

    Raises
    ------
    ValueError
        If ``u`` is not 2-D, if ``structure``'s dimension disagrees with the
        one ``var_types`` or ``u`` implies, or if ``x`` is given.

    See Also
    --------
    pyvinecopulib.core.VinecopBase.fit : Refit an existing vine's pairs.
    pyvinecopulib.core.VinecopBase.select : Reselect its structure and pairs.
    """
    reject_covariates(cls, x)
    resolved: Any = controls
    if resolved is None:
      resolved = FitControlsTorchVinecop()

    eff_dtype = resolved.dtype if resolved.dtype is not None else torch.float64
    eff_device = resolved.device
    u_t = torch.as_tensor(u, dtype=eff_dtype, device=eff_device)
    if u_t.ndim != 2:
      raise ValueError(f"u must be 2-D; got shape {tuple(u_t.shape)}")
    # With left-limit columns present `u` is wider than the vine, so it is
    # `var_types` that fixes the dimension.
    d = len(var_types) if var_types else int(u_t.shape[1])

    # The vine's controls are pair controls: `FitControlsTorchVinecop`
    # derives from `FitControlsTorchBicop`, as its core counterparts do.
    bc_controls = resolved
    cache_integrals = cls._resolve_cache_integrals(resolved.cache_integrals)

    def fit_edge_tll(
      tree: int,
      edge: int,
      u_e: Tensor,
      x_e: Optional[Tensor],
      var_types: Sequence[str] = ("c", "c"),
    ) -> BicopLike[Tensor]:
      # `var_types` here is *this edge's* two types, which the fit engines pass
      # by that keyword; the vine's own list is the enclosing argument.
      # Simplified (unconditional) TLL fit — x_e is None here.
      del tree, edge, x_e
      bc = TorchTllBicop.from_data(
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
    batched_fit = resolved.batched_fit
    if batched_fit is None:
      batched_fit = u_t.device.type == "cuda"

    def fit_level_tll(
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
      return TorchTllBicop.from_data_batched(
        u_level,
        bc_controls,
        cache_integrals=cache_integrals,
        device=u_t.device,
        dtype=eff_dtype,
      )

    # A caller's own level fitter wins, as their `fit_edge` does; otherwise
    # the built-in one, and only where batching is worth its launch overhead.
    level_hook = fit_level or (fit_level_tll if batched_fit else None)
    cond_order: dict[tuple[int, int], tuple[int, ...]] = {}

    if structure is None:
      # Select the structure natively in torch, reusing the pairs fit during
      # selection (reoriented onto their slots via TorchTllBicop.flip) — exactly
      # what Vinecop's selector does, so no re-fit is needed. Kendall's tau
      # via wdm needs a host copy; detach so grad-tracking tensors are
      # accepted.
      structure, pairs, cond_order = cls._select_parts(
        u_t,
        pair_fitter,
        fit_level=level_hook,
        trunc_lvl=resolved.trunc_lvl,
        tree_criterion=resolved.tree_criterion,
        threshold=resolved.threshold,
        tree_algorithm=resolved.tree_algorithm,
        seeds=list(resolved.seeds),
        var_types=list(var_types or []) or None,
        conditioning_set=list(resolved.conditioning_set) or None,
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
        tree_criterion=resolved.tree_criterion,
        threshold=resolved.threshold,
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
        # no-argument `TorchTllBicop` *is* the independence copula -- a 2x2
        # sentinel that short-circuits every method on `is_indep`, exactly
        # rather than to rounding -- so it needs no grid of its own and
        # cannot disagree with its siblings about one. `u_t.device`, as
        # `fit_edge` uses: `resolved.device` is `None` whenever the caller
        # let the data carry the placement.
        TorchTllBicop(device=u_t.device, dtype=eff_dtype)
        if isinstance(p, IndependencePair)
        else continuous_view(p)
        for p in row
      ]
      for row in pairs
    ]
    out = cls(
      pair_copulas=cast("list[list[TorchTllBicop]]", modules),
      structure=structure,
      var_types=list(var_types or []) or None,
    )
    out._set_cond_order(cond_order)
    out.compile_cascades = resolved.compile
    return out

  # --------------------------------------------------------------------- #
  # Helpers                                                                #
  # --------------------------------------------------------------------- #

  def set_pair_copulas(
    self, pair_copulas: list[list[BicopLike[Tensor]]]
  ) -> None:
    """Install fitted pairs, so ``fit`` and ``select`` can return ``self``.

    Parameters
    ----------
    pair_copulas : list of list of BicopLike
        Fitted pairs indexed ``[tree][edge]``, as the fit engines return them.
        A pair carrying discrete variables is stored as its continuous grid
        and re-wrapped on read; an edge left independent by ``threshold`` is
        stored as the independence ``TorchTllBicop``.

    Returns
    -------
    None
    """
    ref = cast("Tensor", self._buffers["_device_ref"])
    # Only real `nn.Module`s go in the `ModuleList`: a discrete edge is
    # re-wrapped on read by `get_pair_copula`, and a thresholded edge arrives
    # as a `core.IndependencePair`, which the no-argument `TorchTllBicop` -- the
    # independence copula exactly rather than to rounding -- stands in for.
    self.pair_copulas = torch.nn.ModuleList(
      [
        torch.nn.ModuleList(
          [
            TorchTllBicop(device=ref.device, dtype=ref.dtype)
            if isinstance(pair, IndependencePair)
            else cast("torch.nn.Module", continuous_view(pair))
            for pair in row
          ]
        )
        for row in pair_copulas
      ]
    )
    # Same reason `load_state_dict` and `_apply` drop them: the stacked state
    # and the compiled cascades hold copies of the grids, not views, so pairs
    # replaced under them leave both answering from the old density.
    self._batched = None
    self._compiled = {}

  def _pair_module(self, tree: int, edge: int) -> TorchTllBicop:
    """The stored (always continuous) pair copula at ``(tree, edge)``."""
    return cast(
      "TorchTllBicop",
      cast("torch.nn.ModuleList", self.pair_copulas[tree])[edge],
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

  def set_extra_state(self, state: object) -> None:
    """Check the non-tensor model identity in a ``state_dict``.

    The pair-copula grids are the only tensors a ``state_dict`` holds, so
    nothing in it identifies the vine they hang on: a checkpoint from a
    differently structured vine carries exactly the same keys, and is refused
    here rather than loaded into a model it does not describe.

    Parameters
    ----------
    state : object
        State returned by ``TorchVinecop.get_extra_state()``; anything else is
        refused.

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
    """A registered buffer to read dtype and device from.

    Returns
    -------
    Tensor
        A tensor carrying the vine's dtype and device.
    """
    # Read the placement off the first pair, whatever kind of pair it is,
    # rather than through the mixin's module-wide search: a pair names it in
    # one hop. A vine truncated at zero has no pair and falls back to the
    # buffer the constructor registers for exactly that case.
    if self.trunc_lvl > 0:
      return _placement_of(self._pair_module(0, 0))
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
    the cap if you actually need many compiled vines at once.
    """
    return self._compile_cascades

  @compile_cascades.setter
  def compile_cascades(self, value: bool) -> None:
    self._compile_cascades = bool(value)

  def _cascade(self, name: str) -> Callable[[Tensor], Tensor]:
    """The named batched cascade, compiled if ``compile_cascades`` is set."""
    base = cast("Callable[[Tensor], Tensor]", getattr(super(), name))
    if not self._compile_cascades:
      return base
    fn = self._compiled.get(name)
    if fn is None:
      fn = self._compile_cascade(base)
      self._compiled[name] = fn
    return fn

  def _compile_cascade(
    self, base: Callable[[Tensor], Tensor]
  ) -> Callable[[Tensor], Tensor]:
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
      return inner(u).clone()

    return graphed

  def _pdf_batched(self, u: Tensor) -> Tensor:
    return self._cascade("_pdf_batched")(u)

  def _rosenblatt_batched(self, u: Tensor) -> Tensor:
    return self._cascade("_rosenblatt_batched")(u)

  def _inverse_rosenblatt_batched(self, u: Tensor) -> Tensor:
    return self._cascade("_inverse_rosenblatt_batched")(u)

  def __getstate__(self) -> dict[str, Any]:
    """The picklable state: everything except the two derived caches.

    Both are pure caches rebuilt on demand, and neither belongs in a pickle.
    Compiled callables cannot be pickled at all. The grid-batched state can,
    which is the trap: it is a *copy* of every pair's grid, so a pickle taken
    after one batched call carried the grids twice -- 2.9x the bytes on a
    3-variable vine -- and restored a state nothing revalidated against the
    pairs it was built from.

    Returns
    -------
    dict
        The instance state, with both caches cleared.
    """
    state = dict(super().__getstate__())
    state["_compiled"] = {}
    state["_batched"] = None
    return state

  def load_state_dict(
    self,
    # Forwarded verbatim to `nn.Module.load_state_dict`, which torch itself
    # declares untyped and has re-signed across releases.
    *args: Any,  # noqa: ANN401
    **kwargs: Any,  # noqa: ANN401
  ) -> "_IncompatibleKeys":
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
    # The stacked state and the compiled cascades are copies of the grids, not
    # views of them, so a load that replaces the grids leaves both answering
    # from the old density. `_apply` drops them for the same reason.
    out = super().load_state_dict(*args, **kwargs)
    self._batched = None
    self._compiled = {}
    return cast("_IncompatibleKeys", out)

  def _apply(
    self,
    fn: Callable[[Tensor], Tensor],
    *args: Any,  # noqa: ANN401 - as `load_state_dict`
    **kwargs: Any,  # noqa: ANN401
  ) -> Self:
    # `.to()`, `.cuda()`, `.cpu()` all route through `_apply`. The
    # BatchedVine container holds buffers — `super()._apply` would move
    # them, but we drop the whole structure so it gets rebuilt from the
    # (already-moved) source pair_copulas on next use; that keeps the
    # wire-up tensors aligned with the destination dtype/device.
    self._batched = None
    # Compiled code is specialized on the tensors it was traced with; the
    # guards would recompile anyway, so drop the stale entries.
    self._compiled = {}
    return cast("Self", super()._apply(fn, *args, **kwargs))

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
    # The same pairs `_build_batched` precomputes, in the same order.
    for tree in range(self.trunc_lvl):
      for edge in range(self.d - tree - 1):
        grid = getattr(self.get_pair_copula(tree, edge), "interp_grid", None)
        if grid is None:
          continue
        out.append(bool(grid.values.requires_grad))
        out.append(bool(grid.grid_points.requires_grad))
    return tuple(out)

  def _pair_revisions(self) -> tuple[int, ...]:
    """How many times each pair has had its grid replaced.

    Kept apart from ``_grad_signature`` because the two answer different
    questions: that one decides whether the state needs the graph, this one
    whether it is state for the right density at all. Refitting a pair the vine
    already holds -- ``vine.get_pair_copula(t, e).fit(u)`` -- replaces its grid
    in place and moves no ``requires_grad`` flag, so nothing else would notice.

    Returns
    -------
    tuple of int
        One entry per pair copula, in tree-then-edge order.
    """
    return tuple(
      int(getattr(self.get_pair_copula(tree, edge), "_revision", 0))
      for tree in range(self.trunc_lvl)
      for edge in range(self.d - tree - 1)
    )

  def _ensure_batched(self) -> "BatchedVine":
    """The batched state, rebuilt when grad tracking has changed under it.

    The state holds a copy of each pair's grid, which goes stale in three ways a
    device move does not cover. ``requires_grad_`` flips a flag in place, so
    state built before it is left behind; state built under ``no_grad``
    -- as ``sample`` / ``cdf`` / ``inverse_rosenblatt`` are evaluated -- holds
    detached copies even where the grids themselves track grad, so it is
    redone once, for the first call that needs the graph; and refitting a pair
    the vine holds replaces its grid, which ``_pair_revisions`` counts.
    Only the last is wrong rather than merely detached, but all three are
    silent where the state is read.

    Returns
    -------
    BatchedVine
        The memoized batched state.
    """
    signature = self._grad_signature()
    revisions = self._pair_revisions()
    wants_graph = torch.is_grad_enabled() and any(signature)
    stamped = getattr(self, "_batched_signature", None)
    if self._batched is not None and (
      stamped is None
      or stamped[0] != signature
      or stamped[2] != revisions
      or (wants_graph and not stamped[1])
    ):
      object.__setattr__(self, "_batched", None)
      # A compiled cascade was traced against the grids the stale state holds.
      object.__setattr__(self, "_compiled", {})
    fresh = self._batched is None
    out = cast("BatchedVine", super()._ensure_batched())
    if fresh:
      object.__setattr__(
        self, "_batched_signature", (signature, wants_graph, revisions)
      )
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

  def _eval_context(self) -> torch.no_grad:
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
  # (stacked, precomputed per-tree-level grids + caches).

  def _build_batched(self) -> "BatchedVine":
    """Precompute the grid-batched state from this vine's ``TorchTllBicop`` pairs.

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
