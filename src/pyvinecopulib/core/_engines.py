"""The two array-agnostic fit engines behind ``VinecopBase``.

``fit_parts`` estimates pair copulas tree by tree along a **fixed** structure;
``select_parts`` chooses the structure from the data as well, an exact port of
``Vinecop.select``'s Dissmann / Wilson search, whose selected matrix it matches
byte for byte. Both **return** the loose parts a caller assembles -- pairs, and
for a selection the structure and each slot's fitted conditioning order --
because a factory needs them before an object exists to install them on.

They live here rather than on the class because they are module functions
already: neither reads ``self`` or ``cls``, and both were ``@staticmethod``s in
a 2677-line class body only for namespacing. ``VinecopBase`` keeps
``_fit_parts`` / ``_select_parts`` as the names external code drives them
through.

Also here, because both engines and the class need them: the ``FitEdge`` /
``FitLevel`` callback aliases, the edge-criterion factory, and the two small
call helpers.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Optional, Sequence, cast

from array_api_compat import array_namespace

from ._covariates import pair_eval, prepare
from ._discrete import (
  check_var_types,
  collapse_data,
  disc_cols,
  edge_columns,
  pair_var_types,
  seed_left_limits,
  stack_edge,
  with_left_limit,
)
from .independence import IndependencePair
from ._reorient import _SlotKey, _slot_key, reorientation
from ._validation import validate_weights
from .bicop_base import flip_of
from .context import ConditioningContext, SimplifiedContext
from .protocols import ArrayT, BicopLike

if TYPE_CHECKING:
  import numpy as np

  from ..pyvinecopulib_ext import RVineStructure

__all__ = ["FitEdge", "FitLevel", "fit_parts", "select_parts"]


def _make_criterion(
  tree_criterion: str,
  n: int,
  weights: Optional[ArrayT] = None,
  criterion_function: Optional[Callable[..., float]] = None,
  x: Optional[ArrayT] = None,
) -> Callable[[ArrayT, ArrayT], float]:
  """Build the edge criterion ``calculate_criterion`` computes.

  Covariates are optional to the criterion, the way they are to a pair copula:
  the built-in dependence measures ignore ``x`` -- they read the
  pseudo-observations and nothing else -- while a caller's own
  ``criterion_function`` may use it to score an edge conditionally. It is
  passed only when there is one, so a criterion function written without it
  keeps working, and one that cannot accept an ``x`` it is handed fails loudly
  rather than scoring the wrong thing.

  Parameters
  ----------
  tree_criterion : str
      Dependence measure, as on ``FitControlsVinecop``.
  n : int
      Number of observations; at or below ten the criterion is zero, as
      upstream's guard has it.
  weights : array, shape (n,), or None, optional
      Observation weights, applied to the criterion itself so a weighted
      selection agrees with :meth:`~pyvinecopulib.core.Vinecop.select`.
  criterion_function : callable, or None, optional
      Required when ``tree_criterion`` is ``"custom"``; maps an ``(n, 2)``
      matrix -- and, when there are covariates, ``x`` by keyword -- to a
      criterion value.
  x : array, shape (n, p), or None, optional
      Exogenous covariates, forwarded to ``criterion_function`` only. The
      built-in measures are unconditional and never see them.

  Returns
  -------
  callable
      Maps an edge's two value columns to a non-negative criterion.
  """
  import numpy as np

  from ..pyvinecopulib_ext import _calculate_tree_criterion

  convert = _to_numpy_default
  w = np.empty(0) if weights is None else np.asarray(convert(weights), float)
  # The binding calls the criterion function with the matrix alone, so the
  # covariates are bound here rather than threaded through C++.
  scorer: Optional[Callable[..., float]] = criterion_function
  if criterion_function is not None and x is not None:
    xa = convert(x)

    def scorer(matrix: np.ndarray) -> float:  # noqa: F811 - the conditional variant
      return float(cast("Any", criterion_function)(matrix, x=xa))

  def criterion(col0: ArrayT, col1: ArrayT) -> float:
    if n <= 10:
      return 0.0
    a, b = convert(col0), convert(col1)
    return float(
      _calculate_tree_criterion(
        np.column_stack((a, b)), tree_criterion, w, scorer
      )
    )

  return criterion


def _to_numpy_default(a: Any) -> np.ndarray:
  """Host NumPy view of ``a``, for an array library that leaves it off the host.

  ``np.asarray`` raises on a tensor that lives on an accelerator, so this
  detaches and transfers before converting.
  """
  import numpy as _np

  detach = getattr(a, "detach", None)
  if detach is not None:
    a = detach()
  cpu = getattr(a, "cpu", None)
  if cpu is not None:
    a = cpu()
  return _np.asarray(a)


#: ``(tree, edge, u_e, x_e) -> BicopLike``, fitting one edge's pair copula: the
#: hook external packages drive conditional fitting through (see
#: :meth:`VinecopBase.fit`). An edge with a discrete argument additionally
#: receives ``var_types=[t1, t2]`` and a four-column ``u_e``, so the alias cannot
#: pin the arity -- a ``Callable`` has no way to express a keyword argument.
FitEdge = Callable[..., BicopLike[Any]]

#: ``(tree, u_level, types) -> list[BicopLike]``, fitting a whole tree level
#: at once: the optional companion to ``FitEdge``, for a lane whose fitter
#: carries a leading pair axis. ``u_level`` stacks the level's edges in
#: ascending edge order and ``types`` gives each edge's pair of variable types,
#: so the callback needs no structural knowledge. A subclass that supplies one
#: gets it preferred over ``fit_edge``; everything else keeps working, since a
#: level is only ever fitted this way when every one of its edges is
#: continuous -- a mixed level cannot stack, its edges having different widths.
# Written out rather than quoted: a string inside a type *alias* is resolved in
# whichever module uses the alias, so quoting it made `Sequence` a name every
# importer had to keep in scope -- and dropping that import broke the docs
# build, which resolves annotations at runtime.
FitLevel = Callable[[int, Any, list[tuple[str, str]]], Sequence[BicopLike[Any]]]


def _fit_edge_call(
  fit_edge: FitEdge,
  tree: int,
  edge: int,
  u_e: ArrayT,
  x_e: Optional[ArrayT],
  var_types: tuple[str, str],
) -> BicopLike[Any]:
  """Call ``fit_edge``, forwarding ``var_types`` only for a discrete edge.

  A callback written for a continuous vine takes four arguments, so a fully
  continuous edge must not hand it a fifth. When the edge does have a discrete
  argument the types go by keyword, which is what makes a callback that cannot
  accept them fail loudly instead of fitting a continuous pair copula to four
  columns of data -- the rule ``pair_eval`` applies to ``x``, for the same
  reason.
  """
  if "d" not in var_types:
    return fit_edge(tree, edge, u_e, x_e)
  return fit_edge(tree, edge, u_e, x_e, var_types=list(var_types))


def fit_parts(
  structure: RVineStructure,
  u: Any,
  fit_edge: FitEdge,
  *,
  context: Optional[ConditioningContext[ArrayT]] = None,
  x: Optional[ArrayT] = None,
  var_types: Optional[list[str]] = None,
  fit_level: Optional[FitLevel] = None,
  tree_criterion: str = "tau",
  threshold: float = 0.0,
  weights: Optional[ArrayT] = None,
  criterion_function: Optional[Callable[[Any], float]] = None,
) -> list[list[BicopLike[Any]]]:
  """Fit pair copulas tree-by-tree along a fixed structure (returns them).

  The engine behind ``fit``, kept separate because a factory needs the pairs
  before an object exists to install them on. The array-agnostic (NumPy or
  PyTorch) analog of :meth:`~pyvinecopulib.core.Vinecop.fit`, with the
  pair-copula fit supplied by the ``fit_edge`` callback, it **returns** the
  fitted pairs as a nested ``[tree][edge]`` list rather than storing them.
  It mirrors the forward pdf traversal with the density
  evaluation replaced by ``fit_edge(tree, edge, u_e, x_e)``; the returned
  pair's ``hfunc1`` / ``hfunc2`` must be valid immediately for tree
  propagation. Conditional fitting is driven through this hook (a
  ``fit_edge`` that fits a conditional pair copula on ``(u_e, x_e)``), with
  ``x_e`` assembled in the same C1 order the cascades use.

  Parameters
  ----------
  structure : RVineStructure
      The (fixed) vine structure to fit along.
  u : array, shape (n, d), (n, d + k) or (n, 2d), dtype float
      Pseudo-observations. With ``k`` discrete variables their left limits
      ``F(x^-)`` are required too; see :meth:`pdf` on the layouts.
  fit_edge : callable
      ``(tree, edge, u_e, x_e) -> BicopLike`` fitting one edge's pair copula.
      An edge with a discrete argument gets a four-column ``u_e`` and the
      additional keyword ``var_types=[t1, t2]``; the pair it returns must read
      that layout, so wrap a continuous one in
      :class:`~pyvinecopulib.core.DiscretePair`.
  context : ConditioningContext, or None, optional
      Conditioning-context policy (default: simplified / unconditional).
  x : array, shape (n, p), or None, optional
      External covariates for conditional fitting, else ``None``.
  var_types : list of str, or None, optional
      Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete), in
      variable order; ``None`` means all continuous.
  fit_level : callable, or None, optional
      ``(tree, u_level, types) -> list[BicopLike]``, fitting a whole tree
      level at once; see ``FitLevel``. Preferred over ``fit_edge``
      for a level whose edges are all continuous and unconditional,
      which is the only shape that stacks -- a discrete edge is four
      columns wide where a continuous one is two. ``None`` fits every
      edge separately, which is what every caller did before the hook
      existed.
  tree_criterion : str, default "tau"
      Dependence measure ``threshold`` compares against, as on
      ``FitControlsVinecop``. Read only when ``threshold`` is positive.
  threshold : float, default 0.0
      Dependence threshold. An edge whose criterion falls below it holds
      :class:`~pyvinecopulib.core.IndependencePair` and is not fitted, as
      it does under selection. At the default nothing is below it.
  weights : array, shape (n,), or None, optional
      Observation weights, applied to the tree criterion so a weighted
      selection agrees with :meth:`~pyvinecopulib.core.Vinecop.select`. They
      reach the pair fits through ``controls``, which a default
      ``bicop_class`` fit reads; a caller's own ``fit_edge`` receives no
      weights and has to apply them itself.
  criterion_function : callable, or None, optional
      Required when ``tree_criterion`` is ``"custom"``; maps an ``(n, 2)``
      matrix -- and, when there are covariates, ``x`` by keyword -- to a
      criterion value.

  Returns
  -------
  list of list of BicopLike
      Fitted pair copulas indexed ``[tree][edge]``.

  Raises
  ------
  ValueError
      If ``var_types`` has the wrong length or an entry outside
      ``{"c", "d"}``, or if ``u``'s column count matches no accepted layout.

  See Also
  --------
  pyvinecopulib.core.VinecopBase.select : Select a structure and fit it.
  pyvinecopulib.core.Vinecop.fit : The reference (in-place) fit.
  """
  if context is None:
    context = SimplifiedContext()
  ua: Any = u
  xp = array_namespace(ua)
  d = int(structure.dim)
  trunc_lvl = int(structure.trunc_lvl)
  order = tuple(int(v) for v in structure.order)
  types = check_var_types(var_types, d)
  pair_types = pair_var_types(structure, types) if "d" in types else None
  ua = collapse_data(ua, d, types, "fit")
  n = ua.shape[0]
  x = prepare(ua, x, int(n))
  # The criterion hands these straight to the binding, so they are checked
  # here rather than there: this is the first point that knows both the
  # weights and the row count they must align with.
  weights = validate_weights(weights, ua[:, 0])
  criterion = _make_criterion(
    tree_criterion,
    int(n),
    weights,
    criterion_function,
    x,
  )
  hfunc1 = xp.zeros((n, d), dtype=ua.dtype, device=ua.device)
  hfunc2 = xp.empty((n, d), dtype=ua.dtype, device=ua.device)
  for j in range(d):
    hfunc2[:, j] = ua[:, order[j] - 1]
  # Parallel left-limit scratch, zero-initialized like the pdf traversal this
  # mirrors; a column is only read where the edge's types say it was written.
  hfunc2_sub: Any = seed_left_limits(ua, d, order, types, disc_cols(types), xp)
  hfunc1_sub: Any = (
    None
    if hfunc2_sub is None
    else xp.zeros((n, d), dtype=ua.dtype, device=ua.device)
  )
  u_nat = (
    xp.asarray(hfunc2, copy=True) if context.assembles_conditioning else None
  )
  cache: dict[tuple[int, int], tuple[int, ...]] = {}

  def edge_context_for(tree: int, edge: int) -> Optional[ArrayT]:
    # Assemble x_e = context(u_D, x); u_D columns in ascending conditioning
    # -tree order (C1), matching VinecopBase._cond_positions / _edge_context.
    if not context.assembles_conditioning and x is None:
      return None
    u_D: Optional[ArrayT] = None
    if context.assembles_conditioning:
      key = (tree, edge)
      if key not in cache:
        cache[key] = tuple(
          int(structure.struct_array(i, edge, natural_order=True)) - 1
          for i in range(tree)
        )
      cols = cache[key]
      if cols and u_nat is not None:
        u_D = u_nat[:, list(cols)]
    return context.edge_context(u_D=u_D, x=x)

  s = structure
  pairs: list[list[BicopLike[Any]]] = []
  for tree in range(trunc_lvl):
    row: list[BicopLike[Any]] = []
    # Every edge of one tree reads only columns finalized by earlier trees
    # -- `min_array(tree, edge) - 1 > edge`, so the column an edge reads
    # second is written later in this same tree -- which is why upstream
    # runs a level on a thread pool. Gathering the level's inputs before
    # fitting any of it is the same reordering, and it is what lets a
    # level be fitted in one call.
    level = [
      edge_columns(
        s, pair_types, tree, edge, hfunc1, hfunc2, hfunc1_sub, hfunc2_sub
      )
      for edge in range(d - tree - 1)
    ]
    inputs = [stack_edge(xp, c0, c1, subs) for c0, c1, subs, _ in level]
    contexts = [edge_context_for(tree, e) for e in range(len(level))]
    # Per-edge type *pairs*, distinct from the per-variable `types` above.
    level_types = [t for _, _, _, t in level]
    # A fixed structure thresholds exactly as selection does: upstream builds
    # the same selector on the given matrix, so `fit_or_reuse_pair_copula`
    # leaves an edge below the threshold holding independence here too.
    skip = [
      threshold > 0.0 and criterion(c0, c1) < threshold
      for c0, c1, _, _ in level
    ]
    to_fit = [e for e, s_e in enumerate(skip) if not s_e]
    fitted: Optional[dict[int, BicopLike[Any]]] = None
    if (
      fit_level is not None
      and to_fit
      and all(contexts[e] is None for e in to_fit)
      and all("d" not in level_types[e] for e in to_fit)
    ):
      got = fit_level(
        tree,
        xp.stack([inputs[e] for e in to_fit], axis=0),
        [level_types[e] for e in to_fit],
      )
      fitted = dict(zip(to_fit, got))
    for edge in range(d - tree - 1):
      _, _, subs, edge_types = level[edge]
      u_e, x_e = inputs[edge], contexts[edge]
      edge_copula: BicopLike[Any]
      if skip[edge]:
        edge_copula = IndependencePair()
      elif fitted is not None:
        edge_copula = fitted[edge]
      else:
        edge_copula = _fit_edge_call(fit_edge, tree, edge, u_e, x_e, edge_types)
      row.append(edge_copula)
      if s.needed_hfunc1(tree, edge):
        hfunc1[:, edge] = pair_eval(edge_copula.hfunc1, u_e, x_e)
        if subs is not None and edge_types[1] == "d":
          hfunc1_sub[:, edge] = pair_eval(
            edge_copula.hfunc1, with_left_limit(u_e, 1), x_e
          )
      if s.needed_hfunc2(tree, edge):
        hfunc2[:, edge] = pair_eval(edge_copula.hfunc2, u_e, x_e)
        if subs is not None and edge_types[0] == "d":
          hfunc2_sub[:, edge] = pair_eval(
            edge_copula.hfunc2, with_left_limit(u_e, 0), x_e
          )
    pairs.append(row)
  return pairs


def select_parts(
  u: Any,
  fit_edge: FitEdge,
  *,
  context: Optional[ConditioningContext[ArrayT]] = None,
  x: Optional[ArrayT] = None,
  fit_level: Optional[FitLevel] = None,
  trunc_lvl: Optional[int] = None,
  tree_criterion: str = "tau",
  threshold: float = 0.0,
  tree_algorithm: str = "mst_prim",
  seeds: Optional[list[int]] = None,
  var_types: Optional[list[str]] = None,
  conditioning_set: Optional[list[int]] = None,
  weights: Optional[ArrayT] = None,
  criterion_function: Optional[Callable[[Any], float]] = None,
) -> tuple[
  RVineStructure,
  list[list[BicopLike[ArrayT]]],
  dict[tuple[int, int], tuple[int, ...]],
]:
  """Select an R-vine structure from data (array-agnostic Dissmann).

  The engine behind ``select`` and ``from_data``, kept separate because a
  factory needs the structure and pairs before an object exists to install
  them on. The array-agnostic (NumPy or PyTorch) analog of
  :meth:`~pyvinecopulib.core.Vinecop.select`, with the pair-copula fit
  supplied by the ``fit_edge`` callback, it **returns** the selected
  structure and pairs rather than storing them. It runs the tree-by-tree
  Dissmann greedy search [1]_.

  Parameters
  ----------
  u : array, shape (n, d), (n, d + k) or (n, 2d), dtype float
      See ``fit_parts``.
  fit_edge : callable
      See ``fit_parts``. The pair must also implement
      :meth:`~pyvinecopulib.core.BicopBase.flip`, which reorients it onto its
      finalized slot.
  context : ConditioningContext, or None, optional
      Conditioning-context policy (default: simplified / unconditional). A
      :class:`~pyvinecopulib.core.NonSimplifiedContext` makes each edge's
      pair copula see its conditioning-set values while it is being fitted,
      not only when it is evaluated.
  x : array, shape (n, p), or None, optional
      External covariates for conditional fitting, else ``None``.
  fit_level : callable, or None, optional
      See ``fit_parts``. Whatever it returns must still be per-slot
      ``flip``-able, since finalization reorients reused pairs.
  trunc_lvl : int, or None, optional
      Maximum number of trees to select (default: ``d - 1``, i.e. untruncated).
  tree_criterion : str, default "tau"
      Dependence measure used for edge weighting: ``"tau"``,
      ``"rho"``, ``"hoeffd"``, ``"mcor"``, ``"cxi"`` or ``"joe"``, matching
      ``FitControlsVinecop``. ``"cxi"`` is Chatterjee's xi, which is
      asymmetric, so the weight is the larger of the two directions.
  threshold : float, default 0.0
      Dependence threshold. It acts twice, as it does in
      ``Vinecop.select``: an edge whose criterion falls below it is
      deprioritized during spanning-tree selection (weight ``1.0``), and if
      it survives anyway it is left holding
      :class:`~pyvinecopulib.core.IndependencePair` rather than being
      fitted. At the default no non-negative criterion is below it.
  tree_algorithm : str, default "mst_prim"
      ``"mst_prim"`` / ``"mst_kruskal"`` (Dissmann) or ``"random_weighted"`` /
      ``"random_unweighted"`` (Wilson); the MST variants maximize dependence.
  seeds : list of int, or None, optional
      RNG seeds for the random tree algorithms (ignored by the MST ones).
  var_types : list of str, or None, optional
      See ``fit_parts``. Given here it also fixes the dimension, so ``u`` may
      carry the extra left-limit columns.
  conditioning_set : list of int, or None, optional
      1-based variables to place at the tail of the selected order, so they can
      be conditioned on with :meth:`sample_conditional`. Every candidate edge
      touching a non-conditioning variable is penalized, which makes the
      conditioning set a self-contained block, and the finalized structure is
      then relabeled onto that tail. Requires an MST ``tree_algorithm``, and
      the pairs must implement
      :meth:`~pyvinecopulib.core.BicopBase.flip`.
  weights : array, shape (n,), or None, optional
      See ``fit_parts``.
  criterion_function : callable, or None, optional
      See ``fit_parts``. Edge weights read the unconditional
      pseudo-observations even under a non-simplified selection, so a
      ``tree_criterion`` that conditions is supplied through this.

  Returns
  -------
  structure : RVineStructure
      The selected vine structure.
  pair_copulas : list of list of BicopLike
      The fitted pair copulas, indexed ``[tree][edge]`` in the structure's
      column order and reoriented onto their slots — ready to host in a vine
      without re-fitting.
  conditioning_order : dict
      ``(tree, edge)`` to the 1-based variable labels of that slot's
      conditioning set, in the C1 order the pair sitting on it was fitted on.
      Install it alongside the pairs with ``_set_cond_order``; it is what a
      non-simplified vine gathers ``u_D`` in.

  See Also
  --------
  pyvinecopulib.core.VinecopBase.fit : Fit pair copulas along a fixed
      structure.
  pyvinecopulib.core.Vinecop.select : The reference (in-place) selector.

  References
  ----------
  .. [1] Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
     *Selecting and estimating regular vine copulae and application to
     financial returns.* Computational Statistics & Data Analysis, 59 (1),
     52-69.
  """
  from ..pyvinecopulib_ext import (
    RVineStructure,
    _select_spanning_tree,
  )

  xp = array_namespace(u)
  n = int(u.shape[0])
  ctx: ConditioningContext[ArrayT] = (
    SimplifiedContext() if context is None else context
  )
  x = prepare(u, x, n)
  # With left-limit columns present, `u` is wider than the vine: `var_types`
  # is what fixes the dimension.
  d = len(var_types) if var_types is not None else int(u.shape[1])
  types = check_var_types(var_types, d)
  offsets = disc_cols(types)
  u = collapse_data(u, d, types, "select")
  seed_list = [int(s) for s in (seeds or [])]
  max_trees = d - 1 if trunc_lvl is None else max(0, min(int(trunc_lvl), d - 1))
  tree_algorithms = (
    "mst_prim",
    "mst_kruskal",
    "random_weighted",
    "random_unweighted",
  )
  if tree_algorithm not in tree_algorithms:
    raise ValueError(
      f"tree_algorithm must be one of {tree_algorithms}; "
      f"got {tree_algorithm!r}."
    )
  cond = [int(v) for v in (conditioning_set or [])]
  in_cond = [False] * d
  if cond:
    # Mirrors `Vinecop::check_conditioning_set`: an MST is what makes the
    # penalty below lay the conditioning block down first.
    if len(cond) >= d:
      raise ValueError("conditioning_set must contain at most d - 1 variables.")
    if any(v < 1 or v > d for v in cond):
      raise ValueError("conditioning_set entries must be in 1, ..., d.")
    if tree_algorithm not in ("mst_prim", "mst_kruskal"):
      raise ValueError(
        "conditioning-aware selection requires an MST tree_algorithm "
        "('mst_prim' or 'mst_kruskal')."
      )
    for v in cond:
      in_cond[v - 1] = True
  # Sentinel "root" shared by every base-tree node, so the first tree's
  # candidate graph is complete (the C++ base tree is a star).
  root = d

  def selection_context(chain: tuple[int, ...]) -> Optional[ArrayT]:
    # x_e for an edge whose conditioning set is `chain`, in the C1 order it
    # is fitted on -- variable indices into `u`, since selection has no
    # structure to read natural-order columns from yet.
    if not ctx.assembles_conditioning and x is None:
      return None
    u_D = u[:, list(chain)] if chain and ctx.assembles_conditioning else None
    return ctx.edge_context(u_D=u_D, x=x)

  # Only the value columns enter, so a discrete vine selects the tree it
  # would select continuous.
  # Checked here for the same reason as in `_fit_parts`: the criterion hands
  # them straight to the binding.
  weights = validate_weights(weights, u[:, 0])
  criterion = _make_criterion(tree_criterion, n, weights, criterion_function, x)

  # A node is one edge of the previous tree (a single variable for the base
  # tree). ``prev`` holds the two previous-tree vertex ids that this edge
  # joined; a shared prev id is the proximity condition and picks which
  # h-function feeds the next tree.
  # A base-tree vertex is a single variable, so both of its slots carry that
  # variable's type and its left limit fills the first slot only -- the base
  # tree is a star, so every edge reads slot 0 (`make_base_tree`).
  nodes: list[dict[str, Any]] = [
    {
      "all_indices": (i,),
      "h1": u[:, i],
      "h2": u[:, i],
      "h1_sub": u[:, d + offsets[i]] if types[i] == "d" else None,
      "h2_sub": None,
      "types": ("d", "d") if types[i] == "d" else ("c", "c"),
      "prev": (root, i),
      # Per endpoint: the conditioning chain an edge built on this node
      # inherits when that endpoint is its diagonal, extended by the
      # variable this node conditioned away. A base-tree node conditions on
      # nothing, so a first-tree edge starts from an empty chain.
      "ext": {i: ()},
    }
    for i in range(d)
  ]

  trees: list[list[tuple[int, int, list[int]]]] = []
  flip_checked = False
  # Per tree: {(conditioned pair, conditioning set) -> (arg1 label, fitted
  # pair, the chain it was fitted on)}, used to place + reorient the pairs
  # onto the finalized slots and to record what each conditions on.
  records: list[dict[_SlotKey, tuple[int, Any, tuple[int, ...]]]] = []
  for _ in range(max_trees):
    m = len(nodes)
    cand: list[tuple[int, int]] = []
    cand_cols: list[tuple[Any, Any]] = []
    cand_subs: list[Optional[tuple[Any, Any]]] = []
    cand_types: list[tuple[str, str]] = []
    cand_crits: list[float] = []
    edge_costs: list[float] = []
    # Candidate enumeration mirrors the C++ selector exactly
    # (tools_select.ipp add_allowed_edges_proximity): the outer loop runs
    # over v0 and the inner over v1 < v0, so an edge's *first* endpoint v0 —
    # which contributes pc_data column 0 and the first conditioned variable —
    # is the larger vertex index, and candidate insertion order is preserved.
    for v0 in range(m):
      prev0 = nodes[v0]["prev"]
      for v1 in range(v0):
        prev1 = nodes[v1]["prev"]
        shared = set(prev0) & set(prev1)
        if not shared:
          continue
        common = min(shared)
        pos0, pos1 = prev0.index(common), prev1.index(common)
        col0 = nodes[v0]["h1"] if pos0 == 0 else nodes[v0]["h2"]
        col1 = nodes[v1]["h1"] if pos1 == 0 else nodes[v1]["h2"]
        # The h-function comes from slot `pos`, the type from the *other*
        # slot: an h-function integrates out its conditioning variable and
        # keeps the other one (`add_pc_info`).
        edge_types = (
          nodes[v0]["types"][1 - pos0],
          nodes[v1]["types"][1 - pos1],
        )
        subs: Optional[tuple[Any, Any]] = None
        if "d" in edge_types:
          # A slot without a left limit is continuous, and its own value is
          # its left limit (`get_hfunc_sub`).
          sub0 = nodes[v0]["h1_sub" if pos0 == 0 else "h2_sub"]
          sub1 = nodes[v1]["h1_sub" if pos1 == 0 else "h2_sub"]
          subs = (
            col0 if sub0 is None else sub0,
            col1 if sub1 is None else sub1,
          )
        # The edge weight reads the value columns only, so the spanning tree
        # a discrete vine selects is the one it would select continuous.
        tau = criterion(col0, col1)
        weight = 1.0 - (tau >= threshold) * tau
        if cond:
          # Base weights lie in [0, 1], so adding `d` keeps them non-negative
          # (Prim requires it) while making every all-conditioning edge
          # strictly cheaper: the minimum spanning tree lays down the
          # conditioning set's own optimal sub-vine first at every tree, which
          # is what makes it a block the relabeling can move to the tail
          # (tools_select.ipp add_allowed_edges_proximity).
          all_cond = all(in_cond[i] for i in nodes[v0]["all_indices"]) and all(
            in_cond[i] for i in nodes[v1]["all_indices"]
          )
          if not all_cond:
            weight += float(d)
        cand.append((v0, v1))
        cand_cols.append((col0, col1))
        cand_subs.append(subs)
        cand_types.append(edge_types)
        cand_crits.append(float(tau))
        edge_costs.append(weight)

    # Ascending candidate index = boost's edge-list (insertion) order, which
    # is the order the C++ selector iterates surviving edges in.
    selected = sorted(
      _select_spanning_tree(m, cand, edge_costs, tree_algorithm, seed_list)
    )

    tree_edges: list[tuple[int, int, list[int]]] = []
    new_nodes: list[dict[str, Any]] = []
    tree_records: dict[_SlotKey, tuple[int, Any, tuple[int, ...]]] = {}
    build_next_level = len(trees) + 1 < max_trees and len(selected) > 1
    # The surviving edges' inputs were all materialized above, off the
    # previous tree's nodes, and each fitted pair is written to a fresh
    # `new_nodes` -- so unlike `fit`, this level needs no reordering to be
    # fitted in one call.
    survivors = [
      stack_edge(xp, cand_cols[e][0], cand_cols[e][1], cand_subs[e])
      for e in selected
    ]
    level_types = [cand_types[e] for e in selected]
    # Each surviving edge's conditioned pair and its conditioning chains,
    # resolved before anything is fitted: a conditional pair's conditioning
    # matrix is part of its input, so it needs the same lookahead
    # `survivors` does. `a_var` is v0's unique variable — the fitted pair's
    # first argument — and `b_var` v1's, in the C++ `set_sym_diff` order.
    # Extending an endpoint node's chain by the variable it conditioned away
    # reads the finalized column downwards, which is what makes the chain
    # the C1 order for the slot that endpoint ends up the diagonal of.
    conditioned: list[tuple[int, int, list[int]]] = []
    chains: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for e in selected:
      v0, v1 = cand[e]
      idx0 = set(nodes[v0]["all_indices"])
      idx1 = set(nodes[v1]["all_indices"])
      a_var = next(iter(idx0 - idx1))
      b_var = next(iter(idx1 - idx0))
      conditioned.append((a_var, b_var, sorted(idx0 & idx1)))
      chains.append((nodes[v0]["ext"][a_var], nodes[v1]["ext"][b_var]))
    # The pair is fitted in the orientation the search built it in, so the
    # `a` chain is the order it is estimated on and the one recorded with it.
    contexts = [selection_context(chain_a) for chain_a, _ in chains]
    # The weight above only decides which edges survive. Whether a surviving
    # edge is *fitted* is a second question with the same answer upstream
    # gives: an edge whose criterion falls below the threshold keeps a
    # default-constructed pair -- independence -- and `select` is never
    # called on it (tools_select.ipp fit_or_reuse_pair_copula). At the
    # default `threshold=0.0` no non-negative criterion is below it, so
    # nothing here is thresholded.
    thresholded = [cand_crits[e] < threshold for e in selected]
    to_fit = [i for i, skip in enumerate(thresholded) if not skip]
    fitted_level: Optional[dict[int, BicopLike[Any]]] = None
    if (
      fit_level is not None
      and to_fit
      and all(contexts[i] is None for i in to_fit)
      and all("d" not in level_types[i] for i in to_fit)
    ):
      got = fit_level(
        len(trees),
        xp.stack([survivors[i] for i in to_fit], axis=0),
        [level_types[i] for i in to_fit],
      )
      fitted_level = dict(zip(to_fit, got))
    for edge_idx, e in enumerate(selected):
      v0, v1 = cand[e]
      subs, edge_types = cand_subs[e], cand_types[e]
      u_e = survivors[edge_idx]
      x_e = contexts[edge_idx]
      pair: BicopLike[ArrayT]
      if thresholded[edge_idx]:
        pair = IndependencePair()
      elif fitted_level is not None:
        pair = fitted_level[edge_idx]
      else:
        pair = _fit_edge_call(
          fit_edge, len(trees), edge_idx, u_e, x_e, edge_types
        )
      if not flip_checked and not thresholded[edge_idx]:
        # `_check_selectable` settles this up front when the vine names a
        # `bicop_class`; behind a caller's own `fit_edge` the class is not
        # knowable until one pair exists, so probe that one rather than
        # discovering it after every edge has been fitted. A thresholded
        # edge is not one of theirs -- `IndependencePair.flip` returns
        # `self` and would pass the probe for them.
        flip_checked = True
        try:
          flip_of(pair)
        except NotImplementedError as err:
          raise NotImplementedError(
            f"{type(pair).__name__} has no `flip`, which structure "
            "selection needs to reorient each pair onto its finalized "
            "slot. Implement it (return the argument-swapped copula), or "
            "supply a structure and fit along it instead."
          ) from err
      a_var, b_var, conditioning = conditioned[edge_idx]
      chain_a, chain_b = chains[edge_idx]
      tree_edges.append((a_var + 1, b_var + 1, [c + 1 for c in conditioning]))
      key = (
        frozenset((a_var + 1, b_var + 1)),
        frozenset(c + 1 for c in conditioning),
      )
      tree_records[key] = (
        a_var + 1,
        pair,
        tuple(c + 1 for c in chain_a),
      )
      if build_next_level:
        new_nodes.append(
          {
            "all_indices": tuple(sorted((a_var, b_var, *conditioning))),
            "h1": pair_eval(pair.hfunc1, u_e, x_e),
            "h2": pair_eval(pair.hfunc2, u_e, x_e),
            # A discrete argument's next tree needs the h-function at the
            # atom's lower end too, exactly as the cascades compute it.
            "h1_sub": (
              pair_eval(pair.hfunc1, with_left_limit(u_e, 1), x_e)
              if edge_types[1] == "d"
              else None
            ),
            "h2_sub": (
              pair_eval(pair.hfunc2, with_left_limit(u_e, 0), x_e)
              if edge_types[0] == "d"
              else None
            ),
            "types": edge_types,
            "prev": (v0, v1),
            "ext": {a_var: chain_a + (b_var,), b_var: chain_b + (a_var,)},
          }
        )
    trees.append(tree_edges)
    records.append(tree_records)
    nodes = new_nodes
    if not build_next_level:
      break

  # Selection finalization via the shared list-of-trees primitive.
  # ``Vinecop.select`` and ``RVineStructure.from_trees`` share one diagonal
  # convention (conditioned[0], flip-free), so this reproduces the compiled
  # selector's matrix exactly.
  structure = RVineStructure.from_trees(d, trees)
  if cond:
    # Place the conditioning set at the tail, as `Vinecop.select` does after
    # finalizing. The placement below keys off each slot's label sets rather
    # than the diagonal policy, so it lands -- and flips -- the selection-time
    # pairs on the relabeled slots with no further bookkeeping.
    structure = reorientation(structure, cond).structure
  # Place each selection-time pair onto its finalized slot, mirroring the
  # peel: the slot at column ``e`` of tree ``t`` hosts the (unique) edge
  # whose conditioned pair is {order[e], struct_array(t, e)} and whose
  # conditioning set is {struct_array(0..t-1, e)}; the pair is flipped iff
  # the diagonal variable differs from its first argument ``a``
  # (rvine_trees.ipp peel).
  pairs: list[list[BicopLike[ArrayT]]] = []
  cond_order: dict[tuple[int, int], tuple[int, ...]] = {}
  for t in range(int(structure.trunc_lvl)):
    row: list[BicopLike[ArrayT]] = []
    for e in range(d - 1 - t):
      diag, key = _slot_key(structure, t, e)
      a_label, pair, chain = records[t][key]
      row.append(pair if a_label == diag else flip_of(pair))
      # `flip` swaps the pair's two arguments; nothing reorders the
      # conditioning columns, so the slot inherits the order its pair was
      # fitted on. On an unswapped slot that is the order the finalized
      # matrix names; on a swapped one it is the other endpoint's chain,
      # which is why it has to be carried rather than re-derived.
      if t:
        cond_order[(t, e)] = chain
    pairs.append(row)
  return structure, pairs, cond_order
