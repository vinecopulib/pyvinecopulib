"""Canonical partial implementation of the vine-copula evaluator contract.

:class:`VinecopBase` is the array-agnostic (NumPy / PyTorch, via
:func:`array_api_compat.array_namespace`) implementation of the vine cascades —
``pdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` / ``sample`` / ``cdf`` — plus
``loglik`` / ``plot`` and the shared sequential-fit engine. It walks the vine
tree by tree, evaluating one pair copula per edge, so a concrete subclass (e.g.
:class:`~pyvinecopulib.core.Vinecop`'s torch counterpart
:class:`~pyvinecopulib.torch.TorchVinecop`) only supplies a small set of hooks
and inherits the whole evaluation surface.

Conditioning is threaded through a pluggable
:class:`~pyvinecopulib.core.ConditioningContext`: each pair-copula call receives an
``x`` matrix (``x_e`` in the cascades) assembled per edge from the edge's
conditioning-set values ``u_D`` and the optional external covariate matrix
``x``. The default :class:`~pyvinecopulib.core.SimplifiedContext` forwards only
``x`` (``None`` in the unconditional case) and skips the ``u_D`` gather,
reproducing the classic simplified cascade at zero extra cost.

The ``pdf`` cascade accumulates the vine density as a **product** of per-edge
copula densities (there is no per-observation ``log_pdf``; ``loglik`` sums the
log-density).

Discrete variables are declared through ``_bind_vine(..., var_types=...)``. A
declared variable with atoms makes the forward cascades carry a **parallel
left-limit scratch** alongside every h-function, so each pair copula whose own
variables include a discrete one receives a four-column
``u_e = [u1, u2, u1^-, u2^-]`` — the layout ``Bicop`` consumes. Which pair sees
which types is fixed by the structure (:meth:`VinecopBase.pair_var_types`), and
data enters in the expanded ``(n, 2d)`` or compact ``(n, d + k)`` layout that
``Vinecop`` accepts. :meth:`VinecopBase.fit` and :meth:`VinecopBase.select` take
``var_types`` too, handing each edge's types to the ``fit_edge`` callback so it
can fit the pair copula the edge actually needs. A pair copula that reads only
two columns is hosted on a discrete edge by wrapping it in
:class:`~pyvinecopulib.core.DiscretePair`, which supplies the difference
quotients from its continuous ``pdf`` / ``cdf`` / ``hfunc1`` / ``hfunc2``.

The only hook a concrete subclass must provide is ``_get_pair_copula``; ``_prep``
(input coercion + unit-box clamp) ships a concrete default, and
``_sample_uniform`` (RNG) is needed only to enable ``sample``. To enable the
grid-batched fast path, override ``_build_batched`` (plus ``_default_batched``).
The batched *cascade loops* are
array-agnostic and live here; only the grid/cache builder returned by
``_build_batched`` is subclass-specific. The concrete subclass calls
``_bind_vine`` once to install the structure and context. Array values are
handled as ``Any`` inside the cascades per the ``pyvinecopulib.core`` typing
policy (the Array API namespace is untyped); the generic ``ArrayT`` lives on the
public signatures.
"""

from __future__ import annotations

import contextlib
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Callable, Optional, cast

from array_api_compat import array_namespace

from ._discrete import (
  check_var_types,
  collapse_data,
  continuous_view,
  disc_cols,
  edge_columns,
  n_discrete,
  pair_var_types,
  seed_left_limits,
  stack_edge,
  with_left_limit,
)
from ._reorient import Reorientation, reorientation
from .bicop_base import _pair_eval
from .context import ConditioningContext, SimplifiedContext
from .._deprecations import _reject_renamed_hook
from .protocols import ArrayT, BicopLike, VinecopLike, _VINECOP_EXAMPLE
from ._trim import trim


def _to_numpy_default(a: Any) -> Any:
  """Host NumPy view of ``a``, for backends that leave it off the host."""
  import numpy as _np

  detach = getattr(a, "detach", None)
  if detach is not None:
    a = detach()
  cpu = getattr(a, "cpu", None)
  if cpu is not None:
    a = cpu()
  return _np.asarray(a)


if TYPE_CHECKING:
  from ..pyvinecopulib_ext import RVineStructure

__all__ = ["VinecopBase"]


#: ``(tree, edge, u_e, x_e) -> BicopLike``, fitting one edge's pair copula: the
#: seam external packages drive conditional fitting through (see
#: :meth:`VinecopBase.fit`). An edge with a discrete argument additionally
#: receives ``var_types=[t1, t2]`` and a four-column ``u_e``, so the alias cannot
#: pin the arity -- a ``Callable`` has no way to express a keyword argument.
FitEdge = Callable[..., BicopLike]


def _fit_edge_call(
  fit_edge: FitEdge,
  tree: int,
  edge: int,
  u_e: Any,
  x_e: Optional[Any],
  var_types: tuple[str, str],
) -> BicopLike:
  """Call ``fit_edge``, forwarding ``var_types`` only for a discrete edge.

  A callback written for a continuous vine takes four arguments, so a fully
  continuous edge must not hand it a fifth. When the edge does have a discrete
  argument the types go by keyword, which is what makes a callback that cannot
  accept them fail loudly instead of fitting a continuous pair copula to four
  columns of data -- the rule :func:`_pair_eval` applies to ``x``, for the same
  reason.
  """
  if "d" not in var_types:
    return fit_edge(tree, edge, u_e, x_e)
  return fit_edge(tree, edge, u_e, x_e, var_types=list(var_types))


class _NotBatchable(Exception):
  """Raised by :meth:`VinecopBase._build_batched` when batching is unavailable.

  The dispatch layer catches it and falls back to the non-batched cascade.
  """


def infer_conditioning_set(
  order: list[int], var_types: list[str], n_cols: int
) -> list[int]:
  """Which variables a ``u_cond`` of this width conditions on.

  The convention, when the caller names no set: the last ``k`` variables of the
  sampling order, with ``k`` recovered from the column count -- a discrete
  conditioner contributes a left-limit column as well. Shared by the copula
  scale and by ``Vinedist``, so a caller sees one rule on both.

  Parameters
  ----------
  order : list of int
      The sampling order, 1-based.
  var_types : list of str
      Per-variable types in variable order.
  n_cols : int
      Width of the conditioning block.

  Returns
  -------
  list of int
      The 1-based conditioning variables, in column order.

  Raises
  ------
  ValueError
      If no order tail has a layout that wide.
  """
  d = len(order)
  for k in range(1, d):
    tail = list(order[d - k :])
    n_disc = sum(1 for v in tail if var_types[v - 1] == "d")
    if n_cols == k + n_disc:
      return tail
  raise ValueError(
    "u_cond has an invalid number of columns; expected k columns for "
    "all-continuous conditioning, or k + k_d columns when k_d conditioning "
    f"variables are discrete, for some k in 1, ..., {d - 1}; got {n_cols}."
  )


class VinecopBase(VinecopLike[ArrayT], ABC):
  """Canonical array-agnostic vine cascades (numpy / torch).

  Concrete subclasses implement ``_get_pair_copula`` (and optionally override
  ``_prep`` / ``_sample_uniform`` / the batched-path hooks) and call
  ``_bind_vine`` once; they
  then inherit the whole evaluator surface — ``pdf`` / ``cdf`` / ``rosenblatt`` /
  ``inverse_rosenblatt`` / ``sample``, ``loglik`` / ``plot`` / ``__repr__``,
  the ``dim`` / ``trunc_lvl`` / ``order`` accessors, and the :meth:`fit` /
  :meth:`select` engines. Not an ``nn.Module``, so it composes with any
  pair-copula implementation (including non-torch ones).

  See Also
  --------
  pyvinecopulib.core.VinecopLike : The contract this implements.
  pyvinecopulib.core.Vinecop : The reference vine.
  pyvinecopulib.core.ConditioningContext : Per-edge conditioning policy.
  """

  # --- layout installed by _bind_vine (hooks / state) ------------------- #
  structure: RVineStructure
  d: int
  trunc_lvl: int
  order: tuple[int, ...]
  inverse_order: tuple[int, ...]

  #: Whether this vine's pairs read the external covariates ``x``. Declared, not
  #: inferred: a vine assembles each edge's conditioning matrix through its
  #: ``ConditioningContext`` and cannot know whether the pairs it hosts accept
  #: one, so a consumer such as ``Vinedist`` reads this flag rather than
  #: forwarding covariates a compiled ``Bicop`` pair would refuse. Set it in a
  #: subclass whose pairs are conditional.
  supports_covariates: bool = False

  def __init_subclass__(cls, **kwargs: Any) -> None:
    """Reject an override of the pre-1.0 hook name.

    Parameters
    ----------
    **kwargs
        Forwarded to ``super().__init_subclass__``.
    """
    super().__init_subclass__(**kwargs)
    _reject_renamed_hook(cls, "_simulate_uniform", "_sample_uniform")

  _context: ConditioningContext
  _cond_pos_cache: dict[tuple[int, int], tuple[int, ...]]
  #: Lazily-built grid-batched state (see :meth:`_build_batched`); ``None`` until
  #: the first batched call. Subclasses invalidate it on device moves.
  _batched: Any
  #: Array namespace of this vine's working arrays; ``None`` until resolved.
  _xp: Any
  _var_types: tuple[str, ...]
  _n_discrete: int
  #: Variable index -> its offset within the compact layout's left-limit block;
  #: meaningful only at discrete variables.
  _disc_cols: tuple[int, ...]
  _pair_types: tuple[tuple[tuple[str, str], ...], ...]

  def _bind_vine(
    self,
    structure: RVineStructure,
    context: Optional[ConditioningContext] = None,
    var_types: Optional[list[str]] = None,
  ) -> None:
    """Install the vine structure + context and derive the order arrays.

    The initialization seam a concrete subclass calls once from its ``__init__``
    (after storing its pair copulas). ``context`` defaults to
    :class:`~pyvinecopulib.core.SimplifiedContext` — the unconditional /
    simplified vine that covers the common case — so most subclasses pass only a
    ``structure``. Advanced subclasses may override this method to install extra
    state, calling ``super()._bind_vine(...)``.

    ``var_types`` is what makes the vine discrete-aware. It is *not* pushed onto
    the pair copulas: :class:`~pyvinecopulib.core.BicopLike` declares no
    attributes and the pairs belong to the subclass, not to ``VinecopBase``. The
    subclass reads :meth:`pair_var_types` to configure each pair it hosts, and
    the cascades hand a pair whose types include ``"d"`` the four-column
    ``[u1, u2, u1^-, u2^-]`` input.

    Parameters
    ----------
    structure : RVineStructure
        The (fixed) vine structure to evaluate along.
    context : ConditioningContext, optional
        Per-edge conditioning-context policy; ``None`` uses
        :class:`~pyvinecopulib.core.SimplifiedContext`.
    var_types : list of str, optional
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete), in
        variable order; ``None`` means all continuous.

    Returns
    -------
    None
        The structure, context, and derived order arrays are stored on ``self``.

    Raises
    ------
    ValueError
        If ``var_types`` has the wrong length or an entry outside
        ``{"c", "d"}``.
    """
    if context is None:
      context = SimplifiedContext()
    self.structure = structure
    self._context = context
    self.d = int(structure.dim)
    self.trunc_lvl = int(structure.trunc_lvl)
    order = tuple(int(v) for v in structure.order)
    self.order = order
    inv = [0] * len(order)
    for j, k in enumerate(order):
      inv[k - 1] = j
    self.inverse_order = tuple(inv)
    self._cond_pos_cache = {}
    self._batched = None
    self._xp = None
    self._bind_var_types(var_types)

  def _bind_var_types(self, var_types: Optional[list[str]]) -> None:
    """Store the variable types and derive the per-edge type table."""
    types = check_var_types(var_types, self.d)
    self._var_types = types
    self._n_discrete = n_discrete(types)
    self._disc_cols = disc_cols(types)
    self._pair_types = pair_var_types(self.structure, types)

  @property
  def var_types(self) -> list[str]:
    """Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete).

    Returns
    -------
    list of str
        One entry per variable, in variable order; all ``"c"`` unless
        ``var_types`` was declared when the vine was bound.
    """
    return list(self._var_types)

  def pair_var_types(self, tree: int, edge: int) -> tuple[str, str]:
    """Variable types of the pair copula hosted at ``(tree, edge)``.

    Derived from :attr:`var_types` and the structure, so it matches what
    ``Vinecop`` assigns to the same slot. A subclass that owns its pair copulas
    uses this to configure them; a pair reporting ``"d"`` for either variable is
    handed the four-column ``[u1, u2, u1^-, u2^-]`` input by the cascades.

    Parameters
    ----------
    tree : int
        Tree index (``0``-based).
    edge : int
        Edge index within the tree (``0``-based).

    Returns
    -------
    tuple of str
        The pair's ``(type of first variable, type of second variable)``.
    """
    return self._pair_types[tree][edge]

  # --- hooks a concrete subclass provides ------------------------------- #
  @abstractmethod
  def _get_pair_copula(self, tree: int, edge: int) -> BicopLike[ArrayT]:
    """Return the pair copula at ``(tree, edge)`` (the one required hook).

    Parameters
    ----------
    tree : int
        Tree index (``0``-based).
    edge : int
        Edge index within the tree (``0``-based).

    Returns
    -------
    BicopLike
        The pair copula hosted at that position.
    """

  def _prep(self, u: ArrayT, name: str, *, values_only: bool = False) -> ArrayT:
    """Coerce ``u`` to the working array, normalize its layout, and clamp it.

    Concrete default: accept any layout the vine's :attr:`var_types` admits,
    reduce it to the compact ``(n, d + k)`` form (``k`` discrete variables), and
    clamp to ``[1e-10, 1 - 1e-10]`` on ``u``'s own array namespace. A subclass
    that needs dtype / device coercion (e.g. accepting NumPy input on a torch
    vine) overrides this.

    An all-continuous vine takes ``(n, d)``, or ``(n, 2d)`` whose left-limit
    block is dropped. With ``k`` discrete variables it takes the expanded
    ``(n, 2d)`` or the compact ``(n, d + k)``; the plain ``(n, d)`` is rejected,
    because silently reusing each value as its own left limit would evaluate a
    continuous density under a discrete model.

    Parameters
    ----------
    u : array, shape (n, d), (n, d + k) or (n, 2d), dtype float
        Pseudo-observations to prepare.
    name : str
        Calling-method name, used only in the shape-error message.
    values_only : bool, default=False
        Return just the ``d`` value columns, and accept a plain ``(n, d)`` input
        even on a discrete vine. Set by the cascades that never read a left
        limit.

    Returns
    -------
    array, shape (n, d + k) or (n, d), dtype float
        ``u`` coerced to the working array, reduced to the compact layout (to
        the ``d`` value columns when ``values_only``), and clamped to
        ``[1e-10, 1 - 1e-10]``.

    Raises
    ------
    ValueError
        If ``u`` is not 2-d or its column count matches no accepted layout.
    """
    ua: Any = u
    xp = self._namespace(ua)
    return cast(ArrayT, trim(xp, self._layout(ua, name, values_only)))

  def _layout(self, ua: Any, name: str, values_only: bool) -> Any:
    """Validate ``ua``'s layout and reduce it to the columns the caller needs."""
    return collapse_data(
      ua, self.d, self._var_types, name, values_only=values_only
    )

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> ArrayT:
    """Draw ``(n, d)`` base uniforms for :meth:`sample` (namespace-dependent RNG).

    Raising default; override it (numpy / torch differ on RNG) to enable
    :meth:`sample`. Named after :func:`pyvinecopulib.utils.sample_uniform`.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    qrng : bool
        Whether to draw a quasi-random (low-discrepancy) sequence.
    seeds : list of int
        RNG seeds.

    Returns
    -------
    array, shape (n, d), dtype float
        Base uniforms in ``[0, 1)``.

    Raises
    ------
    NotImplementedError
        Unless a subclass overrides this hook.
    """
    raise NotImplementedError(
      f"{type(self).__name__} does not implement _sample_uniform; override it "
      "to enable sample()."
    )

  def _default_batched(self) -> bool:
    """Whether ``batched`` defaults to ``True`` (subclass- and device-dependent).

    Returns
    -------
    bool
        The value used when a caller passes ``batched=None`` (``False`` here; a
        grid subclass may key it on the device).
    """
    return False

  def _build_batched(self) -> Any:
    """Build the grid-batched state for the fast path (subclass-specific).

    The default raises :class:`_NotBatchable`, so the dispatch layer falls back
    to the non-batched cascade. A grid subclass overrides this to return an
    object exposing the batched-vine surface the cascades call
    (``level`` / ``grid_points`` / per-level ``gather_inputs`` / ``pdf`` /
    ``hfunc1`` / ``hfunc2`` / ``n_pairs`` / ``needs_h1`` / ``needs_h2``).

    Returns
    -------
    object
        The subclass-specific batched-vine state the batched cascades run on.

    Raises
    ------
    _NotBatchable
        In the default implementation (no grid fast path available).
    """
    raise _NotBatchable(
      f"{type(self).__name__} does not provide a batched fast path."
    )

  def _ensure_batched(self) -> Any:
    """Return the cached batched state, building it once on first use.

    Returns
    -------
    object
        The memoized :meth:`_build_batched` result.
    """
    if self._batched is None:
      # Bypass any framework `__setattr__`: on the torch subclass this value
      # is an `nn.Module`, and a normal assignment would register it as a
      # child, so a derived cache would leak into `state_dict()` and a
      # checkpoint taken after a batched call would not load into a fresh
      # model. It is a memo of the pair copulas, not part of the model.
      object.__setattr__(self, "_batched", self._build_batched())
    return self._batched

  def _eval_context(self):
    """Context manager disabling grad for inverse / sample / cdf.

    Defaults to a no-op; the torch subclass overrides it with
    ``torch.no_grad()``.

    Returns
    -------
    contextlib.AbstractContextManager
        A context manager wrapping the grad-sensitive sections (a
        ``nullcontext`` by default).
    """
    return contextlib.nullcontext()

  # --- conditioning-context assembly ------------------------------------ #
  def _cond_positions(self, tree: int, edge: int) -> tuple[int, ...]:
    """Natural-order column indices of this edge's conditioning set ``D``.

    These are ``struct_array(i, edge, natural_order=True) - 1`` for the
    conditioning trees ``i = 0 .. tree - 1``, i.e. the columns of the
    natural-order observation matrix holding the variables the pair copula
    ``c_{a,b;D}`` conditions on. They are returned in **ascending conditioning
    -tree order** ``i`` — the C1 column order the context assembles ``u_D`` in
    and a conditional pair consumes positionally. All are ``> edge`` (the
    inverse-cascade invariant), so they are finalized before this edge is read.

    Parameters
    ----------
    tree : int
        Tree index (``0``-based).
    edge : int
        Edge index within the tree (``0``-based).

    Returns
    -------
    tuple of int
        Natural-order column indices of the conditioning variables, in ascending
        conditioning-tree order (the C1 order).
    """
    key = (tree, edge)
    cache = self._cond_pos_cache
    if key not in cache:
      s = self.structure
      cache[key] = tuple(
        int(s.struct_array(i, edge, natural_order=True)) - 1
        for i in range(tree)
      )
    return cache[key]

  def _edge_context(
    self,
    tree: int,
    edge: int,
    x: Optional[Any],
    u_nat: Optional[Any],
    hinv2_final: Optional[Any],
  ) -> Optional[Any]:
    """Assemble the per-edge conditioning context ``x_e`` = context(``u_D``, ``x``).

    ``u_D`` is gathered in the C1 column order (see :meth:`_cond_positions`),
    then the context appends the external covariates ``x`` last. The source of
    ``u_D`` differs by direction: the forward cascades pass ``u_nat`` (the
    natural-order observations, columns), the inverse cascade passes
    ``hinv2_final`` (the finalized ``hinv2[0]`` rows, transposed).

    Parameters
    ----------
    tree : int
        Tree index (``0``-based).
    edge : int
        Edge index within the tree (``0``-based).
    x : array, shape (n, p), or None
        External covariates for this call, or ``None``.
    u_nat : array, shape (n, d), or None
        Natural-order observations (forward cascades); ``None`` in the inverse
        direction.
    hinv2_final : array, shape (d, n), or None
        The finalized ``hinv2[0]`` scratch rows (inverse cascade); ``None`` in
        the forward direction.

    Returns
    -------
    array, shape (n, k), or None
        The per-edge conditioning matrix ``x_e``, or ``None`` when the pair is
        unconditional for this edge.
    """
    ctx = self._context
    # Simplified + unconditional: skip the gather entirely (zero cost).
    if not ctx.assembles_conditioning and x is None:
      return None
    u_D: Optional[Any] = None
    if ctx.assembles_conditioning:
      cols = list(self._cond_positions(tree, edge))
      if cols:
        if u_nat is not None:  # forward: read observation columns
          u_D = u_nat[:, cols]
        else:  # inverse: read finalized hinv2[0] rows, transpose to (n, |D|)
          finalized: Any = hinv2_final
          xp = array_namespace(finalized)
          u_D = xp.matrix_transpose(finalized[cols, :])
    return ctx.edge_context(u_D=u_D, x=x)

  # --- discrete left-limit scratch -------------------------------------- #
  def _seed_sub(self, u: Any, xp: Any) -> Optional[Any]:
    """Natural-order left limits from the compact layout, or ``None``.

    ``None`` for an all-continuous vine, which is what switches the whole
    left-limit cascade off. A continuous variable's column holds its own value:
    a pair only ever reads the left-limit column of a variable it declares
    discrete, and this keeps the four-column edge input well defined anyway.
    """
    return seed_left_limits(
      u, self.d, self.order, self._var_types, self._disc_cols, xp
    )

  def _edge_columns(
    self,
    tree: int,
    edge: int,
    hfunc1: Any,
    hfunc2: Any,
    hfunc1_sub: Optional[Any],
    hfunc2_sub: Optional[Any],
  ) -> tuple[Any, Any, Optional[tuple[Any, Any]], tuple[str, str]]:
    """Resolve one edge's pair-copula input columns and its variable types.

    ``m`` is the min-array entry: the natural-order index of the column
    finalized in a previous tree. The second pair input comes from ``hfunc2``
    when ``m`` sits on the natural-order diagonal, else from ``hfunc1``
    (``class.ipp:1026-1034``). The left-limit pair is returned only when the
    edge has a discrete variable, and mirrors ``Bicop::format_data``: a
    continuous variable's left limit is its own value.
    """
    return edge_columns(
      self.structure,
      self._pair_types,
      tree,
      edge,
      hfunc1,
      hfunc2,
      hfunc1_sub,
      hfunc2_sub,
    )

  # --- non-batched cascades (single source of truth) -------------------- #
  def _pdf(self, u: Any, x: Optional[Any]) -> Any:
    """Vine density as a product of per-edge copula densities (``Vinecop::pdf``).

    Parameters
    ----------
    u : array, shape (n, d + k), dtype float
        Prepared pseudo-observations in the compact layout (natural-order
        seeding happens inside).
    x : array, shape (n, p), or None
        External covariates threaded to each pair copula, or ``None``.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values.
    """
    xp = array_namespace(u)
    d, trunc_lvl = self.d, self.trunc_lvl
    n = u.shape[0]
    if trunc_lvl == 0:
      return xp.ones(n, dtype=u.dtype, device=u.device)
    # Dense (n, d) h-function scratch; seed hfunc2 with the observations in
    # natural order (class.ipp:399).
    hfunc1 = xp.zeros((n, d), dtype=u.dtype, device=u.device)
    hfunc2 = xp.empty((n, d), dtype=u.dtype, device=u.device)
    order = self.order
    for j in range(d):
      hfunc2[:, j] = u[:, order[j] - 1]
    # Parallel left-limit scratch, allocated only for a discrete vine. hfunc1_sub
    # needs no seed: tree 0 always reads its second input on the diagonal, so
    # every entry is written before it is read.
    hfunc2_sub: Any = self._seed_sub(u, xp)
    hfunc1_sub: Any = (
      None
      if hfunc2_sub is None
      else xp.zeros((n, d), dtype=u.dtype, device=u.device)
    )
    # Keep an immutable copy of the seeded observations for conditioning-set
    # (u_D) gathers; skipped entirely under a simplified/unconditional vine.
    u_nat = (
      xp.asarray(hfunc2, copy=True)
      if self._context.assembles_conditioning
      else None
    )
    pdf = xp.ones(n, dtype=u.dtype, device=u.device)
    s = self.structure
    for tree in range(trunc_lvl):
      for edge in range(d - tree - 1):
        edge_copula = self._get_pair_copula(tree, edge)
        col0, col1, subs, types = self._edge_columns(
          tree, edge, hfunc1, hfunc2, hfunc1_sub, hfunc2_sub
        )
        u_e = xp.stack(
          [col0, col1] if subs is None else [col0, col1, *subs], axis=-1
        )
        x_e = self._edge_context(tree, edge, x, u_nat, None)
        # Accumulate the density as a product over edges (cwiseProduct,
        # class.ipp:1047).
        pdf = pdf * _pair_eval(edge_copula.pdf, u_e, x_e)
        # h-functions only evaluated if a later tree needs them (class.ipp:1050).
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = _pair_eval(edge_copula.hfunc1, u_e, x_e)
          if subs is not None and types[1] == "d":
            u_h1 = xp.stack([col0, subs[1], *subs], axis=-1)
            hfunc1_sub[:, edge] = _pair_eval(edge_copula.hfunc1, u_h1, x_e)
        if s.needed_hfunc2(tree, edge):
          hfunc2[:, edge] = _pair_eval(edge_copula.hfunc2, u_e, x_e)
          if subs is not None and types[0] == "d":
            u_h2 = xp.stack([subs[0], col1, *subs], axis=-1)
            hfunc2_sub[:, edge] = _pair_eval(edge_copula.hfunc2, u_h2, x_e)
    return pdf

  def _rosenblatt(
    self,
    u: Any,
    x: Optional[Any],
    randomize_discrete: bool = True,
    seeds: Optional[list[int]] = None,
  ) -> Any:
    """Rosenblatt transform (``Vinecop::rosenblatt``).

    Parameters
    ----------
    u : array, shape (n, d + k), dtype float
        Prepared pseudo-observations in the compact layout.
    x : array, shape (n, p), or None
        External covariates threaded to each pair copula, or ``None``.
    randomize_discrete : bool, default=True
        Mix each discrete variable's conditional distribution function with its
        left limit using independent uniforms.
    seeds : list of int, or None, optional
        RNG seeds for that randomization.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms in ``[1e-10, 1 - 1e-10]``.
    """
    xp = array_namespace(u)
    d, trunc_lvl = self.d, self.trunc_lvl
    n = u.shape[0]
    order, inv = self.order, self.inverse_order
    # Seed both h-function scratch matrices with the natural-order observations.
    hfunc2 = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      hfunc2[:, j] = u[:, order[j] - 1]
    hfunc1 = xp.asarray(hfunc2, copy=True)
    # See _pdf on why hfunc1_sub needs no seed.
    hfunc2_sub: Any = self._seed_sub(u, xp)
    hfunc1_sub: Any = (
      None
      if hfunc2_sub is None
      else xp.zeros((n, d), dtype=u.dtype, device=u.device)
    )
    u_nat = (
      xp.asarray(hfunc2, copy=True)
      if self._context.assembles_conditioning
      else None
    )
    s = self.structure
    for tree in range(trunc_lvl):
      for edge in range(d - tree - 1):
        edge_copula = self._get_pair_copula(tree, edge)
        col0, col1, subs, types = self._edge_columns(
          tree, edge, hfunc1, hfunc2, hfunc1_sub, hfunc2_sub
        )
        u_e = xp.stack(
          [col0, col1] if subs is None else [col0, col1, *subs], axis=-1
        )
        x_e = self._edge_context(tree, edge, x, u_nat, None)
        # hfunc1 only if needed downstream; hfunc2 is the running transform.
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = _pair_eval(edge_copula.hfunc1, u_e, x_e)
          if subs is not None and types[1] == "d":
            u_h1 = xp.stack([col0, subs[1], *subs], axis=-1)
            hfunc1_sub[:, edge] = _pair_eval(edge_copula.hfunc1, u_h1, x_e)
        hfunc2[:, edge] = _pair_eval(edge_copula.hfunc2, u_e, x_e)
        if subs is not None and types[0] == "d":
          u_h2 = xp.stack([subs[0], col1, *subs], axis=-1)
          hfunc2_sub[:, edge] = _pair_eval(edge_copula.hfunc2, u_h2, x_e)
    # Scatter the transformed columns back to variable order.
    out = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      out[:, j] = hfunc2[:, inv[j]]
    if randomize_discrete and hfunc2_sub is not None:
      # A discrete variable's conditional distribution function jumps, so the
      # transform is uniform only after mixing the jump's two ends with an
      # independent uniform (Brockwell 2007). Continuous variables have
      # coinciding limits, so the mix leaves them untouched.
      left = xp.empty((n, d), dtype=u.dtype, device=u.device)
      for j in range(d):
        source = hfunc2_sub if self._var_types[j] == "d" else hfunc2
        left[:, j] = source[:, inv[j]]
      r: Any = self._sample_uniform(n, False, list(seeds or []))
      out = out * r + left * (1.0 - r)
    return trim(xp, out)

  def _inverse_rosenblatt(self, u: Any, x: Optional[Any]) -> Any:
    """Inverse Rosenblatt transform (``Vinecop::inverse_rosenblatt``).

    Walks variables from ``d - 2`` down to ``0``; at each ``var`` it fills the
    ``hinv2`` column from the outermost tree inward. The ``(trunc_lvl + 1, d, n)``
    scratch is transposed relative to the forward cascades (variable axis first)
    so a finalized ``hinv2[0, var, :]`` row can seed later inversions.

    There is no left-limit cascade here: the transform produces the values a
    left limit would be taken of, so every pair is evaluated as continuous —
    which is also what makes its output a continuous ``(n, d)`` matrix.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Prepared independent uniforms.
    x : array, shape (n, p), or None
        External covariates threaded to each pair copula, or ``None``.

    Returns
    -------
    array, shape (n, d), dtype float
        Dependent uniforms in ``[1e-10, 1 - 1e-10]``.
    """
    xp = array_namespace(u)
    d, trunc_lvl = self.d, self.trunc_lvl
    n = u.shape[0]
    order, inv = self.order, self.inverse_order
    if trunc_lvl == 0:
      out = xp.empty((n, d), dtype=u.dtype, device=u.device)
      for j in range(d):
        out[:, j] = u[:, order[inv[j]] - 1]
      return out
    hinv2 = xp.empty((trunc_lvl + 1, d, n), dtype=u.dtype, device=u.device)
    hfunc1 = xp.empty_like(hinv2)
    for j in range(d):
      hinv2[min(trunc_lvl, d - j - 1), j, :] = u[:, order[j] - 1]
    hfunc1[0, d - 1, :] = hinv2[0, d - 1, :]
    s = self.structure
    for var in range(d - 2, -1, -1):
      tree_start = min(trunc_lvl - 1, d - var - 2)
      for tree in range(tree_start, -1, -1):
        edge_copula = self._get_pair_copula(tree, var)
        if self._n_discrete and "d" in self._pair_types[tree][var]:
          # The inverse cascade *produces* the values a left limit would be
          # taken of, so it evaluates every pair as continuous -- exactly as
          # ``Vinecop::inverse_rosenblatt`` does.
          edge_copula = continuous_view(edge_copula)
        # Same m / on-diagonal rule as the forward cascades (class.ipp:1026),
        # but the inputs are rows of the transposed hinv2 / hfunc1 scratch.
        m = int(s.min_array(tree, var))
        on_diagonal = m == int(s.struct_array(tree, var, natural_order=True))
        u_e_col0 = hinv2[tree + 1, var, :]
        u_e_col1 = (
          hinv2[tree, m - 1, :] if on_diagonal else hfunc1[tree, m - 1, :]
        )
        u_e = xp.stack([u_e_col0, u_e_col1], axis=-1)
        # Conditioning u_D is read from the finalized hinv2[0] rows (the
        # conditioning variables are finalized before this var by the invariant).
        x_e = self._edge_context(tree, var, x, None, hinv2[0])
        hinv2[tree, var, :] = _pair_eval(edge_copula.hinv2, u_e, x_e)
        # Propagate hfunc1 for the next-inner inversion when needed.
        if var < d - 1 and s.needed_hfunc1(tree, var):
          u_e_after = xp.stack([hinv2[tree, var, :], u_e_col1], axis=-1)
          hfunc1[tree + 1, var, :] = _pair_eval(
            edge_copula.hfunc1, u_e_after, x_e
          )
    out = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      out[:, j] = hinv2[0, inv[j], :]
    return trim(xp, out)

  # --- batched cascades (grid fast path; array-agnostic loops) ---------- #
  #
  # Numerically equivalent to the non-batched cascades on a simplified vine,
  # but each tree level fires one stacked pair-copula call over its edges
  # instead of a Python loop. Only the grid state built by ``_build_batched``
  # is subclass-specific; the loops below are array-agnostic (``xp``). The
  # batched-vine surface used here: ``bv.level(t)`` / ``bv.grid_points`` and,
  # per level, ``gather_inputs`` / ``pdf_h1_h2`` / ``h1_h2`` (each returning
  # ``(N_t, n)`` slices over the ``N_t = n_pairs`` edges, fusing the shared
  # bilinear cell search across pdf + both h-functions) plus the ``needs_h1`` /
  # ``needs_h2`` masks. These receive already-prepped ``u`` (the public methods
  # prep before dispatch).
  def _pdf_batched(self, u: Any) -> Any:
    """Batched vine pdf: product over per-tree-level stacked densities.

    Numerically equivalent to ``_pdf`` on a simplified vine, but each tree
    level fires one stacked (fused) pair-copula call over its edges.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Prepared pseudo-observations.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values.
    """
    xp = self._namespace(u)
    d, trunc_lvl = self.d, self.trunc_lvl
    n = u.shape[0]
    if trunc_lvl == 0:
      return xp.ones(n, dtype=u.dtype, device=u.device)
    bv = self._ensure_batched()
    hfunc1 = xp.zeros((n, d), dtype=u.dtype, device=u.device)
    hfunc2 = xp.empty((n, d), dtype=u.dtype, device=u.device)
    order = self.order
    for j in range(d):
      hfunc2[:, j] = u[:, order[j] - 1]
    pdf = xp.ones(n, dtype=u.dtype, device=u.device)
    for t in range(trunc_lvl):
      lvl = bv.level(t)
      u_e = lvl.gather_inputs(hfunc1, hfunc2)  # (N_t, n, 2)
      # One fused lookup yields pdf + both h-functions (shared cell search).
      pdf_e, h1_e, h2_e = lvl.pdf_h1_h2(bv.grid_points, u_e)
      # Product over the level's edges (axis 0), then into the running product
      # (the batched analog of _pdf's per-edge cwiseProduct).
      pdf = pdf * xp.prod(pdf_e, axis=0)
      # Overwrite the next-tree columns flagged by needs_h{1,2} (mirrors _pdf's
      # gated per-edge writes).
      n_pairs = lvl.n_pairs
      h1_new = xp.matrix_transpose(h1_e)  # (n, N_t)
      h2_new = xp.matrix_transpose(h2_e)
      hfunc1[:, :n_pairs] = xp.where(
        lvl.needs_h1[None, :], h1_new, hfunc1[:, :n_pairs]
      )
      hfunc2[:, :n_pairs] = xp.where(
        lvl.needs_h2[None, :], h2_new, hfunc2[:, :n_pairs]
      )
    return pdf

  def _inverse_rosenblatt_batched(self, u: Any) -> Any:
    """Batched inverse Rosenblatt: one stacked call per dependency wave.

    Bit-identical to :meth:`_inverse_rosenblatt` on a simplified vine: the
    waves reorder the cells without changing what any one of them computes.
    The inverse's dependencies do not reduce to tree levels -- a wave holds
    one cell from almost every tree -- so the grouping is by longest-path
    level of the static ``(var, tree)`` graph, which the subclass's batched
    state levels once at bake time. The scratch is flattened to
    ``((trunc_lvl + 1) * d, n)`` so a cell's slot is one row, and a whole wave
    is one gather per input and one scatter per output.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Prepared independent uniforms.

    Returns
    -------
    array, shape (n, d), dtype float
        Dependent uniforms in ``[1e-10, 1 - 1e-10]``.
    """
    xp = self._namespace(u)
    d, trunc_lvl = self.d, self.trunc_lvl
    n = u.shape[0]
    order, inv = self.order, self.inverse_order
    if trunc_lvl == 0:
      out = xp.empty((n, d), dtype=u.dtype, device=u.device)
      for j in range(d):
        out[:, j] = u[:, order[inv[j]] - 1]
      return out
    bv = self._ensure_batched()
    rows = (trunc_lvl + 1) * d
    hinv2 = xp.empty((rows, n), dtype=u.dtype, device=u.device)
    hfunc1 = xp.empty_like(hinv2)
    for j in range(d):
      hinv2[min(trunc_lvl, d - j - 1) * d + j, :] = u[:, order[j] - 1]
    hfunc1[d - 1, :] = hinv2[d - 1, :]
    for k in range(bv.n_waves):
      bv.wave(k).apply_to(bv.grid_points, hinv2, hfunc1)
    out = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      out[:, j] = hinv2[inv[j], :]
    return trim(xp, out)

  def _rosenblatt_batched(self, u: Any) -> Any:
    """Batched Rosenblatt transform (per-tree-level stacked h-functions).

    Numerically equivalent to ``_rosenblatt`` on a simplified vine.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Prepared pseudo-observations.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms in ``[1e-10, 1 - 1e-10]``.
    """
    xp = self._namespace(u)
    d, trunc_lvl = self.d, self.trunc_lvl
    n = u.shape[0]
    order, inv = self.order, self.inverse_order
    bv = self._ensure_batched()
    hfunc2 = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      hfunc2[:, j] = u[:, order[j] - 1]
    hfunc1 = xp.asarray(hfunc2, copy=True)
    for t in range(trunc_lvl):
      lvl = bv.level(t)
      u_e = lvl.gather_inputs(hfunc1, hfunc2)
      n_pairs = lvl.n_pairs
      # One fused lookup yields both h-functions (shared cell search).
      h1_e, h2_e = lvl.h1_h2(bv.grid_points, u_e)
      h1_new = xp.matrix_transpose(h1_e)
      h2_new = xp.matrix_transpose(h2_e)
      # hfunc2 is overwritten unconditionally at every edge; hfunc1 is gated.
      hfunc2[:, :n_pairs] = h2_new
      hfunc1[:, :n_pairs] = xp.where(
        lvl.needs_h1[None, :], h1_new, hfunc1[:, :n_pairs]
      )
    out = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      out[:, j] = hfunc2[:, inv[j]]
    return trim(xp, out)

  # --- batched dispatch ------------------------------------------------- #
  def _namespace(self, a: Any) -> Any:
    """The array namespace of this vine's working arrays, resolved once.

    :meth:`_prep` coerces every input to the vine's own array type, so the
    namespace is a property of the vine rather than of the call. Resolving it
    per call would put a type-dispatch table walk inside each cascade, which a
    tracing compiler then has to trace through -- ``array_namespace`` is
    memoized on the type, but the memo is itself Python that ends up in the
    graph. Every entry point resolves it through :meth:`_prep`, which runs
    before the cascade, so by the time a cascade asks it is already answered.

    Parameters
    ----------
    a : array
        Any array already coerced to the working type.

    Returns
    -------
    module
        The array-API namespace for ``a``.
    """
    xp = self._xp
    if xp is None:
      xp = array_namespace(a)
      object.__setattr__(self, "_xp", xp)
    return xp

  def __getstate__(self) -> dict:
    """The picklable state: everything but the resolved array namespace.

    :meth:`_namespace` memoizes a *module*, which no pickle can carry, so it
    is dropped and re-resolved on first use. Nothing else about the vine
    depends on it having been resolved.

    Returns
    -------
    dict
        The instance state, with the namespace memo cleared.
    """
    # `object.__getstate__` answers `None` for an empty instance; every
    # concrete vine has state, so this only keeps the type honest.
    raw = cast("dict[str, Any]", super().__getstate__() or {})
    state = dict(raw)
    state["_xp"] = None
    return state

  def _resolve_batched(
    self, requested: Optional[bool], x: Optional[Any]
  ) -> bool:
    """Resolve the ``batched`` flag; force ``False`` for conditional or discrete.

    Declining is the right answer rather than raising: ``batched`` defaults to
    the subclass's :meth:`_default_batched` (device-dependent on the torch
    vine), so a raise would make an ordinary ``pdf(u)`` fail on a discrete vine
    for a reason the caller never asked about. Discreteness is a property of the
    vine, not of the subclass's grid, which is why it is decided here rather
    than through :class:`_NotBatchable`.

    Parameters
    ----------
    requested : bool or None
        The caller's ``batched`` argument; ``None`` defers to
        :meth:`_default_batched`.
    x : array or None
        External covariates for the call; any non-``None`` value forces the
        non-batched cascade.

    Returns
    -------
    bool
        Whether to attempt the batched fast path.
    """
    if self._context.assembles_conditioning or x is not None:
      return False
    if self._n_discrete:
      return False
    if requested is None:
      requested = self._default_batched()
    return bool(requested)

  # --- relabeling onto a chosen sampling-order tail --------------------- #
  _REORIENT_NON_SIMPLIFIED = (
    "conditioning_set is not supported on a non-simplified vine: relabeling "
    "can permute the columns of the conditioning matrix x_e each pair copula "
    "receives, so the reoriented vine is not the same model. Condition on the "
    "variables already at the tail of the order."
  )

  def _reoriented(
    self, conditioning_set: Optional[list[int]]
  ) -> "VinecopBase[ArrayT]":
    """The vine to evaluate for ``conditioning_set`` (``self`` if none needed)."""
    if conditioning_set is None:
      return self
    r = reorientation(self.structure, [int(v) for v in conditioning_set])
    if r.identity:
      return self
    # A relabeling can permute the columns of an edge's conditioning matrix
    # x_e -- their order follows the chain of edges through the diagonal
    # variable, which is what the relabeling changes -- and a conditional pair
    # copula consumes x_e positionally. The relabeled vine would be a different
    # model, not the same one in a different order.
    if self._context.assembles_conditioning:
      raise NotImplementedError(self._REORIENT_NON_SIMPLIFIED)
    return _ReorientedVine(self, r)

  def reorient(
    self, conditioning_set: list[int]
  ) -> tuple["RVineStructure", list[list[BicopLike[ArrayT]]]]:
    """Relabel to an equivalent vine whose order tail is ``conditioning_set``.

    Value-preserving: the relabeled vine has the same density and
    log-likelihood, only a different sampling-order representation, so the given
    variables are drawn first and can be conditioned on with
    :meth:`sample_conditional`. Like :meth:`fit` and :meth:`select` -- and
    unlike :meth:`~pyvinecopulib.core.Vinecop.reorient`, which mutates -- it
    **returns** the relabeled model, since ``VinecopBase`` leaves pair storage
    to the subclass.

    Parameters
    ----------
    conditioning_set : list of int
        1-based variable labels to place at the tail of the order.

    Returns
    -------
    structure : RVineStructure
        The relabeled structure; its order ends with ``conditioning_set``.
    pair_copulas : list of list of BicopLike
        The same pair copulas on their new slots, argument-swapped with
        :meth:`~pyvinecopulib.core.BicopLike.flip` where the slot requires it --
        ready to host in a vine without re-fitting.

    Raises
    ------
    NotImplementedError
        If the vine is non-simplified and a relabeling is actually required, or
        if a pair copula does not implement ``flip``.
    RuntimeError
        If ``conditioning_set`` is empty, holds duplicates or out-of-range
        entries, leaves no variable free, or is not admissible as a
        sampling-order tail. Same messages as
        :meth:`~pyvinecopulib.core.Vinecop.reorient`.

    See Also
    --------
    pyvinecopulib.core.Vinecop.reorient : The reference (in-place) relabeling.
    """
    r = reorientation(self.structure, [int(v) for v in conditioning_set])
    if r.identity:
      return self.structure, [
        [self._get_pair_copula(t, e) for e in range(self.d - 1 - t)]
        for t in range(self.trunc_lvl)
      ]
    if self._context.assembles_conditioning:
      raise NotImplementedError(self._REORIENT_NON_SIMPLIFIED)
    pairs: list[list[BicopLike[ArrayT]]] = []
    for tree in range(self.trunc_lvl):
      row: list[BicopLike[ArrayT]] = []
      for edge in range(self.d - 1 - tree):
        old_edge, flipped = r.locations[(tree, edge)]
        pair = self._get_pair_copula(tree, old_edge)
        row.append(pair.flip() if flipped else pair)
      pairs.append(row)
    return r.structure, pairs

  # --- public evaluator surface (VinecopLike) --------------------------- #
  def pdf(
    self,
    u: ArrayT,
    *,
    num_threads: int = 1,
    x: Optional[ArrayT] = None,
    batched: Optional[bool] = None,
  ) -> ArrayT:
    """Evaluate the vine-copula density ``c(u_1, ..., u_d)``.

    Parameters
    ----------
    u : array, shape (n, d), (n, d + k) or (n, 2d), dtype float
        Pseudo-observations in ``[0, 1]`` (clamped to ``[1e-10, 1 - 1e-10]``).
        With ``k`` discrete variables, the left limits ``F(x^-)`` are required
        too: pass the expanded ``(n, 2d)`` layout, or the compact ``(n, d + k)``
        one that omits the left-limit columns of the continuous variables.
    num_threads : int, default=1
        Accepted for parity with ``Vinecop.pdf()``; ignored.
    x : array, shape (n, p), or None, optional
        External covariates threaded to each pair copula (non-simplified
        vines). ``None`` for the unconditional / simplified case.
    batched : bool or None, optional
        Fire one batched pair-copula call per tree level. ``None`` resolves
        via the subclass default; forced ``False`` when conditioning is active
        or any variable is discrete, and falls back to the non-batched cascade
        if the subclass has no batched fast path (or a pair does not support it).

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values.
    """
    del num_threads
    u_p = self._prep(u, "pdf")
    if self._resolve_batched(batched, x):
      try:
        return cast(ArrayT, self._pdf_batched(u_p))
      except _NotBatchable:
        pass  # no grid fast path available -> non-batched cascade
    return cast(ArrayT, self._pdf(u_p, x))

  def rosenblatt(
    self,
    u: ArrayT,
    *,
    num_threads: int = 1,
    randomize_discrete: bool = True,
    seeds: Optional[list[int]] = None,
    x: Optional[ArrayT] = None,
    batched: Optional[bool] = None,
    conditioning_set: Optional[list[int]] = None,
  ) -> ArrayT:
    """Rosenblatt transform: dependent uniforms to independent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), (n, d + k) or (n, 2d), dtype float
        Pseudo-observations in ``[0, 1]``; see :meth:`pdf` on the layouts.
    num_threads : int, default=1
        Accepted for parity; ignored.
    randomize_discrete : bool, default=True
        For a discrete variable the conditional distribution function jumps, so
        the transform is uniform only after mixing the jump's two ends with an
        independent uniform. Continuous variables are unaffected either way.
    seeds : list of int or None, optional
        RNG seeds for that randomization; forwarded to the subclass's
        base-uniform draw.
    x : array, shape (n, p), or None, optional
        External covariates (non-simplified vines), else ``None``.
    batched : bool or None, optional
        See :meth:`pdf`.
    conditioning_set : list of int or None, optional
        Condition on these 1-based variables instead of the ones at the tail of
        the vine order. It does not subset ``u``, which stays the full matrix;
        what changes is which conditional distributions the output columns
        represent. Passing the order tail itself is the identity.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms in ``[1e-10, 1 - 1e-10]``.

    Raises
    ------
    NotImplementedError
        If ``conditioning_set`` needs a relabeling on a non-simplified vine, or
        if a pair copula does not implement ``flip``.
    RuntimeError
        If ``conditioning_set`` is inadmissible as a sampling-order tail; see
        :meth:`reorient`.
    """
    del num_threads
    view = self._reoriented(conditioning_set)
    if view is not self:
      return view.rosenblatt(
        u,
        randomize_discrete=randomize_discrete,
        seeds=seeds,
        x=x,
        batched=batched,
      )
    u_p = self._prep(u, "rosenblatt")
    if self._resolve_batched(batched, x):
      try:
        return cast(ArrayT, self._rosenblatt_batched(u_p))
      except _NotBatchable:
        pass  # no grid fast path available -> non-batched cascade
    return cast(ArrayT, self._rosenblatt(u_p, x, randomize_discrete, seeds))

  def inverse_rosenblatt(
    self,
    u: ArrayT,
    *,
    num_threads: int = 1,
    x: Optional[ArrayT] = None,
    batched: Optional[bool] = None,
    conditioning_set: Optional[list[int]] = None,
  ) -> ArrayT:
    """Inverse Rosenblatt transform: independent uniforms to dependent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Independent uniforms in ``[0, 1]^d``. A wider layout is accepted and its
        left-limit columns ignored, since the transform needs only the values.
    num_threads : int, default=1
        Accepted for parity; ignored.
    x : array, shape (n, p), or None, optional
        External covariates (non-simplified vines), else ``None``.
    batched : bool or None, optional
        Whether to evaluate whole groups of pair copulas per call rather than
        one at a time; ``None`` takes the subclass default. The inverse
        cascade's dependencies run across tree levels, so the groups here are
        the levels of the dependency graph rather than the trees themselves.
        Ignored when no batched path is available.
    conditioning_set : list of int or None, optional
        See :meth:`rosenblatt`.

    Returns
    -------
    array, shape (n, d), dtype float
        Dependent uniforms in ``[1e-10, 1 - 1e-10]``.

    Raises
    ------
    NotImplementedError
        If ``conditioning_set`` needs a relabeling on a non-simplified vine,
        or if a pair copula does not implement ``flip``.
    RuntimeError
        If ``conditioning_set`` is inadmissible as a sampling-order tail; see
        :meth:`reorient`.
    """
    del num_threads
    view = self._reoriented(conditioning_set)
    if view is not self:
      return view.inverse_rosenblatt(u, x=x, batched=batched)
    u_p = self._prep(u, "inverse_rosenblatt", values_only=True)
    with self._eval_context():
      if self._resolve_batched(batched, x):
        try:
          return cast(ArrayT, self._inverse_rosenblatt_batched(u_p))
        except _NotBatchable:
          pass
      return cast(ArrayT, self._inverse_rosenblatt(u_p, x))

  def sample(
    self,
    n: int,
    *,
    qrng: bool = False,
    num_threads: int = 1,
    seeds: Optional[list[int]] = None,
    x: Optional[ArrayT] = None,
    batched: Optional[bool] = None,
  ) -> ArrayT:
    """Simulate ``n`` samples from the fitted copula.

    Parameters
    ----------
    n : int
        Number of samples (must equal ``rows(x)`` when ``x`` is given).
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random.
    num_threads : int, default=1
        Accepted for parity; ignored.
    seeds : list of int or None, optional
        RNG seeds forwarded to the subclass's base-uniform draw.
    x : array, shape (n, p), or None, optional
        External covariates for a conditional draw (one row per sample).
    batched : bool or None, optional
        Forwarded to :meth:`inverse_rosenblatt`.

    Returns
    -------
    array, shape (n, d), dtype float
        Dependent uniforms in ``[1e-10, 1 - 1e-10]``.
    """
    del num_threads
    seeds = list(seeds) if seeds else []
    if x is not None:
      xa: Any = x
      if xa.shape[0] != n:
        raise ValueError(
          f"sample(n={n}, x=...) requires one covariate row per sample; "
          f"got x with {xa.shape[0]} rows."
        )
    with self._eval_context():
      base_u = self._sample_uniform(n, qrng, seeds)
      return self.inverse_rosenblatt(base_u, x=x, batched=batched)

  def sample_conditional(
    self,
    u_cond: ArrayT,
    *,
    qrng: bool = False,
    num_threads: int = 1,
    seeds: Optional[list[int]] = None,
    conditioning_set: Optional[list[int]] = None,
    x: Optional[ArrayT] = None,
  ) -> ArrayT:
    """Sample from the conditional copula given fixed values of some variables.

    Each row of ``u_cond`` is one conditioning point, and the matching output row
    draws the remaining variables from their distribution conditional on it. To
    draw many samples at one point, pass that point repeated over ``n`` rows.

    Parameters
    ----------
    u_cond : array, shape (n, k) or wider, dtype float
        Conditioning values in ``(0, 1)``, one point per row. Column ``i``
        corresponds to the ``i``-th conditioning variable. A discrete
        conditioning variable needs its left limit ``F(x^-)`` too: pass the
        expanded ``(n, 2k)`` layout, or the compact ``(n, k + k_d)`` one that
        omits the left limits of the continuous conditioners.
    qrng : bool, default=False
        Draw quasi-random base uniforms for the conditioned variables.
    num_threads : int, default=1
        Accepted for parity; ignored.
    seeds : list of int or None, optional
        RNG seeds forwarded to the subclass's base-uniform draw.
    conditioning_set : list of int or None, optional
        The 1-based variables to condition on. ``None`` takes the last
        ``k`` of :attr:`order`, with ``k`` inferred from ``u_cond``'s width.
        **The two forms map the columns differently**: without it, column ``i``
        is the ``i``-th variable of the order tail; with it, column ``i`` is
        ``conditioning_set[i]``.
    x : array, shape (n, p), or None, optional
        External covariates (non-simplified vines), one row per sample.

    Returns
    -------
    array, shape (n, d), dtype float
        Conditional draws. The conditioning variables' columns reproduce
        ``u_cond``; a discrete one is reproduced up to its atom, landing in
        ``[F(x^-), F(x)]``.

    Raises
    ------
    ValueError
        If ``u_cond`` is not 2-d, if its column count matches no admissible
        layout, if a discrete conditioner's left limit exceeds its value, or if
        ``x`` does not have one row per sample.
    NotImplementedError
        If ``conditioning_set`` needs a relabeling on a non-simplified vine.
    RuntimeError
        If ``conditioning_set`` is inadmissible as a sampling-order tail; see
        :meth:`reorient`.

    Notes
    -----
    The free variables are completed with an arbitrary ``0.5`` before the
    forward transform, which cannot affect the result: in natural order, column
    ``j``'s Rosenblatt coordinate reads only columns ``j, ..., d - 1``, so the
    order tail is a self-contained sub-vine. That is also why the conditioning
    set has to *be* the tail -- and why the placeholders are harmless on a
    non-simplified vine too, a tail edge's conditioning columns all lying in the
    tail.

    See Also
    --------
    pyvinecopulib.core.Vinecop.sample_conditional : The reference sampler.
    """
    del num_threads
    seeds = list(seeds) if seeds else []
    ua: Any = u_cond
    if ua.ndim != 2:
      raise ValueError(
        "sample_conditional: u_cond must have shape (n, k); got "
        f"{tuple(ua.shape)}"
      )
    d = self.d
    n, n_cols = int(ua.shape[0]), int(ua.shape[1])
    if x is not None and cast("Any", x).shape[0] != n:
      raise ValueError(
        "sample_conditional: x requires one covariate row per sample; got x "
        f"with {cast('Any', x).shape[0]} rows and u_cond with {n}."
      )
    view = self._reoriented(conditioning_set)
    if conditioning_set is None:
      cond_vars = self._infer_conditioning_set(n_cols)
    else:
      cond_vars = [int(v) for v in conditioning_set]
    k = len(cond_vars)
    n_disc_cond = sum(1 for v in cond_vars if self._var_types[v - 1] == "d")
    # When every conditioner is discrete the two layouts coincide; the expanded
    # one then wins, as it does upstream.
    expanded = n_cols == 2 * k
    if not expanded and n_cols != k + n_disc_cond:
      raise ValueError(
        f"u_cond has wrong number of columns; expected: {2 * k} (n x 2k "
        f"expanded layout) or {k + n_disc_cond} (n x (k + k_d) compact "
        f"layout), actual: {n_cols}."
      )

    with self._eval_context():
      xp = array_namespace(ua)
      u_completed = xp.full(
        (n, d + self._n_discrete), 0.5, dtype=ua.dtype, device=ua.device
      )
      seen = 0
      for i, var in enumerate(cond_vars):
        col = var - 1
        u_completed[:, col] = ua[:, i]
        if self._var_types[col] == "d":
          left = ua[:, k + i if expanded else k + seen]
          if bool(xp.any(left > ua[:, i])):
            raise ValueError(
              "for discrete conditioning variables, the left-limit columns of "
              "u_cond (F(x^-)) must not exceed the value columns (F(x))."
            )
          u_completed[:, d + self._disc_cols[col]] = left
          seen += 1
      # The conditioning variables' own Rosenblatt coordinates; randomized so a
      # discrete conditioner's jump becomes a uniform the inverse can invert.
      w: Any = view.rosenblatt(
        cast(ArrayT, u_completed),
        randomize_discrete=True,
        seeds=seeds,
        x=x,
      )
      base_u: Any = self._sample_uniform(n, qrng, seeds)
      for var in cond_vars:
        base_u[:, var - 1] = w[:, var - 1]
      return view.inverse_rosenblatt(cast(ArrayT, base_u), x=x)

  def _infer_conditioning_set(self, n_cols: int) -> list[int]:
    """The order tail whose layout is ``n_cols`` columns wide."""
    return infer_conditioning_set(
      list(self.order), list(self._var_types), n_cols
    )

  def cdf(
    self,
    u: ArrayT,
    *,
    N: int = 10000,
    qrng: bool = True,
    num_threads: int = 1,
    seeds: Optional[list[int]] = None,
    x: Optional[ArrayT] = None,
    block_size: int = 4096,
    batched: Optional[bool] = None,
  ) -> ArrayT:
    """Evaluate the joint CDF at each query row via quasi-Monte-Carlo.

    Parameters
    ----------
    u : array, shape (m, d), (m, d + k) or (m, 2d), dtype float
        Query points in ``[0, 1]``; see :meth:`pdf` on the layouts. Only the
        ``d`` value columns enter the estimate — ``C(u)`` is a right limit — but
        a discrete vine still requires the wider layout, as ``Vinecop.cdf()``
        does.
    N : int, default=10000
        Number of Monte-Carlo samples.
    qrng : bool, default=True
        Draw quasi-random samples (matches ``Vinecop.cdf()``).
    num_threads : int, default=1
        Accepted for parity; ignored.
    seeds : list of int or None, optional
        RNG seeds forwarded to :meth:`sample`.
    x : array or None, optional
        Must be ``None`` — a per-row conditional CDF is not supported (the
        Monte-Carlo dominance estimate cannot condition each query row on a
        different covariate without per-``x`` resampling).
    block_size : int, default=4096
        Query rows processed per iteration (peak-memory control).
    batched : bool or None, optional
        Forwarded to :meth:`sample`, which draws the Monte-Carlo sample.

    Returns
    -------
    array, shape (m,), dtype float
        CDF values in ``[0, 1]``.

    Raises
    ------
    NotImplementedError
        If ``x`` is not ``None``.
    """
    del num_threads
    if x is not None:
      raise NotImplementedError(
        "Conditional cdf (x is not None) is not supported: the Monte-Carlo "
        "dominance estimate would need O(unique_x * N) resampling."
      )
    seeds = list(seeds) if seeds else []
    with self._eval_context():
      prepped: Any = self._prep(u, "cdf")
      # Only the value block enters the dominance count: C(u) is a right limit.
      u_t: Any = prepped[:, : self.d]
      samples: Any = self.sample(N, qrng=qrng, seeds=seeds, batched=batched)
      xp = array_namespace(u_t)
      m = u_t.shape[0]
      out = xp.empty(m, dtype=u_t.dtype, device=u_t.device)
      for start in range(0, m, block_size):
        end = min(start + block_size, m)
        dominated = xp.all(
          samples[None, :, :] <= u_t[start:end][:, None, :], axis=-1
        )
        out[start:end] = xp.mean(xp.astype(dominated, u_t.dtype), axis=1)
      return cast(ArrayT, out)

  # --- convenience surface (loglik / plot / accessors) ------------------ #
  def loglik(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Total log-likelihood ``sum(log c(u))`` of the vine at ``u``.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d``.
    x : array, shape (n, p), or None, optional
        External covariates for a non-simplified vine, else ``None``.

    Returns
    -------
    array, shape (), dtype float
        The summed log-density (a differentiable scalar under autograd, e.g.
        PyTorch); the per-observation vine density is floored at ``1e-20`` before
        the log.
    """
    dens: Any = self.pdf(u, x=x)
    xp = array_namespace(dens)
    return cast(ArrayT, xp.sum(xp.log(xp.clip(dens, 1e-20))))

  @property
  def dim(self) -> int:
    """Number of variables in the vine.

    Returns
    -------
    int
        The vine dimension ``d``.
    """
    return self.d

  @property
  def matrix(self) -> Any:
    """R-vine structure matrix (from :attr:`structure`).

    Returns
    -------
    ndarray, shape (d, d), dtype int
        The structure matrix; used by :meth:`plot`.
    """
    return self.structure.matrix

  def get_pair_copula(self, tree: int, edge: int) -> BicopLike[ArrayT]:
    """Return the pair copula at ``(tree, edge)``.

    Parameters
    ----------
    tree : int
        Tree index (``0``-based).
    edge : int
        Edge index within the tree (``0``-based).

    Returns
    -------
    BicopLike
        The pair copula hosted at that position.
    """
    return self._get_pair_copula(tree, edge)

  def plot(
    self,
    tree: Optional[list[int]] = None,
    add_edge_labels: bool = True,
    layout: str = "graphviz",
    vars_names: Optional[list[str]] = None,
  ) -> None:
    """Plot the vine tree structure with networkx.

    Draws one panel per requested tree from :attr:`structure`, mirroring
    :meth:`pyvinecopulib.core.Vinecop.plot`.

    Parameters
    ----------
    tree : list of int, or None, optional
        Tree indices to plot; all trees when ``None``.
    add_edge_labels : bool, default=True
        Annotate edges with their conditioned / conditioning sets.
    layout : str, default="graphviz"
        ``"graphviz"`` (needs pydot + graphviz) or ``"spring_layout"``.
    vars_names : list of str, or None, optional
        Variable names; the integer indices are used when ``None``.

    Returns
    -------
    None
        The figure is drawn with matplotlib.
    """
    from .._python_helpers.vinecop import vinecop_plot

    vinecop_plot(self, tree, add_edge_labels, layout, vars_names)

  def __repr__(self) -> str:
    return (
      f"{type(self).__name__}(dim={self.d}, trunc_lvl={self.trunc_lvl}, "
      f"order={list(self.order)})"
    )

  # --- shared sequential-fit engine ------------------------------------- #
  @staticmethod
  def fit(
    structure: Any,
    u: Any,
    fit_edge: FitEdge,
    *,
    context: Optional[ConditioningContext] = None,
    x: Optional[Any] = None,
    var_types: Optional[list[str]] = None,
  ) -> list[list[BicopLike]]:
    """Fit pair copulas tree-by-tree along a fixed structure (returns them).

    The array-agnostic (NumPy or PyTorch) analog of
    :meth:`~pyvinecopulib.core.Vinecop.fit`, with the pair-copula fit supplied
    by the ``fit_edge`` callback. It differs in one way: rather than mutating a
    vine in place (as :meth:`~pyvinecopulib.core.Vinecop.fit` does), it
    **returns** the fitted pairs as a nested ``[tree][edge]`` list — since
    ``VinecopBase`` leaves pair storage to the subclass, there is no single
    object to mutate. It mirrors the forward pdf traversal with the density
    evaluation replaced by ``fit_edge(tree, edge, u_e, x_e)``; the returned
    pair's ``hfunc1`` / ``hfunc2`` must be valid immediately for tree
    propagation. Conditional fitting is driven through this seam (a
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
    context : ConditioningContext, optional
        Conditioning-context policy (default: simplified / unconditional).
    x : array, shape (n, p), or None, optional
        External covariates for conditional fitting, else ``None``.
    var_types : list of str, optional
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete), in
        variable order; ``None`` means all continuous.

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
    pair_types = pair_var_types(structure, types) if n_discrete(types) else None
    ua = collapse_data(ua, d, types, "fit")
    n = ua.shape[0]
    hfunc1 = xp.zeros((n, d), dtype=ua.dtype, device=ua.device)
    hfunc2 = xp.empty((n, d), dtype=ua.dtype, device=ua.device)
    for j in range(d):
      hfunc2[:, j] = ua[:, order[j] - 1]
    # Parallel left-limit scratch, zero-initialized like the pdf traversal this
    # mirrors; a column is only read where the edge's types say it was written.
    hfunc2_sub: Any = seed_left_limits(
      ua, d, order, types, disc_cols(types), xp
    )
    hfunc1_sub: Any = (
      None
      if hfunc2_sub is None
      else xp.zeros((n, d), dtype=ua.dtype, device=ua.device)
    )
    u_nat = (
      xp.asarray(hfunc2, copy=True) if context.assembles_conditioning else None
    )
    cache: dict[tuple[int, int], tuple[int, ...]] = {}

    def edge_context_for(tree: int, edge: int) -> Optional[Any]:
      # Assemble x_e = context(u_D, x); u_D columns in ascending conditioning
      # -tree order (C1), matching VinecopBase._cond_positions / _edge_context.
      if not context.assembles_conditioning and x is None:
        return None
      u_D: Optional[Any] = None
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
    pairs: list[list[BicopLike]] = []
    for tree in range(trunc_lvl):
      row: list[BicopLike] = []
      for edge in range(d - tree - 1):
        col0, col1, subs, edge_types = edge_columns(
          s, pair_types, tree, edge, hfunc1, hfunc2, hfunc1_sub, hfunc2_sub
        )
        u_e = stack_edge(xp, col0, col1, subs)
        x_e = edge_context_for(tree, edge)
        edge_copula = _fit_edge_call(fit_edge, tree, edge, u_e, x_e, edge_types)
        row.append(edge_copula)
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = _pair_eval(edge_copula.hfunc1, u_e, x_e)
          if subs is not None and edge_types[1] == "d":
            hfunc1_sub[:, edge] = _pair_eval(
              edge_copula.hfunc1, with_left_limit(u_e, 1), x_e
            )
        if s.needed_hfunc2(tree, edge):
          hfunc2[:, edge] = _pair_eval(edge_copula.hfunc2, u_e, x_e)
          if subs is not None and edge_types[0] == "d":
            hfunc2_sub[:, edge] = _pair_eval(
              edge_copula.hfunc2, with_left_limit(u_e, 0), x_e
            )
      pairs.append(row)
    return pairs

  @staticmethod
  def select(
    u: Any,
    fit_edge: FitEdge,
    *,
    trunc_lvl: Optional[int] = None,
    tree_criterion: str = "tau",
    threshold: float = 0.0,
    tree_algorithm: str = "mst_prim",
    seeds: Optional[list[int]] = None,
    to_numpy: Optional[Callable[[Any], Any]] = None,
    var_types: Optional[list[str]] = None,
    conditioning_set: Optional[list[int]] = None,
  ) -> tuple[RVineStructure, list[list[BicopLike[ArrayT]]]]:
    """Select an R-vine structure from data (array-agnostic Dissmann).

    The array-agnostic (NumPy or PyTorch) analog of
    :meth:`~pyvinecopulib.core.Vinecop.select`, with the pair-copula fit
    supplied by the ``fit_edge`` callback. It differs in one way: rather than
    mutating a vine in place (as :meth:`~pyvinecopulib.core.Vinecop.select`
    does), it **returns** the selected structure and pairs — ``VinecopBase``
    leaves pair storage to the subclass, so there is no single object to
    mutate. It runs the tree-by-tree Dissmann greedy search
    [1]_: for each tree it builds a candidate graph honoring the proximity
    condition, weights every candidate edge by ``1 - |tau|`` (``tau`` is the
    dependence measure named by ``tree_criterion``, Kendall's tau by default,
    via ``wdm``), and keeps a spanning tree (``tree_algorithm``) —
    maximum-dependence for the MST variants, Wilson-weighted random for the
    random ones. Each surviving edge's pair copula is fit by the ``fit_edge``
    callback, whose h-functions feed the next tree.

    The fitted pairs are returned reused, never re-fit: each is placed on its
    slot in the finalized structure and reoriented with its
    :meth:`~pyvinecopulib.core.BicopLike.flip` where the slot's orientation
    requires it. Selection is for a simplified vine — edge weights use the
    unconditional pseudo-observations, so ``fit_edge`` receives ``x_e = None``.

    Parameters
    ----------
    u : array, shape (n, d), (n, d + k) or (n, 2d), dtype float
        Pseudo-observations on any array-API namespace (NumPy or PyTorch). With
        ``k`` discrete variables their left limits ``F(x^-)`` are required too;
        see :meth:`pdf` on the layouts.
    fit_edge : callable
        ``(tree, edge, u_e, x_e) -> BicopLike`` fitting one edge's pair copula;
        its ``hfunc1`` / ``hfunc2`` must be valid immediately, and it must
        implement :meth:`~pyvinecopulib.core.BicopLike.flip` (used to
        reorient reused pairs onto
        their finalized slots). ``x_e`` is always ``None`` here. An edge with a
        discrete argument gets a four-column ``u_e`` and the additional keyword
        ``var_types=[t1, t2]``; the pair it returns must read that layout, so
        wrap a continuous one in
        :class:`~pyvinecopulib.core.DiscretePair`.
    trunc_lvl : int, optional
        Maximum number of trees to select (default: ``d - 1``, i.e. untruncated).
    tree_criterion : str, default "tau"
        Dependence measure passed to ``wdm`` for edge weighting: ``"tau"``,
        ``"rho"``, ``"hoeffd"``, ``"mcor"``, ``"cxi"`` or ``"joe"``, matching
        ``FitControlsVinecop``. ``"cxi"`` is Chatterjee's xi, which is
        asymmetric, so the weight is the larger of the two directions.
    threshold : float, default 0.0
        Dependence threshold: edges with criterion below it are deprioritized
        (weight ``1.0``) during spanning-tree selection.
    tree_algorithm : str, default "mst_prim"
        ``"mst_prim"`` / ``"mst_kruskal"`` (Dissmann) or ``"random_weighted"`` /
        ``"random_unweighted"`` (Wilson).
    seeds : list of int, optional
        RNG seeds for the random tree algorithms (ignored by the MST ones).
    to_numpy : callable, optional
        Maps a 1-d array to a NumPy array for the ``wdm`` call. Defaults to
        :func:`numpy.asarray`; PyTorch callers pass one that detaches and moves
        to host (e.g. ``lambda t: t.detach().cpu().numpy()``).
    var_types : list of str, optional
        Per-variable types, ``"c"`` (continuous) or ``"d"`` (discrete), in
        variable order; ``None`` means all continuous. Given, it also fixes the
        dimension, so ``u`` may carry the extra left-limit columns.
    conditioning_set : list of int or None, optional
        1-based variables to place at the tail of the selected order, so they can
        be conditioned on with :meth:`sample_conditional`. Every candidate edge
        touching a non-conditioning variable is penalized, which makes the
        conditioning set a self-contained block, and the finalized structure is
        then relabeled onto that tail. Requires an MST ``tree_algorithm``, and
        the pairs must implement
        :meth:`~pyvinecopulib.core.BicopLike.flip`.

    Returns
    -------
    structure : RVineStructure
        The selected vine structure.
    pair_copulas : list of list of BicopLike
        The fitted pair copulas, indexed ``[tree][edge]`` in the structure's
        column order and reoriented onto their slots — ready to host in a vine
        without re-fitting.

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
    import numpy as np

    from ..pyvinecopulib_ext import (
      RVineStructure,
      _select_spanning_tree,
      wdm,
    )

    xp = array_namespace(u)
    n = int(u.shape[0])
    # With left-limit columns present, `u` is wider than the vine: `var_types`
    # is what fixes the dimension.
    d = len(var_types) if var_types is not None else int(u.shape[1])
    types = check_var_types(var_types, d)
    offsets = disc_cols(types)
    u = collapse_data(u, d, types, "select")
    seed_list = [int(s) for s in (seeds or [])]
    # ``np.asarray`` raises on a GPU tensor, so the default routes any
    # non-NumPy array through the array API's own host transfer. A caller
    # that knows its backend can still pass a cheaper ``to_numpy``.
    convert = _to_numpy_default if to_numpy is None else to_numpy
    max_trees = (
      d - 1 if trunc_lvl is None else max(0, min(int(trunc_lvl), d - 1))
    )
    cond = [int(v) for v in (conditioning_set or [])]
    in_cond = [False] * d
    if cond:
      # Mirrors `Vinecop::check_conditioning_set`: an MST is what makes the
      # penalty below lay the conditioning block down first.
      if len(cond) >= d:
        raise ValueError(
          "conditioning_set must contain at most d - 1 variables."
        )
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

    def criterion(col0: Any, col1: Any) -> float:
      # |wdm| with the selector's n > 10 guard and NaN -> 0. Only the value
      # columns enter, so a discrete vine selects the tree it would select
      # continuous; the weighted / NaN-compaction corrections are skipped, this
      # path being unweighted.
      if n <= 10:
        return 0.0
      a, b = convert(col0), convert(col1)
      if tree_criterion == "cxi":
        # Chatterjee's xi measures how far one argument is a function of the
        # other, so it is not symmetric -- and a spanning-tree edge weight has
        # to be. `pairwise_cxi` takes the larger of the two directions; taking
        # one of them would pick a different tree.
        value = max(
          float(wdm(a, b, tree_criterion)), float(wdm(b, a, tree_criterion))
        )
      else:
        value = float(abs(wdm(a, b, tree_criterion)))
      return 0.0 if not np.isfinite(value) else value

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
      }
      for i in range(d)
    ]

    trees: list[list[tuple[int, int, list[int]]]] = []
    # Per tree: {(conditioned pair, conditioning set) -> (arg1 label, fitted
    # pair)}, used to place + reorient the pairs onto the finalized slots.
    records: list[dict[Any, tuple[int, Any]]] = []
    for _ in range(max_trees):
      m = len(nodes)
      cand: list[tuple[int, int]] = []
      cand_cols: list[tuple[Any, Any]] = []
      cand_subs: list[Optional[tuple[Any, Any]]] = []
      cand_types: list[tuple[str, str]] = []
      weights: list[float] = []
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
            all_cond = all(
              in_cond[i] for i in nodes[v0]["all_indices"]
            ) and all(in_cond[i] for i in nodes[v1]["all_indices"])
            if not all_cond:
              weight += float(d)
          cand.append((v0, v1))
          cand_cols.append((col0, col1))
          cand_subs.append(subs)
          cand_types.append(edge_types)
          weights.append(weight)

      # Ascending candidate index = boost's edge-list (insertion) order, which
      # is the order the C++ selector iterates surviving edges in.
      selected = sorted(
        _select_spanning_tree(m, cand, weights, tree_algorithm, seed_list)
      )

      tree_edges: list[tuple[int, int, list[int]]] = []
      new_nodes: list[dict[str, Any]] = []
      tree_records: dict[Any, tuple[int, Any]] = {}
      for edge_idx, e in enumerate(selected):
        v0, v1 = cand[e]
        col0, col1 = cand_cols[e]
        subs, edge_types = cand_subs[e], cand_types[e]
        u_e = stack_edge(xp, col0, col1, subs)
        pair = _fit_edge_call(
          fit_edge, len(trees), edge_idx, u_e, None, edge_types
        )
        indices0 = set(nodes[v0]["all_indices"])
        indices1 = set(nodes[v1]["all_indices"])
        # Conditioned pair in the C++ set_sym_diff order: v0's unique variable
        # first (the fitted pair's first argument), then v1's.
        a_var = next(iter(indices0 - indices1))
        b_var = next(iter(indices1 - indices0))
        conditioning = sorted(indices0 & indices1)
        tree_edges.append((a_var + 1, b_var + 1, [c + 1 for c in conditioning]))
        key = (
          frozenset((a_var + 1, b_var + 1)),
          frozenset(c + 1 for c in conditioning),
        )
        tree_records[key] = (a_var + 1, pair)
        new_nodes.append(
          {
            "all_indices": tuple(sorted((a_var, b_var, *conditioning))),
            "h1": pair.hfunc1(u_e),
            "h2": pair.hfunc2(u_e),
            # A discrete argument's next tree needs the h-function at the atom's
            # lower end too, exactly as the cascades compute it.
            "h1_sub": (
              pair.hfunc1(with_left_limit(u_e, 1))
              if edge_types[1] == "d"
              else None
            ),
            "h2_sub": (
              pair.hfunc2(with_left_limit(u_e, 0))
              if edge_types[0] == "d"
              else None
            ),
            "types": edge_types,
            "prev": (v0, v1),
          }
        )
      trees.append(tree_edges)
      records.append(tree_records)
      nodes = new_nodes
      if len(nodes) <= 1:
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
    order = [int(v) for v in structure.order]
    pairs: list[list[BicopLike[ArrayT]]] = []
    for t in range(int(structure.trunc_lvl)):
      row: list[BicopLike[ArrayT]] = []
      for e in range(d - 1 - t):
        diag = order[e]
        partner = int(structure.struct_array(t, e, natural_order=False))
        conditioning_key = frozenset(
          int(structure.struct_array(i, e, natural_order=False))
          for i in range(t)
        )
        a_label, pair = records[t][
          (frozenset((diag, partner)), conditioning_key)
        ]
        row.append(pair if a_label == diag else pair.flip())
      pairs.append(row)
    return structure, pairs


class _ReorientedVine(VinecopBase[ArrayT]):
  """A vine evaluated in a relabeled sampling order, without copying it.

  The counterpart of the compiled ``VinecopView``: it binds the relabeled
  structure and resolves every pair copula back to the viewed vine's slot,
  swapping the pair's arguments where the relabeling requires it. Every cascade
  comes from :class:`VinecopBase`; only the pair lookup is redirected.

  Parameters
  ----------
  base : VinecopBase
      The vine being viewed; its pair copulas, RNG, dtype / device coercion and
      grad control are all reused.
  relabeling : Reorientation
      The relabeled structure and the slot map onto ``base``'s slots.
  """

  def __init__(
    self, base: VinecopBase[ArrayT], relabeling: Reorientation
  ) -> None:
    self._base = base
    self._locations = relabeling.locations
    self._flipped: dict[tuple[int, int], BicopLike[ArrayT]] = {}
    # A declared capability is not inherited by a view of a vine that has it:
    # `supports_covariates` is a class attribute, so without this the view
    # reports `False` and `Vinedist` would drop every covariate it forwards.
    self.supports_covariates = base.supports_covariates
    self._bind_vine(
      relabeling.structure, base._context, var_types=base.var_types
    )

  def _get_pair_copula(self, tree: int, edge: int) -> BicopLike[ArrayT]:
    old_edge, flipped = self._locations[(tree, edge)]
    if not flipped:
      return self._base._get_pair_copula(tree, old_edge)
    key = (tree, old_edge)
    pair = self._flipped.get(key)
    if pair is None:
      # `flip` can be costly -- a grid pair rebuilds itself and re-bakes its
      # integral caches -- and each slot is read once per cascade pass, of which
      # conditional sampling makes two. Callers relabeling repeatedly should use
      # `reorient()` once and host the pairs it returns.
      pair = self._base._get_pair_copula(tree, old_edge).flip()
      self._flipped[key] = pair
    return pair

  # dtype / device coercion, RNG placement and grad control belong to the vine
  # being viewed. `_default_batched` / `_build_batched` are deliberately *not*
  # delegated: the base's batched state is baked against the base's structure and
  # edge order, so the view stays on the non-batched cascade.
  def _prep(self, u: ArrayT, name: str, *, values_only: bool = False) -> ArrayT:
    return self._base._prep(u, name, values_only=values_only)

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> ArrayT:
    return self._base._sample_uniform(n, qrng, seeds)

  def _eval_context(self):
    return self._base._eval_context()


# Shares the worked example with :class:`~pyvinecopulib.core.VinecopLike` (see
# protocols.py) so the contract and its canonical base never drift apart.
VinecopBase.__doc__ = (VinecopBase.__doc__ or "") + _VINECOP_EXAMPLE
