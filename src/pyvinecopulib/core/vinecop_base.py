"""Canonical partial implementation of the vine-copula evaluator contract.

:class:`VinecopBase` is the array-agnostic (NumPy / PyTorch, via
:func:`array_api_compat.array_namespace`) implementation of the vine cascades —
``pdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` / ``simulate`` / ``cdf`` — plus
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

The only hook a concrete subclass must provide is ``_get_pair_copula``; ``_prep``
(input coercion + unit-box clamp) ships a concrete default, and
``_simulate_uniform`` (RNG) is needed only to enable ``simulate``. To enable the
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

from .context import ConditioningContext, SimplifiedContext
from .protocols import ArrayT, BicopLike, VinecopLike, _VINECOP_EXAMPLE

if TYPE_CHECKING:
  from ..pyvinecopulib_ext import RVineStructure

__all__ = ["VinecopBase"]

_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10

#: (tree, edge, u_e, x_e) -> fitted pair copula. The seam external packages
#: drive conditional fitting through (see :meth:`VinecopBase.fit`).
FitEdge = Callable[[int, int, Any, Optional[Any]], BicopLike]


class _NotBatchable(Exception):
  """Raised by :meth:`VinecopBase._build_batched` when batching is unavailable.

  The dispatch layer catches it and falls back to the non-batched cascade.
  """


class VinecopBase(VinecopLike[ArrayT], ABC):
  """Canonical array-agnostic vine cascades (numpy / torch).

  Concrete subclasses implement ``_get_pair_copula`` (and optionally override
  ``_prep`` / ``_simulate_uniform`` / the batched-path hooks) and call
  ``_bind_vine`` once; they
  then inherit the whole evaluator surface — ``pdf`` / ``cdf`` / ``rosenblatt`` /
  ``inverse_rosenblatt`` / ``simulate``, ``loglik`` / ``plot`` / ``__repr__``,
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
  _context: ConditioningContext
  _cond_pos_cache: dict[tuple[int, int], tuple[int, ...]]
  #: Lazily-built grid-batched state (see :meth:`_build_batched`); ``None`` until
  #: the first batched call. Subclasses invalidate it on device moves.
  _batched: Any

  def _bind_vine(
    self,
    structure: RVineStructure,
    context: Optional[ConditioningContext] = None,
  ) -> None:
    """Install the vine structure + context and derive the order arrays.

    The initialization seam a concrete subclass calls once from its ``__init__``
    (after storing its pair copulas). ``context`` defaults to
    :class:`~pyvinecopulib.core.SimplifiedContext` — the unconditional /
    simplified vine that covers the common case — so most subclasses pass only a
    ``structure``. Advanced subclasses may override this method to install extra
    state, calling ``super()._bind_vine(...)``.

    Parameters
    ----------
    structure : RVineStructure
        The (fixed) vine structure to evaluate along.
    context : ConditioningContext, optional
        Per-edge conditioning-context policy; ``None`` uses
        :class:`~pyvinecopulib.core.SimplifiedContext`.

    Returns
    -------
    None
        The structure, context, and derived order arrays are stored on ``self``.
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

  def _prep(self, u: ArrayT, name: str) -> ArrayT:
    """Coerce ``u`` to the working array and clamp to the unit box.

    Concrete default: shape-check, then clamp to ``[1e-10, 1 - 1e-10]`` on
    ``u``'s own array namespace. A subclass that needs dtype / device coercion
    (e.g. accepting NumPy input on a torch vine) overrides this.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations to prepare.
    name : str
        Calling-method name, used only in the shape-error message.

    Returns
    -------
    array, shape (n, d), dtype float
        ``u`` coerced to the working array and clamped to ``[1e-10, 1 - 1e-10]``.
    """
    ua: Any = u
    xp = array_namespace(ua)
    if ua.ndim != 2 or ua.shape[1] != self.d:
      raise ValueError(
        f"{name}: u must have shape (n, {self.d}); got {tuple(ua.shape)}"
      )
    return cast(ArrayT, xp.clip(ua, _TRIM_LO, _TRIM_HI))

  def _simulate_uniform(self, n: int, qrng: bool, seeds: list[int]) -> ArrayT:
    """Draw ``(n, d)`` base uniforms for :meth:`simulate` (namespace-dependent RNG).

    Raising default; override it (numpy / torch differ on RNG) to enable
    :meth:`simulate`. Named after :func:`pyvinecopulib.utils.simulate_uniform`.

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
      f"{type(self).__name__} does not implement _simulate_uniform; override it "
      "to enable simulate()."
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
      self._batched = self._build_batched()
    return self._batched

  def _eval_context(self):
    """Context manager disabling grad for inverse / simulate / cdf.

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

  # --- non-batched cascades (single source of truth) -------------------- #
  def _pdf(self, u: Any, x: Optional[Any]) -> Any:
    """Vine density as a product of per-edge copula densities (``Vinecop::pdf``).

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Prepared pseudo-observations (natural-order seeding happens inside).
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
        # m: min-array — natural-order index of the column finalized in a
        # previous tree; the second pair input comes from hfunc2 when m sits on
        # the natural-order diagonal, else from hfunc1 (class.ipp:1026-1034).
        m = int(s.min_array(tree, edge))
        on_diagonal = m == int(s.struct_array(tree, edge, natural_order=True))
        u_e_col0 = hfunc2[:, edge]
        u_e_col1 = hfunc2[:, m - 1] if on_diagonal else hfunc1[:, m - 1]
        u_e = xp.stack([u_e_col0, u_e_col1], axis=-1)
        x_e = self._edge_context(tree, edge, x, u_nat, None)
        # Accumulate the density as a product over edges (cwiseProduct,
        # class.ipp:1047).
        pdf = pdf * edge_copula.pdf(u_e, x_e)
        # h-functions only evaluated if a later tree needs them (class.ipp:1050).
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = edge_copula.hfunc1(u_e, x_e)
        if s.needed_hfunc2(tree, edge):
          hfunc2[:, edge] = edge_copula.hfunc2(u_e, x_e)
    return pdf

  def _rosenblatt(self, u: Any, x: Optional[Any]) -> Any:
    """Rosenblatt transform (``Vinecop::rosenblatt``).

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Prepared pseudo-observations.
    x : array, shape (n, p), or None
        External covariates threaded to each pair copula, or ``None``.

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
    u_nat = (
      xp.asarray(hfunc2, copy=True)
      if self._context.assembles_conditioning
      else None
    )
    s = self.structure
    for tree in range(trunc_lvl):
      for edge in range(d - tree - 1):
        edge_copula = self._get_pair_copula(tree, edge)
        # See _pdf for the m / on-diagonal second-input rule (class.ipp:1026).
        m = int(s.min_array(tree, edge))
        on_diagonal = m == int(s.struct_array(tree, edge, natural_order=True))
        u_e_col0 = hfunc2[:, edge]
        u_e_col1 = hfunc2[:, m - 1] if on_diagonal else hfunc1[:, m - 1]
        u_e = xp.stack([u_e_col0, u_e_col1], axis=-1)
        x_e = self._edge_context(tree, edge, x, u_nat, None)
        # hfunc1 only if needed downstream; hfunc2 is the running transform.
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = edge_copula.hfunc1(u_e, x_e)
        hfunc2[:, edge] = edge_copula.hfunc2(u_e, x_e)
    # Scatter the transformed columns back to variable order.
    out = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      out[:, j] = hfunc2[:, inv[j]]
    return xp.clip(out, _TRIM_LO, _TRIM_HI)

  def _inverse_rosenblatt(self, u: Any, x: Optional[Any]) -> Any:
    """Inverse Rosenblatt transform (``Vinecop::inverse_rosenblatt``).

    Walks variables from ``d - 2`` down to ``0``; at each ``var`` it fills the
    ``hinv2`` column from the outermost tree inward. The ``(trunc_lvl + 1, d, n)``
    scratch is transposed relative to the forward cascades (variable axis first)
    so a finalized ``hinv2[0, var, :]`` row can seed later inversions.

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
        hinv2[tree, var, :] = edge_copula.hinv2(u_e, x_e)
        # Propagate hfunc1 for the next-inner inversion when needed.
        if var < d - 1 and s.needed_hfunc1(tree, var):
          u_e_after = xp.stack([hinv2[tree, var, :], u_e_col1], axis=-1)
          hfunc1[tree + 1, var, :] = edge_copula.hfunc1(u_e_after, x_e)
    out = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      out[:, j] = hinv2[0, inv[j], :]
    return xp.clip(out, _TRIM_LO, _TRIM_HI)

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
    xp = array_namespace(u)
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
    xp = array_namespace(u)
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
    return xp.clip(out, _TRIM_LO, _TRIM_HI)

  # --- batched dispatch ------------------------------------------------- #
  def _resolve_batched(
    self, requested: Optional[bool], x: Optional[Any]
  ) -> bool:
    """Resolve the ``batched`` flag; force ``False`` for conditional calls.

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
    if requested is None:
      requested = self._default_batched()
    return bool(requested)

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
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d`` (clamped to
        ``[1e-10, 1 - 1e-10]``).
    num_threads : int, default=1
        Accepted for parity with ``Vinecop.pdf()``; ignored.
    x : array, shape (n, p), or None, optional
        External covariates threaded to each pair copula (non-simplified
        vines). ``None`` for the unconditional / simplified case.
    batched : bool or None, optional
        Fire one batched pair-copula call per tree level. ``None`` resolves
        via the subclass default; forced ``False`` when conditioning is active,
        and falls back to the non-batched cascade if the subclass has no batched
        fast path (or a pair does not support it).

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
    x: Optional[ArrayT] = None,
    batched: Optional[bool] = None,
  ) -> ArrayT:
    """Rosenblatt transform: dependent uniforms to independent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d``.
    num_threads : int, default=1
        Accepted for parity; ignored.
    x : array, shape (n, p), or None, optional
        External covariates (non-simplified vines), else ``None``.
    batched : bool or None, optional
        See :meth:`pdf`.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms in ``[1e-10, 1 - 1e-10]``.
    """
    del num_threads
    u_p = self._prep(u, "rosenblatt")
    if self._resolve_batched(batched, x):
      try:
        return cast(ArrayT, self._rosenblatt_batched(u_p))
      except _NotBatchable:
        pass  # no grid fast path available -> non-batched cascade
    return cast(ArrayT, self._rosenblatt(u_p, x))

  def inverse_rosenblatt(
    self,
    u: ArrayT,
    *,
    num_threads: int = 1,
    x: Optional[ArrayT] = None,
    batched: Optional[bool] = None,
  ) -> ArrayT:
    """Inverse Rosenblatt transform: independent uniforms to dependent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Independent uniforms in ``[0, 1]^d``.
    num_threads : int, default=1
        Accepted for parity; ignored.
    x : array, shape (n, p), or None, optional
        External covariates (non-simplified vines), else ``None``.
    batched : bool or None, optional
        Must be ``False`` / ``None``. The inverse cascade has cross-tree
        dependencies that the per-tree-level wavefront cannot satisfy, so
        ``batched=True`` raises :exc:`NotImplementedError`.

    Returns
    -------
    array, shape (n, d), dtype float
        Dependent uniforms in ``[1e-10, 1 - 1e-10]``.

    Raises
    ------
    NotImplementedError
        If ``batched`` is ``True``.
    """
    del num_threads
    if batched:
      raise NotImplementedError(
        "batched=True is not implemented for inverse_rosenblatt. The inverse "
        "cascade has a 2-D dependency graph that does not reduce to "
        "per-tree-level waves. Pass batched=False instead."
      )
    u_p = self._prep(u, "inverse_rosenblatt")
    with self._eval_context():
      return cast(ArrayT, self._inverse_rosenblatt(u_p, x))

  def simulate(
    self,
    n: int,
    *,
    qrng: bool = False,
    num_threads: int = 1,
    seeds: Optional[list[int]] = None,
    x: Optional[ArrayT] = None,
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
          f"simulate(n={n}, x=...) requires one covariate row per sample; "
          f"got x with {xa.shape[0]} rows."
        )
    with self._eval_context():
      base_u = self._simulate_uniform(n, qrng, seeds)
      return self.inverse_rosenblatt(base_u, x=x)

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
  ) -> ArrayT:
    """Evaluate the joint CDF at each query row via quasi-Monte-Carlo.

    Parameters
    ----------
    u : array, shape (m, d), dtype float
        Query points in ``[0, 1]^d``.
    N : int, default=10000
        Number of Monte-Carlo samples.
    qrng : bool, default=True
        Draw quasi-random samples (matches ``Vinecop.cdf()``).
    num_threads : int, default=1
        Accepted for parity; ignored.
    seeds : list of int or None, optional
        RNG seeds forwarded to :meth:`simulate`.
    x : array or None, optional
        Must be ``None`` — a per-row conditional CDF is not supported (the
        Monte-Carlo dominance estimate cannot condition each query row on a
        different covariate without per-``x`` resampling).
    block_size : int, default=4096
        Query rows processed per iteration (peak-memory control).

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
      u_t: Any = self._prep(u, "cdf")
      samples: Any = self.simulate(N, qrng=qrng, seeds=seeds)
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
    u : array, shape (n, d), dtype float
        Pseudo-observations.
    fit_edge : callable
        ``(tree, edge, u_e, x_e) -> BicopLike`` fitting one edge's pair copula.
    context : ConditioningContext, optional
        Conditioning-context policy (default: simplified / unconditional).
    x : array, shape (n, p), or None, optional
        External covariates for conditional fitting, else ``None``.

    Returns
    -------
    list of list of BicopLike
        Fitted pair copulas indexed ``[tree][edge]``.

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
    order = [int(v) for v in structure.order]
    n = ua.shape[0]
    hfunc1 = xp.zeros((n, d), dtype=ua.dtype, device=ua.device)
    hfunc2 = xp.empty((n, d), dtype=ua.dtype, device=ua.device)
    for j in range(d):
      hfunc2[:, j] = ua[:, order[j] - 1]
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
        # Same m / on-diagonal second-input rule as _pdf (class.ipp:1026).
        m = int(s.min_array(tree, edge))
        on_diagonal = m == int(s.struct_array(tree, edge, natural_order=True))
        u_e_col0 = hfunc2[:, edge]
        u_e_col1 = hfunc2[:, m - 1] if on_diagonal else hfunc1[:, m - 1]
        u_e = xp.stack([u_e_col0, u_e_col1], axis=-1)
        x_e = edge_context_for(tree, edge)
        edge_copula = fit_edge(tree, edge, u_e, x_e)
        row.append(edge_copula)
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = edge_copula.hfunc1(u_e, x_e)
        if s.needed_hfunc2(tree, edge):
          hfunc2[:, edge] = edge_copula.hfunc2(u_e, x_e)
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
    u : array, shape (n, d), dtype float
        Pseudo-observations on any array-API namespace (NumPy or PyTorch).
    fit_edge : callable
        ``(tree, edge, u_e, x_e) -> BicopLike`` fitting one edge's pair copula;
        its ``hfunc1`` / ``hfunc2`` must be valid immediately, and it must
        implement :meth:`~pyvinecopulib.core.BicopLike.flip` (used to
        reorient reused pairs onto
        their finalized slots). ``x_e`` is always ``None`` here.
    trunc_lvl : int, optional
        Maximum number of trees to select (default: ``d - 1``, i.e. untruncated).
    tree_criterion : str, default "tau"
        Dependence measure passed to ``wdm`` for edge weighting (e.g. ``"tau"``,
        ``"rho"``, ``"hoeffd"``).
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
    d = int(u.shape[1])
    seed_list = [int(s) for s in (seeds or [])]
    convert = np.asarray if to_numpy is None else to_numpy
    max_trees = (
      d - 1 if trunc_lvl is None else max(0, min(int(trunc_lvl), d - 1))
    )
    # Sentinel "root" shared by every base-tree node, so the first tree's
    # candidate graph is complete (the C++ base tree is a star).
    root = d

    def criterion(col0: Any, col1: Any) -> float:
      # |wdm| with the selector's n > 10 guard and NaN -> 0. The weighted /
      # tie corrections are skipped: this path is continuous and unweighted.
      if n <= 10:
        return 0.0
      value = float(abs(wdm(convert(col0), convert(col1), tree_criterion)))
      return 0.0 if not np.isfinite(value) else value

    # A node is one edge of the previous tree (a single variable for the base
    # tree). ``prev`` holds the two previous-tree vertex ids that this edge
    # joined; a shared prev id is the proximity condition and picks which
    # h-function feeds the next tree.
    nodes: list[dict[str, Any]] = [
      {
        "all_indices": (i,),
        "h1": u[:, i],
        "h2": u[:, i],
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
          col0 = (
            nodes[v0]["h1"] if prev0.index(common) == 0 else nodes[v0]["h2"]
          )
          col1 = (
            nodes[v1]["h1"] if prev1.index(common) == 0 else nodes[v1]["h2"]
          )
          tau = criterion(col0, col1)
          cand.append((v0, v1))
          cand_cols.append((col0, col1))
          weights.append(1.0 - (tau >= threshold) * tau)

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
        u_e = xp.stack([col0, col1], axis=-1)
        pair = fit_edge(len(trees), edge_idx, u_e, None)
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
            "h1": pair.hfunc1(u_e, None),
            "h2": pair.hfunc2(u_e, None),
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


# Shares the worked example with :class:`~pyvinecopulib.core.VinecopLike` (see
# protocols.py) so the contract and its canonical base never drift apart.
VinecopBase.__doc__ = (VinecopBase.__doc__ or "") + _VINECOP_EXAMPLE
