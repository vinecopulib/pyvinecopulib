"""Canonical partial implementation of the vine-copula evaluator contract.

:class:`VinecopBase` is the array-backend-agnostic (NumPy / PyTorch, via
:func:`array_api_compat.array_namespace`) implementation of the vine cascades —
``pdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` / ``simulate`` / ``cdf`` — plus
the shared sequential-fit engine. It is a direct port of the C++
``Vinecop::pdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` tree-by-tree
h-function cascade, so a concrete backend (e.g.
:class:`pyvinecopulib.torch.TorchVinecop`) only supplies a small set of hooks and
inherits the whole evaluation surface.

Conditioning is threaded through a pluggable
:class:`~pyvinecopulib.core.ContextPolicy`: each pair-copula call receives a
``cond`` matrix assembled per edge from the edge's conditioning-set values
``u_D`` and an optional external covariate matrix ``x``. The default
:class:`~pyvinecopulib.core.SimplifiedContext` passes ``cond=None`` everywhere
and skips the ``u_D`` gather, reproducing the classic simplified cascade at zero
extra cost.

Hooks a concrete backend must provide (see the abstract members): ``_pair``,
``_prep``, ``_draw_base_u``; and, to enable the batched fast path,
``_pdf_batched`` / ``_rosenblatt_batched`` / ``_default_batched``. The
concrete backend calls :meth:`_bind_vine` once to install the structure and
context. Array values are handled as ``Any`` inside the cascades per the
``pyvinecopulib.core`` typing policy (the Array API namespace is untyped); the
generic ``ArrayT`` lives on the public signatures.
"""

from __future__ import annotations

import contextlib
from abc import ABC, abstractmethod
from typing import Any, Callable, Optional, cast

from array_api_compat import array_namespace

from ._context import ContextPolicy, SimplifiedContext
from ._protocols import ArrayT, BicopLike, VinecopLike

__all__ = ["VinecopBase"]

_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10

#: (tree, edge, u_e, cond) -> fitted pair copula. The seam external packages
#: drive conditional fitting through (see :meth:`VinecopBase.sequential_fit`).
FitEdge = Callable[[int, int, Any, Optional[Any]], BicopLike]


class _NotBatchable(Exception):
  """Raised by :meth:`VinecopBase._build_batched` when batching is unavailable.

  The dispatch layer catches it and falls back to the non-batched cascade.
  """


class VinecopBase(VinecopLike[ArrayT], ABC):
  """Canonical array-agnostic vine cascades (numpy / torch).

  Concrete subclasses implement ``_pair`` / ``_prep`` / ``_draw_base_u`` (and,
  optionally, the batched-path hooks) and call :meth:`_bind_vine` once; they
  then inherit ``pdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` / ``simulate``
  / ``cdf`` and the :meth:`sequential_fit` engine. Not an ``nn.Module``.
  """

  # --- layout installed by _bind_vine (hooks / state) ------------------- #
  structure: Any
  d: int
  trunc_lvl: int
  order: tuple[int, ...]
  inverse_order: tuple[int, ...]
  _context: ContextPolicy
  _cond_pos_cache: dict[tuple[int, int], tuple[int, ...]]

  def _bind_vine(self, structure: Any, context: ContextPolicy) -> None:
    """Install the vine structure + context and derive the order arrays."""
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

  # --- hooks a concrete backend provides -------------------------------- #
  @abstractmethod
  def _pair(self, tree: int, edge: int) -> BicopLike[ArrayT]:
    """Return the pair copula at ``(tree, edge)``."""

  @abstractmethod
  def _prep(self, u: ArrayT, name: str) -> ArrayT:
    """Coerce ``u`` to the working dtype/device and clamp to the unit box."""

  @abstractmethod
  def _draw_base_u(self, n: int, qrng: bool, seeds: list[int]) -> ArrayT:
    """Draw ``(n, d)`` base uniforms for :meth:`simulate` (RNG is backend-specific)."""

  def _default_batched(self) -> bool:
    """Whether ``batched`` defaults to ``True`` (backend/device dependent)."""
    return False

  def _pdf_batched(self, u: Any) -> Any:
    raise NotImplementedError(
      f"{type(self).__name__} does not provide a batched pdf path."
    )

  def _rosenblatt_batched(self, u: Any) -> Any:
    raise NotImplementedError(
      f"{type(self).__name__} does not provide a batched rosenblatt path."
    )

  def _eval_context(self):
    """Context manager disabling grad for inverse / simulate / cdf.

    Defaults to a no-op; the torch backend overrides it with
    ``torch.no_grad()``.
    """
    return contextlib.nullcontext()

  # --- conditioning-context assembly ------------------------------------ #
  def _cond_positions(self, tree: int, edge: int) -> tuple[int, ...]:
    key = (tree, edge)
    cache = self._cond_pos_cache
    if key not in cache:
      s = self.structure
      cache[key] = tuple(
        int(s.struct_array(i, edge, natural_order=True)) - 1
        for i in range(tree)
      )
    return cache[key]

  def _edge_cond(
    self,
    tree: int,
    edge: int,
    x: Optional[Any],
    u_nat: Optional[Any],
    hinv2_final: Optional[Any],
  ) -> Optional[Any]:
    """Assemble ``cond`` for one edge (forward: ``u_nat``; inverse: ``hinv2_final``)."""
    ctx = self._context
    if not ctx.assembles_conditioning and x is None:
      return None
    u_D: Optional[Any] = None
    if ctx.assembles_conditioning:
      dpos = self._cond_positions(tree, edge)
      if dpos:
        cols = list(dpos)
        if u_nat is not None:
          u_D = u_nat[:, cols]
        else:
          hf: Any = hinv2_final
          xp = array_namespace(hf)
          u_D = xp.matrix_transpose(hf[cols, :])
    return ctx.edge_context(u_D=u_D, x=x)

  # --- non-batched cascades (single source of truth) -------------------- #
  def _pdf(self, u: Any, x: Optional[Any]) -> Any:
    xp = array_namespace(u)
    d, trunc_lvl = self.d, self.trunc_lvl
    n = u.shape[0]
    if trunc_lvl == 0:
      return xp.ones(n, dtype=u.dtype, device=u.device)
    hfunc1 = xp.zeros((n, d), dtype=u.dtype, device=u.device)
    hfunc2 = xp.empty((n, d), dtype=u.dtype, device=u.device)
    order = self.order
    for j in range(d):
      hfunc2[:, j] = u[:, order[j] - 1]
    u_nat = (
      xp.asarray(hfunc2, copy=True)
      if self._context.assembles_conditioning
      else None
    )
    log_pdf = xp.zeros(n, dtype=u.dtype, device=u.device)
    s = self.structure
    for tree in range(trunc_lvl):
      for edge in range(d - tree - 1):
        cop = self._pair(tree, edge)
        m = int(s.min_array(tree, edge))
        sarr = int(s.struct_array(tree, edge, natural_order=True))
        col0 = hfunc2[:, edge]
        col1 = hfunc2[:, m - 1] if m == sarr else hfunc1[:, m - 1]
        u_e = xp.stack([col0, col1], axis=-1)
        cond = self._edge_cond(tree, edge, x, u_nat, None)
        log_pdf = log_pdf + cop.log_pdf(u_e, cond)
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = cop.hfunc1(u_e, cond)
        if s.needed_hfunc2(tree, edge):
          hfunc2[:, edge] = cop.hfunc2(u_e, cond)
    return xp.exp(log_pdf)

  def _rosenblatt(self, u: Any, x: Optional[Any]) -> Any:
    xp = array_namespace(u)
    d, trunc_lvl = self.d, self.trunc_lvl
    n = u.shape[0]
    order, inv = self.order, self.inverse_order
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
        cop = self._pair(tree, edge)
        m = int(s.min_array(tree, edge))
        sarr = int(s.struct_array(tree, edge, natural_order=True))
        col0 = hfunc2[:, edge]
        col1 = hfunc2[:, m - 1] if m == sarr else hfunc1[:, m - 1]
        u_e = xp.stack([col0, col1], axis=-1)
        cond = self._edge_cond(tree, edge, x, u_nat, None)
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = cop.hfunc1(u_e, cond)
        hfunc2[:, edge] = cop.hfunc2(u_e, cond)
    out = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      out[:, j] = hfunc2[:, inv[j]]
    return xp.clip(out, _TRIM_LO, _TRIM_HI)

  def _inverse_rosenblatt(self, u: Any, x: Optional[Any]) -> Any:
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
        cop = self._pair(tree, var)
        m = int(s.min_array(tree, var))
        sarr = int(s.struct_array(tree, var, natural_order=True))
        col0 = hinv2[tree + 1, var, :]
        col1 = hinv2[tree, m - 1, :] if m == sarr else hfunc1[tree, m - 1, :]
        u_e = xp.stack([col0, col1], axis=-1)
        cond = self._edge_cond(tree, var, x, None, hinv2[0])
        hinv2[tree, var, :] = cop.hinv2(u_e, cond)
        if var < d - 1 and s.needed_hfunc1(tree, var):
          u_e_after = xp.stack([hinv2[tree, var, :], col1], axis=-1)
          hfunc1[tree + 1, var, :] = cop.hfunc1(u_e_after, cond)
    out = xp.empty((n, d), dtype=u.dtype, device=u.device)
    for j in range(d):
      out[:, j] = hinv2[0, inv[j], :]
    return xp.clip(out, _TRIM_LO, _TRIM_HI)

  # --- batched dispatch ------------------------------------------------- #
  def _resolve_batched(
    self, requested: Optional[bool], x: Optional[Any]
  ) -> bool:
    """Resolve the ``batched`` flag; force ``False`` for conditional calls."""
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
        Accepted for parity with ``pv.Vinecop.pdf``; ignored.
    x : array, shape (n, p), or None, optional
        External covariates threaded to each pair copula (non-simplified
        vines). ``None`` for the unconditional / simplified case.
    batched : bool or None, optional
        Fire one batched pair-copula call per tree level. ``None`` resolves
        via the backend default; forced ``False`` when conditioning is active.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values.
    """
    del num_threads
    u_p = self._prep(u, "pdf")
    if self._resolve_batched(batched, x):
      return cast(ArrayT, self._pdf_batched(u_p))
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
      return cast(ArrayT, self._rosenblatt_batched(u_p))
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
        RNG seeds forwarded to the backend's base-uniform draw.
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
      base_u = self._draw_base_u(n, qrng, seeds)
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
        Draw quasi-random samples (matches ``pv.Vinecop.cdf``).
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

  # --- shared sequential-fit engine ------------------------------------- #
  @staticmethod
  def sequential_fit(
    structure: Any,
    u: Any,
    fit_edge: FitEdge,
    *,
    context: Optional[ContextPolicy] = None,
    x: Optional[Any] = None,
  ) -> list[list[BicopLike]]:
    """Fit pair copulas tree-by-tree along the vine, returning the nested list.

    Mirrors the forward pdf traversal with the density evaluation replaced by
    ``fit_edge(tree, edge, u_e, cond)`` — the Python analogue of the C++
    ``Vinecop::fit`` ``fit_edge`` lambda. The returned pair's ``hfunc1`` /
    ``hfunc2`` must be valid immediately for tree propagation. External
    packages drive conditional fitting through this seam (a ``fit_edge`` that
    fits a conditional pair copula on ``(u_e, cond)``).

    Parameters
    ----------
    structure : RVineStructure
        The (fixed) vine structure to fit along.
    u : array, shape (n, d), dtype float
        Pseudo-observations.
    fit_edge : callable
        ``(tree, edge, u_e, cond) -> BicopLike`` fitting one edge's pair copula.
    context : ContextPolicy, optional
        Conditioning-context policy (default: simplified / unconditional).
    x : array, shape (n, p), or None, optional
        External covariates for conditional fitting, else ``None``.

    Returns
    -------
    list of list of BicopLike
        Fitted pair copulas indexed ``[tree][edge]``.
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

    def cond_for(tree: int, edge: int) -> Optional[Any]:
      ctx = context
      if not ctx.assembles_conditioning and x is None:
        return None
      u_D: Optional[Any] = None
      if ctx.assembles_conditioning:
        key = (tree, edge)
        if key not in cache:
          cache[key] = tuple(
            int(structure.struct_array(i, edge, natural_order=True)) - 1
            for i in range(tree)
          )
        dpos = cache[key]
        if dpos and u_nat is not None:
          u_D = u_nat[:, list(dpos)]
      return ctx.edge_context(u_D=u_D, x=x)

    s = structure
    pairs: list[list[BicopLike]] = []
    for tree in range(trunc_lvl):
      row: list[BicopLike] = []
      for edge in range(d - tree - 1):
        m = int(s.min_array(tree, edge))
        sarr = int(s.struct_array(tree, edge, natural_order=True))
        col0 = hfunc2[:, edge]
        col1 = hfunc2[:, m - 1] if m == sarr else hfunc1[:, m - 1]
        u_e = xp.stack([col0, col1], axis=-1)
        cond = cond_for(tree, edge)
        cop = fit_edge(tree, edge, u_e, cond)
        row.append(cop)
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = cop.hfunc1(u_e, cond)
        if s.needed_hfunc2(tree, edge):
          hfunc2[:, edge] = cop.hfunc2(u_e, cond)
      pairs.append(row)
    return pairs
