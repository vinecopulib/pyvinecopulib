"""PyTorch evaluator for a fitted R-vine copula.

Wraps a fitted :class:`pyvinecopulib.Vinecop` (with *Transformed Local
Likelihood* pair copulas — :data:`pyvinecopulib.families.tll`) and
exposes ``pdf`` / ``cdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` /
``simulate`` on top of :class:`TorchBicop` for every pair copula. The
whole evaluation chain stays in PyTorch, so the vine can move to GPU
with ``.to("cuda")`` and be composed with autograd-aware downstream
code.

The two routes to a fitted vine are :meth:`TorchVinecop.from_vinecop`
(lift a C++-fitted :class:`pyvinecopulib.Vinecop`) and
:meth:`TorchVinecop.from_data` (fit directly from pseudo-observations
in pure PyTorch given a fixed structure). Fit-time controls live on
:class:`FitControlsTorchVinecop`, which mirrors
:class:`pyvinecopulib.FitControlsVinecop`.

The cascade is a direct port of the C++
:class:`pyvinecopulib.Vinecop` tree-by-tree h-function chain in
``Vinecop::pdf`` / ``Vinecop::rosenblatt`` /
``Vinecop::inverse_rosenblatt``: dense ``(n, d)`` scratch matrices,
fixed natural-order traversal, byte-for-byte agreement with the C++
evaluator on the same fit. ``pdf`` / ``rosenblatt`` additionally
accept ``batched=True`` to fire one stacked bicop call per tree
level.

See Also
--------
pyvinecopulib.Vinecop : The C++ counterpart.
FitControlsTorchVinecop : Fit-time controls.
"""

from __future__ import annotations

from typing import Optional, cast

import numpy as np
import torch
from torch import Tensor

from ..core import ContextPolicy, SimplifiedContext, VinecopBase
from ..pyvinecopulib_ext import (
  RVineStructure,
  Vinecop,
  indep as _INDEP_FAMILY,
  tll as _TLL_FAMILY,
)
from ..utils import simulate_uniform as _simulate_uniform
from ._batched import BatchedVine
from ._controls import FitControlsTorchVinecop
from ._interp import _TRIM_HI, _TRIM_LO
from .bicop import TorchBicop


class TorchVinecop(VinecopBase[torch.Tensor], torch.nn.Module):
  """PyTorch R-vine copula evaluator built on `TorchBicop`.

  Mirrors the public surface of `pyvinecopulib.Vinecop` — ``pdf`` /
  ``rosenblatt`` / ``inverse_rosenblatt`` / ``simulate`` — but keeps
  the entire evaluation chain in PyTorch so the vine can move to
  GPU with ``.to(device)`` and compose with autograd-aware
  downstream code. Continuous variables only, single batch, no
  threading. Build a `TorchVinecop` via `from_data` (pure-torch TLL
  fit) or `from_vinecop` (lifts a C++-fitted `Vinecop`).

  The cascade is a dense ``(n, d)``-scratch port of the C++
  evaluator, byte-for-byte parity with `Vinecop`. ``pdf`` /
  ``rosenblatt`` accept ``batched=True`` to fire a single batched
  bicop call per tree level. ``inverse_rosenblatt(batched=True)``
  raises (cross-tree deps in the inverse cascade). The default
  ``batched=None`` resolves to ``True`` on CUDA (3–7x faster) and
  ``False`` on CPU.

  Parameters
  ----------
  pair_copulas : list of list of TorchBicop
      Indexed ``[tree][edge]`` and shaped like
      ``Vinecop::pair_copulas_``: tree 0 has ``d - 1`` edges, tree 1
      has ``d - 2``, etc., up to ``trunc_lvl``. May also be passed
      as a `torch.nn.ModuleList` of `torch.nn.ModuleList`.
  structure : RVineStructure
      Vine structure whose accessors (``min_array``,
      ``struct_array``, ``needed_hfunc1`` / ``needed_hfunc2``)
      describe how to walk the trees.
  context : ContextPolicy or None, default=None
      Conditioning-context policy that assembles each pair copula's
      ``cond`` per edge. ``None`` uses
      :class:`~pyvinecopulib.core.SimplifiedContext` (an unconditional /
      simplified vine).
  """

  d: int
  trunc_lvl: int
  pair_copulas: torch.nn.ModuleList

  def __init__(
    self,
    pair_copulas: list[list[TorchBicop]],
    structure,
    *,
    context: Optional[ContextPolicy] = None,
  ) -> None:
    # Initialise nn.Module explicitly: TorchVinecop also subclasses VinecopBase
    # (a Protocol-derived ABC), whose __init__ chain would otherwise shadow
    # nn.Module's under super().
    torch.nn.Module.__init__(self)
    if context is None:
      context = SimplifiedContext()
    # Install structure + context + derived order arrays (VinecopBase hook).
    self._bind_vine(structure, context)

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

    # Lazy-built stacked-per-tree-level state (see `batched=` on the public
    # methods). Constructed on first batched call; cleared in `_apply` so
    # device moves invalidate it.
    self._batched: Optional["BatchedVine"] = None

  # --------------------------------------------------------------------- #
  # Constructor                                                            #
  # --------------------------------------------------------------------- #

  @classmethod
  def from_vinecop(
    cls,
    cop: Vinecop,
    cache_integrals: bool = True,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchVinecop":
    """Lifts a fitted `Vinecop` into the torch backend.

    Each pair copula is lifted via `TorchBicop.from_bicop`, so the
    resulting cascade matches the C++ evaluator to within
    bilinear-interpolation precision on the shared grid.

    Parameters
    ----------
    cop : Vinecop
        A fitted `pyvinecopulib.Vinecop` whose pair copulas are all
        TLL or independence families. Discrete variables are not
        supported.
    cache_integrals : bool, default=True
        Forwarded to `TorchBicop.from_bicop` for every pair copula.
    device : torch.device or None, default=None
        Placement of the underlying tensors.
    dtype : torch.dtype, default=torch.float64
        Precision of the underlying tensors.

    Returns
    -------
    TorchVinecop
        A `TorchVinecop` mirroring ``cop``.
    """
    var_types = list(cop.var_types)
    if any(v != "c" for v in var_types):
      raise ValueError(
        "TorchVinecop is continuous-only; got var_types="
        f"{var_types!r}. Discrete variables are not yet supported."
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

    return cls(pair_copulas=pair_copulas_torch, structure=cop.structure)

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
    """Builds a `TorchVinecop` from a structure (or matrix) and pair-copulas.

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
        Per-variable types. Only ``"c"`` (continuous) is supported.
    device : torch.device or None, default=None
        Placement of the independence pair copulas that fill
        missing edges.
    dtype : torch.dtype, default=torch.float64
        Precision of the underlying tensors.

    Returns
    -------
    TorchVinecop
        A `TorchVinecop`.
    """
    if (structure is None) == (matrix is None):
      raise ValueError("Provide exactly one of `structure` or `matrix`.")
    if structure is None:
      structure = RVineStructure.from_matrix(np.asarray(matrix))

    d = int(structure.dim)
    if var_types and len(var_types) != d:
      raise ValueError(f"var_types has {len(var_types)} entries, expected {d}")
    if any(v != "c" for v in var_types):
      raise NotImplementedError(
        "TorchVinecop is continuous-only; got var_types="
        f"{var_types!r}. Discrete variables are not yet supported."
      )

    trunc_lvl = int(structure.trunc_lvl)
    if not pair_copulas:
      # Independence vine: TorchBicop() defaults to the independence copula.
      pair_copulas = [
        [TorchBicop(device=device, dtype=dtype) for _ in range(d - 1 - t)]
        for t in range(trunc_lvl)
      ]
    return cls(pair_copulas=pair_copulas, structure=structure)

  @classmethod
  def from_data(
    cls,
    u,
    structure,
    controls: Optional[FitControlsTorchVinecop] = None,
  ) -> "TorchVinecop":
    """Fits a pure-PyTorch TLL vine on ``u`` given a fixed structure.

    Tree-by-tree cascade mirroring ``Vinecop::select_families`` with
    a user-provided structure: at each ``(tree, edge)`` collects the
    pair of pseudo-obs columns by the same rule as the legacy
    pdf/rosenblatt cascade, fits a `TorchBicop` on them via
    `TorchBicop.from_data`, then propagates ``hfunc1`` / ``hfunc2``
    forward. The structure is *fixed* (not selected); the cascade
    matches the C++ TLL fit to machine precision when
    ``cache_integrals=False``.

    Parameters
    ----------
    u : ndarray or Tensor, shape (n, d), dtype float
        Pseudo-observations.
    structure : RVineStructure
        Vine skeleton.
    controls : FitControlsTorchVinecop or None, default=None
        Pair-copula fit controls bundled with vine-level placement /
        cascade knobs. `None` defaults to TLL on a 30x30
        normal-spaced grid, float64.

    Returns
    -------
    TorchVinecop
        A fitted `TorchVinecop`.
    """
    if controls is None:
      controls = FitControlsTorchVinecop()

    eff_dtype = controls.dtype if controls.dtype is not None else torch.float64
    eff_device = controls.device
    u_t = torch.as_tensor(u, dtype=eff_dtype, device=eff_device)
    if u_t.ndim != 2:
      raise ValueError(f"u must be 2-D; got shape {tuple(u_t.shape)}")
    n, d = u_t.shape
    if int(structure.dim) != d:
      raise ValueError(
        f"structure.dim={structure.dim} does not match u.shape[1]={d}"
      )

    bc_controls = controls.bicop_controls

    def fit_edge(
      tree: int, edge: int, u_e: Tensor, cond: Optional[Tensor]
    ) -> TorchBicop:
      # Simplified (unconditional) TLL fit — cond is None here.
      del tree, edge, cond
      return TorchBicop.from_data(
        u_e,
        bc_controls,
        cache_integrals=controls.cache_integrals,
        device=u_t.device,
        dtype=eff_dtype,
      )

    # Shared tree-by-tree fit engine (SimplifiedContext -> cond=None). The
    # fit_edge closure returns TorchBicop, so the nested list is concrete.
    pairs = cls.sequential_fit(structure, u_t, fit_edge)
    return cls(
      pair_copulas=cast("list[list[TorchBicop]]", pairs), structure=structure
    )

  # --------------------------------------------------------------------- #
  # Helpers                                                                #
  # --------------------------------------------------------------------- #

  def _pair(self, tree: int, edge: int) -> TorchBicop:
    """Indexed access to a pair copula.

    The two-level :class:`torch.nn.ModuleList` is typed as ``Module`` after
    indexing, so :func:`cast` is used to recover the concrete element type
    for the static checker.
    """
    return cast(
      TorchBicop, cast(torch.nn.ModuleList, self.pair_copulas[tree])[edge]
    )

  def _ref_tensor(self) -> Tensor:
    """A registered buffer we can crib dtype/device from."""
    # Every TorchBicop registers its interpolation grid; reuse the first.
    return self._pair(0, 0).interp_grid.values

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

  def _prep(self, u: Tensor, name: str) -> Tensor:
    ref = self._ref_tensor()
    u = torch.as_tensor(u, dtype=ref.dtype, device=ref.device)
    if u.ndim != 2 or u.shape[1] != self.d:
      raise ValueError(
        f"{name}: u must have shape (n, {self.d}); got {tuple(u.shape)}"
      )
    return u.clamp(_TRIM_LO, _TRIM_HI)

  def _apply(self, fn, *args, **kwargs):
    # `.to()`, `.cuda()`, `.cpu()` all route through `_apply`. The
    # BatchedVine container holds buffers — `super()._apply` would move
    # them, but we drop the whole structure so it gets re-baked from the
    # (already-moved) source pair_copulas on next use; that keeps the
    # wire-up tensors aligned with the destination dtype/device.
    self._batched = None
    return super()._apply(fn, *args, **kwargs)

  # --------------------------------------------------------------------- #
  # VinecopBase hooks: RNG for simulate + grad control                     #
  # --------------------------------------------------------------------- #

  def _draw_base_u(self, n: int, qrng: bool, seeds: list[int]) -> Tensor:
    """Draw ``(n, d)`` base uniforms on the fitted grid's dtype/device.

    Pseudo-random via ``torch.rand`` (the first seed seeds a fresh
    ``torch.Generator``), or quasi-random via
    ``pyvinecopulib.utils.simulate_uniform`` when ``qrng=True``.

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
      u_np = _simulate_uniform(n, self.d, qrng=True, seeds=list(seeds))
      return torch.as_tensor(u_np, dtype=dtype, device=device)
    gen: Optional[torch.Generator] = None
    if seeds:
      gen = torch.Generator(device=device).manual_seed(int(seeds[0]))
    return torch.rand(n, self.d, generator=gen, dtype=dtype, device=device)

  def _eval_context(self):
    """Disable autograd for ``inverse_rosenblatt`` / ``simulate`` / ``cdf``.

    Returns
    -------
    torch.autograd.grad_mode.no_grad
        A ``torch.no_grad()`` context manager.
    """
    return torch.no_grad()

  # ====================================================================== #
  # Batched cascades (`batched=True`)                                        #
  # ====================================================================== #
  #
  # Lazy-built BatchedVine (stacked / pre-baked per-tree-level state). Same
  # cascade math as the unbatched paths, but each tree level fires a single
  # batched bicop call instead of a Python loop over edges.

  def _ensure_batched(self) -> "BatchedVine":
    if self._batched is None:
      self._batched = BatchedVine.from_torch_vinecop(self)
    return self._batched

  # ---- pdf ------------------------------------------------------------- #

  def _pdf_batched(self, u: Tensor) -> Tensor:
    """Batched pdf with dense ``(n, d)`` scratch.

    Mirrors :meth:`_pdf` exactly except the inner per-edge loop is
    collapsed into one :meth:`BatchedTreeLevel.pdf` call per tree level.
    """
    u = self._prep(u, "pdf")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d, trunc_lvl = self.d, self.trunc_lvl

    if trunc_lvl == 0:
      return torch.ones(n, dtype=dtype, device=device)

    bv = self._ensure_batched()
    hfunc1 = torch.zeros(n, d, dtype=dtype, device=device)
    hfunc2 = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      hfunc2[:, j] = u[:, self.order[j] - 1]

    log_pdf = torch.zeros(n, dtype=dtype, device=device)
    for t in range(trunc_lvl):
      lvl = bv.level(t)
      u_e = lvl.gather_inputs(hfunc1, hfunc2)  # (N_t, n, 2)
      log_pdf = log_pdf + lvl.log_pdf(bv.grid_points, u_e).sum(dim=0)
      # Update next-tree inputs: compute h1/h2 for every pair, then
      # selectively overwrite columns whose needs_h{1,2} flag is set.
      N_t = lvl.n_pairs
      h1_new = lvl.hfunc1(bv.grid_points, u_e).t()  # (n, N_t)
      h2_new = lvl.hfunc2(bv.grid_points, u_e).t()
      hfunc1[:, :N_t] = torch.where(
        lvl.needs_h1[None, :], h1_new, hfunc1[:, :N_t]
      )
      hfunc2[:, :N_t] = torch.where(
        lvl.needs_h2[None, :], h2_new, hfunc2[:, :N_t]
      )
    return log_pdf.exp()

  # ---- rosenblatt ------------------------------------------------------ #

  def _rosenblatt_batched(self, u: Tensor) -> Tensor:
    """Batched rosenblatt with dense ``(n, d)`` scratch."""
    u = self._prep(u, "rosenblatt")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d, trunc_lvl = self.d, self.trunc_lvl

    bv = self._ensure_batched()
    hfunc2 = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      hfunc2[:, j] = u[:, self.order[j] - 1]
    hfunc1 = hfunc2.clone()

    for t in range(trunc_lvl):
      lvl = bv.level(t)
      u_e = lvl.gather_inputs(hfunc1, hfunc2)
      N_t = lvl.n_pairs
      h1_new = lvl.hfunc1(bv.grid_points, u_e).t()
      h2_new = lvl.hfunc2(bv.grid_points, u_e).t()
      # hfunc2 is unconditionally overwritten at every edge in the
      # cascade (`hfunc2[:, edge] = cop.hfunc2(u_e)`); hfunc1 is
      # gated by needs_h1.
      hfunc2[:, :N_t] = h2_new
      hfunc1[:, :N_t] = torch.where(
        lvl.needs_h1[None, :], h1_new, hfunc1[:, :N_t]
      )

    out = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      out[:, j] = hfunc2[:, self.inverse_order[j]]
    return out.clamp(_TRIM_LO, _TRIM_HI)
