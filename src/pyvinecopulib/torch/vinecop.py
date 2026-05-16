"""PyTorch ``TorchVinecop`` — port of the multivariate vinecop evaluation chain.

Wraps a fitted :class:`pyvinecopulib.Vinecop` (with TLL pair copulas) and
exposes :meth:`pdf`, :meth:`rosenblatt`, and :meth:`inverse_rosenblatt`
on top of :class:`TorchBicop` for every pair copula. The whole evaluation
chain stays in PyTorch, so the vine can move to GPU and be composed with
autograd-aware downstream code.

Each entry point ships with two equivalent implementations, switchable
via the ``impl=`` kwarg:

- ``impl="legacy"`` (default) — direct port of the C++ Vinecop's
  tree-by-tree h-function cascade in :func:`Vinecop::pdf` /
  :func:`Vinecop::rosenblatt` / :func:`Vinecop::inverse_rosenblatt`.
  Dense ``(n, d)`` scratch matrices, fixed pv-natural traversal order,
  byte-for-byte agreement with ``pv.Vinecop``.
- ``impl="lazy"`` — dict-based bookkeeping inspired by
  ``torchvinecopulib.VineCop``. Pseudo-obs are keyed by
  ``(v, *cond_ing)`` and materialized on first access; per-level
  garbage collection keeps live memory smaller, with no compute change.
  ``inverse_rosenblatt`` additionally uses a reference-counted walk to
  free intermediate pseudo-obs as soon as they're no longer needed.

Fitting is delegated to the C++ library; use
:meth:`TorchVinecop.from_vinecop` after fitting with
``pv.Vinecop.from_data(..., controls=FitControlsVinecop(family_set=[pv.tll]))``.
"""

from __future__ import annotations

from collections import Counter
from typing import Optional, cast

import torch
from torch import Tensor

from ..pyvinecopulib_ext import (
  Vinecop,
  indep as _INDEP_FAMILY,
  tll as _TLL_FAMILY,
)
from ._batched import BatchedVine
from ._interp import _TRIM_HI, _TRIM_LO
from .bicop import TorchBicop

_VALID_IMPLS = ("legacy", "lazy")


def _check_impl(impl: str) -> None:
  if impl not in _VALID_IMPLS:
    raise ValueError(f"impl must be one of {_VALID_IMPLS}; got {impl!r}")


class TorchVinecop(torch.nn.Module):
  """PyTorch evaluator over an R-vine of :class:`TorchBicop` pair copulas.

  Port of the C++ ``Vinecop`` class's ``pdf`` / ``rosenblatt`` /
  ``inverse_rosenblatt`` entry points. Continuous-only; single batch; no
  threading. The fit lives on the C++ side — this class just consumes the
  resulting per-edge density grids and the vine structure.

  Parameters
  ----------
  pair_copulas:
    Nested list (or :class:`torch.nn.ModuleList` of
    :class:`torch.nn.ModuleList`) of :class:`TorchBicop`, indexed
    ``[tree][edge]`` and shaped like ``Vinecop::pair_copulas_``: tree 0
    has ``d - 1`` edges, tree 1 has ``d - 2``, etc., up to ``trunc_lvl``.
  structure:
    The R-vine structure (a :class:`pyvinecopulib.RVineStructure`) whose
    accessors describe how to walk the trees. Held as-is so the inner
    loops query it directly.
  """

  d: int
  trunc_lvl: int
  pair_copulas: torch.nn.ModuleList

  def __init__(
    self,
    pair_copulas: list[list[TorchBicop]],
    structure,
  ) -> None:
    super().__init__()
    self.structure = structure
    self.d = int(structure.dim)
    self.trunc_lvl = int(structure.trunc_lvl)

    self.order = [int(x) for x in structure.order]  # 1-indexed
    self.inverse_order = [0] * self.d
    for j, k in enumerate(self.order):
      self.inverse_order[k - 1] = j

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

    # Lazy-impl metadata (built once, used by all `_*_lazy` methods).
    self._build_lazy_structure()

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
    cache_integrals: bool = False,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchVinecop":
    """Build a ``TorchVinecop`` from a fitted :class:`pyvinecopulib.Vinecop`.

    Each pair copula in ``cop.pair_copulas`` must be a ``tll`` or
    independence family; anything else raises. Discrete variables are
    not yet supported.
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
  def from_data(
    cls,
    u,
    structure,
    *,
    grid_size: int = 30,
    mult: float = 1.0,
    cache_integrals: bool = False,
    grid_type: str = "normal",
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchVinecop":
    """Fit a pure-PyTorch TLL vine on ``u`` given a fixed ``structure``.

    Tree-by-tree cascade mirroring ``Vinecop::select_families`` with a
    user-provided structure: at each ``(tree, edge)`` we collect the pair
    of pseudo-obs columns (``col0``, ``col1``) by the same rule as the
    legacy pdf/rosenblatt cascade, fit a :class:`TorchBicop` on them via
    :meth:`TorchBicop.from_data`, then propagate ``hfunc1`` / ``hfunc2``
    forward according to ``needed_hfunc1`` / ``needed_hfunc2``.

    The structure is *fixed* (not selected); this measures pair-copula
    estimation error only, matching ``pv.Vinecop.from_data(u,
    controls=FitControlsVinecop(family_set=[tll], structure=...))`` to
    machine precision (per-pair agreement is ~1e-11; the cascade
    preserves this when ``cache_integrals=False``).

    Args:
      u: ``(n, d)`` pseudo-observations; np.ndarray or Tensor.
      structure: :class:`pyvinecopulib.RVineStructure` describing the vine
        skeleton.
      grid_size, mult, grid_type: forwarded to :meth:`TorchBicop.from_data`.
      cache_integrals, device, dtype: forwarded to :meth:`TorchBicop.from_data`.
    """
    u_t = torch.as_tensor(u, dtype=dtype, device=device)
    if u_t.ndim != 2:
      raise ValueError(f"u must be 2-D; got shape {tuple(u_t.shape)}")
    n, d = u_t.shape
    if int(structure.dim) != d:
      raise ValueError(
        f"structure.dim={structure.dim} does not match u.shape[1]={d}"
      )

    order = [int(s) for s in structure.order]
    hfunc1 = torch.zeros(n, d, dtype=dtype, device=u_t.device)
    hfunc2 = torch.empty(n, d, dtype=dtype, device=u_t.device)
    for j in range(d):
      hfunc2[:, j] = u_t[:, order[j] - 1]

    trunc_lvl = int(structure.trunc_lvl)
    pair_copulas: list[list[TorchBicop]] = []
    s = structure
    for tree in range(trunc_lvl):
      tree_pairs: list[TorchBicop] = []
      for edge in range(d - tree - 1):
        m = int(s.min_array(tree, edge))
        sarr = int(s.struct_array(tree, edge, natural_order=True))
        col0 = hfunc2[:, edge]
        col1 = hfunc2[:, m - 1] if m == sarr else hfunc1[:, m - 1]
        u_e = torch.stack([col0, col1], dim=-1)
        bc = TorchBicop.from_data(
          u_e,
          grid_size=grid_size,
          mult=mult,
          cache_integrals=cache_integrals,
          grid_type=grid_type,
          device=u_t.device,
          dtype=dtype,
        )
        tree_pairs.append(bc)
        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = bc.hfunc1(u_e)
        if s.needed_hfunc2(tree, edge):
          hfunc2[:, edge] = bc.hfunc2(u_e)
      pair_copulas.append(tree_pairs)

    return cls(pair_copulas=pair_copulas, structure=structure)

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

  def _prep(self, u: Tensor, name: str) -> Tensor:
    ref = self._ref_tensor()
    u = torch.as_tensor(u, dtype=ref.dtype, device=ref.device)
    if u.ndim != 2 or u.shape[1] != self.d:
      raise ValueError(
        f"{name}: u must have shape (n, {self.d}); got {tuple(u.shape)}"
      )
    return u.clamp(_TRIM_LO, _TRIM_HI)

  # --------------------------------------------------------------------- #
  # Lazy-impl metadata                                                     #
  # --------------------------------------------------------------------- #

  def _build_lazy_structure(self) -> None:
    """Pre-compute the dict-based view of the tree cascade.

    Walks the C++ Vinecop's h-function cascade once symbolically,
    tracking which ``(v, *cond_ing)`` pseudo-obs sits in each
    natural-order column at each tree level. Produces:

    - ``_struct_obs``: per-level ``{(v, *cond_ing) -> bicop_name}``.
    - ``_struct_bcp``: ``{bicop_name -> {cond_ed, cond_ing, tree, edge, is_indep}}``.
    - ``_tree_bidep``: per-level ordered list of ``(v_l, v_r, *cond_ing)``.
    - ``_rosenblatt_keys``: the pseudo-obs key whose tensor value is
      the natural-order column ``j`` of the Rosenblatt output (matches
      the legacy impl's final ``hfunc2[:, j]``).
    """
    d = self.d
    s = self.structure

    hfunc1_id: list[tuple[int, ...]] = [(-1,)] * d
    hfunc2_id: list[tuple[int, ...]] = [(self.order[j] - 1,) for j in range(d)]

    self._struct_obs: list[dict[tuple[int, ...], str]] = [{} for _ in range(d)]
    for key in hfunc2_id:
      self._struct_obs[0][key] = ""

    self._struct_bcp: dict[str, dict] = {}
    self._tree_bidep: list[list[tuple[int, ...]]] = [
      [] for _ in range(max(0, d - 1))
    ]

    for tree in range(self.trunc_lvl):
      for edge in range(d - tree - 1):
        m = int(s.min_array(tree, edge))
        sarr = int(s.struct_array(tree, edge, natural_order=True))
        col0 = hfunc2_id[edge]
        col1 = hfunc2_id[m - 1] if m == sarr else hfunc1_id[m - 1]
        v_a, *s_a = col0
        v_b, *s_b = col1
        assert tuple(s_a) == tuple(s_b), (
          f"proximity condition violated at (tree={tree}, edge={edge}): "
          f"cond_ing left={s_a!r}, right={s_b!r}"
        )
        cond_ing = tuple(s_a)
        # Preserve the legacy impl's column order in the bicop name
        # and the `cond_ed` tuple: the C++ cascade calls each bicop
        # with col0 being the "left" variable (the one at natural-order
        # column `edge`) and col1 being the "right" variable (at column
        # m-1). The TLL density grid is asymmetric, so swapping these
        # would change the result.
        name = f"{v_a},{v_b}"

        self._tree_bidep[tree].append((v_a, v_b, *cond_ing))
        self._struct_bcp[name] = {
          "cond_ed": (v_a, v_b),
          "cond_ing": cond_ing,
          "tree": tree,
          "edge": edge,
          "is_indep": self._pair(tree, edge).is_indep,
        }

        # After this bicop:
        #   hfunc1[:, edge] -> pseudo-obs for v_b given (v_a, *cond_ing)
        #   hfunc2[:, edge] -> pseudo-obs for v_a given (v_b, *cond_ing)
        ci_for_va = tuple(sorted((*cond_ing, v_b)))
        ci_for_vb = tuple(sorted((*cond_ing, v_a)))
        hfunc1_id[edge] = (v_b, *ci_for_vb)
        hfunc2_id[edge] = (v_a, *ci_for_va)
        self._struct_obs[tree + 1][hfunc1_id[edge]] = name
        self._struct_obs[tree + 1][hfunc2_id[edge]] = name

    # The Rosenblatt output at natural-order column j is hfunc2[:, j] at
    # the end of the cascade. Below the truncation level the column is
    # untouched (the C++ leaves hfunc2[:, j] equal to the input column).
    self._rosenblatt_keys: list[tuple[int, ...]] = list(hfunc2_id)

  # --------------------------------------------------------------------- #
  # Public entry points (dispatch on `impl`)                               #
  # --------------------------------------------------------------------- #

  def pdf(
    self,
    u: Tensor,
    *,
    impl: str = "legacy",
    batched: bool = False,
  ) -> Tensor:
    """Bivariate-cascaded copula density at ``u``.

    See :meth:`_pdf_legacy` and :meth:`_pdf_lazy` for the two
    implementations; both produce numerically equivalent output.

    ``batched=True`` evaluates every pair-copula at one tree level in a
    single batched bicop call (one ``(N_t, m, m)`` interp per tree level
    instead of N_t scalar interps). Orthogonal to ``impl`` — the dense
    (legacy) and dict-based (lazy) storage strategies still apply; only
    the per-level inner loop is collapsed.
    """
    _check_impl(impl)
    fn = self._pick_pdf(impl, batched)
    return fn(u)

  def rosenblatt(
    self,
    u: Tensor,
    *,
    impl: str = "legacy",
    batched: bool = False,
  ) -> Tensor:
    """Rosenblatt transform: dependent uniforms ``u`` → independent ``w``.

    See :meth:`pdf` for the semantics of ``batched``.
    """
    _check_impl(impl)
    fn = self._pick_rosenblatt(impl, batched)
    return fn(u)

  @torch.no_grad()
  def inverse_rosenblatt(
    self,
    u: Tensor,
    *,
    impl: str = "legacy",
    batched: bool = False,
  ) -> Tensor:
    """Inverse Rosenblatt transform: independent uniforms → dependent.

    See :meth:`pdf` for the semantics of ``batched``.
    """
    _check_impl(impl)
    fn = self._pick_inverse_rosenblatt(impl, batched)
    return fn(u)

  # --------------------------------------------------------------------- #
  # Dispatch helpers                                                       #
  # --------------------------------------------------------------------- #

  def _pick_pdf(self, impl: str, batched: bool):
    if batched:
      return (
        self._pdf_batched_legacy if impl == "legacy" else self._pdf_batched_lazy
      )
    return self._pdf_legacy if impl == "legacy" else self._pdf_lazy

  def _pick_rosenblatt(self, impl: str, batched: bool):
    if batched:
      return (
        self._rosenblatt_batched_legacy
        if impl == "legacy"
        else self._rosenblatt_batched_lazy
      )
    return (
      self._rosenblatt_legacy if impl == "legacy" else self._rosenblatt_lazy
    )

  def _pick_inverse_rosenblatt(self, impl: str, batched: bool):
    if batched:
      return (
        self._inverse_rosenblatt_batched_legacy
        if impl == "legacy"
        else self._inverse_rosenblatt_batched_lazy
      )
    return (
      self._inverse_rosenblatt_legacy
      if impl == "legacy"
      else self._inverse_rosenblatt_lazy
    )

  def _apply(self, fn, *args, **kwargs):
    # `.to()`, `.cuda()`, `.cpu()` all route through `_apply`. The
    # BatchedVine container holds buffers — `super()._apply` would move
    # them, but we drop the whole structure so it gets re-baked from the
    # (already-moved) source pair_copulas on next use; that keeps the
    # wire-up tensors aligned with the destination dtype/device.
    self._batched = None
    return super()._apply(fn, *args, **kwargs)

  # --------------------------------------------------------------------- #
  # pdf — legacy impl                                                     #
  # --------------------------------------------------------------------- #

  def _pdf_legacy(self, u: Tensor) -> Tensor:
    """Direct port of ``Vinecop::pdf`` (class.ipp:916–1010)."""
    u = self._prep(u, "pdf")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d, trunc_lvl = self.d, self.trunc_lvl

    if trunc_lvl == 0:
      return torch.ones(n, dtype=dtype, device=device)

    hfunc1 = torch.zeros(n, d, dtype=dtype, device=device)
    hfunc2 = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      hfunc2[:, j] = u[:, self.order[j] - 1]

    log_pdf = torch.zeros(n, dtype=dtype, device=device)
    s = self.structure
    for tree in range(trunc_lvl):
      for edge in range(d - tree - 1):
        cop = self._pair(tree, edge)
        m = int(s.min_array(tree, edge))
        sarr = int(s.struct_array(tree, edge, natural_order=True))
        col0 = hfunc2[:, edge]
        col1 = hfunc2[:, m - 1] if m == sarr else hfunc1[:, m - 1]
        u_e = torch.stack([col0, col1], dim=-1)

        log_pdf = log_pdf + cop.log_pdf(u_e)

        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = cop.hfunc1(u_e)
        if s.needed_hfunc2(tree, edge):
          hfunc2[:, edge] = cop.hfunc2(u_e)

    return log_pdf.exp()

  # --------------------------------------------------------------------- #
  # rosenblatt — legacy impl                                              #
  # --------------------------------------------------------------------- #

  def _rosenblatt_legacy(self, u: Tensor) -> Tensor:
    """Direct port of ``Vinecop::rosenblatt`` (class.ipp:1527–1643)."""
    u = self._prep(u, "rosenblatt")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d, trunc_lvl = self.d, self.trunc_lvl

    hfunc2 = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      hfunc2[:, j] = u[:, self.order[j] - 1]
    hfunc1 = hfunc2.clone()

    s = self.structure
    for tree in range(trunc_lvl):
      for edge in range(d - tree - 1):
        cop = self._pair(tree, edge)
        m = int(s.min_array(tree, edge))
        sarr = int(s.struct_array(tree, edge, natural_order=True))
        col0 = hfunc2[:, edge]
        col1 = hfunc2[:, m - 1] if m == sarr else hfunc1[:, m - 1]
        u_e = torch.stack([col0, col1], dim=-1)

        if s.needed_hfunc1(tree, edge):
          hfunc1[:, edge] = cop.hfunc1(u_e)
        hfunc2[:, edge] = cop.hfunc2(u_e)

    U = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      U[:, j] = hfunc2[:, self.inverse_order[j]]
    return U.clamp(_TRIM_LO, _TRIM_HI)

  # --------------------------------------------------------------------- #
  # inverse_rosenblatt — legacy impl                                      #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def _inverse_rosenblatt_legacy(self, u: Tensor) -> Tensor:
    """Direct port of ``Vinecop::inverse_rosenblatt`` (class.ipp:1682–1770)."""
    u = self._prep(u, "inverse_rosenblatt")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d, trunc_lvl = self.d, self.trunc_lvl

    if trunc_lvl == 0:
      out = torch.empty(n, d, dtype=dtype, device=device)
      for j in range(d):
        out[:, j] = u[:, self.order[self.inverse_order[j]] - 1]
      return out

    hinv2 = torch.empty(trunc_lvl + 1, d, n, dtype=dtype, device=device)
    hfunc1 = torch.empty_like(hinv2)
    for j in range(d):
      hinv2[min(trunc_lvl, d - j - 1), j, :] = u[:, self.order[j] - 1]
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
        U_e = torch.stack([col0, col1], dim=-1)
        hinv2[tree, var, :] = cop.hinv2(U_e)
        if var < d - 1 and s.needed_hfunc1(tree, var):
          U_e_after = torch.stack([hinv2[tree, var, :], col1], dim=-1)
          hfunc1[tree + 1, var, :] = cop.hfunc1(U_e_after)

    U_vine = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      U_vine[:, j] = hinv2[0, self.inverse_order[j], :]
    return U_vine.clamp(_TRIM_LO, _TRIM_HI)

  # --------------------------------------------------------------------- #
  # pdf — lazy impl                                                        #
  # --------------------------------------------------------------------- #

  def _pdf_lazy(self, u: Tensor) -> Tensor:
    """Dict-based dual of :meth:`_pdf_legacy`.

    Materializes each ``(v, *cond_ing)`` pseudo-obs on first access via
    a recursive ``_ensure`` helper, accumulates ``log_pdf`` over edges
    of each tree level, and clears each level's dict after consumption.
    """
    u = self._prep(u, "pdf")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d = self.d

    if self.trunc_lvl == 0:
      return torch.ones(n, dtype=dtype, device=device)

    dct_obs: list[dict[tuple[int, ...], Tensor]] = [{} for _ in range(d)]
    for v in range(d):
      dct_obs[0][(v,)] = u[:, v]

    def _ensure(lv: int, key: tuple[int, ...]) -> None:
      if key in dct_obs[lv]:
        return
      name = self._struct_obs[lv][key]
      info = self._struct_bcp[name]
      cond_ing = info["cond_ing"]
      # cond_ed is in cascade column order: (col0_var, col1_var).
      v_a, v_b = info["cond_ed"]
      key_a = (v_a, *cond_ing)
      key_b = (v_b, *cond_ing)
      for kk in (key_a, key_b):
        if kk not in dct_obs[lv - 1]:
          _ensure(lv - 1, kk)
      if info["is_indep"]:
        same_key = key_a if key[0] == v_a else key_b
        dct_obs[lv][key] = dct_obs[lv - 1][same_key]
        return
      cop = self._pair(info["tree"], info["edge"])
      u_e = torch.stack(
        [dct_obs[lv - 1][key_a], dct_obs[lv - 1][key_b]], dim=-1
      )
      # key[0] == v_a: pseudo-obs for the left variable v_a given v_b + s
      # -> hfunc2 (returns U1 conditional on U2).
      # key[0] == v_b: pseudo-obs for v_b given v_a + s -> hfunc1.
      dct_obs[lv][key] = cop.hfunc2(u_e) if key[0] == v_a else cop.hfunc1(u_e)

    log_pdf = torch.zeros(n, dtype=dtype, device=device)
    for lv in range(self.trunc_lvl):
      for v_a, v_b, *cond_ing_list in self._tree_bidep[lv]:
        cond_ing = tuple(cond_ing_list)
        key_a = (v_a, *cond_ing)
        key_b = (v_b, *cond_ing)
        for kk in (key_a, key_b):
          if kk not in dct_obs[lv]:
            _ensure(lv, kk)
        info = self._struct_bcp[f"{v_a},{v_b}"]
        cop = self._pair(info["tree"], info["edge"])
        u_e = torch.stack([dct_obs[lv][key_a], dct_obs[lv][key_b]], dim=-1)
        log_pdf = log_pdf + cop.log_pdf(u_e)
      if lv > 0:
        dct_obs[lv - 1].clear()

    return log_pdf.exp()

  # --------------------------------------------------------------------- #
  # rosenblatt — lazy impl                                                 #
  # --------------------------------------------------------------------- #

  def _rosenblatt_lazy(self, u: Tensor) -> Tensor:
    """Dict-based dual of :meth:`_rosenblatt_legacy`.

    Materializes every key in ``self._rosenblatt_keys`` (one per
    natural-order column) via the same lazy walk as :meth:`_pdf_lazy`,
    then reorders to original variable indices via ``inverse_order``.
    """
    u = self._prep(u, "rosenblatt")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d = self.d

    if self.trunc_lvl == 0:
      return u.clone()

    dct_obs: list[dict[tuple[int, ...], Tensor]] = [{} for _ in range(d)]
    for v in range(d):
      dct_obs[0][(v,)] = u[:, v]

    # `key_level[key]` is the level at which `key` lives; used by the
    # _ensure helper to avoid scanning all levels per access. We can
    # precompute it from `_struct_obs` once per call.
    key_level: dict[tuple[int, ...], int] = {}
    for lv in range(d):
      for k in self._struct_obs[lv]:
        key_level[k] = lv

    def _ensure(key: tuple[int, ...]) -> None:
      lv = key_level[key]
      if key in dct_obs[lv]:
        return
      name = self._struct_obs[lv][key]
      info = self._struct_bcp[name]
      cond_ing = info["cond_ing"]
      v_a, v_b = info["cond_ed"]
      key_a = (v_a, *cond_ing)
      key_b = (v_b, *cond_ing)
      for kk in (key_a, key_b):
        _ensure(kk)
      if info["is_indep"]:
        same_key = key_a if key[0] == v_a else key_b
        dct_obs[lv][key] = dct_obs[key_level[same_key]][same_key]
        return
      cop = self._pair(info["tree"], info["edge"])
      u_e = torch.stack(
        [
          dct_obs[key_level[key_a]][key_a],
          dct_obs[key_level[key_b]][key_b],
        ],
        dim=-1,
      )
      dct_obs[lv][key] = cop.hfunc2(u_e) if key[0] == v_a else cop.hfunc1(u_e)

    # Materialize every Rosenblatt-output key.
    for key in self._rosenblatt_keys:
      _ensure(key)

    U = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      key = self._rosenblatt_keys[self.inverse_order[j]]
      U[:, j] = dct_obs[key_level[key]][key]
    return U.clamp(_TRIM_LO, _TRIM_HI)

  # --------------------------------------------------------------------- #
  # inverse_rosenblatt — lazy impl                                         #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def _inverse_rosenblatt_lazy(self, u: Tensor) -> Tensor:
    """Reference-counted ``hinv`` walk, dual of :meth:`_inverse_rosenblatt_legacy`.

    Source pseudo-obs (the deepest-cond_ing pseudo-obs along each
    natural-order column) are initialized from ``u``; each is then
    climbed upward via ``hinv1`` / ``hinv2`` calls until its
    conditioning set is empty, freeing intermediate values via ref
    counting once they're no longer needed.
    """
    u = self._prep(u, "inverse_rosenblatt")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d, trunc_lvl = self.d, self.trunc_lvl

    if trunc_lvl == 0:
      out = torch.empty(n, d, dtype=dtype, device=device)
      for j in range(d):
        out[:, j] = u[:, self.order[self.inverse_order[j]] - 1]
      return out

    # The C++ cascade puts u[:, order[j]-1] into hinv2[min(trunc_lvl, d-j-1), j],
    # which corresponds (in the lazy dict-keyed view) to the
    # natural-order-column-j pseudo-obs at the *deepest* level reached
    # by that column. Find that key:
    #   - for j < d - 1: at level min(trunc_lvl, d - j - 1) the
    #     natural-order-column j has been touched (since the cascade
    #     visits column j up to tree d - j - 2 inclusive). The key is
    #     the one we recorded as `hfunc2_id[j]` after running the
    #     symbolic cascade in `_build_lazy_structure`. That's exactly
    #     `_rosenblatt_keys[j]`.
    #   - for j == d - 1: never touched in the cascade — stays
    #     `(self.order[d-1] - 1,)`.
    source_keys: list[tuple[int, ...]] = list(self._rosenblatt_keys)
    # _rosenblatt_keys[d - 1] is already `(self.order[d-1] - 1,)` since
    # column d-1 is never updated, so no special case needed.

    ref_count, source_list = self._ref_count_hinv(source_keys)
    dct_obs: dict[tuple[int, ...], Tensor] = {}
    for j, key in enumerate(source_keys):
      dct_obs[key] = u[:, self.order[j] - 1]

    def _decrement(key: tuple[int, ...]) -> None:
      ref_count[key] -= 1
      if ref_count[key] <= 0 and len(key) > 1:
        # Only pseudo-obs with non-empty cond_ing are GC'd; the
        # ``(v,)`` outputs are kept around for the final stack.
        dct_obs.pop(key, None)

    def _visit(key: tuple[int, ...]) -> tuple[int, ...]:
      """Climb ``key`` one tree level toward shallower cond_ing.

      Returns the parent key ``(v_down, *cond_ing_up)``. Recursively
      fills any missing same-level "sibling" pseudo-obs via downward
      hfunc evaluations.
      """
      v_down, *cond_ing_down_list = key
      lv = len(cond_ing_down_list)
      name = self._struct_obs[lv][key]
      info = self._struct_bcp[name]
      # cond_ed is in cascade column order: (col0_var, col1_var).
      v_a, v_b = info["cond_ed"]
      cond_ing_up = info["cond_ing"]
      is_down_b = v_down == v_b
      sibling_key = (v_a, *cond_ing_up) if is_down_b else (v_b, *cond_ing_up)
      if sibling_key not in dct_obs:
        _materialize_down(sibling_key)

      key_up = (v_down, *cond_ing_up)
      if info["is_indep"]:
        dct_obs[key_up] = dct_obs[key]
      else:
        cop = self._pair(info["tree"], info["edge"])
        p = dct_obs[key]
        if is_down_b:
          # v_down corresponds to v_b (col 1). The forward map is
          #   p = hfunc1((u_a, u_b)) = H(u_b | u_a) [returns u_b col].
          # Invert with hinv1((u_a, p)) → returns u_b at cond_ing_up.
          u_a_val = dct_obs[sibling_key]
          inp = torch.stack([u_a_val, p], dim=-1)
          dct_obs[key_up] = cop.hinv1(inp)
        else:
          # v_down corresponds to v_a (col 0). Forward map is
          #   p = hfunc2((u_a, u_b)) = H(u_a | u_b) [returns u_a col].
          # Invert with hinv2((p, u_b)) → returns u_a at cond_ing_up.
          u_b_val = dct_obs[sibling_key]
          inp = torch.stack([p, u_b_val], dim=-1)
          dct_obs[key_up] = cop.hinv2(inp)

      # GC: every step consumes three pseudo-obs refs (matching the
      # increment pattern in `_ref_count_hinv`): the two top-level
      # cousins (one of which is the produced `key_up`) and `key`
      # itself. The next iteration will re-increment `key_up`'s ref
      # at its consumption site, so this decrement is safe.
      for k_done in (key, sibling_key, key_up):
        _decrement(k_done)
      return key_up

    def _materialize_down(key: tuple[int, ...]) -> None:
      """Compute ``key`` by descending hfunc calls.

      Used when the sibling needed at a `_visit` call hasn't been
      produced yet (an upstream cond_ing pseudo-obs that wasn't on the
      direct path from a source).
      """
      if key in dct_obs:
        return
      lv = len(key) - 1
      assert lv > 0, f"missing source key {key!r}"
      name = self._struct_obs[lv][key]
      info = self._struct_bcp[name]
      v_a, v_b = info["cond_ed"]
      cond_ing_up = info["cond_ing"]
      key_a = (v_a, *cond_ing_up)
      key_b = (v_b, *cond_ing_up)
      for kk in (key_a, key_b):
        if kk not in dct_obs:
          _materialize_down(kk)
      if info["is_indep"]:
        same_key = key_a if key[0] == v_a else key_b
        dct_obs[key] = dct_obs[same_key]
      else:
        cop = self._pair(info["tree"], info["edge"])
        u_e = torch.stack([dct_obs[key_a], dct_obs[key_b]], dim=-1)
        dct_obs[key] = cop.hfunc2(u_e) if key[0] == v_a else cop.hfunc1(u_e)

    # Climb each source to its empty-cond_ing root.
    for src in source_list:
      cur = src
      while len(cur) > 1:
        nxt = _visit(cur)
        cur = nxt

    U_vine = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      U_vine[:, j] = dct_obs[(j,)]
    return U_vine.clamp(_TRIM_LO, _TRIM_HI)

  # --------------------------------------------------------------------- #
  # Reference-counting walk for `_inverse_rosenblatt_lazy`                 #
  # --------------------------------------------------------------------- #

  def _ref_count_hinv(
    self, source_keys: list[tuple[int, ...]]
  ) -> tuple[Counter, list[tuple[int, ...]]]:
    """Count consumers of each pseudo-obs along the upward hinv walk.

    Mirrors ``torchvinecopulib.VineCop.ref_count_hfunc`` but specialized
    to the natural-order traversal. Returns a ``Counter`` keyed on
    pseudo-obs keys plus the ordered ``source_list`` (sorted from
    shallowest to deepest cond_ing — the order in which we climb).
    """
    ref_count: Counter = Counter()

    def _visit_count(key: tuple[int, ...], is_hinv: bool) -> tuple[int, ...]:
      if len(key) == 1:
        ref_count[key] += 1
        return key
      v_down, *cond_ing_down_list = key
      lv = len(cond_ing_down_list)
      name = self._struct_obs[lv][key]
      info = self._struct_bcp[name]
      v_a, v_b = info["cond_ed"]
      cond_ing_up = info["cond_ing"]
      is_down_b = v_down == v_b
      if is_hinv:
        sibling_key = (v_a, *cond_ing_up) if is_down_b else (v_b, *cond_ing_up)
        frontier = [sibling_key]
      else:
        frontier = [(v_a, *cond_ing_up), (v_b, *cond_ing_up)]
      for fk in frontier:
        if ref_count[fk] == 0:
          _visit_count(fk, is_hinv=False)
      # Increment refs on the three keys touched at this step.
      key_up = (v_down, *cond_ing_up)
      for fk in [(v_a, *cond_ing_up), (v_b, *cond_ing_up), key]:
        ref_count[fk] += 1
      if is_hinv:
        return key_up
      return key

    # Sort sources by ascending cond_ing length so shallower sources
    # are climbed first (matches the alt-design `lst_source` order).
    source_list = sorted(source_keys, key=lambda k: len(k))
    for src in source_list:
      cur = src
      if len(cur) == 1:
        ref_count[cur] += 1
        continue
      while len(cur) > 1:
        cur = _visit_count(cur, is_hinv=True)
    return ref_count, source_list

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

  def _pdf_batched_legacy(self, u: Tensor) -> Tensor:
    """Batched pdf with dense ``(n, d)`` scratch (legacy storage).

    Mirrors :meth:`_pdf_legacy` exactly except the inner per-edge loop is
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

  def _pdf_batched_lazy(self, u: Tensor) -> Tensor:
    """Same as :meth:`_pdf_batched_legacy`: once a tree level fires one
    batched bicop call, the per-edge dict-vs-dense distinction collapses
    (both impls store one ``(n, d)`` scratch pair across the cascade).
    Kept under a separate name so callers can parametrize on ``impl``
    consistently with the non-batched paths.
    """
    return self._pdf_batched_legacy(u)

  # ---- rosenblatt ------------------------------------------------------ #

  def _rosenblatt_batched_legacy(self, u: Tensor) -> Tensor:
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
      # legacy cascade (`hfunc2[:, edge] = cop.hfunc2(u_e)`); hfunc1 is
      # gated by needs_h1.
      hfunc2[:, :N_t] = h2_new
      hfunc1[:, :N_t] = torch.where(
        lvl.needs_h1[None, :], h1_new, hfunc1[:, :N_t]
      )

    out = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      out[:, j] = hfunc2[:, self.inverse_order[j]]
    return out.clamp(_TRIM_LO, _TRIM_HI)

  def _rosenblatt_batched_lazy(self, u: Tensor) -> Tensor:
    # Same rationale as ``_pdf_batched_lazy``: batched-per-level fan-out
    # makes dict-vs-dense moot.
    return self._rosenblatt_batched_legacy(u)

  # ---- inverse_rosenblatt --------------------------------------------- #

  @torch.no_grad()
  def _inverse_rosenblatt_batched_legacy(self, u: Tensor) -> Tensor:
    """``batched=True`` fallback for inverse_rosenblatt.

    The inverse cascade's dependency graph isn't tree-level batchable:
    iteration ``(var, tree)`` reads ``hfunc1[tree, m-1, :]`` which was
    written at ``(m-1, tree-1)`` — a different tree level. A tree-outer
    reordering therefore needs a full DAG topological sort across the
    ``(var, tree)`` lattice, not just within-tree waves. That's a v2
    rewrite; for v1 we route ``batched=True`` to the per-pair legacy
    cascade (still bit-equal to C++ via :func:`solve_bisection`).
    """
    return self._inverse_rosenblatt_legacy(u)

  @torch.no_grad()
  def _inverse_rosenblatt_batched_lazy(self, u: Tensor) -> Tensor:
    """Same fallback as :meth:`_inverse_rosenblatt_batched_legacy`."""
    return self._inverse_rosenblatt_lazy(u)
