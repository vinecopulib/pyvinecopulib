"""PyTorch ``TorchVinecop`` — port of the multivariate vinecop evaluation chain.

Wraps a fitted :class:`pyvinecopulib.Vinecop` (with TLL pair copulas) and
exposes :meth:`pdf`, :meth:`rosenblatt`, and :meth:`inverse_rosenblatt`
on top of :class:`TorchBicop` for every pair copula. The whole evaluation
chain stays in PyTorch, so the vine can move to GPU and be composed with
autograd-aware downstream code.

Fitting is delegated to the C++ library; use
:meth:`TorchVinecop.from_vinecop` after fitting with
``pv.Vinecop.from_data(..., controls=FitControlsVinecop(family_set=[pv.tll]))``.
"""

from __future__ import annotations

from typing import Optional, cast

import torch
from torch import Tensor

from ..pyvinecopulib_ext import (
  Vinecop,
  indep as _INDEP_FAMILY,
  tll as _TLL_FAMILY,
)
from ._interp import _TRIM_HI, _TRIM_LO
from .bicop import TorchBicop


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
  # pdf                                                                    #
  # --------------------------------------------------------------------- #

  def pdf(self, u: Tensor) -> Tensor:
    """Bivariate-cascaded copula density at ``u``.

    Direct port of ``Vinecop::pdf`` (class.ipp:916–1010), continuous-only
    and single-batch.
    """
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
  # rosenblatt                                                             #
  # --------------------------------------------------------------------- #

  def rosenblatt(self, u: Tensor) -> Tensor:
    """Rosenblatt transform: dependent uniforms ``u`` → independent ``w``.

    Direct port of ``Vinecop::rosenblatt`` (class.ipp:1527–1643),
    continuous-only and single-batch (no discrete-randomization branch).
    """
    u = self._prep(u, "rosenblatt")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d, trunc_lvl = self.d, self.trunc_lvl

    hfunc2 = torch.empty(n, d, dtype=dtype, device=device)
    for j in range(d):
      hfunc2[:, j] = u[:, self.order[j] - 1]
    hfunc1 = hfunc2.clone()  # mirrors C++ "hfunc1 = hfunc2" init

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
  # inverse_rosenblatt                                                     #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def inverse_rosenblatt(self, u: Tensor) -> Tensor:
    """Inverse Rosenblatt transform: independent uniforms → dependent.

    Direct port of ``Vinecop::inverse_rosenblatt`` (class.ipp:1682–1770),
    continuous-only and single-batch (no memory-split heuristic).
    """
    u = self._prep(u, "inverse_rosenblatt")
    n = u.shape[0]
    ref = self._ref_tensor()
    dtype, device = ref.dtype, ref.device
    d, trunc_lvl = self.d, self.trunc_lvl

    if trunc_lvl == 0:
      # Independent vine: just permute by the structure order then back.
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
