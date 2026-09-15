"""The data layouts a vine's discrete variables imply.

A discrete variable is described by two numbers, ``F(x)`` and its left limit
``F(x^-)``, so a vine with ``k`` of them reads an ``(n, d + k)`` argument rather
than an ``(n, d)`` one. This module owns that bookkeeping: normalizing a caller's
layout to the compact one, deriving the types each pair-copula slot sees from
the structure alone, and gathering an edge's columns out of the cascade's
h-function buffers.

The mixed-discrete evaluation those layouts feed is
:class:`~pyvinecopulib.core.DiscreteBicop`, in ``core/bicop_discrete.py``.

Internal: the vine cascades and the fit engines call these helpers.
"""

from __future__ import annotations

from types import ModuleType
from typing import Any, Optional, cast

from array_api_compat import array_namespace

from ..pyvinecopulib_ext import RVineStructure
from .protocols import ArrayT

__all__ = ["collapse_data"]


def check_var_types(var_types: Optional[list[str]], d: int) -> tuple[str, ...]:
  """Normalize and validate a vine's per-variable types.

  Parameters
  ----------
  var_types : list of str, or None, optional
      Per-variable types, ``"c"`` or ``"d"``; ``None`` means all continuous.
  d : int
      Dimension the types must cover.

  Returns
  -------
  tuple of str
      The validated types, one per variable.

  Raises
  ------
  ValueError
      If the length is not ``d`` or an entry is outside ``{"c", "d"}``.
  """
  types = ("c",) * d if var_types is None else tuple(var_types)
  if len(types) != d:
    raise ValueError(f"var_types has {len(types)} entries, expected {d}")
  bad = [t for t in types if t not in ("c", "d")]
  if bad:
    raise ValueError(f"var_types entries must be 'c' or 'd'; got {bad[0]!r}")
  return types


def disc_cols(var_types: tuple[str, ...]) -> tuple[int, ...]:
  """Offsets of the left-limit columns within the compact layout's second block.

  Parameters
  ----------
  var_types : tuple of str
      Per-variable types.

  Returns
  -------
  tuple of int
      Variable ``i``'s left limit sits at column ``d + result[i]``; the entry is
      meaningless (``0``) at a continuous variable, which has no such column.
  """
  offsets, seen = [0] * len(var_types), 0
  for i, t in enumerate(var_types):
    if t == "d":
      offsets[i] = seen
      seen += 1
  return tuple(offsets)


def collapse_data(
  u: ArrayT,
  d: int,
  var_types: tuple[str, ...],
  name: str,
  *,
  values_only: bool = False,
) -> ArrayT:
  """Validate ``u``'s column layout and reduce it to the columns needed.

  Accepts the layouts ``Vinecop`` accepts: ``(n, d)`` for an all-continuous
  model, and the compact ``(n, d + k)`` or expanded ``(n, 2d)`` form when ``k``
  of the ``d`` variables are discrete. The plain ``(n, d)`` is rejected on a
  discrete model, because silently reusing each value as its own left limit
  would evaluate a continuous density.

  Parameters
  ----------
  u : array, shape (n, d), (n, d + k) or (n, 2d), dtype float
      The matrix to validate.
  d : int
      Dimension of the model.
  var_types : tuple of str
      Per-variable types.
  name : str
      Calling-method name, used only in the error message.
  values_only : bool, default=False
      Return just the ``d`` value columns, and accept a plain ``(n, d)`` input
      even on a discrete model. Set by the callers that never read a left limit.

  Returns
  -------
  array, shape (n, d + k) or (n, d), dtype float
      ``u`` in the compact layout, or its value block when ``values_only``.

  Raises
  ------
  ValueError
      If ``u`` is not 2-d or its column count matches no accepted layout.
  """
  ua: Any = u
  k = var_types.count("d")
  accepted = {d + k, 2 * d} | ({d} if values_only else set())
  if ua.ndim != 2 or int(ua.shape[1]) not in accepted:
    shapes = ", ".join(f"(n, {c})" for c in sorted(accepted))
    raise ValueError(
      f"{name}: u must have shape {shapes} for a vine with var_types="
      f"{list(var_types)}; got {tuple(ua.shape)}"
    )
  if values_only or k == 0:
    return cast("ArrayT", ua[:, :d])
  if int(ua.shape[1]) == d + k:
    return u
  # Expanded (n, 2d) -> compact (n, d + k): keep the left-limit columns of the
  # discrete variables only, in variable order.
  xp = array_namespace(ua)
  keep = [d + i for i, t in enumerate(var_types) if t == "d"]
  return cast("ArrayT", xp.concat([ua[:, :d], ua[:, keep]], axis=1))


def pair_var_types(
  structure: RVineStructure, var_types: tuple[str, ...]
) -> tuple[tuple[tuple[str, str], ...], ...]:
  """Per-edge variable types implied by a structure and its variable types.

  Port of ``Vinecop::set_var_types_internal``: tree 0 reads the natural-order
  variable types, and every later tree inherits from the pair that produced its
  input column -- ``hfunc2`` carries the first variable's type, ``hfunc1`` the
  second's. The types are therefore a function of the structure and never have
  to be stored on a pair copula.

  Parameters
  ----------
  structure : RVineStructure
      The vine structure.
  var_types : tuple of str
      Per-variable types, in variable order.

  Returns
  -------
  tuple of tuple of tuple of str
      ``result[tree][edge]`` is that edge's ``(type1, type2)``.
  """
  d = int(structure.dim)
  trunc_lvl = int(structure.trunc_lvl)
  order = [int(v) for v in structure.order]
  natural = [var_types[v - 1] for v in order]
  table: list[tuple[tuple[str, str], ...]] = []
  for tree in range(trunc_lvl):
    row: list[tuple[str, str]] = []
    for edge in range(d - tree - 1):
      m = int(structure.min_array(tree, edge))
      on_diagonal = m == int(structure.struct_array(tree, edge, True))
      if tree == 0:
        row.append((natural[edge], natural[m - 1]))
      else:
        prev = table[tree - 1]
        row.append((prev[edge][0], prev[m - 1][0 if on_diagonal else 1]))
    table.append(tuple(row))
  return tuple(table)


def seed_left_limits(
  u: ArrayT,
  d: int,
  order: tuple[int, ...],
  var_types: tuple[str, ...],
  offsets: tuple[int, ...],
  xp: ModuleType,
) -> Optional[ArrayT]:
  """Natural-order left limits read off the compact layout, or ``None``.

  ``None`` for an all-continuous model, which is what switches the whole
  left-limit cascade off. A continuous variable's column holds its own value: a
  pair only ever reads the left-limit column of a variable it declares discrete,
  and this keeps the four-column edge input well defined anyway.

  Parameters
  ----------
  u : array, shape (n, d + k), dtype float
      Prepared observations in the compact layout.
  d : int
      Dimension of the model.
  order : tuple of int
      The structure's variable order (1-based).
  var_types : tuple of str
      Per-variable types.
  offsets : tuple of int
      Left-limit column offsets, from :func:`disc_cols`.
  xp : module
      The array namespace of ``u``.

  Returns
  -------
  array, shape (n, d), dtype float, or None
      Left limits in natural order, or ``None`` when nothing is discrete.
  """
  if "d" not in var_types:
    return None
  ua: Any = u
  sub = xp.empty((ua.shape[0], d), dtype=ua.dtype, device=ua.device)
  for j in range(d):
    v = order[j] - 1
    sub[:, j] = ua[:, d + offsets[v] if var_types[v] == "d" else v]
  return cast("ArrayT", sub)


def edge_columns(
  structure: RVineStructure,
  pair_types: Optional[tuple[tuple[tuple[str, str], ...], ...]],
  tree: int,
  edge: int,
  hfunc1: ArrayT,
  hfunc2: ArrayT,
  hfunc1_sub: Optional[ArrayT],
  hfunc2_sub: Optional[ArrayT],
) -> tuple[ArrayT, ArrayT, Optional[tuple[ArrayT, ArrayT]], tuple[str, str]]:
  """Resolve one edge's pair-copula input columns and its variable types.

  ``m`` is the min-array entry: the natural-order index of the column finalized
  in a previous tree. The second pair input comes from ``hfunc2`` when ``m`` sits
  on the natural-order diagonal, else from ``hfunc1``
  (``class.ipp:1026-1034``). The left-limit pair is returned only when the edge
  has a discrete variable, and mirrors ``Bicop::format_data``: a continuous
  variable's left limit is its own value.

  Parameters
  ----------
  structure : RVineStructure
      The vine structure being walked.
  pair_types : tuple of tuple of tuple of str, or None, optional
      Per-edge types from :func:`pair_var_types`; ``None`` when all continuous.
  tree : int
      Tree index (``0``-based).
  edge : int
      Edge index within the tree (``0``-based).
  hfunc1, hfunc2 : array, shape (n, d), dtype float
      The h-function scratch matrices.
  hfunc1_sub, hfunc2_sub : array, shape (n, d), dtype float, or None, optional
      The left-limit scratch matrices; ``None`` when all continuous.

  Returns
  -------
  col0, col1 : array, shape (n,), dtype float
      The pair's two value inputs.
  subs : tuple of array, or None
      Their left limits, or ``None`` when the edge is fully continuous.
  types : tuple of str
      The edge's ``(type1, type2)``.
  """
  h1: Any = hfunc1
  h2: Any = hfunc2
  h2_sub: Any = hfunc2_sub
  m = int(structure.min_array(tree, edge))
  on_diagonal = m == int(structure.struct_array(tree, edge, True))
  col0 = cast("ArrayT", h2[:, edge])
  col1 = cast("ArrayT", h2[:, m - 1] if on_diagonal else h1[:, m - 1])
  types = ("c", "c") if pair_types is None else pair_types[tree][edge]
  if hfunc2_sub is None or "d" not in types:
    return col0, col1, None, types
  sub0 = cast("ArrayT", h2_sub[:, edge]) if types[0] == "d" else col0
  if types[1] != "d":
    sub1 = col1
  elif on_diagonal:
    sub1 = cast("ArrayT", h2_sub[:, m - 1])
  else:
    sub1 = cast("ArrayT", cast("Any", hfunc1_sub)[:, m - 1])
  return col0, col1, (sub0, sub1), types


def stack_edge(
  xp: ModuleType,
  col0: ArrayT,
  col1: ArrayT,
  subs: Optional[tuple[ArrayT, ArrayT]],
) -> ArrayT:
  """Assemble a pair-copula argument from its value columns and left limits.

  Parameters
  ----------
  xp : module
      The array namespace to build on.
  col0, col1 : array, shape (n,), dtype float
      The pair's two value inputs.
  subs : tuple of array, or None, optional
      Their left limits, or ``None`` for a fully continuous edge.

  Returns
  -------
  array, shape (n, 2) or (n, 4), dtype float
      ``[u1, u2]``, or ``[u1, u2, u1^-, u2^-]`` when left limits are given.
  """
  cols = [col0, col1] if subs is None else [col0, col1, *subs]
  return cast("ArrayT", xp.stack(cols, axis=-1))


def with_left_limit(u_e: ArrayT, arg: int) -> ArrayT:
  """A four-column edge input with one argument replaced by its left limit.

  The input the cascade needs for a left-limit h-function: conditioning on (or
  evaluating at) the lower end of a discrete argument's atom.

  Parameters
  ----------
  u_e : array, shape (n, 4), dtype float
      An edge input ``[u1, u2, u1^-, u2^-]``.
  arg : int
      Which argument to replace, ``0`` or ``1``.

  Returns
  -------
  array, shape (n, 4), dtype float
      ``u_e`` with column ``arg`` replaced by column ``2 + arg``.
  """
  edge: Any = u_e
  xp = array_namespace(edge)
  cols = [edge[:, 0], edge[:, 1], edge[:, 2], edge[:, 3]]
  cols[arg] = edge[:, 2 + arg]
  return cast("ArrayT", xp.stack(cols, axis=-1))
