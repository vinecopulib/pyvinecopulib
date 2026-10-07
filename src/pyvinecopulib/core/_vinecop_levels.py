"""A vine's tree levels and inverse waves, wired for one stacked call each.

The forward cascades walk a vine one tree level at a time, and the inverse one
dependency wave at a time. This module owns what that needs from the structure
alone -- which scratch column feeds each edge, which h-functions the next tree
reads, and how the inverse's cells group into waves -- and the default way to
evaluate such a group: one pair at a time, through each pair's own dispatchers.

A pair class that can evaluate a group in one call supplies ``_stack_pairs``,
returning an object with the same three methods as :class:`LoopedPairs`;
``VinecopBase._build_batched`` asks for it.

Internal: the cascades in ``vinecop_base`` call these helpers.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from ..pyvinecopulib_ext import RVineStructure
from ._covariates import pair_eval
from .protocols import ArrayT, BicopLike, Namespace

#: The per-edge masks of one tree level, kept as Python lists: a pair class's
#: stack and the loop over pairs read them, never an array op.
_LEVEL_MASKS = ("needs_h1", "needs_h2", "disc1", "disc2")


def _where(flags: list[bool]) -> list[int]:
  return [k for k, f in enumerate(flags) if f]


def level_wiring(
  structure: RVineStructure,
  pair_types: tuple[tuple[tuple[str, str], ...], ...],
  tree: int,
) -> dict[str, list[Any]]:
  """One tree level's per-edge tables, derived from the structure alone.

  ``m`` is the min-array entry: the natural-order index of the column finalized
  in a previous tree. An edge's second input comes from ``hfunc2`` when ``m``
  sits on the natural-order diagonal, else from ``hfunc1``
  (``class.ipp:1026-1034``).

  Parameters
  ----------
  structure : RVineStructure
      The vine structure.
  pair_types : tuple of tuple of tuple of str
      Per-edge variable types, ``pair_types[tree][edge]``.
  tree : int
      Tree index (``0``-based).

  Returns
  -------
  dict
      Per edge, its masks (Python lists of bool); and as positions, its input
      columns (``col0_src`` / ``col1_src``, with ``h1_pos`` / ``h1_src`` the
      edges whose second input reads ``hfunc1``), its continuous arguments
      (``cont1_pos`` / ``cont2_pos``) and the outputs a later tree reads
      (``w_*``).
  """
  n_pairs = int(structure.dim) - tree - 1
  col1_src: list[int] = []
  reads_h1: list[bool] = []
  out: dict[str, list[Any]] = {k: [] for k in _LEVEL_MASKS}
  for edge in range(n_pairs):
    m = int(structure.min_array(tree, edge))
    types = pair_types[tree][edge]
    col1_src.append(m - 1)
    reads_h1.append(m != int(structure.struct_array(tree, edge, True)))
    out["needs_h1"].append(bool(structure.needed_hfunc1(tree, edge)))
    out["needs_h2"].append(bool(structure.needed_hfunc2(tree, edge)))
    out["disc1"].append(types[0] == "d")
    out["disc2"].append(types[1] == "d")
  h1_pos = _where(reads_h1)
  # The index tables, as positions so a write or a gather is plain indexing.
  out.update(
    col0_src=list(range(n_pairs)),
    col1_src=col1_src,
    h1_pos=h1_pos,
    h1_src=[col1_src[k] for k in h1_pos],
    cont1_pos=_where([not f for f in out["disc1"]]),
    cont2_pos=_where([not f for f in out["disc2"]]),
    w_h1=_where(out["needs_h1"]),
    w_h2=_where(out["needs_h2"]),
    w_h1s=_where(
      [a and b for a, b in zip(out["needs_h1"], out["disc2"], strict=True)]
    ),
    w_h2s=_where(
      [a and b for a, b in zip(out["needs_h2"], out["disc1"], strict=True)]
    ),
    w_h2s_all=_where(out["disc1"]),
  )
  return out


def inverse_waves(
  structure: RVineStructure,
  trunc_lvl: int,
  cond_positions: Callable[[int, int], Sequence[int]] | None = None,
) -> list[list[tuple[int, int]]]:
  """Group the inverse cascade's ``(var, tree)`` cells into parallel waves.

  The inverse walks variables outward and, within each, trees inward. Cell
  ``(var, tree)`` reads ``hinv2[tree + 1, var]`` -- the same variable one tree
  further out -- and, at the same tree, either ``hinv2[tree, m - 1]`` or
  ``hfunc1[tree, m - 1]``, the latter written by cell ``(m - 1, tree - 1)``.
  Both predecessors are fixed by the structure, so the dependency graph is
  static and can be levelled once.

  The grouping is *not* the tree level -- each wave holds one cell from almost
  every tree -- and it is not the anti-diagonal either: with ``m - 1 == var +
  1`` off the diagonal, which is the generic D-vine cell, ``(var + 1, tree -
  1)`` lands on the same anti-diagonal as ``(var, tree)``. Levelling the actual
  graph is both correct and tighter than any fixed key.

  Parameters
  ----------
  structure : RVineStructure
      The vine structure, read for ``min_array`` / ``struct_array``.
  trunc_lvl : int
      Truncation level.
  cond_positions : callable, or None, optional
      ``(tree, var) -> columns`` of the cell's conditioning set, for a vine
      whose pairs read it: the cell then also waits for each of those
      variables to be finalized, i.e. for cell ``(column, 0)``.

  Returns
  -------
  list of list of tuple
      Cells per wave, in execution order. Every cell in a wave is
      independent of the others, so a wave is one stacked call.
  """
  d = int(structure.dim)
  deps: dict[tuple[int, int], set[tuple[int, int]]] = {}
  for var in range(d - 2, -1, -1):
    for tree in range(min(trunc_lvl - 1, d - var - 2), -1, -1):
      pred: set[tuple[int, int]] = set()
      if tree + 1 <= min(trunc_lvl - 1, d - var - 2):
        pred.add((var, tree + 1))
      m = int(structure.min_array(tree, var))
      if m == int(structure.struct_array(tree, var, True)):
        pred.add((m - 1, tree))
      elif tree - 1 >= 0:
        pred.add((m - 1, tree - 1))
      if cond_positions is not None:
        pred.update((c, 0) for c in cond_positions(tree, var))
      deps[var, tree] = pred
  for cell in deps:
    deps[cell] &= deps.keys()

  # Longest-path level of each cell; the descending sweep is already a
  # topological order, so one pass settles it.
  depth: dict[tuple[int, int], int] = {}
  for cell in sorted(deps, key=lambda c: (-c[0], -c[1])):
    depth[cell] = 0 if not deps[cell] else 1 + max(depth[p] for p in deps[cell])
  n_waves = max(depth.values()) + 1 if depth else 0
  return [sorted(c for c in deps if depth[c] == k) for k in range(n_waves)]


def wave_wiring(
  structure: RVineStructure, cells: list[tuple[int, int]]
) -> dict[str, list[Any]]:
  """One inverse wave's per-cell tables, as rows of the flattened scratch.

  The inverse scratch is ``(trunc_lvl + 1, d, n)``, flattened to
  ``((trunc_lvl + 1) * d, n)``, so a cell's slot is one row and a whole wave is
  one gather per input and one scatter per output.

  Parameters
  ----------
  structure : RVineStructure
      The vine structure.
  cells : list of tuple of int
      The wave's ``(var, tree)`` cells.

  Returns
  -------
  dict
      Positions: each cell's input rows (``col0_src`` / ``col1_src``, with
      ``h1_pos`` / ``h1_src`` the cells whose second input reads ``hfunc1``),
      its output row, and -- only for the cells whose ``hfunc1`` the
      next-inner inversion reads -- ``h1_rows`` and ``out_hfunc1``.
  """
  d = int(structure.dim)
  keys = ("col0_src", "col1_src", "h1_pos", "h1_src", "out_hinv2")
  out: dict[str, list[Any]] = {k: [] for k in (*keys, "h1_rows", "out_hfunc1")}
  for slot, (var, tree) in enumerate(cells):
    m = int(structure.min_array(tree, var))
    out["col0_src"].append((tree + 1) * d + var)
    out["col1_src"].append(tree * d + (m - 1))
    if m != int(structure.struct_array(tree, var, True)):
      out["h1_pos"].append(slot)
      out["h1_src"].append(tree * d + (m - 1))
    out["out_hinv2"].append(tree * d + var)
    if var < d - 1 and bool(structure.needed_hfunc1(tree, var)):
      out["h1_rows"].append(slot)
      out["out_hfunc1"].append((tree + 1) * d + var)
  return out


def as_arrays(
  tables: dict[str, list[Any]],
  xp: Namespace[ArrayT],
  device: Any,  # noqa: ANN401 - a namespace's device object
) -> dict[str, Any]:
  """The integer and boolean tables as arrays on ``device``.

  Parameters
  ----------
  tables : dict
      What :func:`level_wiring` or :func:`wave_wiring` returns.
  xp : module
      The array namespace to build on.
  device : object
      The device the cascade's scratch lives on.

  Returns
  -------
  dict
      The position tables, each an ``int64`` array; the masks stay lists.
  """
  # The dtype is a module attribute the array API standard names, which the
  # `Namespace` protocol leaves out, as it leaves out `bool`, whose name would
  # shadow the builtin in every annotation of the class.
  ns: Any = xp
  return {
    key: xp.asarray(values, dtype=ns.int64, device=device)
    for key, values in tables.items()
    if key not in _LEVEL_MASKS
  }


def gather(
  xp: Namespace[ArrayT],
  w: dict[str, Any],
  hfunc1: ArrayT,
  hfunc2: ArrayT,
  hfunc1_sub: ArrayT | None,
  hfunc2_sub: ArrayT | None,
) -> ArrayT:
  """One tree level's inputs, gathered from the h-function scratch.

  Given the left-limit scratch as well, the result is the four-column
  ``[u1, u2, u1^-, u2^-]`` a discrete slot reads, gathered through the same
  wiring; a continuous argument's left limit is its own value, so a stale entry
  of the scratch is never read.

  Parameters
  ----------
  xp : module
      The array namespace of the scratch.
  w : dict
      The level's position tables, as :func:`as_arrays` returns them.
  hfunc1, hfunc2 : array, shape (n, d), dtype float
      The h-function scratch.
  hfunc1_sub, hfunc2_sub : array, shape (n, d), dtype float, or None
      The left-limit scratch, for a vine with discrete variables.

  Returns
  -------
  array, shape (P, n, 2) or (P, n, 4), dtype float
      One input per edge.
  """

  # Plain indexing rather than namespace calls: it is what both libraries do
  # natively, and this runs once per level of every cascade.
  def pick(h1: Any, h2: Any) -> tuple[Any, Any]:  # noqa: ANN401 - scratch arrays
    col0 = h2[:, w["col0_src"]]
    col1 = h2[:, w["col1_src"]]
    if int(w["h1_pos"].shape[0]):
      col1[:, w["h1_pos"]] = h1[:, w["h1_src"]]
    return col0, col1

  col0, col1 = pick(hfunc1, hfunc2)
  if hfunc2_sub is None:
    return xp.stack([col0.T, col1.T], axis=-1)
  sub0, sub1 = pick(hfunc1_sub, hfunc2_sub)
  if int(w["cont1_pos"].shape[0]):
    sub0[:, w["cont1_pos"]] = col0[:, w["cont1_pos"]]
  if int(w["cont2_pos"].shape[0]):
    sub1[:, w["cont2_pos"]] = col1[:, w["cont2_pos"]]
  return xp.stack([col0.T, col1.T, sub0.T, sub1.T], axis=-1)


def apply_wave(
  xp: Namespace[ArrayT],
  w: dict[str, Any],
  stacked: Any,  # noqa: ANN401 - LoopedPairs or a pair class's own stack
  hinv2: ArrayT,
  hfunc1: ArrayT,
) -> None:
  """Invert one wave's pairs, in place on the flattened scratch.

  Parameters
  ----------
  xp : module
      The array namespace of the scratch.
  w : dict
      The wave's position tables, as :func:`as_arrays` returns them.
  stacked : object
      The wave's pairs, evaluated through ``hinv2`` and ``hfunc1``.
  hinv2, hfunc1 : array, shape ((trunc_lvl + 1) * d, n), dtype float
      The flattened inverse scratch, written in place.

  Returns
  -------
  None
  """
  h2: Any = hinv2
  h1: Any = hfunc1
  col0 = h2[w["col0_src"]]
  col1 = h2[w["col1_src"]]
  if int(w["h1_pos"].shape[0]):
    col1[w["h1_pos"]] = h1[w["h1_src"]]
  inv = stacked.hinv2(xp.stack([col0, col1], axis=-1))
  h2[w["out_hinv2"]] = inv
  rows = w["h1_rows"]
  if int(rows.shape[0]) == 0:
    return
  u_after = xp.stack([inv[rows], col1[rows]], axis=-1)
  h1[w["out_hfunc1"]] = stacked.hfunc1(u_after, rows)


def stack_of(
  pairs: Sequence[BicopLike[Any]], tables: dict[str, list[Any]] | None
) -> Any:  # noqa: ANN401 - a pair class's own stack
  """The pairs' own stack, or ``None`` where their class supplies none.

  Parameters
  ----------
  pairs : sequence of BicopLike
      A tree level's pairs, or an inverse wave's continuous ones.
  tables : dict, or None
      The level's :func:`level_wiring` tables, which say which pairs are
      discrete and which outputs the next tree reads; ``None`` for a wave.

  Returns
  -------
  object, or None
      The stack, evaluated as :class:`LoopedPairs` is; ``None`` when the pairs
      are of several classes or their class has no ``_stack_pairs``.
  """
  kinds = {type(p) for p in pairs}
  if len(kinds) != 1:
    return None
  hook = getattr(kinds.pop(), "_stack_pairs", None)
  if hook is None:
    return None
  if tables is None:
    off = [False] * len(pairs)
    return hook(pairs, disc1=off, disc2=off, needs_h1=off, needs_h2=off)
  return hook(
    pairs,
    disc1=tables["disc1"],
    disc2=tables["disc2"],
    needs_h1=tables["needs_h1"],
    needs_h2=tables["needs_h2"],
  )


class LoopedPairs:
  """A group of pair copulas evaluated one at a time, as a stacked level.

  The default a level or wave falls back to wherever its pairs supply no stack
  of their own: each pair is called through its own dispatchers, on exactly the
  input the per-edge walk would hand it, so the result is what that walk
  computes, bit for bit.

  Parameters
  ----------
  xp : module
      The array namespace of the inputs.
  pairs : sequence of BicopLike
      The group's pairs, in edge (or cell) order.
  tables : dict
      The group's tables as Python lists: :func:`level_wiring`'s for a tree
      level, :func:`wave_wiring`'s for an inverse wave.
  x : list, or None, optional
      Each pair's conditioning matrix, or ``None`` throughout.
  """

  def __init__(
    self,
    xp: Namespace[ArrayT],
    pairs: Sequence[BicopLike[Any]],
    tables: dict[str, list[Any]],
    x: Sequence[Any] | None = None,
  ) -> None:
    self._xp = xp
    self._pairs = list(pairs)
    self._tables = tables
    self._x = list(x) if x is not None else [None] * len(self._pairs)

  def _blank(self, like: Any) -> Any:  # noqa: ANN401 - arrays
    """An uninitialized output shaped like ``like``, filled row by row."""
    return self._xp.empty(like.shape, dtype=like.dtype, device=like.device)

  def evaluate(
    self, u: ArrayT, *, with_pdf: bool, every_h2: bool
  ) -> tuple[Any, Any, Any, Any, Any]:
    """``(pdf, hfunc1, hfunc2, hfunc1^-, hfunc2^-)`` for one tree level.

    Only what the next tree reads is evaluated: ``hfunc1`` (and on a discrete
    second argument its left limit) where the edge's ``needs_h1``, ``hfunc2``
    (and on a discrete first argument its left limit) where ``needs_h2`` or,
    with ``every_h2``, at every edge. The other rows are left unset, and the
    cascade reads none of them.

    Parameters
    ----------
    u : array, shape (P, n, 2) or (P, n, 4), dtype float
        The level's inputs, from :func:`gather`.
    with_pdf : bool
        Whether to evaluate the densities; ``None`` in their place otherwise.
    every_h2 : bool
        Evaluate ``hfunc2`` at every edge, the running transform of a
        Rosenblatt cascade.

    Returns
    -------
    tuple of array
        Five ``(P, n)`` arrays, the first ``None`` unless ``with_pdf``; on a
        two-column input the left limits are the values themselves.
    """
    xp = self._xp
    ua: Any = u
    t = self._tables
    four = int(ua.shape[2]) == 4
    like = ua[:, :, 0]
    pdf = self._blank(like) if with_pdf else None
    h1, h2 = self._blank(like), self._blank(like)
    h1s, h2s = (self._blank(like), self._blank(like)) if four else (h1, h2)
    for p, pair in enumerate(self._pairs):
      d1, d2 = t["disc1"][p], t["disc2"][p]
      edge = ua[p]
      u_e = edge if (d1 or d2) else edge[:, :2]
      x_e = self._x[p]
      if pdf is not None:
        pdf[p] = pair_eval(pair.pdf, u_e, x=x_e)
      want_h1 = t["needs_h1"][p]
      want_h2 = every_h2 or t["needs_h2"][p]
      if want_h1:
        h1[p] = pair_eval(pair.hfunc1, u_e, x=x_e)
      if want_h2:
        h2[p] = pair_eval(pair.hfunc2, u_e, x=x_e)
      if four and want_h1 and d2:
        u_h1 = xp.stack(
          [edge[:, 0], edge[:, 3], edge[:, 2], edge[:, 3]], axis=-1
        )
        h1s[p] = pair_eval(pair.hfunc1, u_h1, x=x_e)
      if four and want_h2 and d1:
        u_h2 = xp.stack(
          [edge[:, 2], edge[:, 1], edge[:, 2], edge[:, 3]], axis=-1
        )
        h2s[p] = pair_eval(pair.hfunc2, u_h2, x=x_e)
    return pdf, h1, h2, h1s, h2s

  def hinv2(self, u: ArrayT) -> Any:  # noqa: ANN401 - an array
    """Each cell's ``hinv2``, on its continuous pair.

    Parameters
    ----------
    u : array, shape (K, n, 2), dtype float
        The wave's inputs.

    Returns
    -------
    array, shape (K, n), dtype float
        The inverted values.
    """
    ua: Any = u
    out = self._blank(ua[:, :, 0])
    for k, pair in enumerate(self._pairs):
      out[k] = pair_eval(pair.hinv2, ua[k], x=self._x[k])
    return out

  def hfunc1(self, u: ArrayT, rows: Any) -> Any:  # noqa: ANN401 - arrays
    """``hfunc1`` at some of the wave's cells.

    Parameters
    ----------
    u : array, shape (R, n, 2), dtype float
        The inputs, one per listed cell.
    rows : array, shape (R,), dtype int
        Which of the wave's cells they belong to; read from the wave's own
        tables here, so no device array is walked.

    Returns
    -------
    array, shape (R, n), dtype float
        The h-function values.
    """
    del rows
    ua: Any = u
    out = self._blank(ua[:, :, 0])
    for i, k in enumerate(self._tables["h1_rows"]):
      out[i] = pair_eval(self._pairs[k].hfunc1, ua[i], x=self._x[k])
    return out
