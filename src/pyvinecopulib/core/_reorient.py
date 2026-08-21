"""Relabel a vine structure onto a chosen sampling-order tail.

Internal. Backs :meth:`~pyvinecopulib.core.VinecopBase.reorient`, the
``conditioning_set`` argument of the Rosenblatt transforms and of conditional
sampling, and conditioning-aware structure selection.

Given a structure and a set of 1-based variables, it produces an equivalent
structure whose order *ends* with that set, together with the map sending each
``(tree, edge)`` slot of the new structure to the slot of the old one holding the
same pair copula, and whether that pair's arguments have to be swapped. The
relabeling is value-preserving: the density and log-likelihood are unchanged, and
only which conditional distributions the sampling order walks changes.

The peel that steers the diagonal towards a chosen tail lives upstream and is
reachable only through the compiled ``Vinecop.reorient``, so that is what runs --
on a throwaway vine of independence copulas, whose pair copulas are discarded.
Borrowing it keeps this exactly as admissible as ``Vinecop`` is, down to the
error messages, and is the same trade ``VinecopBase.select`` makes when it
borrows ``_select_spanning_tree`` and ``RVineStructure.from_trees``.
"""

from __future__ import annotations

from typing import Any, NamedTuple

__all__ = ["Reorientation", "reorientation"]


class Reorientation(NamedTuple):
  """A structure relabeling and the pair-copula slot map it induces.

  Attributes
  ----------
  structure : RVineStructure
      The relabeled structure, whose order ends with the requested variables.
  locations : dict
      ``(tree, edge)`` of the relabeled structure -> ``(edge, flipped)`` of the
      original one. The tree index is not part of the value: a relabeling
      permutes edges within a tree and may swap a pair copula's arguments, but
      never moves a pair between trees.
  identity : bool
      Whether the requested variables already formed the order tail, in which
      case ``structure`` is the original object and ``locations`` is empty.
  """

  structure: Any
  locations: dict[tuple[int, int], tuple[int, bool]]
  identity: bool


def _slot_key(structure: Any, tree: int, edge: int) -> tuple[int, Any]:
  """The diagonal variable and label sets identifying one matrix slot."""
  # A slot hosts the edge with conditioned pair {order[edge],
  # struct_array(tree, edge)} and conditioning set {struct_array(i, edge) :
  # i < tree} -- the same keying `VinecopBase.select` places its pairs by.
  diag = int(structure.order[edge])
  partner = int(structure.struct_array(tree, edge, natural_order=False))
  conditioning = frozenset(
    int(structure.struct_array(i, edge, natural_order=False))
    for i in range(tree)
  )
  return diag, (frozenset((diag, partner)), conditioning)


def reorientation(structure: Any, conditioning_set: list[int]) -> Reorientation:
  """Relabel ``structure`` so that its order ends with ``conditioning_set``.

  Parameters
  ----------
  structure : RVineStructure
      The structure to relabel. A truncated one relabels like any other: the
      trees above the truncation are independence, so the peel has nothing to
      move there and the slot map covers only the trees that exist.
  conditioning_set : list of int
      1-based variable labels to place at the tail of the order.

  Returns
  -------
  Reorientation
      The relabeled structure, the slot map onto ``structure``'s slots, and
      whether the relabeling is the identity.

  Raises
  ------
  RuntimeError
      If ``conditioning_set`` is empty, holds duplicates or entries outside
      ``1, ..., d``, does not leave a variable free, or is not admissible as a
      sampling-order tail. The messages are the compiled relabeling's, so a
      caller can catch and match the same thing whichever class it went
      through.
  """
  from ..pyvinecopulib_ext import Vinecop

  d = int(structure.dim)
  cond = [int(v) for v in conditioning_set]
  order = [int(v) for v in structure.order]

  # Validated ahead of the identity check, as the compiled relabeling does:
  # otherwise the whole variable set would look like a tail and pass, and
  # `sample_conditional` would be left with nothing to draw.
  if not 1 <= len(cond) < d:
    raise RuntimeError(
      "conditioning_set must contain between 1 and d - 1 variables."
    )
  if len(set(cond)) != len(cond):
    raise RuntimeError("conditioning_set must not contain duplicates.")
  if any(v < 1 or v > d for v in cond):
    raise RuntimeError("conditioning_set entries must be in 1, ..., d.")
  # Preserve the exact representation when the requested variables already form
  # the tail: evaluating through it is then bit-identical to evaluating the vine
  # directly, and no pair copula has to be flipped.
  if set(cond) == set(order[d - len(cond) :]):
    return Reorientation(structure, {}, True)

  # An empty store means independence everywhere, which is all the peel needs:
  # it reads the structure and never the pairs.
  placeholder = Vinecop.from_structure(structure=structure)
  placeholder.reorient(cond)
  relabeled = placeholder.structure

  # A truncated structure holds no entries above its truncation level, so the
  # map covers the trees it stores; the trees above are independence, and a
  # relabeling of independence is independence.
  trunc = int(structure.trunc_lvl)
  index: list[dict[Any, tuple[int, int]]] = []
  for tree in range(trunc):
    by_key: dict[Any, tuple[int, int]] = {}
    for edge in range(d - 1 - tree):
      diag, key = _slot_key(structure, tree, edge)
      by_key[key] = (edge, diag)
    index.append(by_key)

  locations: dict[tuple[int, int], tuple[int, bool]] = {}
  for tree in range(trunc):
    for edge in range(d - 1 - tree):
      diag, key = _slot_key(relabeled, tree, edge)
      old_edge, old_diag = index[tree][key]
      # Flipped exactly when the diagonal variable is no longer the pair's
      # first argument -- the rule the selector's finalization also applies.
      locations[(tree, edge)] = (old_edge, old_diag != diag)
  return Reorientation(relabeled, locations, False)
