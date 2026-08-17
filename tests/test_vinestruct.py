import os

import numpy as np
import pytest

import pyvinecopulib as pv

from .helpers import compare_rvinestructure


def _trees_from_structure(
  s: pv.RVineStructure,
) -> list[list[tuple[int, int, list[int]]]]:
  """Invert an RVineStructure into per-tree ``(a, b, C)`` edge lists.

  Mirrors the upstream ``RVineTrees(order, struct_array)`` reader: edge
  ``(t, e)`` is ``(order[e], struct_array(t, e), {struct_array(0..t-1, e)})``,
  with labels read in labelled (non-natural) order.
  """
  d = int(s.dim)
  trunc_lvl = int(s.trunc_lvl)
  order = [int(v) for v in s.order]
  sa = s.get_struct_array(natural_order=False)
  trees = []
  for t in range(trunc_lvl):
    tree = [
      (order[e], int(sa[t][e]), [int(sa[i][e]) for i in range(t)])
      for e in range(d - 1 - t)
    ]
    trees.append(tree)
  return trees


def _edge_sets(
  trees: list[list[tuple[int, int, list[int]]]],
) -> list[set[tuple[frozenset[int], frozenset[int]]]]:
  """Per-tree set of unordered edges ``({a, b}, C)`` (encoding-independent)."""
  return [
    {(frozenset((a, b)), frozenset(c)) for (a, b, c) in tree} for tree in trees
  ]


def _dvine_trees(order: list[int]) -> list[list[tuple[int, int, list[int]]]]:
  """The canonical D-vine edge lists for a given variable order."""
  d = len(order)
  trees = []
  for t in range(d - 1):
    tree = [
      (order[e], order[e + t + 1], [order[e + 1 + k] for k in range(t)])
      for e in range(d - 1 - t)
    ]
    trees.append(tree)
  return trees


def test_rvinestructure_from_trees_dvine() -> None:
  # A hand-built D-vine decomposition assembles into a valid structure that
  # represents the same vine (same per-tree edge sets), for several orders.
  for order in ([1, 2, 3, 4, 5], [3, 1, 4, 2, 5], [2, 5, 1, 4, 3, 6]):
    d = len(order)
    trees = _dvine_trees(order)
    s = pv.RVineStructure.from_trees(d, trees)
    assert isinstance(s, pv.RVineStructure)
    assert s.dim == d
    assert s.trunc_lvl == d - 1
    assert _edge_sets(_trees_from_structure(s)) == _edge_sets(trees)
    # The assembled matrix is a valid R-vine: it round-trips through
    # from_struct_array unchanged.
    rebuilt = pv.RVineStructure.from_struct_array(
      s.order, s.get_struct_array(natural_order=False)
    )
    assert np.array_equal(np.asarray(rebuilt.matrix), np.asarray(s.matrix))


def test_rvinestructure_from_trees_roundtrips_selected_vine() -> None:
  # For a vine chosen by the C++ selector, inverting its structure into trees
  # and re-assembling yields the same vine (same per-tree edge sets), and the
  # result is deterministic and insensitive to within-tree edge order.
  for seed in range(5):
    rng = np.random.default_rng(seed)
    d = 6
    u = pv.to_pseudo_obs(
      rng.standard_normal((500, d)) @ rng.standard_normal((d, d))
    )
    vc = pv.Vinecop.from_data(u)
    trees = _trees_from_structure(vc.structure)
    s = pv.RVineStructure.from_trees(d, trees)
    assert _edge_sets(_trees_from_structure(s)) == _edge_sets(trees)
    # Determinism + within-tree edge-order invariance.
    s_again = pv.RVineStructure.from_trees(d, trees)
    shuffled = [list(reversed(tree)) for tree in trees]
    s_shuf = pv.RVineStructure.from_trees(d, shuffled)
    assert np.array_equal(np.asarray(s.matrix), np.asarray(s_again.matrix))
    assert np.array_equal(np.asarray(s.matrix), np.asarray(s_shuf.matrix))


def test_rvinestructure_from_trees_truncated_and_empty() -> None:
  # A truncated decomposition (fewer trees than d-1) assembles at the right
  # truncation level; an empty decomposition is a fully truncated structure.
  order = [1, 2, 3, 4, 5]
  full = _dvine_trees(order)
  truncated = full[:2]
  s = pv.RVineStructure.from_trees(5, truncated)
  assert s.dim == 5
  assert s.trunc_lvl == 2
  assert _edge_sets(_trees_from_structure(s)) == _edge_sets(truncated)

  s_empty = pv.RVineStructure.from_trees(5, [])
  assert s_empty.dim == 5
  assert s_empty.trunc_lvl == 0


def test_rvinestructure_from_trees_rejects_bad_edge_count() -> None:
  # Tree t must hold exactly d - 1 - t edges; otherwise assembly fails.
  bad = [[(1, 2, []), (2, 3, [])]]  # tree 0 of a d=5 vine needs 4 edges
  with pytest.raises(RuntimeError):
    pv.RVineStructure.from_trees(5, bad)


def test_rvinestructure(unique_json_path: str) -> None:
  d = 5
  # Test RVineStructure class
  rvine = pv.RVineStructure(d)
  assert isinstance(rvine, pv.RVineStructure)
  assert rvine.dim == d
  assert rvine.trunc_lvl == d - 1
  assert rvine.order == list(range(1, d + 1))
  matrix = rvine.matrix
  assert isinstance(matrix, np.ndarray)
  assert matrix.shape == (d, d)
  assert matrix.dtype == np.uint64
  assert np.all(np.logical_and(matrix >= 0, matrix <= d))

  # Should be the same as the previous one
  dvine = pv.DVineStructure(rvine.order)
  compare_rvinestructure(dvine, rvine, True)

  # Test to_json and from_json
  new_rvine = pv.RVineStructure.from_json(rvine.to_json())
  compare_rvinestructure(rvine, new_rvine)
  filename = os.fspath(unique_json_path)
  rvine.to_file(filename)
  new_rvine = pv.RVineStructure.from_file(filename)
  compare_rvinestructure(rvine, new_rvine)

  # CBOR round-trip: a ``.cbor`` filename selects the binary format
  # (vinecopulib#684).
  cbor_filename = filename.removesuffix(".json") + ".cbor"
  rvine.to_file(cbor_filename)
  compare_rvinestructure(rvine, pv.RVineStructure.from_file(cbor_filename))
  with open(cbor_filename, "rb") as f:
    assert f.read(1) != b"{"


def test_rvinestructure_get_trees_faithful_roundtrip() -> None:
  # get_trees() is the faithful inverse of from_trees(): re-assembling the
  # decomposition reproduces the exact R-vine matrix, not just an equivalent
  # vine (vinecopulib#698).
  for seed in range(5):
    rng = np.random.default_rng(seed)
    d = 6
    u = pv.to_pseudo_obs(
      rng.standard_normal((500, d)) @ rng.standard_normal((d, d))
    )
    s = pv.Vinecop.from_data(u).structure
    st = s.get_trees()
    assert _edge_sets(st) == _edge_sets(_trees_from_structure(s))
    s2 = pv.RVineStructure.from_trees(s.dim, st)
    assert s2 == s
    assert np.array_equal(np.asarray(s2.matrix), np.asarray(s.matrix))
    assert list(s2.order) == list(s.order)


def test_rvinestructure_bulk_triangular_getters_match_scalar_accessors() -> (
  None
):
  s = pv.RVineStructure.simulate(d=6, seeds=[1, 2, 3])
  min_array = s.get_min_array()
  hfunc1, hfunc2 = s.get_needed_hfunc1(), s.get_needed_hfunc2()

  assert len(min_array) == len(hfunc1) == len(hfunc2) == int(s.trunc_lvl)
  for tree, row in enumerate(min_array):
    assert len(row) == int(s.dim) - 1 - tree
    for edge, value in enumerate(row):
      assert value == s.min_array(tree, edge)
      assert hfunc1[tree][edge] == s.needed_hfunc1(tree, edge)
      assert hfunc2[tree][edge] == s.needed_hfunc2(tree, edge)
