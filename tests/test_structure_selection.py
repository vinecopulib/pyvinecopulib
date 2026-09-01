"""Tests for the array-agnostic vine structure-selection engine.

Covers the boost-based spanning-tree primitive (``_select_spanning_tree``)
and ``VinecopBase.select``, whose contract is exact parity with ``Vinecop``:
same selected structure (identical R-vine matrix encoding) and same reused
pair copulas (identical density, no re-fit).
"""

from typing import Any, Optional, cast

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import BicopLike, VinecopBase

# Internal C++ primitive backing Python structure selection (boost prim /
# kruskal / Wilson). Imported from the extension directly as it has no public
# wrapper; it is an implementation detail of ``VinecopBase.select``.
from pyvinecopulib.pyvinecopulib_ext import _select_spanning_tree

from .conftest import HostedVinecop


class _CppBicopLike:
  """Adapt a ``Bicop`` to the ``(u, x)`` :class:`BicopLike` signature.

  ``VinecopBase.select`` / ``fit`` call ``hfunc1(u, x_e)``; the
  ``Bicop`` methods take ``(u, num_threads)`` instead, so wrap it (delegating
  the full ``BicopLike`` surface, including ``flip``). Used only to drive the
  array-agnostic selector from a NumPy backend in these tests.
  """

  def __init__(self, bicop: pv.Bicop) -> None:
    self._b = bicop

  def pdf(self, u: Any, x: Any = None) -> Any:
    return self._b.pdf(np.asarray(u))

  def cdf(self, u: Any, x: Any = None) -> Any:
    return self._b.cdf(np.asarray(u))

  def hfunc1(self, u: Any, x: Any = None) -> Any:
    return self._b.hfunc1(np.asarray(u))

  def hfunc2(self, u: Any, x: Any = None) -> Any:
    return self._b.hfunc2(np.asarray(u))

  def hinv1(self, u: Any, x: Any = None) -> Any:
    return self._b.hinv1(np.asarray(u))

  def hinv2(self, u: Any, x: Any = None) -> Any:
    return self._b.hinv2(np.asarray(u))

  def sample(
    self,
    n: int,
    *,
    x: Any = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> Any:
    return self._b.sample(n, qrng=qrng, seeds=seeds or [])

  def flip(self) -> "_CppBicopLike":
    return _CppBicopLike(self._b.flip())


_GAUSSIAN = pv.FitControlsBicop(family_set=[pv.families.gaussian])
_TLL = pv.FitControlsBicop(family_set=[pv.families.tll])


def _gaussian_fit_edge(
  tree: int, edge: int, u_e: object, x_e: object
) -> _CppBicopLike:
  return _CppBicopLike(pv.Bicop.from_data(np.asarray(u_e), controls=_GAUSSIAN))


def _tll_fit_edge(
  tree: int, edge: int, u_e: object, x_e: object
) -> _CppBicopLike:
  return _CppBicopLike(pv.Bicop.from_data(np.asarray(u_e), controls=_TLL))


def _vine_controls(
  family: pv.BicopFamily, trunc_lvl: int = 20
) -> pv.FitControlsVinecop:
  return pv.FitControlsVinecop(
    family_set=[family], trunc_lvl=trunc_lvl, num_threads=1
  )


def _correlated_pseudo_obs(seed: int, d: int, n: int = 400) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return pv.to_pseudo_obs(
    rng.standard_normal((n, d)) @ rng.standard_normal((d, d))
  )


# ---------------------------------------------------------------------------
# Spanning-tree primitive
# ---------------------------------------------------------------------------


def test_select_spanning_tree_mst_is_unique_and_correct() -> None:
  # Diamond graph on 4 vertices; weight = 1 - |tau|, so lower weight is the
  # stronger edge that an MST should prefer.
  n = 4
  edges = [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3), (0, 3)]
  weights = [0.1, 0.2, 0.15, 0.9, 0.85, 0.95]
  prim = sorted(_select_spanning_tree(n, edges, weights, "mst_prim"))
  kruskal = sorted(_select_spanning_tree(n, edges, weights, "mst_kruskal"))
  # The three light edges form the unique MST.
  assert prim == [0, 1, 2]
  # On distinct weights the MST is unique, so Prim and Kruskal agree.
  assert prim == kruskal


def test_select_spanning_tree_early_exit_returns_all() -> None:
  # A candidate graph that is already a spanning tree is returned as-is.
  n = 4
  tree_edges = [(0, 1), (1, 2), (2, 3)]
  sel = _select_spanning_tree(n, tree_edges, [0.5, 0.5, 0.5], "mst_prim")
  assert sorted(sel) == [0, 1, 2]


def test_select_spanning_tree_random_is_reproducible_and_valid() -> None:
  n = 5
  rng = np.random.default_rng(0)
  edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
  weights = list(rng.uniform(0.0, 1.0, size=len(edges)))
  for algo in ("random_weighted", "random_unweighted"):
    t1 = sorted(_select_spanning_tree(n, edges, weights, algo, [42]))
    t2 = sorted(_select_spanning_tree(n, edges, weights, algo, [42]))
    # Reproducible given the same seed, and a valid spanning tree.
    assert t1 == t2
    assert len(t1) == n - 1


def test_select_spanning_tree_weighted_random_accepts_zero_weights() -> None:
  n = 5
  edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
  selected = _select_spanning_tree(
    n, edges, [1.0] * len(edges), "random_weighted", [42]
  )
  assert len(selected) == n - 1


# ---------------------------------------------------------------------------
# VinecopBase.select — exact parity with Vinecop
# ---------------------------------------------------------------------------


def test_select_matches_vinecop_matrix() -> None:
  # The selector reproduces Vinecop's R-vine matrix *exactly* (not just an
  # equivalent encoding), when both use identical pair fits — including a
  # truncated case, where the finalization is also sensitive to the
  # within-tree edge order.
  for seed in range(3):
    for d, trunc in ((4, 20), (5, 20), (7, 20), (8, 3)):
      u = _correlated_pseudo_obs(seed, d)
      mine, _ = VinecopBase.select(u, _gaussian_fit_edge, trunc_lvl=trunc)
      cpp = pv.Vinecop.from_data(
        u, controls=_vine_controls(pv.families.gaussian, trunc)
      ).structure
      assert np.array_equal(np.asarray(mine.matrix), np.asarray(cpp.matrix))
      assert list(mine.order) == list(cpp.order)


def test_select_reused_pairs_match_vinecop_exactly() -> None:
  # THE flip/placement regression test. With TLL pairs — whose fits are
  # orientation-sensitive, so any error in the pair placement, the flip rule,
  # the pc_data column roles, or the matrix encoding shows up as an O(1)
  # density error — the vine assembled from select's reused pairs must match
  # ``Vinecop.from_data`` to floating-point noise.
  for seed in (0, 5):
    d = 5
    u = _correlated_pseudo_obs(seed, d, n=800)
    structure, pairs = VinecopBase.select(u, _tll_fit_edge)
    auto = pv.Vinecop.from_data(u, controls=_vine_controls(pv.families.tll))
    assert np.array_equal(
      np.asarray(structure.matrix), np.asarray(auto.structure.matrix)
    )
    rng = np.random.default_rng(seed)
    grid = rng.uniform(0.02, 0.98, size=(300, d))
    np.testing.assert_allclose(
      HostedVinecop(pairs, structure).pdf(grid),
      auto.pdf(grid),
      rtol=1e-12,
      atol=1e-12,
    )


def test_select_refit_pdf_parity_with_vinecop() -> None:
  # Refitting on the selected structure also reproduces the density of the
  # auto-selected vine for orientation-symmetric (Gaussian) pairs.
  rng = np.random.default_rng(123)
  for seed in range(3):
    d = 5
    u = _correlated_pseudo_obs(seed, d)
    mine, _ = VinecopBase.select(u, _gaussian_fit_edge)
    controls = _vine_controls(pv.families.gaussian)
    auto = pv.Vinecop.from_data(u, controls=controls)
    refit = pv.Vinecop.from_data(u, structure=mine, controls=controls)
    grid = rng.uniform(0.02, 0.98, size=(200, d))
    # 1.4e-10 measured. It was 1e-8 while the margin renormalization was a
    # fixed three-sweep budget, which left `fit(a, b)` and `fit(b, a)`
    # different models; it converges now, so this can be pinned two decades
    # tighter.
    np.testing.assert_allclose(
      refit.pdf(grid), auto.pdf(grid), rtol=1e-9, atol=1e-9
    )


def test_select_respects_truncation() -> None:
  d = 6
  u = _correlated_pseudo_obs(0, d)
  for trunc in (1, 2, 3):
    s, pairs = VinecopBase.select(u, _gaussian_fit_edge, trunc_lvl=trunc)
    assert s.dim == d
    assert s.trunc_lvl == trunc
    assert [len(row) for row in pairs] == [d - 1 - t for t in range(trunc)]


def test_select_random_weighted_is_reproducible() -> None:
  d = 5
  u = _correlated_pseudo_obs(1, d)
  s1, _ = VinecopBase.select(
    u, _gaussian_fit_edge, tree_algorithm="random_weighted", seeds=[1, 2, 3]
  )
  s2, _ = VinecopBase.select(
    u, _gaussian_fit_edge, tree_algorithm="random_weighted", seeds=[1, 2, 3]
  )
  # Same seeds -> identical structure; and a valid full R-vine.
  assert np.array_equal(np.asarray(s1.matrix), np.asarray(s2.matrix))
  assert s1.trunc_lvl == d - 1


def test_select_random_weighted_accepts_small_samples() -> None:
  u = _correlated_pseudo_obs(1, 5, n=10)
  structure, _ = VinecopBase.select(
    u,
    _gaussian_fit_edge,
    tree_algorithm="random_weighted",
    seeds=[1, 2, 3],
  )
  assert structure.trunc_lvl == 4


def test_select_rejects_unknown_tree_algorithm() -> None:
  u = _correlated_pseudo_obs(1, 4, n=50)
  with pytest.raises(ValueError, match="tree_algorithm must be one of"):
    VinecopBase.select(
      u, _gaussian_fit_edge, tree_algorithm="definitely_not_an_algorithm"
    )


def test_select_skips_hfunctions_after_final_tree() -> None:
  calls = {"hfunc1": 0, "hfunc2": 0}

  class CountingPair(_CppBicopLike):
    def hfunc1(self, u: Any, x: Any = None) -> Any:
      calls["hfunc1"] += 1
      return super().hfunc1(u, x)

    def hfunc2(self, u: Any, x: Any = None) -> Any:
      calls["hfunc2"] += 1
      return super().hfunc2(u, x)

    def flip(self) -> "CountingPair":
      return self

  def fit_edge(tree: int, edge: int, u_e: object, x_e: object) -> CountingPair:
    return CountingPair(pv.Bicop())

  VinecopBase.select(_correlated_pseudo_obs(1, 6, n=50), fit_edge, trunc_lvl=1)
  assert calls == {"hfunc1": 0, "hfunc2": 0}


def test_bicop_base_flip_default_raises() -> None:
  # BicopBase provides a raising flip default: implementing it is only needed
  # to host a custom pair in structure selection.
  from .conftest import GaussianBicop

  with pytest.raises(NotImplementedError, match="flip"):
    GaussianBicop(base_rho=0.5).flip()


def test_compiled_bicop_hosted_unwrapped_matches_vinecop() -> None:
  # A simplified vine may host the compiled ``Bicop`` directly: it satisfies
  # ``BicopLike`` nominally and takes no conditioning argument, so the cascade
  # must not hand it one. Every other pair here goes through ``_CppBicopLike``,
  # which accepts (and drops) ``x``; this pins the unwrapped case.
  d = 4
  pairs = [
    [
      pv.Bicop(family=pv.families.gaussian, parameters=np.array([[r]]))
      for r in row
    ]
    for row in ([0.5, 0.4, 0.3], [0.25, 0.2], [0.15])
  ]
  structure = pv.RVineStructure.from_order(list(range(1, d + 1)))
  ref = pv.Vinecop.from_structure(structure=structure, pair_copulas=pairs)
  # Cast: the conformance is nominal -- ``Bicop.pdf`` takes per-row
  # ``parameters`` where ``BicopLike.pdf`` takes keyword-only ``x``.
  mine = HostedVinecop(cast("list[list[BicopLike[Any]]]", pairs), structure)

  rng = np.random.default_rng(0)
  u = rng.uniform(0.02, 0.98, size=(256, d))
  np.testing.assert_allclose(mine.pdf(u), ref.pdf(u), rtol=1e-10, atol=1e-10)
  np.testing.assert_allclose(
    mine.rosenblatt(u), ref.rosenblatt(u), rtol=1e-10, atol=1e-10
  )
  np.testing.assert_allclose(
    mine.inverse_rosenblatt(u),
    ref.inverse_rosenblatt(u),
    rtol=1e-10,
    atol=1e-10,
  )


@pytest.mark.parametrize(
  "tree_criterion", ["tau", "rho", "hoeffd", "mcor", "cxi", "joe"]
)
def test_select_matches_vinecop_for_every_tree_criterion(
  tree_criterion: str,
) -> None:
  """The port computes the edge weight itself, so it has to accept the same
  measures the compiled selector does -- including Chatterjee's xi, which is
  the only asymmetric one and therefore the one a symmetric shortcut would
  silently get wrong.

  Note the spelling: ``wdm`` takes ``"chatterjee"`` / ``"cxi"`` / ``"xi"``,
  where ``FitControlsVinecop.tree_criterion`` takes only ``"cxi"``.
  """
  d = 5
  u = _correlated_pseudo_obs(1, d)
  mine, _ = VinecopBase.select(
    u, _gaussian_fit_edge, tree_criterion=tree_criterion
  )
  controls = pv.FitControlsVinecop(
    family_set=[pv.families.gaussian],
    tree_criterion=tree_criterion,
    num_threads=1,
  )
  theirs = pv.Vinecop.from_data(u, controls=controls)
  np.testing.assert_array_equal(
    np.asarray(mine.matrix), np.asarray(theirs.structure.matrix)
  )
