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
from pyvinecopulib.core import BicopBase, BicopLike, NonSimplifiedContext

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
class _Controls:
  """A minimal ``ControlsLike``: the engines read settings through ``to_dict``.

  Used rather than ``FitControlsVinecop`` so a test can hand the engine a
  value the compiled controls would reject at construction, and so each test
  sets only the knobs it cares about.
  """

  def __init__(self, **settings: Any) -> None:
    self._settings = settings

  def to_dict(self) -> dict:
    return dict(self._settings)


def _base_controls(**settings: Any) -> _Controls:
  return _Controls(**settings)


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
      mine = HostedVinecop.from_data(
        u, fit_edge=_gaussian_fit_edge, controls=_base_controls(trunc_lvl=trunc)
      ).structure
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
    mine = HostedVinecop.from_data(u, fit_edge=_tll_fit_edge)
    auto = pv.Vinecop.from_data(u, controls=_vine_controls(pv.families.tll))
    assert np.array_equal(
      np.asarray(mine.structure.matrix), np.asarray(auto.structure.matrix)
    )
    rng = np.random.default_rng(seed)
    grid = rng.uniform(0.02, 0.98, size=(300, d))
    np.testing.assert_allclose(
      mine.pdf(grid), auto.pdf(grid), rtol=1e-12, atol=1e-12
    )


def test_select_refit_pdf_parity_with_vinecop() -> None:
  # Refitting on the selected structure also reproduces the density of the
  # auto-selected vine for orientation-symmetric (Gaussian) pairs.
  rng = np.random.default_rng(123)
  for seed in range(3):
    d = 5
    u = _correlated_pseudo_obs(seed, d)
    mine = HostedVinecop.from_data(u, fit_edge=_gaussian_fit_edge).structure
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
    vine = HostedVinecop.from_data(
      u, fit_edge=_gaussian_fit_edge, controls=_base_controls(trunc_lvl=trunc)
    )
    assert vine.structure.dim == d
    assert vine.structure.trunc_lvl == trunc
    # Every slot the truncated structure declares is populated, and none above.
    assert [
      len([e for e in range(d - 1 - t) if vine.get_pair_copula(t, e)])
      for t in range(trunc)
    ] == [d - 1 - t for t in range(trunc)]


def test_select_random_weighted_is_reproducible() -> None:
  d = 5
  u = _correlated_pseudo_obs(1, d)
  controls = _base_controls(tree_algorithm="random_weighted", seeds=[1, 2, 3])
  s1 = HostedVinecop.from_data(
    u, fit_edge=_gaussian_fit_edge, controls=controls
  ).structure
  s2 = HostedVinecop.from_data(
    u, fit_edge=_gaussian_fit_edge, controls=controls
  ).structure
  # Same seeds -> identical structure; and a valid full R-vine.
  assert np.array_equal(np.asarray(s1.matrix), np.asarray(s2.matrix))
  assert s1.trunc_lvl == d - 1


def test_select_random_weighted_accepts_small_samples() -> None:
  u = _correlated_pseudo_obs(1, 5, n=10)
  structure = HostedVinecop.from_data(
    u,
    fit_edge=_gaussian_fit_edge,
    controls=_base_controls(tree_algorithm="random_weighted", seeds=[1, 2, 3]),
  ).structure
  assert structure.trunc_lvl == 4


def test_select_rejects_unknown_tree_algorithm() -> None:
  u = _correlated_pseudo_obs(1, 4, n=50)
  with pytest.raises(ValueError, match="tree_algorithm must be one of"):
    HostedVinecop.from_data(
      u,
      fit_edge=_gaussian_fit_edge,
      controls=_base_controls(tree_algorithm="definitely_not_an_algorithm"),
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

  HostedVinecop.from_data(
    _correlated_pseudo_obs(1, 6, n=50),
    fit_edge=fit_edge,
    controls=_base_controls(trunc_lvl=1),
  )
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
  mine = HostedVinecop.from_data(
    u,
    fit_edge=_gaussian_fit_edge,
    controls=_base_controls(tree_criterion=tree_criterion),
  ).structure
  controls = pv.FitControlsVinecop(
    family_set=[pv.families.gaussian],
    tree_criterion=tree_criterion,
    num_threads=1,
  )
  theirs = pv.Vinecop.from_data(u, controls=controls)
  np.testing.assert_array_equal(
    np.asarray(mine.matrix), np.asarray(theirs.structure.matrix)
  )


# ---------------------------------------------------------------------------
# Declared parts: bicop_class, the up-front flip gate, and weighted selection
# ---------------------------------------------------------------------------


class _BicopVine(HostedVinecop):
  """A vine that names the pair copula it fits, so `from_data` needs no
  callback."""

  bicop_class = pv.Bicop


def test_bicop_class_fits_without_a_fit_edge_callback() -> None:
  # Naming the pair class is what replaces the callback: a pair class is
  # itself a fitter, so the vine can fit its own edges.
  u = _correlated_pseudo_obs(0, 5)
  vine = _BicopVine.from_data(u)
  assert vine.dim == 5
  assert np.all(np.isfinite(vine.pdf(u)))
  assert isinstance(vine.get_pair_copula(0, 0), pv.Bicop)


def test_a_named_pair_class_without_flip_is_refused_before_any_fitting() -> (
  None
):
  # Selection finalizes by reorienting pairs, so a pair without `flip` cannot
  # be selected with. Naming the class makes that knowable before the data is
  # touched -- the point of `bicop_class` over an opaque callback.
  calls = {"n": 0}

  class _NoFlipPair(pv.core.BicopBase[Any]):
    def __init__(self) -> None:
      calls["n"] += 1

    def pdf(self, u: Any, *, x: Any = None) -> Any:
      return np.ones(u.shape[0])

    def hfunc1(self, u: Any, *, x: Any = None) -> Any:
      return u[:, 1]

    def hfunc2(self, u: Any, *, x: Any = None) -> Any:
      return u[:, 0]

  class _NoFlipVine(HostedVinecop):
    bicop_class = _NoFlipPair

  with pytest.raises(NotImplementedError, match="has no `flip`"):
    _NoFlipVine.from_data(_correlated_pseudo_obs(0, 5))
  # Nothing was constructed, so nothing was fitted.
  assert calls["n"] == 0


def _weak_link_chain(seed: int, n: int = 1500) -> np.ndarray:
  """A Markov chain whose weakest link joins the two lowest-numbered columns.

  The spanning tree has to select that edge to stay connected, and it sorts
  first by candidate index -- so a threshold above it makes the *first* edge of
  tree 0 the thresholded one.
  """
  rng = np.random.default_rng(seed)
  cols = [rng.normal(size=n)]
  for r in (0.25, 0.95, 0.95):
    cols.append(r * cols[-1] + np.sqrt(1 - r * r) * rng.normal(size=n))
  return pv.to_pseudo_obs(np.column_stack(cols))


def test_a_thresholded_edge_does_not_consume_the_flip_probe() -> None:
  """A pair without `flip` must be caught behind a caller's own `fit_edge` too.

  The probe fires once. A thresholded edge holds an `IndependencePair`, whose
  `flip` returns `self`, so letting one consume the probe passes it on behalf
  of a pair that cannot flip -- which then fails from the finalizing
  reorientation, after every other edge has been fitted.
  """
  u = _weak_link_chain(0)
  controls = _vine_controls(pv.families.gaussian)
  controls.threshold = 0.22

  # The setup this test needs, asserted so it cannot quietly go vacuous: in the
  # selection loop's own order the first edge of tree 0 is thresholded, so the
  # first pair the probe can legitimately see is a later, fitted one.
  seen: list[int] = []

  def recording_fit_edge(tree: int, edge: int, u_e: Any, x_e: Any) -> Any:
    if tree == 0:
      seen.append(edge)
    return pv.Bicop.from_data(
      u_e, pv.FitControlsBicop(family_set=[pv.families.gaussian])
    )

  HostedVinecop.from_data(u, controls=controls, fit_edge=recording_fit_edge)
  assert seen and min(seen) > 0

  class _NoFlipPair(pv.core.BicopBase[Any]):
    def pdf(self, u: Any, *, x: Any = None) -> Any:
      return np.ones(u.shape[0])

    def hfunc1(self, u: Any, *, x: Any = None) -> Any:
      return u[:, 1]

    def hfunc2(self, u: Any, *, x: Any = None) -> Any:
      return u[:, 0]

  def fit_edge(tree: int, edge: int, u_e: Any, x_e: Any) -> Any:
    return _NoFlipPair()

  with pytest.raises(NotImplementedError, match="has no `flip`"):
    HostedVinecop.from_data(u, controls=controls, fit_edge=fit_edge)


@pytest.mark.parametrize(
  "setting,value",
  [
    ("select_trunc_lvl", True),
    ("select_threshold", True),
    ("select_families", False),
    ("show_trace", True),
  ],
)
def test_vine_level_settings_it_cannot_honor_are_refused(
  setting: str, value: bool
) -> None:
  """Refused, not dropped -- the `ControlsLike` contract.

  These four are on `FitControlsVinecop` only, so no pair-copula fit can read
  them either. Dropping one returns a model the caller's controls do not
  describe, and silently disagrees with `Vinecop.select` on the same object.
  """
  controls = _vine_controls(pv.families.gaussian)
  setattr(controls, setting, value)
  with pytest.raises(ValueError, match=setting):
    _BicopVine.from_data(_correlated_pseudo_obs(0, 4), controls=controls)


def test_the_defaults_of_those_settings_still_select() -> None:
  """The refusal keys off a non-default value, not off the field existing."""
  vine = _BicopVine.from_data(
    _correlated_pseudo_obs(0, 4),
    controls=_vine_controls(pv.families.gaussian),
  )
  assert vine.dim == 4


@pytest.mark.parametrize("seed,d", [(0, 5), (1, 6), (2, 7)])
def test_weighted_selection_matches_vinecop(seed: int, d: int) -> None:
  # The tree criterion is weighted now, so a weighted array-agnostic selection
  # reproduces `Vinecop.select`'s matrix exactly. Before the weights reached
  # the criterion this silently disagreed.
  rng = np.random.default_rng(seed)
  u = _correlated_pseudo_obs(seed, d, n=800)
  controls = _vine_controls(pv.families.gaussian)
  controls.weights = rng.uniform(0.2, 2.0, u.shape[0])

  mine = HostedVinecop.from_data(
    u, fit_edge=_gaussian_fit_edge, controls=controls
  ).structure
  theirs = pv.Vinecop.from_data(u, controls=controls).structure
  assert np.array_equal(np.asarray(mine.matrix), np.asarray(theirs.matrix))


def test_unweighted_selection_is_unchanged_by_the_weights_plumbing() -> None:
  # The criterion's new weights argument defaults to empty, which is what the
  # call site passed before, so every unweighted selection is untouched.
  u = _correlated_pseudo_obs(3, 6)
  plain = HostedVinecop.from_data(u, fit_edge=_gaussian_fit_edge).structure
  unit = _vine_controls(pv.families.gaussian)
  unit.weights = np.ones(u.shape[0])
  weighted = HostedVinecop.from_data(
    u, fit_edge=_gaussian_fit_edge, controls=unit
  ).structure
  assert np.array_equal(np.asarray(plain.matrix), np.asarray(weighted.matrix))


@pytest.mark.parametrize(
  ("bad", "match"),
  [
    ("all_nan", "at least one observation standing"),
    ("all_zero", "at least one observation standing"),
    ("inf", "must not contain infinite"),
    ("negative", "nonnegative"),
    ("wrong_length", "one weight per observation"),
  ],
)
def test_selection_refuses_unusable_weights(bad: str, match: str) -> None:
  # The criterion hands the weights straight to the binding, so the engines
  # validate them with the same validator every other layer uses -- they used
  # to go through unchecked.
  u = _correlated_pseudo_obs(0, 4, n=200)
  n = u.shape[0]
  weights = {
    "all_nan": np.full(n, np.nan),
    "all_zero": np.zeros(n),
    "inf": np.where(np.arange(n) == 3, np.inf, 1.0),
    "negative": np.where(np.arange(n) == 3, -1.0, 1.0),
    "wrong_length": np.ones(n - 1),
  }[bad]
  controls = _Controls(tree_criterion="tau", weights=weights)
  with pytest.raises(ValueError, match=match):
    HostedVinecop.from_data(u, fit_edge=_gaussian_fit_edge, controls=controls)


def test_selection_accepts_the_drop_markers_wdm_honors() -> None:
  # `wdm` removes missing rows by default and gives a NaN weight and a zero
  # weight the same answer, so the validator must let both through rather than
  # refusing a weighting the criterion below it handles.
  u = _correlated_pseudo_obs(0, 4, n=200)
  n = u.shape[0]
  as_nan = _Controls(
    tree_criterion="tau",
    weights=np.where(np.arange(n) % 2 == 0, 1.0, np.nan),
  )
  as_zero = _Controls(
    tree_criterion="tau",
    weights=np.where(np.arange(n) % 2 == 0, 1.0, 0.0),
  )
  dropped = HostedVinecop.from_data(
    u, fit_edge=_gaussian_fit_edge, controls=as_nan
  ).structure
  zeroed = HostedVinecop.from_data(
    u, fit_edge=_gaussian_fit_edge, controls=as_zero
  ).structure
  assert np.array_equal(np.asarray(dropped.matrix), np.asarray(zeroed.matrix))


def test_custom_tree_criterion_is_callable_instead_of_raising() -> None:
  # `tree_criterion="custom"` had nothing to call, so it raised from C++.
  u = _correlated_pseudo_obs(1, 5)
  seen = {"n": 0}

  def criterion(pair: Any) -> float:
    seen["n"] += 1
    return float(abs(np.corrcoef(pair, rowvar=False)[0, 1]))

  controls = _base_controls(
    tree_criterion="custom", criterion_function=criterion
  )
  vine = HostedVinecop.from_data(
    u, fit_edge=_gaussian_fit_edge, controls=controls
  )
  assert seen["n"] > 0
  assert vine.structure.dim == 5


# --------------------------------------------------------------------------- #
# Conditional pairs: the null hypothesis, where the x-dependence is switched   #
# off and a non-simplified vine must reduce to the simplified one exactly.     #
# --------------------------------------------------------------------------- #


class _ConditionalGaussian(BicopBase[np.ndarray]):
  """A fitted Gaussian pair wrapped in a conditional shell.

  The correlation is pushed through ``tanh`` in *z* space around whatever the
  ordinary Gaussian fit estimated::

      rho(x) = rho_max * tanh(atanh(rho_hat / rho_max) + slope * g(x))

  At ``slope == 0`` that is ``rho_hat`` and the shell evaluates through the
  scalar-parameter path, so it *is* the compiled Gaussian pair. That is what
  makes it a null hypothesis: a non-simplified vine built from these must
  reduce exactly to the simplified one, and structure selection must recover
  what the compiled selector recovers. With ``slope != 0`` the same class is
  genuinely conditional, so any difference is attributable to the covariates
  and not to the wrapper.

  The per-row correlation reaches the density through ``Bicop``'s own per-row
  ``parameters`` argument, rather than a hand-written Gaussian, so the two legs
  differ in exactly one thing.
  """

  supports_covariates = True

  def __init__(
    self,
    *,
    slope: float = 0.0,
    rho_max: float = 0.95,
    bicop: Optional[pv.Bicop] = None,
  ) -> None:
    self._slope = float(slope)
    self._rho_max = float(rho_max)
    self._bicop = bicop

  def fit(
    self,
    u: np.ndarray,
    /,
    controls: Any = None,
    var_types: Optional[list[str]] = None,
    *,
    x: Optional[np.ndarray] = None,
  ) -> "_ConditionalGaussian":
    del controls, var_types, x
    self._bicop = pv.Bicop.from_data(np.asarray(u), controls=_GAUSSIAN)
    return self

  def _per_row(self, x: Optional[np.ndarray]) -> Optional[np.ndarray]:
    """Per-row correlations, or ``None`` to take the scalar path."""
    if x is None or self._slope == 0.0:
      return None
    assert self._bicop is not None
    rho = float(np.asarray(self._bicop.parameters).ravel()[0])
    xa = np.asarray(x, dtype=float)
    weights = np.arange(1, xa.shape[1] + 1, dtype=float)
    shift = self._slope * (xa * weights).sum(axis=-1) / xa.shape[1]
    # A fitted correlation can sit outside `rho_max`, which would make the
    # inverse tanh undefined; the anchor only has to be *some* z whose tanh
    # recovers `rho` when the shift is zero.
    anchor = np.arctanh(np.clip(rho / self._rho_max, -1 + 1e-12, 1 - 1e-12))
    z = anchor + shift
    return (self._rho_max * np.tanh(z)).reshape(-1, 1)

  def _call(self, name: str, u: np.ndarray, x: Optional[np.ndarray]) -> Any:
    assert self._bicop is not None, "fit the pair first"
    method = getattr(self._bicop, name)
    per_row = self._per_row(x)
    return method(u) if per_row is None else method(u, per_row)

  def pdf(self, u: np.ndarray, *, x: Optional[np.ndarray] = None) -> Any:
    return self._call("pdf", u, x)

  def hfunc1(self, u: np.ndarray, *, x: Optional[np.ndarray] = None) -> Any:
    return self._call("hfunc1", u, x)

  def hfunc2(self, u: np.ndarray, *, x: Optional[np.ndarray] = None) -> Any:
    return self._call("hfunc2", u, x)

  def flip(self) -> "_ConditionalGaussian":
    assert self._bicop is not None
    return _ConditionalGaussian(
      slope=self._slope, rho_max=self._rho_max, bicop=self._bicop.flip()
    )


def _conditional_pairs(d: int, slope: float) -> list[list[Any]]:
  return [
    [_ConditionalGaussian(slope=slope) for _ in range(d - 1 - t)]
    for t in range(d - 1)
  ]


@pytest.mark.parametrize(("seed", "d"), [(0, 4), (1, 5), (2, 6)])
def test_a_conditional_pair_with_no_x_dependence_recovers_the_cpp_structure(
  seed: int, d: int
) -> None:
  """The null hypothesis: switch the covariate term off and nothing changes.

  With ``slope == 0`` the shell is the compiled Gaussian pair, so selecting
  through a *non-simplified* vine must recover byte-for-byte what
  ``Vinecop.select`` recovers, and evaluate to the same density. This is the
  leg that isolates the wrapper, the placement and the finalizing ``flip``
  from the conditioning itself -- it is insensitive to the C1 column order,
  because the pair never reads a column.
  """
  u = _correlated_pseudo_obs(seed, d, n=600)

  mine = HostedVinecop(
    _conditional_pairs(d, slope=0.0),
    pv.RVineStructure.from_order(list(range(1, d + 1))),
    context=NonSimplifiedContext(),
  )
  mine.select(
    u,
    fit_edge=lambda t, e, u_e, x_e, var_types=("c", "c"): (
      _ConditionalGaussian(slope=0.0).fit(u_e)
    ),
  )
  theirs = pv.Vinecop.from_data(
    u, controls=_vine_controls(pv.families.gaussian)
  )

  assert np.array_equal(
    np.asarray(mine.structure.matrix), np.asarray(theirs.structure.matrix)
  )
  assert list(mine.structure.order) == list(theirs.structure.order)

  grid = np.random.default_rng(seed + 99).uniform(0.02, 0.98, size=(300, d))
  np.testing.assert_allclose(
    mine.pdf(grid), theirs.pdf(grid), rtol=1e-12, atol=1e-12
  )


def test_the_same_pair_with_x_dependence_on_is_a_different_model() -> None:
  """The control leg: with ``slope != 0`` the covariates must actually bite.

  Otherwise the null-hypothesis test above would pass for the trivial reason
  that the shell ignores ``x`` in both legs.
  """
  d, seed = 5, 3
  u = _correlated_pseudo_obs(seed, d, n=600)
  structure = pv.RVineStructure.from_order(list(range(1, d + 1)))
  grid = np.random.default_rng(7).uniform(0.02, 0.98, size=(200, d))

  fitted = [
    HostedVinecop(
      _conditional_pairs(d, slope=slope),
      structure,
      context=NonSimplifiedContext(),
    )
    for slope in (0.0, 1.5)
  ]
  for vine, slope in zip(fitted, (0.0, 1.5)):
    vine.fit(
      u,
      fit_edge=lambda t, e, u_e, x_e, var_types=("c", "c"), s=slope: (
        _ConditionalGaussian(slope=s).fit(u_e)
      ),
    )

  flat, conditional = (vine.pdf(grid) for vine in fitted)
  assert not np.allclose(flat, conditional, rtol=1e-6, atol=1e-6)


@pytest.mark.xfail(
  strict=True,
  reason="`select` fits the pairs with x_e=None and then re-installs the "
  "conditioning context, so the vine evaluates as a model it was never "
  "fitted as. Fixing it needs the C1 column order, which is only fixed once "
  "the structure is -- resolved at finalization from the same `diag` that "
  "already decides the `flip`. Remove this marker when that lands.",
)
def test_select_fits_the_model_it_then_evaluates() -> None:
  """A selected non-simplified vine must be fitted as what it evaluates as.

  With the covariate term on, `select` hands every `fit_edge` `x_e = None` and
  then re-installs the context, so the pairs are estimated as a simplified vine
  and read as a conditional one. The log-likelihood the object reports is then
  not the one its own fit maximized.
  """
  d, seed, slope = 4, 0, 1.5
  u = _correlated_pseudo_obs(seed, d, n=600)
  structure = pv.RVineStructure.from_order(list(range(1, d + 1)))

  def fit_edge(t, e, u_e, x_e, var_types=("c", "c")):
    return _ConditionalGaussian(slope=slope).fit(u_e)

  selected = HostedVinecop(
    _conditional_pairs(d, slope=slope),
    structure,
    context=NonSimplifiedContext(),
  )
  selected.select(u, fit_edge=fit_edge)

  # The same pairs, read as the simplified vine they were actually fitted as.
  as_fitted = HostedVinecop(
    [list(row) for row in selected._pairs], selected.structure
  )
  np.testing.assert_allclose(
    float(selected.loglik(u)), float(as_fitted.loglik(u)), rtol=1e-9
  )
