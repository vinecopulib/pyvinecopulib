"""Tests for conditional sampling in the array-agnostic vine cascade.

``VinecopBase.sample_conditional`` holds a subset of the variables fixed and
draws the rest from their conditional distribution, and ``reorient`` relabels a
vine so that a chosen subset *is* the sampling-order tail — which is what makes
conditioning on it possible. The contract these tests pin is parity with
``Vinecop``: hosting its pair copulas and its RNG, the array-agnostic sampler
reproduces the compiled one exactly, and the relabeling reproduces the compiled
``reorient`` slot for slot.
"""

import copy
import itertools
import math
from typing import Any

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import NonSimplifiedContext, VinecopBase

from .conftest import GaussianBicop, HostedVinecop, host_vinecop

_D = 5
_N = 500
#: Rotation-sensitive families, so a mistaken flip or a mistaken slot in the
#: relabeling is an O(1) density error rather than a no-op.
_ROTATING = pv.FitControlsVinecop(
  family_set=[pv.families.clayton, pv.families.gumbel], num_threads=1
)
_GAUSSIAN_PAIR = pv.FitControlsBicop(family_set=[pv.families.gaussian])
#: Cross-implementation parity bound. Hosting the compiled pair copulas and the
#: compiled RNG, both samplers walk identical numbers, so on one toolchain they
#: agree to the last bit -- but not across them: x86_64 contracts the
#: multiply-adds that arm64 leaves separate, which CI measured at 4e-15 on the
#: sibling discrete-cascade tests. Orders of magnitude below any real error,
#: which mis-indexes a column and lands at O(0.1).
_PARITY_RTOL = 1e-12
_PARITY_ATOL = 1e-12


def _data(seed: int, d: int = _D, n: int = _N) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return pv.to_pseudo_obs(
    rng.standard_normal((n, d)) @ rng.standard_normal((d, d))
  )


def _fitted(seed: int = 2) -> tuple[pv.Vinecop, np.ndarray]:
  u = _data(seed)
  return pv.Vinecop.from_data(u, controls=_ROTATING), u


def _gaussian_fit_edge(
  tree: int, edge: int, u_e: Any, x_e: Any, var_types: Any = ("c", "c")
) -> Any:
  del tree, edge, x_e
  return pv.Bicop.from_data(
    np.asarray(u_e), controls=_GAUSSIAN_PAIR, var_types=list(var_types)
  )


def _toy_vine(context: Any = None) -> HostedVinecop:
  """A 4-d D-vine of toy conditional Gaussian pairs (no ``flip``)."""
  structure = pv.RVineStructure.from_order([1, 2, 3, 4])
  pairs = [
    [GaussianBicop(base_rho=r) for r in row]
    for row in ([0.5, 0.4, 0.3], [0.25, 0.2], [0.15])
  ]
  return HostedVinecop(pairs, structure, context=context)


# ---------------------------------------------------------------------------
# sample_conditional
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k", [1, 2])
def test_sample_conditional_matches_vinecop(k: int) -> None:
  # The headline: same pair copulas, same base-uniform draw, so the two samplers
  # walk identical numbers and agree to the last bit on one toolchain (see
  # `_PARITY` on why the bound is not exactly zero).
  cop, _ = _fitted()
  mine = host_vinecop(cop)
  u_cond = np.random.RandomState(k).uniform(0.05, 0.95, size=(40, k))
  np.testing.assert_allclose(
    mine.sample_conditional(u_cond, seeds=[1, 2, 3]),
    cop.sample_conditional(u_cond, seeds=[1, 2, 3]),
    rtol=_PARITY_RTOL,
    atol=_PARITY_ATOL,
  )


@pytest.mark.parametrize("k", [1, 2])
def test_sample_conditional_reproduces_the_conditioning_columns(k: int) -> None:
  cop, _ = _fitted()
  mine = host_vinecop(cop)
  order = [int(v) for v in cop.structure.order]
  u_cond = np.random.RandomState(k + 7).uniform(0.05, 0.95, size=(50, k))
  sim = mine.sample_conditional(u_cond, seeds=[1, 2, 3])
  assert sim.shape == (50, _D)
  cond_cols = [v - 1 for v in order[-k:]]
  # Exact for k = 1: natural column d - 1 is never an edge, so its coordinate
  # passes through untouched. The other tail columns round-trip through matched
  # hfunc2 / hinv2 pairs and land at machine precision.
  np.testing.assert_allclose(sim[:, cond_cols], u_cond, rtol=1e-12, atol=1e-12)
  assert np.all((sim >= 0.0) & (sim <= 1.0))
  np.testing.assert_array_equal(
    sim, mine.sample_conditional(u_cond, seeds=[1, 2, 3])
  )


def test_sample_conditional_draws_the_free_variables_conditionally() -> None:
  # The exactness tests above only pin the *fixed* columns; this is the one that
  # would catch a wrong conditional law. Compare the conditional mean against an
  # unconditional sample filtered to the same neighborhood.
  cop, _ = _fitted()
  mine = host_vinecop(cop)
  order = [int(v) for v in cop.structure.order]
  fixed, free = order[-1] - 1, order[0] - 1
  point = 0.8
  conditional = mine.sample_conditional(np.full((4000, 1), point), seeds=[5, 6])
  pool = mine.sample(60000, seeds=[8, 9])
  near = np.abs(pool[:, fixed] - point) < 0.02
  assert near.sum() > 500
  assert abs(conditional[:, free].mean() - pool[near, free].mean()) < 0.05


def test_sample_conditional_rejects_a_non_2d_u_cond() -> None:
  cop, _ = _fitted()
  with pytest.raises(ValueError, match=r"shape \(n, k\)"):
    host_vinecop(cop).sample_conditional(np.full((5,), 0.5))


def test_sample_conditional_rejects_an_unusable_column_count() -> None:
  # k must leave at least one variable free, so d columns match no layout.
  cop, _ = _fitted()
  with pytest.raises(ValueError, match="invalid number of columns"):
    host_vinecop(cop).sample_conditional(np.full((5, _D), 0.5))


def test_sample_conditional_works_on_a_non_simplified_vine() -> None:
  # The order tail is a self-contained sub-vine whatever the pair copulas
  # condition on, so the 0.5 placeholders cannot leak in here either.
  vine = _toy_vine(NonSimplifiedContext())
  u_cond = np.random.RandomState(3).uniform(0.05, 0.95, size=(25, 1))
  sim = vine.sample_conditional(u_cond)
  assert sim.shape == (25, 4)
  np.testing.assert_allclose(
    sim[:, int(vine.order[-1]) - 1], u_cond[:, 0], rtol=1e-12, atol=1e-12
  )


def test_vine_entry_points_require_row_aligned_external_covariates() -> None:
  """The cascade never broadcasts one conditional-design row over a batch."""
  vine = _toy_vine()
  u = np.full((5, 4), 0.5)
  for x in (np.zeros(5), np.zeros((1, 1))):
    calls = (
      lambda: vine.pdf(u, x=x),
      lambda: vine.rosenblatt(u, x=x),
      lambda: vine.inverse_rosenblatt(u, x=x),
      lambda: vine.loglik(u, x=x),
      lambda: vine.sample(5, x=x),
      lambda: vine.sample_conditional(np.full((5, 1), 0.5), x=x),
    )
    for call in calls:
      with pytest.raises(ValueError, match="one row per observation|shape"):
        call()


def test_fixed_structure_fit_requires_row_aligned_covariates() -> None:
  """Conditional pair fitting receives the same validated X rows as evaluation."""
  structure = pv.RVineStructure.from_order([1, 2, 3])
  u = np.full((5, 3), 0.5)
  with pytest.raises(ValueError, match="one row per observation"):
    VinecopBase.fit(
      structure,
      u,
      _gaussian_fit_edge,
      x=np.zeros((1, 1)),
    )


def test_vinecopbase_loglik_preserves_extreme_tail_density() -> None:
  """The hosted likelihood matches the compiled oracle below the old floor."""
  pair = pv.Bicop.from_family(
    pv.families.gaussian, parameters=np.array([[0.99]])
  )
  structure = pv.RVineStructure.from_order([1, 2])
  ref = pv.Vinecop.from_structure(structure=structure, pair_copulas=[[pair]])
  mine = HostedVinecop([[pair]], structure)
  u = np.array([[1e-6, 1 - 1e-6], [1e-5, 1 - 1e-5]])
  np.testing.assert_allclose(mine.loglik(u), ref.loglik(u), rtol=1e-14)


# ---------------------------------------------------------------------------
# Discrete conditioners
# ---------------------------------------------------------------------------


def _discrete_tail_vine() -> tuple[pv.Vinecop, np.ndarray]:
  """A vine whose order tail is a ``Binomial(4, 0.5)`` count variable."""
  rs = np.random.RandomState(11)
  n = 800
  counts = rs.binomial(4, 0.5, size=n)
  cdf = np.cumsum([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0
  cont = pv.to_pseudo_obs(rs.normal(size=(n, 2)))
  data = np.column_stack(
    [cdf[counts], cont, np.where(counts > 0, cdf[counts - 1], 0.0)]
  )
  controls = pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  controls.conditioning_set = [1]
  cop = pv.Vinecop.from_data(data, var_types=["d", "c", "c"], controls=controls)
  assert [int(v) for v in cop.structure.order][-1] == 1
  return cop, cdf


def test_sample_conditional_on_a_discrete_variable() -> None:
  cop, cdf = _discrete_tail_vine()
  mine = host_vinecop(cop, ["d", "c", "c"])
  # Condition variable 1 on the atom x = 2: its value and its left limit.
  u_cond = np.tile([[cdf[2], cdf[1]]], (20, 1))
  sim = mine.sample_conditional(u_cond, seeds=[1, 2, 3])
  np.testing.assert_allclose(
    sim,
    cop.sample_conditional(u_cond, seeds=[1, 2, 3]),
    rtol=_PARITY_RTOL,
    atol=_PARITY_ATOL,
  )
  # A discrete conditioner is reproduced up to its atom, not exactly.
  assert np.all(sim[:, 0] >= cdf[1] - 1e-9)
  assert np.all(sim[:, 0] <= cdf[2] + 1e-9)


def test_sample_conditional_rejects_a_left_limit_above_its_value() -> None:
  cop, cdf = _discrete_tail_vine()
  mine = host_vinecop(cop, ["d", "c", "c"])
  bad = np.tile([[cdf[1], cdf[2]]], (5, 1))
  with pytest.raises(ValueError, match="must not exceed the value columns"):
    mine.sample_conditional(bad)


def test_sample_conditional_needs_the_left_limit_column() -> None:
  # Dropping it leaves the model under-specified: one column matches no
  # admissible layout for this vine.
  cop, cdf = _discrete_tail_vine()
  mine = host_vinecop(cop, ["d", "c", "c"])
  with pytest.raises(ValueError, match="invalid number of columns"):
    mine.sample_conditional(np.full((5, 1), cdf[2]))


def test_sample_conditional_with_multiple_discrete_conditioners_matches_cpp() -> (
  None
):
  """The accepted multi-atom behavior follows the authoritative sampler."""
  rs = np.random.RandomState(31)
  n = 700
  counts = [rs.binomial(4, p, size=n) for p in (0.35, 0.65)]

  def binomial_cdf(p: float) -> np.ndarray:
    masses = np.array(
      [math.comb(4, k) * p**k * (1.0 - p) ** (4 - k) for k in range(5)]
    )
    return np.cumsum(masses)

  cdfs = [binomial_cdf(p) for p in (0.35, 0.65)]
  values = [cdf[col] for cdf, col in zip(cdfs, counts, strict=True)]
  left = [
    np.where(col > 0, cdf[np.maximum(col - 1, 0)], 0.0)
    for cdf, col in zip(cdfs, counts, strict=True)
  ]
  cont = pv.to_pseudo_obs(rs.normal(size=(n, 2)))
  data = np.column_stack([*values, cont, *left])
  controls = pv.FitControlsVinecop(
    family_set=[pv.families.gaussian], num_threads=1
  )
  controls.conditioning_set = [1, 2]
  cop = pv.Vinecop.from_data(
    data, var_types=["d", "d", "c", "c"], controls=controls
  )
  mine = host_vinecop(cop, ["d", "d", "c", "c"])
  tail = [int(v) for v in cop.structure.order[-2:]]
  atoms = {1: 2, 2: 3}
  row_values = [cdfs[v - 1][atoms[v]] for v in tail]
  row_left = [cdfs[v - 1][atoms[v] - 1] for v in tail]
  u_cond = np.tile([*row_values, *row_left], (30, 1))

  got = mine.sample_conditional(u_cond, conditioning_set=tail, seeds=[2, 5, 7])
  expected = cop.sample_conditional(
    u_cond, conditioning_set=tail, seeds=[2, 5, 7]
  )
  np.testing.assert_allclose(
    got, expected, rtol=_PARITY_RTOL, atol=_PARITY_ATOL
  )
  # Do not impose a stronger atom-reproduction invariant than the reference:
  # with several discrete conditioners later-drawn coordinates may sit just
  # outside their atom, while the free-variable law remains correct.


# ---------------------------------------------------------------------------
# reorient and the conditioning_set views
# ---------------------------------------------------------------------------


def test_reorient_is_value_preserving_and_matches_vinecop() -> None:
  # Exhaustive over 1- and 2-variable sets: every admissible one must be
  # value-preserving and reproduce the compiled relabeling's matrix, and at least
  # one must be rejected as inadmissible.
  cop, u = _fitted()
  mine = host_vinecop(cop)
  grid = np.random.default_rng(4).uniform(0.02, 0.98, size=(200, _D))
  pdf0, ll0 = cop.pdf(grid), cop.loglik(u)
  reoriented = rejected = 0
  for k in (1, 2):
    for cand in itertools.combinations(range(1, _D + 1), k):
      cs = list(cand)
      try:
        structure, pairs = mine.reorient(cs)
      except RuntimeError:
        rejected += 1
        continue
      order = [int(v) for v in structure.order]
      assert set(order[-k:]) == set(cs)
      relabeled = HostedVinecop(pairs, structure)
      np.testing.assert_allclose(
        relabeled.pdf(grid), pdf0, rtol=1e-11, atol=1e-11
      )
      np.testing.assert_allclose(
        relabeled.loglik(u), ll0, rtol=1e-11, atol=1e-11
      )
      reference = copy.deepcopy(cop)
      reference.reorient(cs)
      assert np.array_equal(
        np.asarray(structure.matrix),
        np.asarray(reference.structure.matrix),
      )
      reoriented += 1
  assert reoriented, "expected at least one admissible non-suffix tail"
  assert rejected, "expected at least one inadmissible set to be rejected"


def test_conditioning_set_evaluates_through_the_relabeled_vine() -> None:
  # The keyword form must equal a vine rebuilt from `reorient`, equal the
  # compiled view, and leave the model alone.
  cop, u = _fitted()
  mine = host_vinecop(cop)
  order_before = [int(v) for v in mine.structure.order]
  checked = 0
  for k in (1, 2):
    for cand in itertools.combinations(range(1, _D + 1), k):
      cs = list(cand)
      try:
        structure, pairs = mine.reorient(cs)
      except RuntimeError:
        continue
      w = mine.rosenblatt(u, conditioning_set=cs)
      np.testing.assert_allclose(
        w,
        HostedVinecop(pairs, structure).rosenblatt(u),
        rtol=1e-12,
        atol=1e-12,
      )
      np.testing.assert_allclose(
        w, cop.rosenblatt(u, conditioning_set=cs), rtol=1e-10, atol=1e-10
      )
      back = mine.inverse_rosenblatt(w, conditioning_set=cs)
      # The clamp at the unit-box boundary leaves a few points the inverse
      # cannot resolve, so measure the 99th percentile rather than the max.
      assert np.percentile(np.abs(back - u), 99) < 1e-8
      checked += 1
  assert checked
  assert [int(v) for v in mine.structure.order] == order_before


def test_conditioning_set_of_the_order_tail_is_the_identity() -> None:
  cop, u = _fitted()
  mine = host_vinecop(cop)
  tail = [int(v) for v in cop.structure.order][-2:]
  np.testing.assert_array_equal(
    mine.rosenblatt(u, conditioning_set=tail), mine.rosenblatt(u)
  )
  w = mine.rosenblatt(u)
  np.testing.assert_array_equal(
    mine.inverse_rosenblatt(w, conditioning_set=tail),
    mine.inverse_rosenblatt(w),
  )
  u_cond = np.random.RandomState(2).uniform(0.05, 0.95, size=(30, 2))
  np.testing.assert_array_equal(
    mine.sample_conditional(u_cond, conditioning_set=tail, seeds=[4]),
    mine.sample_conditional(u_cond, seeds=[4]),
  )


def test_conditioning_set_column_mapping_follows_the_set() -> None:
  # With the argument, column i is conditioning_set[i] -- not the i-th variable
  # of the order tail. Reversing the set permutes what is consumed.
  cop, _ = _fitted()
  mine = host_vinecop(cop)
  tail = [int(v) for v in cop.structure.order][-2:]
  u_cond = np.array([[0.2, 0.8]] * 12)
  forward = mine.sample_conditional(u_cond, conditioning_set=tail, seeds=[3])
  backward = mine.sample_conditional(
    u_cond, conditioning_set=tail[::-1], seeds=[3]
  )
  np.testing.assert_allclose(
    forward[:, [v - 1 for v in tail]], u_cond, rtol=1e-10, atol=1e-10
  )
  np.testing.assert_allclose(
    backward[:, [v - 1 for v in tail[::-1]]], u_cond, rtol=1e-10, atol=1e-10
  )


@pytest.mark.parametrize("trunc_lvl", [0, 1, 2])
def test_reorient_handles_a_truncated_vine(trunc_lvl: int) -> None:
  """A truncated vine relabels like any other.

  The trees above the truncation are independence, so the peel has nothing to
  move there and the slot map covers only the trees the structure stores. The
  reference is the compiled in-place `reorient`, evaluated directly.
  """
  controls = pv.FitControlsVinecop(
    family_set=[pv.families.gaussian], trunc_lvl=trunc_lvl, num_threads=1
  )
  u = _data(6)
  cop = pv.Vinecop.from_data(u, controls=controls)
  # The tail of the fitted order is admissible by construction, whatever the
  # truncation, so the case does not depend on which structure was selected.
  cond = [int(v) for v in cop.structure.order[-2:]][::-1]

  structure, pairs = host_vinecop(cop).reorient(cond)
  assert len(pairs) == int(structure.trunc_lvl) == trunc_lvl

  ref = pv.Vinecop.from_json(cop.to_json())
  ref.reorient(cond)
  assert list(structure.order) == list(ref.structure.order)
  np.testing.assert_allclose(
    host_vinecop(cop).rosenblatt(u, conditioning_set=cond),
    ref.rosenblatt(u),
    rtol=1e-10,
    atol=1e-10,
  )


@pytest.mark.parametrize(
  "bad", [[], [1, 1], [0], [99], [1, 2, 3, 4, 5]], ids=str
)
def test_reorient_rejects_inadmissible_sets(bad: list[int]) -> None:
  cop, _ = _fitted()
  with pytest.raises(RuntimeError):
    host_vinecop(cop).reorient(bad)


def test_reorient_requires_flip_on_the_pairs() -> None:
  # `GaussianBicop` (conftest) leaves `flip` at the raising default: it is only
  # needed to reuse a pair in selection, or to move it in a relabeling.
  with pytest.raises(NotImplementedError, match="flip"):
    _toy_vine().reorient([1])


def test_conditioning_set_is_refused_on_a_non_simplified_vine() -> None:
  vine = _toy_vine(NonSimplifiedContext())
  u = np.random.default_rng(5).uniform(0.05, 0.95, size=(20, 4))
  with pytest.raises(NotImplementedError, match="non-simplified vine"):
    vine.rosenblatt(u, conditioning_set=[1])
  # The tail is still fine: it needs no relabeling at all.
  tail = [int(vine.order[-1])]
  np.testing.assert_array_equal(
    vine.rosenblatt(u, conditioning_set=tail), vine.rosenblatt(u)
  )


# ---------------------------------------------------------------------------
# Conditioning-aware selection
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cs", [[2], [2, 3], [1, 4], [1, 3, 5]], ids=str)
def test_select_conditioning_set_matches_vinecop(cs: list[int]) -> None:
  # The penalty and the relabeling together must reproduce the compiled
  # conditioning-aware selector's matrix exactly, and the reused pairs must land
  # on the same slots -- an O(1) density error otherwise.
  u = _data(11)
  structure, pairs = VinecopBase.select(
    u, _gaussian_fit_edge, conditioning_set=cs
  )
  controls = pv.FitControlsVinecop(
    family_set=[pv.families.gaussian], num_threads=1
  )
  controls.conditioning_set = cs
  reference = pv.Vinecop.from_data(u, controls=controls)
  order = [int(v) for v in structure.order]
  assert set(order[-len(cs) :]) == set(cs)
  assert np.array_equal(
    np.asarray(structure.matrix), np.asarray(reference.structure.matrix)
  )
  grid = np.random.default_rng(12).uniform(0.02, 0.98, size=(200, _D))
  np.testing.assert_allclose(
    HostedVinecop(pairs, structure).pdf(grid),
    reference.pdf(grid),
    rtol=1e-12,
    atol=1e-12,
  )


@pytest.mark.parametrize(
  "kwargs,match",
  [
    ({"conditioning_set": [1, 2, 3, 4, 5]}, "at most d - 1"),
    ({"conditioning_set": [0]}, r"in 1, \.\.\., d"),
    ({"conditioning_set": [99]}, r"in 1, \.\.\., d"),
    (
      {"conditioning_set": [1], "tree_algorithm": "random_weighted"},
      "requires an MST",
    ),
  ],
  ids=["too-many", "zero", "out-of-range", "non-mst"],
)
def test_select_validates_the_conditioning_set(
  kwargs: dict[str, Any], match: str
) -> None:
  with pytest.raises(ValueError, match=match):
    VinecopBase.select(_data(13), _gaussian_fit_edge, **kwargs)


@pytest.mark.parametrize("trunc_lvl", [1, 2])
def test_select_combines_a_conditioning_set_with_truncation(
  trunc_lvl: int,
) -> None:
  """Conditioning-aware selection accepts a truncated model.

  The penalty lays the conditioning block down first in whatever trees are
  selected; the relabeling then only has to reach those, so truncation is not
  an obstacle. The selected order must end with the requested set, and the
  pair store must hold exactly the selected trees.
  """
  cond = [2, 5]
  structure, pairs = VinecopBase.select(
    _data(13),
    _gaussian_fit_edge,
    conditioning_set=cond,
    trunc_lvl=trunc_lvl,
  )
  assert int(structure.trunc_lvl) == trunc_lvl
  assert len(pairs) == trunc_lvl
  assert set(int(v) for v in structure.order[-2:]) == set(cond)
