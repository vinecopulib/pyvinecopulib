"""Tests for discrete data layouts in the array-agnostic vine cascade.

``VinecopBase`` carries a parallel left-limit cascade so a vine with discrete
variables evaluates exactly as ``Vinecop`` does. The contract these tests pin is
parity of ``pdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` on
mixed continuous / discrete data, in both the expanded ``(n, 2d)`` and the
compact ``(n, d + k)`` layout — the same pair copulas run through both
evaluators, so any difference is the cascade's.
"""

from typing import Any, Optional, cast

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import BicopLike, Vinedist, VinecopBase
from pyvinecopulib.core import Kde1d

from .conftest import GaussianBicop

#: Cross-language parity tolerance. The cascade and the compiled `Vinecop`
#: multiply the same per-edge densities in the same order, so on one toolchain
#: they agree bit-for-bit -- but not across toolchains: FMA contraction and libm
#: differ, and CI measured up to 2.7e-15 absolute (1.2e-15 relative) on
#: manylinux, musllinux, macOS arm64 and Windows. This bound is ~400x that and
#: still orders of magnitude below any real cascade error, which mis-indexes a
#: column or an h-function and lands at O(0.1).
_PARITY_RTOL = 1e-12
_PARITY_ATOL = 1e-12


def _assert_parity(actual: Any, desired: Any) -> None:
  """Assert agreement with the compiled reference to the parity bounds.

  Parameters
  ----------
  actual : array
      The array-agnostic cascade's output.
  desired : array
      The compiled ``Vinecop``'s output.

  Returns
  -------
  None
  """
  np.testing.assert_allclose(
    actual, desired, rtol=_PARITY_RTOL, atol=_PARITY_ATOL
  )


_D = 4
# Correlations for the D-vine's three trees, deliberately mixed in sign so a
# swapped pair argument shows up as an O(1) density error.
_RHOS = ([0.6, -0.4, 0.3], [0.25, 0.2], [0.15])
#: Cumulative-mass boundaries of the ordered-categorical variables used for
#: ``"d"``: five atoms, and deliberately strictly inside the unit interval so
#: that ``VinecopBase``'s input clamp is a no-op and the comparison against
#: ``Vinecop`` -- which instead trims inside each pair copula -- isolates the
#: cascade.
_BOUNDS = np.array([0.04, 0.14, 0.39, 0.69, 0.89, 0.96])

# All four regimes the cascade has to cover: nothing discrete (which must stay
# byte-identical to the pre-existing continuous path), one of four, an
# interleaved pair, and everything discrete.
_VAR_TYPES = [
  ["c", "c", "c", "c"],
  ["d", "c", "c", "c"],
  ["c", "d", "c", "d"],
  ["d", "d", "d", "d"],
]
_MIXED = _VAR_TYPES[1:]
_LAYOUTS = ["expanded", "compact"]


class _ListVinecop(VinecopBase[Any]):
  """Minimal ``VinecopBase`` hosting a nested list of pair copulas (NumPy).

  Declares ``var_types`` and then stamps each hosted pair with the types the
  base class derives for its slot, which is what a subclass owning its pairs is
  expected to do. ``_sample_uniform`` matches the C++ draw so the randomized
  Rosenblatt path is comparable.
  """

  def __init__(
    self,
    pairs: list[list[pv.Bicop]],
    structure: pv.RVineStructure,
    var_types: Optional[list[str]] = None,
  ) -> None:
    self._pairs = pairs
    self._bind_vine(structure, var_types=var_types)
    for tree, row in enumerate(pairs):
      for edge, pair in enumerate(row):
        pair.var_types = list(self.pair_var_types(tree, edge))

  def _get_pair_copula(self, tree: int, edge: int) -> BicopLike[Any]:
    # Cast: the conformance is nominal -- `Bicop.pdf` takes per-row
    # `parameters` where `BicopLike.pdf` takes keyword-only `x`.
    return cast("BicopLike[Any]", self._pairs[tree][edge])

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> np.ndarray:
    return pv.utils.sample_uniform(n, self.d, qrng, seeds)


class _NeverBatchedVinecop(_ListVinecop):
  """A vine that advertises the batched fast path but must never be asked."""

  def _default_batched(self) -> bool:
    return True

  def _build_batched(self) -> Any:
    raise AssertionError("the batched fast path was entered")


def _gaussian_pairs() -> list[list[pv.Bicop]]:
  return [
    [
      pv.Bicop(family=pv.families.gaussian, parameters=np.array([[r]]))
      for r in row
    ]
    for row in _RHOS
  ]


def _expanded_data(var_types: list[str], seed: int, n: int = 500) -> np.ndarray:
  """An ``(n, 2d)`` sample: ordered categories at ``"d"``, uniforms at ``"c"``.

  A discrete variable contributes the pair ``(F(x), F(x^-))`` read off the
  cumulative atom masses; a continuous one has coinciding limits.
  """
  rng = np.random.default_rng(seed)
  values, limits = [], {}
  for j, t in enumerate(var_types):
    if t == "d":
      level = rng.integers(0, len(_BOUNDS) - 1, n)
      values.append(_BOUNDS[level + 1])
      limits[j] = _BOUNDS[level]
    else:
      values.append(rng.uniform(0.01, 0.99, n))
  second = [limits.get(j, values[j]) for j in range(len(var_types))]
  return np.column_stack(values + second)


def _to_compact(u: np.ndarray, var_types: list[str]) -> np.ndarray:
  d = len(var_types)
  keep = [d + j for j, t in enumerate(var_types) if t == "d"]
  return np.column_stack([u[:, :d]] + ([u[:, keep]] if keep else []))


def _both(
  var_types: list[str], cls: type[_ListVinecop] = _ListVinecop
) -> tuple[_ListVinecop, pv.Vinecop]:
  """One structure and one set of pair copulas, hosted in both evaluators."""
  structure = pv.RVineStructure.from_order(list(range(1, len(var_types) + 1)))
  ref = pv.Vinecop.from_structure(
    structure=structure, pair_copulas=_gaussian_pairs(), var_types=var_types
  )
  return cls(_gaussian_pairs(), structure, var_types=var_types), ref


# ---------------------------------------------------------------------------
# The per-edge type table
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("var_types", _VAR_TYPES)
def test_pair_var_types_match_vinecop(var_types: list[str]) -> None:
  # The table is derived from the structure, so it must agree edge-for-edge with
  # what Vinecop assigns to the same slots -- everything downstream keys off it.
  mine, ref = _both(var_types)
  assert mine.var_types == var_types
  for tree in range(mine.trunc_lvl):
    for edge in range(mine.d - tree - 1):
      assert list(mine.pair_var_types(tree, edge)) == list(
        ref.get_pair_copula(tree, edge).var_types
      )


def test_bind_vine_validates_var_types() -> None:
  structure = pv.RVineStructure.from_order(list(range(1, _D + 1)))
  with pytest.raises(ValueError, match="2 entries, expected 4"):
    _ListVinecop(_gaussian_pairs(), structure, var_types=["c", "d"])
  with pytest.raises(ValueError, match="'c' or 'd'"):
    _ListVinecop(_gaussian_pairs(), structure, var_types=["c", "z", "c", "c"])


# ---------------------------------------------------------------------------
# Cascade parity with Vinecop
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("layout", _LAYOUTS)
@pytest.mark.parametrize("var_types", _VAR_TYPES)
def test_pdf_matches_vinecop(var_types: list[str], layout: str) -> None:
  # Both evaluators call the same compiled Bicop on the same inputs in the same
  # order, so the products agree bit-for-bit; anything else means the cascade
  # assembled a different edge input.
  mine, ref = _both(var_types)
  u = _expanded_data(var_types, seed=1)
  arg = u if layout == "expanded" else _to_compact(u, var_types)
  _assert_parity(mine.pdf(arg), ref.pdf(arg))


@pytest.mark.parametrize("layout", _LAYOUTS)
@pytest.mark.parametrize("var_types", _VAR_TYPES)
def test_rosenblatt_matches_vinecop(var_types: list[str], layout: str) -> None:
  mine, ref = _both(var_types)
  u = _expanded_data(var_types, seed=2)
  arg = u if layout == "expanded" else _to_compact(u, var_types)
  _assert_parity(
    mine.rosenblatt(arg, randomize_discrete=False),
    ref.rosenblatt(arg, randomize_discrete=False),
  )


@pytest.mark.parametrize("var_types", _VAR_TYPES)
def test_randomized_rosenblatt_matches_vinecop(
  var_types: list[str],
) -> None:
  # The default path. Same seeds through the same `sample_uniform`, so the
  # mixing weights are identical and the outputs must be too.
  mine, ref = _both(var_types)
  u = _expanded_data(var_types, seed=3)
  _assert_parity(
    mine.rosenblatt(u, seeds=[1, 2, 3]), ref.rosenblatt(u, seeds=[1, 2, 3])
  )


@pytest.mark.parametrize("var_types", _VAR_TYPES)
def test_randomized_rosenblatt_moves_only_discrete_columns(
  var_types: list[str],
) -> None:
  # Randomization mixes a jump's two ends, and a continuous variable's ends
  # coincide -- so it must leave exactly the continuous columns where they were,
  # up to the rounding of `v * r + v * (1 - r)`.
  mine, _ = _both(var_types)
  u = _expanded_data(var_types, seed=4)
  plain = mine.rosenblatt(u, randomize_discrete=False)
  mixed = mine.rosenblatt(u, seeds=[7])
  for j, t in enumerate(var_types):
    shift = float(np.abs(plain[:, j] - mixed[:, j]).max())
    if t == "c":
      assert shift < 1e-15, (j, shift)
    else:
      assert shift > 1e-3, (j, shift)


@pytest.mark.parametrize("var_types", _VAR_TYPES)
def test_inverse_rosenblatt_matches_vinecop(
  var_types: list[str],
) -> None:
  # The inverse cascade has no left limits to work with, so it evaluates every
  # pair as continuous -- as Vinecop does -- and its output is continuous.
  mine, ref = _both(var_types)
  w = np.random.default_rng(5).uniform(0.02, 0.98, size=(200, _D))
  _assert_parity(mine.inverse_rosenblatt(w), ref.inverse_rosenblatt(w))


@pytest.mark.parametrize("var_types", _VAR_TYPES)
def test_cdf_matches_vinecop(var_types: list[str]) -> None:
  mine, ref = _both(var_types)
  u = _expanded_data(var_types, seed=9, n=50)
  arg = u if any(t == "d" for t in var_types) else u[:, :_D]
  np.testing.assert_allclose(
    mine.cdf(arg, N=2000, seeds=[3, 4]),
    ref.cdf(arg, N=2000, seeds=[3, 4]),
    rtol=0,
    # Both count dominated points among quasi-random draws made with the same
    # seeds, and the two inverse-Rosenblatt cascades agree exactly, so the
    # sample of draws is identical. Any residual is a tie at the comparison
    # boundary: one draw out of N.
    atol=1.0 / 2000,
  )


# ---------------------------------------------------------------------------
# Layouts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("var_types", _MIXED)
def test_expanded_and_compact_layouts_agree(var_types: list[str]) -> None:
  # The compact layout is the expanded one with the continuous variables'
  # left-limit columns dropped, so it carries the same information.
  mine, _ = _both(var_types)
  u = _expanded_data(var_types, seed=6)
  compact = _to_compact(u, var_types)
  np.testing.assert_array_equal(mine.pdf(u), mine.pdf(compact))
  np.testing.assert_array_equal(
    mine.rosenblatt(u, randomize_discrete=False),
    mine.rosenblatt(compact, randomize_discrete=False),
  )


def test_continuous_vine_also_accepts_the_expanded_layout() -> None:
  # An all-continuous vine takes (n, 2d) and drops the second block, matching
  # Vinecop; the left limits are redundant there.
  var_types = ["c"] * _D
  mine, _ = _both(var_types)
  u = _expanded_data(var_types, seed=7)
  np.testing.assert_array_equal(mine.pdf(u), mine.pdf(u[:, :_D]))


@pytest.mark.parametrize("var_types", _MIXED)
def test_discrete_vine_rejects_the_bare_value_layout(
  var_types: list[str],
) -> None:
  # Reusing each value as its own left limit would silently evaluate a
  # continuous density under a discrete model, so `(n, d)` has to raise.
  mine, _ = _both(var_types)
  u = _expanded_data(var_types, seed=8)[:, :_D]
  with pytest.raises(ValueError, match="must have shape"):
    mine.pdf(u)
  with pytest.raises(ValueError, match="must have shape"):
    mine.rosenblatt(u)
  # ... but the inverse transform reads no left limit, so it accepts it.
  assert mine.inverse_rosenblatt(u).shape == u.shape


# ---------------------------------------------------------------------------
# The batched fast path declines
# ---------------------------------------------------------------------------


def test_batched_declines_on_a_discrete_vine() -> None:
  # The batched wavefront has no left-limit lane, so the dispatcher resolves
  # `batched` to False rather than silently evaluating a continuous density. A
  # raise is not an option: `batched=None` resolves to a device-dependent
  # subclass default, so an ordinary pdf(u) call would start failing.
  mine, ref = _both(["d", "c", "c", "c"], cls=_NeverBatchedVinecop)
  u = _expanded_data(["d", "c", "c", "c"], seed=10)
  _assert_parity(mine.pdf(u, batched=True), ref.pdf(u))
  _assert_parity(
    mine.rosenblatt(u, batched=True, randomize_discrete=False),
    ref.rosenblatt(u, randomize_discrete=False),
  )
  # The same vine without discrete variables does reach the fast path.
  cont, _ = _both(["c"] * _D, cls=_NeverBatchedVinecop)
  with pytest.raises(AssertionError, match="batched fast path"):
    cont.pdf(u[:, :_D], batched=True)


def test_inverse_rosenblatt_hosts_a_pair_without_as_continuous() -> None:
  # `as_continuous` is an optional capability: a custom pair that does not
  # advertise one is continuous-only already, so the inverse cascade hands it the
  # two-column input and the result cannot depend on the declared var_types.
  class _CustomPairVinecop(_ListVinecop):
    def __init__(
      self, structure: pv.RVineStructure, var_types: list[str]
    ) -> None:
      self._custom = [[GaussianBicop(base_rho=r) for r in row] for row in _RHOS]
      self._bind_vine(structure, var_types=var_types)

    def _get_pair_copula(self, tree: int, edge: int) -> BicopLike[Any]:
      return self._custom[tree][edge]

  structure = pv.RVineStructure.from_order(list(range(1, _D + 1)))
  assert not hasattr(GaussianBicop(base_rho=0.5), "as_continuous")
  w = np.random.default_rng(12).uniform(0.02, 0.98, size=(64, _D))
  np.testing.assert_array_equal(
    _CustomPairVinecop(structure, ["d", "c", "d", "c"]).inverse_rosenblatt(w),
    _CustomPairVinecop(structure, ["c"] * _D).inverse_rosenblatt(w),
  )


# ---------------------------------------------------------------------------
# End to end through Vinedist
# ---------------------------------------------------------------------------


def test_vinedist_with_a_discrete_margin_uses_the_cascade() -> None:
  # A Vinedist over a VinecopBase copula and one discrete margin must agree with
  # the same distribution built on the compiled Vinecop, and with the
  # hand-rolled `log c(u) + sum_j log f_j(x_j)` factorization.
  var_types = ["d", "c", "c", "c"]
  mine, ref = _both(var_types)
  rng = np.random.default_rng(11)
  n = 300
  counts = rng.integers(0, 5, n).astype(float)
  x = np.column_stack([counts] + [rng.normal(size=n) for _ in range(_D - 1)])
  margins = [Kde1d(type="discrete", xmin=0.0).fit(counts)] + [
    Kde1d().fit(x[:, j]) for j in range(1, _D)
  ]
  assert [m.var_type for m in margins] == var_types
  dist_mine = Vinedist(mine, margins)
  _assert_parity(dist_mine.logpdf(x), Vinedist(ref, margins).logpdf(x))

  # The factorization, spelled out: copula density on the compact layout times
  # the marginal masses / densities.
  u = np.clip(
    np.column_stack(
      [m.cdf(x[:, j]) for j, m in enumerate(margins)]
      + [margins[0].cdf_left(x[:, 0])]
    ),
    1e-10,
    1 - 1e-10,
  )
  expected = np.log(mine.pdf(u))
  for j, m in enumerate(margins):
    expected = expected + np.log(m.pdf(x[:, j]))
  np.testing.assert_allclose(dist_mine.logpdf(x), expected, rtol=1e-12)

  # And the joint cdf reaches the copula in a layout it accepts at all -- it
  # needs no left limits, but a discrete copula rejects the bare `(n, d)` one.
  assert dist_mine.cdf(x[:20], N=500, seeds=[2]).shape == (20,)
