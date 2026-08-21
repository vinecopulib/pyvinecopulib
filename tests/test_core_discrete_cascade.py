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
from pyvinecopulib.core import (
  BicopLike,
  DiscretePair,
  VinecopBase,
)
from pyvinecopulib.core._discrete import continuous_view

from .conftest import GaussianBicop, HostedVinecop

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


class _ListVinecop(HostedVinecop):
  """A hosted vine that stamps its pairs with the types derived for their slot.

  The one thing it adds over :class:`HostedVinecop`: a subclass that *owns* its
  pair copulas is expected to declare each pair's types, and the compiled
  ``Bicop`` needs them set for its discrete evaluation to be exercised at all.
  """

  def __init__(
    self,
    pairs: list[list[pv.Bicop]],
    structure: pv.RVineStructure,
    var_types: Optional[list[str]] = None,
  ) -> None:
    super().__init__(pairs, structure, var_types)
    for tree, row in enumerate(pairs):
      for edge, pair in enumerate(row):
        pair.var_types = list(self.pair_var_types(tree, edge))


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
# DiscretePair: a continuous pair copula hosted on a discrete edge
# ---------------------------------------------------------------------------

#: Fit controls restricted to one family, so a pair fitted here and one fitted
#: by ``Vinecop`` differ only by the optimizer, not by family selection.
_GAUSSIAN_PAIR = pv.FitControlsBicop(family_set=[pv.families.gaussian])
_GAUSSIAN_VINE = pv.FitControlsVinecop(
  family_set=[pv.families.gaussian], num_threads=1
)
#: Rotations of an asymmetric family: `DiscretePair` differences the *public*
#: (rotated) distribution function where the compiled pair differences the
#: unrotated one, so a rotation is what exercises the equivalence of the two.
_ROTATIONS = [0, 90, 180, 270]


def _pair_edge_data(seed: int, n: int = 400) -> np.ndarray:
  """An ``(n, 4)`` pair-copula input with atoms of realistic width."""
  rng = np.random.default_rng(seed)
  values = rng.uniform(0.2, 0.95, size=(n, 2))
  return np.column_stack([values, values - rng.uniform(0.02, 0.15, (n, 2))])


@pytest.mark.parametrize("var_types", [["d", "c"], ["c", "d"], ["d", "d"]])
def test_discrete_pair_matches_bicop_unrotated(var_types: list[str]) -> None:
  # Rotation 0 sends identical doubles into the identical compiled leaf, so on
  # one toolchain the difference quotients agree bit for bit. Not across them:
  # CI measured 4.0e-15 absolute (7.6e-15 relative) on x86_64, where FMA
  # contraction reassociates the quotient, against exact equality on arm64.
  par = np.array([[0.5]])
  wrapped = pv.Bicop.from_family(pv.families.gaussian, parameters=par)
  ref = pv.Bicop.from_family(
    pv.families.gaussian, parameters=par, var_types=var_types
  )
  pair = DiscretePair(wrapped, (var_types[0], var_types[1]))
  u = _pair_edge_data(seed=7)
  for method in ("pdf", "cdf", "hfunc1", "hfunc2"):
    _assert_parity(
      np.asarray(getattr(pair, method)(u)),
      np.asarray(getattr(ref, method)(u)),
    )


@pytest.mark.parametrize("rotation", _ROTATIONS)
@pytest.mark.parametrize("var_types", [["d", "c"], ["c", "d"], ["d", "d"]])
def test_discrete_pair_matches_bicop_rotated(
  rotation: int, var_types: list[str]
) -> None:
  # A rotated distribution function carries an affine part that cancels in every
  # quotient only up to rounding, and the two-discrete density divides that
  # cancellation by the product of the atom widths. With atoms of realistic
  # width (>= 0.02) the measured worst case is 4e-13; regenerate this bound
  # rather than widening it.
  par = np.array([[2.0]])
  wrapped = pv.Bicop.from_family(
    pv.families.clayton, rotation=rotation, parameters=par
  )
  ref = pv.Bicop.from_family(
    pv.families.clayton,
    rotation=rotation,
    parameters=par,
    var_types=var_types,
  )
  pair = DiscretePair(wrapped, (var_types[0], var_types[1]))
  u = _pair_edge_data(seed=11)
  for method in ("pdf", "hfunc1", "hfunc2"):
    np.testing.assert_allclose(
      np.asarray(getattr(pair, method)(u)),
      np.asarray(getattr(ref, method)(u)),
      rtol=1e-11,
      atol=1e-11,
      err_msg=method,
    )


#: Every family ``Bicop`` can fit. The pair-level parity tests ran over
#: ``{gaussian, clayton}`` only, which is why a family whose discrete surface
#: disagreed with all the others went unnoticed
#: (`vinecopulib#739 <https://github.com/vinecopulib/vinecopulib/pull/739>`_).
_FAMILIES = [
  pv.families.tll,
  pv.families.gaussian,
  pv.families.student,
  pv.families.clayton,
  pv.families.gumbel,
  pv.families.frank,
  pv.families.joe,
  pv.families.bb1,
  pv.families.bb6,
  pv.families.bb7,
  pv.families.bb8,
  pv.families.tawn,
]


def _fitted_pair(
  family: pv.BicopFamily, var_types: list[str], seed: int = 3
) -> pv.Bicop:
  """One family fitted on a discrete or mixed edge, atoms a twelfth wide."""
  rng = np.random.default_rng(seed)
  z = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.7], [0.7, 1.0]], size=2000)
  u = pv.utils.to_pseudo_obs(z)
  sub = np.clip(u - 1.0 / 12, 1e-10, None)
  data = np.column_stack(
    [u] + [sub[:, j] for j, ty in enumerate(var_types) if ty == "d"]
  )
  return pv.Bicop.from_data(
    data,
    controls=pv.FitControlsBicop(family_set=[family]),
    var_types=var_types,
  )


@pytest.mark.parametrize("family", _FAMILIES)
@pytest.mark.parametrize("var_types", [["d", "c"], ["c", "d"], ["d", "d"]])
def test_discrete_pair_matches_every_fitted_family(
  family: pv.BicopFamily, var_types: list[str]
) -> None:
  """The same fitted pair, differenced here and by ``Bicop``, on every family.

  ``tll`` is the one that matters: it is the default family, so this is the
  default discrete path, and it is the family whose compiled surface was wrong
  until the pin bump. ``as_continuous`` strips the atom declaration without
  touching the parameters, so the two objects are the same copula.
  """
  ref = _fitted_pair(family, var_types)
  pair = DiscretePair(ref.as_continuous(), (var_types[0], var_types[1]))
  u = _pair_edge_data(seed=13)
  for method in ("pdf", "cdf", "hfunc1", "hfunc2"):
    np.testing.assert_allclose(
      np.asarray(getattr(pair, method)(u)),
      np.asarray(getattr(ref, method)(u)),
      rtol=1e-11,
      atol=1e-11,
      err_msg=f"{family!r} {var_types} {method}",
    )


@pytest.mark.parametrize("levels", [2, 8, 32])
@pytest.mark.parametrize("family", [pv.families.tll, pv.families.gaussian])
def test_the_atom_masses_of_a_discrete_edge_sum_to_one(
  family: pv.BicopFamily, levels: int
) -> None:
  """``sum_atoms c(u1, u2) * (u1 - u1^-) == 1`` for every ``u2``.

  This needs no reference implementation and no tolerance argument: the
  difference quotients telescope to ``h2(1, u2) - h2(0, u2)``, which is one
  exactly. It is the assertion that decides which of two candidate discrete
  surfaces is right, and it is what the midpoint rule the compiled ``tll``
  pair used to apply violated -- by 2% at two atoms, and not converging away.
  """
  ref = _fitted_pair(family, ["d", "c"])
  pair = DiscretePair(ref.as_continuous(), ("d", "c"))
  edges = np.arange(levels + 1) / levels
  hi, lo = edges[1:], edges[:-1]
  for u2 in (0.1, 0.5, 0.9):
    col = np.full(levels, u2)
    u = np.column_stack([hi, col, lo, col])
    for name, dens in (
      ("DiscretePair", np.asarray(pair.pdf(u))),
      ("Bicop", np.asarray(ref.pdf(u))),
    ):
      np.testing.assert_allclose(
        float(np.sum(dens * (hi - lo))),
        1.0,
        rtol=1e-9,
        atol=0.0,
        err_msg=f"{name} {family!r} {levels} atoms at u2={u2}",
      )


def test_discrete_pair_narrow_atoms_fall_back_to_the_derivative() -> None:
  # Below the 5e-5 threshold the quotient is replaced by the derivative, which
  # is the continuous quantity: an atom that narrow is a continuous variable in
  # all but name.
  par = np.array([[0.5]])
  cop = pv.Bicop.from_family(pv.families.gaussian, parameters=par)
  pair = DiscretePair(cop, ("d", "d"))
  values = np.array([[0.4, 0.7], [0.25, 0.55]])
  u = np.column_stack([values, values - 1e-9])
  np.testing.assert_allclose(
    np.asarray(pair.pdf(u)), np.asarray(cop.pdf(values)), rtol=1e-9, atol=1e-9
  )
  np.testing.assert_allclose(
    np.asarray(pair.hfunc1(u)),
    np.asarray(cop.hfunc1(values)),
    rtol=1e-9,
    atol=1e-9,
  )


def test_discrete_pair_hinv_inverts_the_discrete_hfunc() -> None:
  par = np.array([[0.6]])
  cop = pv.Bicop.from_family(pv.families.gaussian, parameters=par)
  pair = DiscretePair(cop, ("d", "d"))
  u = _pair_edge_data(seed=5, n=200)
  # hinv1 solves hfunc1 in the second argument, so column 1 is the level.
  level = u[:, 1]
  back = np.asarray(pair.hinv1(u))
  recovered = np.asarray(
    pair.hfunc1(np.column_stack([u[:, 0], back, u[:, 2], back]))
  )
  np.testing.assert_allclose(recovered, level, rtol=1e-8, atol=1e-8)
  # Discreteness is not silently dropped: inverting the atom's difference
  # quotient is a different problem from inverting the derivative.
  continuous = np.asarray(cop.hinv1(u[:, :2]))
  assert np.abs(back - continuous).max() > 1e-3


def test_discrete_pair_flip_swaps_the_variable_types() -> None:
  cop = pv.Bicop.from_family(
    pv.families.clayton, rotation=90, parameters=np.array([[2.0]])
  )
  pair = DiscretePair(cop, ("d", "c"))
  flipped = pair.flip()
  assert flipped.var_types == ["c", "d"]
  u = _pair_edge_data(seed=3, n=100)
  swapped = u[:, [1, 0, 3, 2]]
  np.testing.assert_allclose(
    np.asarray(flipped.pdf(swapped)), np.asarray(pair.pdf(u)), rtol=1e-12
  )


def test_discrete_pair_requires_a_cdf() -> None:
  # The h-functions of a discrete argument are difference quotients of the
  # copula's distribution function, and the cascades evaluate an h-function at
  # every edge -- so a pair without a `cdf` cannot sit on a discrete edge. Its
  # mixed *density* is a quotient of h-functions and needs none, which is why the
  # failure surfaces here rather than on `pdf`.
  u = _pair_edge_data(seed=1, n=10)
  pair = DiscretePair(GaussianBicop(base_rho=0.5), ("d", "c"))
  with pytest.raises(NotImplementedError, match="cdf"):
    pair.hfunc1(u)
  with pytest.raises(NotImplementedError, match="cdf"):
    DiscretePair(GaussianBicop(base_rho=0.5), ("d", "d")).pdf(u)


def test_discrete_pair_rejects_a_two_column_input() -> None:
  pair = DiscretePair(
    pv.Bicop.from_family(pv.families.gaussian, parameters=np.array([[0.5]])),
    ("d", "c"),
  )
  with pytest.raises(ValueError, match=r"shape \(n, 4\)"):
    pair.pdf(np.array([[0.5, 0.5]]))


class _WrappingVinecop(_ListVinecop):
  """A vine whose pairs learn about discreteness only through the wrapper.

  The counterpart of :class:`_ListVinecop`, which stamps ``var_types`` onto the
  compiled pairs it hosts: here the hosted pairs stay continuous and
  ``DiscretePair`` supplies the mixed-discrete surface, which is the route a
  custom pair copula takes.
  """

  def __init__(
    self,
    pairs: list[list[pv.Bicop]],
    structure: pv.RVineStructure,
    var_types: Optional[list[str]] = None,
  ) -> None:
    self._pairs = pairs
    self._bind_vine(structure, var_types=var_types)

  def _get_pair_copula(self, tree: int, edge: int) -> BicopLike[Any]:
    pair = cast("BicopLike[Any]", self._pairs[tree][edge])
    types = self.pair_var_types(tree, edge)
    if "d" not in types:
      return pair
    return DiscretePair(pair, types)


@pytest.mark.parametrize("var_types", _MIXED)
def test_wrapped_pairs_match_vinecop(var_types: list[str]) -> None:
  # The design invariant: the pair copulas need not know they sit on a discrete
  # edge. Gaussian pairs at rotation 0, so this is exact.
  mine, ref = _both(var_types, cls=_WrappingVinecop)
  u = _to_compact(_expanded_data(var_types, seed=21), var_types)
  _assert_parity(mine.pdf(u), ref.pdf(u))
  _assert_parity(
    mine.rosenblatt(u, randomize_discrete=False),
    ref.rosenblatt(u, randomize_discrete=False),
  )
  _assert_parity(mine.loglik(u), ref.loglik(u))


# ---------------------------------------------------------------------------
# The fit engines on discrete data
# ---------------------------------------------------------------------------


def _discrete_fit_edge(
  tree: int,
  edge: int,
  u_e: Any,
  x_e: Any,
  var_types: Any = ("c", "c"),
) -> BicopLike[Any]:
  """Fit one Gaussian pair copula, declaring the edge's variable types."""
  del tree, edge, x_e
  # Cast: the conformance is nominal -- `Bicop.cdf` takes per-row `parameters`
  # where `BicopLike.cdf` takes keyword-only `x`.
  return cast(
    "BicopLike[Any]",
    pv.Bicop.from_data(
      np.asarray(u_e), controls=_GAUSSIAN_PAIR, var_types=list(var_types)
    ),
  )


def _as_bicops(pairs: list[list[Any]]) -> list[list[pv.Bicop]]:
  """The compiled pairs a fit engine returns, whose static type it widens.

  ``fit`` and ``select`` are typed against ``BicopLike``; the callback above
  builds ``Bicop`` objects, which is what the hosts here need (they read
  ``var_types``).
  """
  return cast("list[list[pv.Bicop]]", pairs)


def _dependent_expanded(
  var_types: list[str], seed: int, n: int = 600
) -> np.ndarray:
  """An ``(n, 2d)`` sample with dependence, so a selected tree is meaningful.

  A one-factor latent Gaussian pushed through the rank transform, with the
  declared variables collapsed onto the ordered categories of ``_BOUNDS``.
  """
  rng = np.random.default_rng(seed)
  d = len(var_types)
  base = rng.standard_normal((n, 1))
  z = 0.7 * base + 0.7 * rng.standard_normal((n, d))
  values, second = [], []
  for j, t in enumerate(var_types):
    p = pv.to_pseudo_obs(z[:, [j]]).ravel()
    if t == "d":
      level = np.searchsorted(_BOUNDS[1:-1], p)
      values.append(_BOUNDS[level + 1])
      second.append(_BOUNDS[level])
    else:
      values.append(p)
      second.append(p)
  return np.column_stack(values + second)


def _order_structure(d: int) -> pv.RVineStructure:
  return pv.RVineStructure.from_order(list(range(1, d + 1)))


@pytest.mark.parametrize("var_types", _MIXED)
def test_fit_matches_vinecop(var_types: list[str]) -> None:
  # Both engines assemble each edge's four columns from the same left-limit
  # cascade and hand them to the same compiled pair fit, so the fitted vines
  # agree to the parity bounds; a mis-assembled column lands at O(0.1).
  structure = _order_structure(len(var_types))
  u = _to_compact(_dependent_expanded(var_types, seed=4), var_types)
  pairs = VinecopBase.fit(structure, u, _discrete_fit_edge, var_types=var_types)
  mine = _ListVinecop(_as_bicops(pairs), structure, var_types=var_types)
  ref = pv.Vinecop.from_data(
    u, structure=structure, var_types=var_types, controls=_GAUSSIAN_VINE
  )
  _assert_parity(mine.pdf(u), ref.pdf(u))
  _assert_parity(mine.loglik(u), ref.loglik(u))


@pytest.mark.parametrize("var_types", _MIXED)
def test_select_matches_vinecop(var_types: list[str]) -> None:
  # The edge weights read only the value columns, so the selected R-vine matrix
  # must match the compiled selector's exactly -- byte for byte, as it does for
  # continuous data.
  u = _to_compact(_dependent_expanded(var_types, seed=9), var_types)
  structure, pairs = VinecopBase.select(
    u, _discrete_fit_edge, var_types=var_types
  )
  auto = pv.Vinecop.from_data(u, var_types=var_types, controls=_GAUSSIAN_VINE)
  assert np.array_equal(
    np.asarray(structure.matrix), np.asarray(auto.structure.matrix)
  )
  assert list(structure.order) == list(auto.structure.order)
  mine = _ListVinecop(_as_bicops(pairs), structure, var_types=var_types)
  _assert_parity(mine.pdf(u), auto.pdf(u))


@pytest.mark.parametrize("var_types", _MIXED)
def test_selected_pairs_carry_the_derived_variable_types(
  var_types: list[str],
) -> None:
  # The design rests on the per-edge types being a function of the structure:
  # `select` fits with the types its graph derives, then places and flips the
  # pairs onto a structure whose types are re-derived. The two must agree, or
  # every hosting vine would silently disagree with the fit.
  u = _to_compact(_dependent_expanded(var_types, seed=13), var_types)
  structure, pairs = VinecopBase.select(
    u, _discrete_fit_edge, var_types=var_types
  )
  # `_WrappingVinecop` is the host that does *not* stamp its pairs, so the types
  # compared here are the ones `select` fitted with.
  bicops = _as_bicops(pairs)
  host = _WrappingVinecop(bicops, structure, var_types=var_types)
  for tree, row in enumerate(bicops):
    for edge, pair in enumerate(row):
      assert list(pair.var_types) == list(host.pair_var_types(tree, edge)), (
        tree,
        edge,
      )


def test_fit_edge_receives_the_edge_types_and_four_columns() -> None:
  var_types = ["d", "c", "d", "c"]
  seen: dict[tuple[int, int], tuple[int, tuple[str, ...]]] = {}

  def recording(
    tree: int, edge: int, u_e: Any, x_e: Any, var_types: Any = ("c", "c")
  ) -> BicopLike[Any]:
    seen[(tree, edge)] = (int(np.asarray(u_e).shape[1]), tuple(var_types))
    return _discrete_fit_edge(tree, edge, u_e, x_e, var_types)

  structure = _order_structure(len(var_types))
  u = _to_compact(_dependent_expanded(var_types, seed=2), var_types)
  pairs = VinecopBase.fit(structure, u, recording, var_types=var_types)
  host = _WrappingVinecop(_as_bicops(pairs), structure, var_types=var_types)
  assert len(seen) == sum(len(row) for row in pairs)
  for (tree, edge), (n_cols, types) in seen.items():
    expected = host.pair_var_types(tree, edge)
    assert types == expected, (tree, edge)
    # The left limits are only handed over where the edge needs them.
    assert n_cols == (4 if "d" in expected else 2), (tree, edge)


def test_a_continuous_fit_edge_fails_loudly_on_a_discrete_edge() -> None:
  # A callback written for a continuous vine cannot fit four columns of data;
  # it must say so rather than quietly fitting the value columns.
  def continuous_only(
    tree: int, edge: int, u_e: Any, x_e: Any
  ) -> BicopLike[Any]:
    return _discrete_fit_edge(tree, edge, u_e, x_e)

  var_types = ["d", "c", "c", "c"]
  u = _to_compact(_dependent_expanded(var_types, seed=6), var_types)
  with pytest.raises(TypeError, match="var_types"):
    VinecopBase.fit(
      _order_structure(len(var_types)),
      u,
      continuous_only,
      var_types=var_types,
    )
  with pytest.raises(TypeError, match="var_types"):
    VinecopBase.select(u, continuous_only, var_types=var_types)


@pytest.mark.parametrize("engine", ["fit", "select"])
def test_fit_engines_reject_a_missing_left_limit_block(engine: str) -> None:
  var_types = ["d", "c", "c", "c"]
  d = len(var_types)
  u = _to_compact(_dependent_expanded(var_types, seed=8), var_types)
  with pytest.raises(ValueError, match=f"{engine}: u must have shape"):
    if engine == "fit":
      VinecopBase.fit(
        _order_structure(d), u[:, :d], _discrete_fit_edge, var_types=var_types
      )
    else:
      VinecopBase.select(u[:, :d], _discrete_fit_edge, var_types=var_types)


@pytest.mark.parametrize("engine", ["fit", "select"])
def test_fit_engines_reject_an_unknown_variable_type(engine: str) -> None:
  var_types = ["d", "x", "c", "c"]
  u = _dependent_expanded(["d", "c", "c", "c"], seed=8)
  with pytest.raises(ValueError, match="var_types entries must be 'c' or 'd'"):
    if engine == "fit":
      VinecopBase.fit(
        _order_structure(4), u, _discrete_fit_edge, var_types=var_types
      )
    else:
      VinecopBase.select(u, _discrete_fit_edge, var_types=var_types)


def test_fit_checks_var_types_against_the_structure() -> None:
  # `fit` walks a given structure, so that structure fixes the dimension.
  u = _dependent_expanded(["d", "c", "c", "c"], seed=8)
  with pytest.raises(ValueError, match="var_types has 3 entries, expected 4"):
    VinecopBase.fit(
      _order_structure(4), u, _discrete_fit_edge, var_types=["c", "c", "c"]
    )


def test_select_takes_its_dimension_from_var_types() -> None:
  # `select` has no structure to read the dimension from, and `u` is wider than
  # the vine once left limits are present -- so `var_types` is what fixes it,
  # and a length that does not match the data surfaces as a layout error.
  var_types = ["d", "c", "c", "c"]
  u = _to_compact(_dependent_expanded(var_types, seed=8), var_types)
  assert u.shape[1] == 5
  with pytest.raises(ValueError, match="select: u must have shape"):
    VinecopBase.select(u, _discrete_fit_edge, var_types=var_types[:3])


def _count_expanded(
  var_types: list[str], seed: int, n: int = 800
) -> np.ndarray:
  """An ``(n, 2d)`` sample whose atoms are a genuine count variable's.

  Unlike ``_expanded_data``, the discrete margins here are a ``Binomial(4, 0.5)``
  distribution function: the top atom has ``F(x) = 1`` and the bottom
  ``F(x^-) = 0``, the two places a pair copula's own trimming has to agree with
  the cascade's.
  """
  rng = np.random.default_rng(seed)
  cdf = np.cumsum([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0
  values, second = [], []
  for t in var_types:
    if t == "d":
      counts = rng.binomial(4, 0.5, n)
      values.append(cdf[counts])
      second.append(np.where(counts > 0, cdf[counts - 1], 0.0))
    else:
      v = rng.uniform(0.02, 0.98, n)
      values.append(v)
      second.append(v)
  return np.column_stack(values + second)


@pytest.mark.parametrize("var_types", _MIXED)
def test_wrapped_pairs_match_vinecop_on_count_margins(
  var_types: list[str],
) -> None:
  # Tail-dependent pairs on a real count margin, whose atoms reach both ends of
  # the unit interval. `DiscretePair` has to trim its four columns before
  # differencing them, exactly as the compiled pair copula does -- the width of
  # the *trimmed* atom is the denominator. Trimming afterwards leaves a relative
  # error of order 1e-10 in the quotient, which a tail-dependent density of
  # order 1e3 turns into an absolute 4e-7.
  structure = pv.RVineStructure.from_order(list(range(1, len(var_types) + 1)))
  thetas = ([2.0, 1.5, 1.8][: len(var_types) - 1], [1.3, 1.4], [1.2])
  pairs = [
    [
      pv.Bicop.from_family(pv.families.gumbel, parameters=np.array([[t]]))
      for t in row
    ]
    for row in thetas[: len(var_types) - 1]
  ]
  u = _to_compact(_count_expanded(var_types, seed=17), var_types)
  mine = _WrappingVinecop(pairs, structure, var_types=var_types)
  ref = pv.Vinecop.from_structure(
    structure=structure,
    pair_copulas=[[p for p in row] for row in pairs],
    var_types=var_types,
  )
  # The parity bounds, not exact equality: the same cross-toolchain last-bit
  # spread applies, and it is nine orders of magnitude below the 4e-7 the
  # untrimmed quotient produced.
  _assert_parity(mine.pdf(u), ref.pdf(u))


def test_a_discrete_pair_without_a_continuous_view_is_rejected() -> None:
  # The inverse cascade evaluates every pair as continuous, so a pair copula that
  # declares atoms has to offer a continuous view. Saying so beats the column-
  # count error the pair itself would raise several frames down.
  class _NoView:
    var_types = ["d", "c"]

    def pdf(self, u: Any) -> Any:
      raise AssertionError("not reached")

  with pytest.raises(ValueError, match="no as_continuous"):
    continuous_view(_NoView())
