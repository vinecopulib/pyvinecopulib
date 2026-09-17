"""Tests for the marginal half of the sklearn estimators.

What this file pins is the delegation. That an estimator *is* a `Vinedist` --
the same layout, the same log-density, the same draws -- rather than a second
implementation of Sklar's theorem that can drift from the first. That the
default `distribution=` still means what it meant before margins were
configurable: a `Kde1d` per column, to the last bit. That a `margin_controls=`
the caller gives is honored per column and never mutated, so `clone`
reproduces the estimator. That a family search runs per column and reads what
`schema_` declared about that column rather than re-inferring it from the
sample. And that the ``{0, 1}`` dummies of an expanded unordered categorical
are fitted on the support they actually have, rather than on a padded grid that
puts mass on values that cannot occur.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("pandas")

import pandas as pd
from sklearn.base import clone

import pyvinecopulib as pv
from pyvinecopulib.core import FitControlsMargin, Kde1d, Vinedist
from pyvinecopulib.margins import SciPyMargin
from pyvinecopulib.sklearn import VineDensity, VineRegressor

from .helpers import AtomicMargin, FlatMargin


class ParametricVineDensity(Vinedist):
  """A `Vinedist` whose margins are SciPy families: the `margin_class` route."""

  margin_class = SciPyMargin


class AtomicVinedist(Vinedist):
  """A `Vinedist` whose margins report mass rather than a density."""

  margin_class = AtomicMargin


class _DiscreteMargin(FlatMargin):
  """A margin declaring atoms before it is fitted, which no shipped one does."""

  @property
  def var_type(self) -> str:
    return "d"


@pytest.fixture
def cat_df() -> pd.DataFrame:
  """A frame whose ordered categorical has non-integer levels."""
  rs = np.random.RandomState(0)
  return pd.DataFrame(
    {
      "a": rs.normal(size=200),
      "grade": pd.Categorical(
        rs.choice([1.5, 2.5, 3.5], 200),
        categories=[1.5, 2.5, 3.5],
        ordered=True,
      ),
    }
  )


# --- the default is the old pipeline ---------------------------------------- #


def test_default_margins_reproduce_the_kde1d_pipeline(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """`margins=None` is the pre-margins pipeline, to floating-point rounding."""
  X, _, _ = sample_array_data
  eps = 1e-10
  kdes = []
  for j in range(X.shape[1]):
    kde = pv.core.Kde1d()
    kde.fit(X[:, j])
    kdes.append(kde)

  def transform(z: np.ndarray) -> np.ndarray:
    return np.column_stack(
      [np.clip(k.cdf(z[:, j]), eps, 1 - eps) for j, k in enumerate(kdes)]
    )

  vine = pv.Vinecop.from_data(
    data=transform(X),
    var_types=["c"] * X.shape[1],
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
    ),
  )
  head = X[:20]
  expected = np.log(np.asarray(vine.pdf(transform(head)))) + sum(
    np.log(np.asarray(kdes[j].pdf(head[:, j]))) for j in range(X.shape[1])
  )
  got = VineDensity().fit(X).score_samples(head)
  np.testing.assert_allclose(got, expected, rtol=1e-12, atol=1e-14)


# --- the fitted distribution ------------------------------------------------ #


def test_fit_publishes_the_distribution(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """`distribution_` is the fitted model, and agrees with the estimator."""
  X, _, _ = sample_array_data
  est = VineDensity().fit(X)
  assert isinstance(est.distribution_, pv.Vinedist)
  assert est.distribution_.dim == est.n_features_in_
  np.testing.assert_allclose(est.distribution_.pdf(X[:10]), est.pdf(X[:10]))
  np.testing.assert_allclose(
    est.distribution_.sample(5, seeds=[1, 2, 3]).shape, (5, 2)
  )


def test_the_held_copula_still_answers_as_the_vine(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """``distribution_`` holds the fitted vine itself, surface and all."""
  X, _, _ = sample_array_data
  est = VineDensity().fit(X)
  copula: Any = est.distribution_.vinecop
  assert copula.dim == est.n_features_in_
  assert copula.var_types == ["c", "c"]
  np.testing.assert_array_equal(
    np.asarray(copula.structure.matrix), np.asarray(est.structure_.matrix)
  )
  np.testing.assert_allclose(
    np.asarray(copula.pdf(est._to_u_scale(X[:10]))),
    est.pdf(X[:10], copula_only=True),
  )


def test_joint_distribution_of_the_regressor_leads_with_the_response(
  regression_data: tuple[np.ndarray, np.ndarray, np.ndarray, float],
) -> None:
  """The regressor's distribution is over ``(Y, X)``, the response first."""
  X, y, _, _ = regression_data
  est = VineRegressor().fit(X, y)
  assert est.distribution_.dim == X.shape[1] + 1
  np.testing.assert_allclose(
    est.distribution_.logpdf(np.column_stack([y[:10], X[:10]])),
    est._pdf_samples(X[:10], y=y[:10], log=True),
  )


# --- bounds of an expanded dummy -------------------------------------------- #


def test_expanded_dummies_are_fitted_on_their_own_support(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """A ``{0, 1}`` dummy is bounded, so no mass lands on an impossible value."""
  X_df, expanded = sample_dataframe_data
  est = VineDensity().fit(X_df)
  dummies = [j for j, name in enumerate(expanded) if name.startswith("cat1_")]
  assert dummies
  for j in dummies:
    # `support` / `var_type` / `is_fitted` are optional capabilities, not
    # members of `MarginLike`, so they are read off a loosely typed binding.
    margin: Any = est.distribution_.margins[j]
    assert margin.support == (0.0, 1.0)
    assert margin.pdf(np.array([2.0])).item() == 0.0
    assert margin.cdf(np.array([-1.0])).item() == 0.0
  # A continuous column states no support, so it stays unbounded.
  unbounded: Any = est.distribution_.margins[0]
  assert unbounded.support == (-np.inf, np.inf)


def test_ordered_categorical_is_fitted_on_its_declared_levels() -> None:
  """A categorical declares its levels, so counts get no density below zero."""
  rng = np.random.default_rng(0)
  counts = rng.poisson(3.0, 400)
  X_df = pd.DataFrame(
    {
      "cnt": pd.Categorical(counts, ordered=True),
      "z": rng.normal(size=400),
    }
  )
  est = VineDensity(random_state=0).fit(X_df)
  assert est.schema_["supports"][0] == (
    float(np.min(counts)),
    float(np.max(counts)),
  )
  margin: Any = est.distribution_.margins[0]
  assert margin.var_type == "d"
  assert margin.support == (float(np.min(counts)), float(np.max(counts)))
  # Padding the grid below the smallest level is what put mass on impossible
  # counts; the declared support removes it exactly rather than approximately.
  assert margin.cdf(np.array([-1.0])).item() == 0.0
  assert margin.pdf(np.array([-1.0])).item() == 0.0


def test_the_schema_reaches_a_column_the_specification_addresses() -> None:
  """`margin_controls` does not cost a margin what the input declared.

  The declaration and the controls are separate arguments for this reason: a
  broadcast controls object addressing every column would otherwise have to
  carry `var_type=None, support=None` for the ones it says nothing about, and
  every margin would re-infer both from the sample -- strictly less than
  `schema_` already knew. What the declaration buys is
  visible in the family that wins: a float column whose observations happen to
  be whole numbers reads as counts on its own, and is selected among the
  continuous families only because the input said the column was continuous.
  """
  pytest.importorskip("scipy")
  rng = np.random.default_rng(0)
  whole = rng.integers(0, 50, size=400) * 1.0
  X_df = pd.DataFrame(
    {
      "grade": pd.Categorical(
        rng.integers(0, 5, size=400), categories=range(5), ordered=True
      ),
      "whole": whole,
    }
  )
  est = VineDensity(distribution=ParametricVineDensity, random_state=0).fit(
    X_df
  )
  ordered: Any = est.distribution_.margins[0]
  assert ordered.var_type == "d"
  assert ordered.family_name in ("poisson", "nbinom", "geom")
  continuous: Any = est.distribution_.margins[1]
  assert continuous.var_type == "c"
  # Which is not what the column says about itself, left to a bare search.
  assert SciPyMargin().select(whole).var_type == "d"


def test_a_declared_support_reaches_a_selected_family() -> None:
  """Bounds a specification-addressed column carries steer the search too.

  The bounded candidates are the ones that can pin their endpoints, and they
  are reachable in no other way: nothing in a sample on ``(0, 1)`` says the
  variable is a proportion rather than a positive variable that stayed small.
  """
  pytest.importorskip("scipy")
  rng = np.random.default_rng(0)
  X = np.column_stack([rng.uniform(0.0, 1.0, 300), rng.normal(size=300)])
  est = VineDensity(distribution=ParametricVineDensity, random_state=0)
  est.schema_ = {
    "var_types": ["c", "c"],
    "supports": [(0.0, 1.0), None],
  }
  est.fit(X)
  bounded: Any = est.distribution_.margins[0]
  assert bounded.family_name in ("uniform", "beta")
  assert bounded.support == (0.0, 1.0)
  # The column that declared nothing is not bounded to its own sample.
  free: Any = est.distribution_.margins[1]
  assert free.support == (-np.inf, np.inf)


def test_margin_summary_is_published(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """One row per fitted variable, in the joint model's order."""
  X_df, _ = sample_dataframe_data
  est = VineDensity(random_state=0).fit(X_df)
  rows = est.margin_summary_
  assert len(rows) == est.n_model_features_
  assert [row["var_type"] for row in rows] == est.distribution_.var_types
  assert all(row["family"] == "kde1d" for row in rows)


def test_cdf_works_with_columns_that_have_atoms(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """A copula with atoms validates the whole layout, and gets it."""
  X_df, _ = sample_dataframe_data
  est = VineDensity(random_state=0).fit(X_df)
  values = est.cdf(X_df.iloc[:10])
  assert values.shape == (10,)
  assert np.all((values >= 0.0) & (values <= 1.0))


def test_sample_returns_the_expanded_feature_space(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """Dummies come back as columns: collapsing them is not an inverse."""
  X_df, expanded = sample_dataframe_data
  est = VineDensity(random_state=0).fit(X_df)
  drawn = est.sample(20, random_state=1)
  assert drawn.shape == (20, len(expanded))
  j = expanded.index("cat1_B")
  assert set(np.unique(drawn[:, j])) <= {0.0, 1.0}


# --- specifications --------------------------------------------------------- #


def test_clone_round_trips_the_distribution_parameter(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """`clone` reproduces an estimator that fits the same margins."""
  X, _, _ = sample_array_data
  est = VineDensity(distribution=ParametricVineDensity).fit(X)
  cloned = clone(est)
  assert cloned.distribution is ParametricVineDensity
  np.testing.assert_allclose(
    cloned.fit(X).score_samples(X[:10]), est.score_samples(X[:10])
  )


def test_a_density_less_margin_reports_itself(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """A margin with no density says so, from the margin, when one is asked for.

  `MarginLike` requires `pdf`, so a margin that declines it is outside the
  contract rather than a case the estimators screen for: the refusal is the
  margin's own and reaches whoever asked. `VineRegressor` never asks -- it
  reads the copula density and the response's `icdf` only -- so the same
  margin serves it.
  """
  X, _, _ = sample_array_data
  est = VineDensity(distribution=AtomicVinedist).fit(X)
  with pytest.raises(NotImplementedError, match="no density"):
    est.score_samples(X[:10])
  VineRegressor(distribution=AtomicVinedist, use_grid=False).fit(X, X[:, 0])


def test_a_wrong_length_sequence_is_refused(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """A per-column sequence must carry one entry per column."""
  X, _, _ = sample_array_data
  with pytest.raises(ValueError, match="length 1, but there are 2"):
    VineDensity(margin_controls=[FitControlsMargin()]).fit(X)


def test_a_mapping_leaves_the_other_columns_unconfigured(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """Addressing one column must not configure the rest.

  The keys are the *expanded* feature names, so an unordered categorical's
  dummies are addressable one at a time -- and a column the mapping does not
  name searches the curated set, not the one entry the mapping carries.
  """
  pytest.importorskip("scipy")
  X_df, expanded = sample_dataframe_data
  est = VineDensity(
    distribution=ParametricVineDensity,
    margin_controls={"cont1": FitControlsMargin(family_set=["laplace"])},
    random_state=0,
  ).fit(X_df)
  margins = est.distribution_.margins
  addressed: Any = margins[0]
  assert addressed.family_name == "laplace"
  # The dummy searched the discrete candidates its own declaration admits.
  dummy: Any = margins[expanded.index("cat1_B")]
  assert dummy.var_type == "d"
  assert dummy.family_name != "laplace"


def test_a_preset_schema_supplies_what_an_array_cannot_carry() -> None:
  """`schema_` declares per-column variable types and bounds before `fit`."""
  rng = np.random.default_rng(0)
  X = np.column_stack([rng.poisson(3.0, size=200) * 1.0, rng.normal(size=200)])
  est = VineDensity()
  est.schema_ = {
    "var_types": ["d", "c"],
    "supports": [(0.0, 20.0), None],
  }
  est.fit(X)
  counts: Any = est.distribution_.margins[0]
  assert counts.var_type == "d"
  assert counts.support == (0.0, 20.0)
  assert est.distribution_.var_types == ["d", "c"]


# --- family selection ------------------------------------------------------- #


def test_parametric_margins_select_a_family_per_column(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """`distribution=ParametricVineDensity` chooses a family per column, and publishes it.

  The margin *is* the winner, so what it selected is read off the margin and
  reported by `margin_summary_` alongside every other column's.
  """
  pytest.importorskip("scipy")
  X, _, _ = sample_array_data
  est = VineDensity(distribution=ParametricVineDensity).fit(X)
  margins: Sequence[Any] = est.distribution_.margins
  assert all(isinstance(margin, SciPyMargin) for margin in margins)
  assert all(margin.is_fitted for margin in margins)
  rows = est.margin_summary_
  assert [row["family"] for row in rows] == [m.family_name for m in margins]
  assert all(row["margin"] == "SciPyMargin" for row in rows)
  assert all(row["loglik"] is not None for row in rows)


def test_the_estimator_honors_a_named_family_and_chooses_an_unnamed_one() -> (
  None
):
  """The estimator asks each margin to `select`, which respects a named family.

  `select` is the weaker requirement: a margin with nothing to choose reduces
  it to its own `fit`, and naming a family *is* the choice -- so the
  specification stands. Leaving the family out is what asks for the search.
  """
  pytest.importorskip("scipy")
  rs = np.random.RandomState(3)
  X = np.column_stack([rs.gamma(2.0, 1.0, 200), rs.gamma(3.0, 1.0, 200)])

  named = VineDensity(
    distribution=ParametricVineDensity,
    margin_controls=FitControlsMargin(family_set=["norm"]),
    random_state=0,
  ).fit(X)
  named_margins: Sequence[Any] = named.distribution_.margins
  assert [m.family_name for m in named_margins] == [
    "norm",
    "norm",
  ]

  chosen = VineDensity(distribution=ParametricVineDensity, random_state=0).fit(
    X
  )
  chosen_margins: Sequence[Any] = chosen.distribution_.margins
  assert "norm" not in [m.family_name for m in chosen_margins]


def test_selection_runs_per_expanded_column(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """Every expanded column is selected on its own, dummies included.

  A ``{0, 1}`` dummy is a count and the continuous columns are not, so the two
  get candidate sets of their own: a probability mass and a Lebesgue density
  are not comparable on one information criterion.
  """
  pytest.importorskip("scipy")
  X_df, expanded = sample_dataframe_data
  est = VineDensity(distribution=ParametricVineDensity, random_state=0).fit(
    X_df
  )
  rows = est.margin_summary_
  assert len(rows) == len(expanded)
  by_name = dict(zip(expanded, rows, strict=False))
  counts = ("poisson", "nbinom", "geom")
  assert by_name["cat1_B"]["var_type"] == "d"
  assert by_name["cat1_B"]["family"] in counts
  assert by_name["cont1"]["var_type"] == "c"
  assert by_name["cont1"]["family"] not in counts


def test_a_failed_family_search_names_its_column(cat_df: pd.DataFrame) -> None:
  """A search that refuses every candidate names the column, as a fit does.

  The levels of an ordered categorical are the variable's support, so the
  column is declared discrete and searched over the count families -- which
  none of these non-integer levels can be.
  """
  pytest.importorskip("scipy")
  with pytest.raises(ValueError, match=r"margin for 'grade': no parametric"):
    VineDensity(distribution=ParametricVineDensity).fit(cat_df)


# --- the response margin ---------------------------------------------------- #


def test_the_response_margin_comes_from_the_margin_class(
  regression_data: tuple[np.ndarray, np.ndarray, np.ndarray, float],
) -> None:
  """The response is variable zero, so it gets the same `margin_class`."""
  pytest.importorskip("scipy")
  X, y, _, _ = regression_data
  est = VineRegressor(distribution=ParametricVineDensity).fit(X, y)
  assert all(isinstance(m, SciPyMargin) for m in est.distribution_.margins)

  # And the default lane gives every variable, response included, a `Kde1d`.
  default = VineRegressor().fit(X, y)
  assert all(isinstance(m, Kde1d) for m in default.distribution_.margins)


def test_a_discrete_response_margin_is_refused(
  regression_data: tuple[np.ndarray, np.ndarray, np.ndarray, float],
) -> None:
  """The joint layout leads with the response and gives it no left limit.

  Checked twice, because a margin that chooses its own family has no variable
  type to declare until it has chosen one: the refusal on the specification
  cannot see a count family coming, and the one on the fitted margin can.
  """
  pytest.importorskip("scipy")
  X, y, _, _ = regression_data
  counts = np.round(np.abs(y) * 8.0)
  with pytest.raises(ValueError, match="response margin must be continuous"):
    VineRegressor(
      distribution=ParametricVineDensity,
      margin_controls={0: FitControlsMargin(family_set=["poisson"])},
      use_grid=False,
    ).fit(X, counts)

  class _Discrete(Vinedist):
    """A distribution whose margins declare atoms before anything is fitted."""

    margin_class = _DiscreteMargin

  with pytest.raises(ValueError, match="response margin must be continuous"):
    VineRegressor(distribution=_Discrete, use_grid=False).fit(X, counts)


def test_a_failing_margin_names_its_column(cat_df: pd.DataFrame) -> None:
  """A margin sees one array and cannot say which; the estimator can.

  `Kde1d` models a discrete variable on the integer lattice, so an ordered
  categorical whose levels are not integers cannot be one. Both the bound and
  the data are checked, and they fail at different points -- one while the
  declaration is validated, one while the margin is fitted -- so both have to
  name the column.
  """
  with pytest.raises(ValueError, match=r"margin for 'grade': discrete bounds"):
    VineDensity().fit(cat_df)

  # The second is reachable only from a declaration the input did not make,
  # which for an array is a pre-set `schema_`; the position is then the name.
  rs = np.random.RandomState(1)
  plain = np.column_stack([rs.normal(size=200), rs.normal(size=200)])
  est = VineDensity()
  est.schema_ = {
    "var_types": ["d", "c"],
    "supports": [None, None],
  }
  with pytest.raises(
    ValueError, match=r"margin for 'variable 0': discrete data"
  ):
    est.fit(plain)


def test_an_integer_categorical_is_fitted_on_its_declared_support() -> None:
  """The levels are the support, and the grid runs half a unit past each end.

  A jittered observation fills the cell around its level, so the grid has to
  cover `[min - 0.5, max + 0.5]`; snapping it to the levels themselves -- which
  is what happened before kde1d#37 -- truncates the outermost half-cells.
  """
  rs = np.random.RandomState(1)
  df = pd.DataFrame(
    {
      "a": rs.normal(size=300),
      "k": pd.Categorical(
        rs.choice([0, 1, 2, 3], 300), categories=[0, 1, 2, 3], ordered=True
      ),
    }
  )
  est = VineDensity().fit(df)
  margin: Any = est.distribution_.margins[1]
  assert est.schema_["supports"][1] == (0.0, 3.0)
  grid = np.asarray(margin.grid_points)
  assert grid[0] == pytest.approx(-0.5)
  assert grid[-1] == pytest.approx(3.5)
  # The masses still live on the levels, and nowhere else.
  assert margin.pdf(np.arange(4.0)).sum() == pytest.approx(1.0)
  assert margin.pdf(np.array([-1.0, 4.0])).tolist() == [0.0, 0.0]
