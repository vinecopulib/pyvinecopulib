"""Tests for the ``margins=`` half of the sklearn estimators.

What this file pins is the delegation. That an estimator *is* a `Vinedist` --
the same layout, the same log-density, the same draws -- rather than a second
implementation of Sklar's theorem that can drift from the first. That
`margins=None` still means what it meant before margins were configurable: a
`Kde1d` per column, to the last bit. That a specification the caller gives is
honored per column and never mutated, so `clone` reproduces the estimator. And
that the ``{0, 1}`` dummies of an expanded unordered categorical are fitted on
the support they actually have, rather than on a padded grid that puts mass on
values that cannot occur.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("pandas")

import pandas as pd  # noqa: E402  (import order: after the availability check)

from sklearn.base import clone  # noqa: E402

import pyvinecopulib as pv  # noqa: E402
from pyvinecopulib.core import Kde1d  # noqa: E402
from pyvinecopulib.margins import (  # noqa: E402
  MarginSelector,
  ParametricMargin,
)
from pyvinecopulib.sklearn import VineDensity, VineRegressor  # noqa: E402

from .helpers import AtomicMargin  # noqa: E402


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
  """Wrapping the vine for the backend does not hide the vine's own surface."""
  X, _, _ = sample_array_data
  est = VineDensity().fit(X)
  copula = est.distribution_.copula
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
  assert est.schema_["bounds"][0] == (
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
  """A `margins=` argument does not cost the margin what the input declared.

  The per-variable default only covers the columns a specification leaves
  unaddressed, so a broadcast alias used to hand every selector
  `var_type=None, bounds=None` and have it re-infer both from the sample --
  strictly less than `schema_` already knew.
  """
  rng = np.random.default_rng(0)
  counts = rng.integers(0, 5, size=400)
  X_df = pd.DataFrame(
    {
      "grade": pd.Categorical(counts, categories=range(5), ordered=True),
      "z": rng.normal(size=400),
    }
  )
  est = VineDensity(margins="parametric", random_state=0).fit(X_df)
  selector: Any = est.distribution_.margins[0]
  assert selector.bounds == (0.0, 4.0)
  assert selector.var_type == "d"
  # The unbounded continuous column is left alone.
  other: Any = est.distribution_.margins[1]
  assert other.bounds is None


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


def test_fit_does_not_mutate_the_margins_argument(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """Specifications are fitted on copies, so the parameters survive `fit`."""
  X, _, _ = sample_array_data
  spec = Kde1d()
  est = VineDensity(margins=spec).fit(X)
  assert est.get_params()["margins"] is spec
  assert not spec.is_fitted
  fitted: tuple[Any, ...] = est.distribution_.margins
  assert all(margin.is_fitted for margin in fitted)
  # Refitting therefore re-estimates rather than reusing the first fit.
  est.fit(X[:100])
  assert not spec.is_fitted


def test_clone_round_trips_the_margins_parameter(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """`clone` reproduces an estimator that fits the same margins."""
  X, _, _ = sample_array_data
  est = VineDensity(margins="parametric").fit(X)
  cloned = clone(est)
  assert cloned.margins == "parametric"
  np.testing.assert_allclose(
    cloned.fit(X).score_samples(X[:10]), est.score_samples(X[:10])
  )


def test_a_density_less_margin_is_refused_at_fit_time(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """A margin that reports no density cannot serve an estimator that does.

  Caught at `fit`, not at the first `score_samples`, so the estimator never
  reaches a state where its own contract cannot be met. Every shipped margin
  has a density, so the guard is exercised through a double -- a third-party
  margin may well be atomic, and this is what turns that into an error here
  instead of a wrong number later.
  """
  X, _, _ = sample_array_data
  with pytest.raises(ValueError, match="needs a density"):
    VineDensity(margins=AtomicMargin()).fit(X)
  with pytest.raises(ValueError, match='margins="kde"'):
    VineDensity(margins=AtomicMargin()).fit(X)
  # `VineRegressor` reads the copula density and the response `icdf` only, so
  # the same margin is fine there.
  VineRegressor(margins=AtomicMargin(), use_grid=False).fit(X, X[:, 0])


def test_a_wrong_length_sequence_is_refused(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """A per-column sequence must carry one entry per column."""
  X, _, _ = sample_array_data
  with pytest.raises(ValueError, match="length 1, but there are 2"):
    VineDensity(margins=[Kde1d()]).fit(X)


def test_a_mapping_leaves_the_other_columns_on_the_inferred_default(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """Addressing one column must not silently retype the rest."""
  X_df, expanded = sample_dataframe_data
  est = VineDensity(margins={"cont1": ParametricMargin("norm")}).fit(X_df)
  margins = est.distribution_.margins
  assert isinstance(margins[0], ParametricMargin)
  dummy: Any = margins[expanded.index("cat1_B")]
  assert isinstance(dummy, Kde1d)
  assert dummy.var_type == "d"
  assert dummy.support == (0.0, 1.0)


def test_a_preset_schema_supplies_what_an_array_cannot_carry() -> None:
  """`schema_` declares per-column variable types and bounds before `fit`."""
  rng = np.random.default_rng(0)
  X = np.column_stack([rng.poisson(3.0, size=200) * 1.0, rng.normal(size=200)])
  est = VineDensity()
  est.schema_ = {
    "kde1d_types": ["discrete", "continuous"],
    "bounds": [(0.0, 20.0), None],
  }
  est.fit(X)
  counts: Any = est.distribution_.margins[0]
  assert counts.var_type == "d"
  assert counts.support == (0.0, 20.0)
  assert est.distribution_.var_types == ["d", "c"]


def test_a_fixed_foreign_margin_is_used_as_given(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """An already-fitted margin is not re-estimated."""
  stats = pytest.importorskip("scipy.stats")
  X, _, _ = sample_array_data
  est = VineDensity(margins=[stats.norm(0.0, 1.0), stats.norm(0.0, 1.0)]).fit(X)
  np.testing.assert_allclose(
    np.asarray(est.distribution_.margins[0].cdf(np.array([0.0]))), 0.5
  )


# --- family selection ------------------------------------------------------- #


def test_selection_report_is_empty_without_a_selector(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """Nothing selected a family, so the table is empty rather than absent."""
  X, _, _ = sample_array_data
  assert VineDensity().fit(X).selection_report_ == []


def test_parametric_margins_report_their_selection(
  sample_array_data: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
  """`margins="parametric"` selects per column and names the column it fitted."""
  pytest.importorskip("scipy")
  X, _, _ = sample_array_data
  est = VineDensity(margins="parametric").fit(X)
  assert all(isinstance(m, MarginSelector) for m in est.distribution_.margins)
  assert {row["column"] for row in est.selection_report_} == {"x0", "x1"}
  assert "selected" in {row["status"] for row in est.selection_report_}


def test_selection_report_names_dataframe_columns(
  sample_dataframe_data: tuple[pd.DataFrame, list[str]],
) -> None:
  """With named features, the report is keyed by the expanded column name."""
  pytest.importorskip("scipy")
  X_df, expanded = sample_dataframe_data
  est = VineDensity(margins="parametric").fit(X_df)
  assert {row["column"] for row in est.selection_report_} == set(expanded)


# --- the response margin ---------------------------------------------------- #


def test_the_response_margin_follows_a_broadcast_specification(
  regression_data: tuple[np.ndarray, np.ndarray, np.ndarray, float],
) -> None:
  """An alias covers the response too; a per-column sequence does not."""
  X, y, _, _ = regression_data
  broadcast = VineRegressor(margins="parametric").fit(X, y)
  assert isinstance(broadcast.distribution_.margins[0], MarginSelector)

  per_column = VineRegressor(margins=[MarginSelector(), MarginSelector()]).fit(
    X, y
  )
  assert isinstance(per_column.distribution_.margins[0], Kde1d)
  assert isinstance(per_column.distribution_.margins[1], MarginSelector)


def test_a_discrete_response_margin_is_refused(
  regression_data: tuple[np.ndarray, np.ndarray, np.ndarray, float],
) -> None:
  """The joint layout leads with the response and gives it no left limit."""
  X, y, _, _ = regression_data
  with pytest.raises(ValueError, match="response margin must be continuous"):
    VineRegressor(margins=Kde1d(type="discrete"), use_grid=False).fit(
      X, np.round(y)
    )


def test_a_failing_margin_names_its_column() -> None:
  """A margin sees one array and cannot say which; the estimator can.

  `Kde1d` models a discrete variable on the integer lattice, so an ordered
  categorical whose levels are not integers cannot be one. Both the bound and
  the data are checked, and they fail at different points -- one while the
  specification is built, one while it is fitted -- so both have to name the
  column.
  """
  rs = np.random.RandomState(0)
  df = pd.DataFrame(
    {
      "a": rs.normal(size=200),
      "grade": pd.Categorical(
        rs.choice([1.5, 2.5, 3.5], 200),
        categories=[1.5, 2.5, 3.5],
        ordered=True,
      ),
    }
  )
  with pytest.raises(ValueError, match=r"margin for 'grade': discrete bounds"):
    VineDensity().fit(df)

  plain = pd.DataFrame({"a": rs.normal(size=200), "b": rs.normal(size=200)})
  with pytest.raises(ValueError, match=r"margin for 'a': discrete data"):
    VineDensity(margins=Kde1d(type="discrete")).fit(plain)


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
  margin = est.distribution_.margins[1]
  assert est.schema_["bounds"][1] == (0.0, 3.0)
  grid = np.asarray(margin.grid_points)
  assert grid[0] == pytest.approx(-0.5)
  assert grid[-1] == pytest.approx(3.5)
  # The masses still live on the levels, and nowhere else.
  assert margin.pdf(np.arange(4.0)).sum() == pytest.approx(1.0)
  assert margin.pdf(np.array([-1.0, 4.0])).tolist() == [0.0, 0.0]


def test_a_named_margin_survives_a_default_that_would_refuse_the_column() -> (
  None
):
  """The default is only built where it is needed.

  A column can be one the default margin refuses -- an ordered categorical with
  non-integer levels cannot be a discrete `Kde1d`. A caller who names a margin
  for every column has answered that already, and should not be stopped by a
  default their specification never uses.
  """
  rs = np.random.RandomState(2)
  df = pd.DataFrame(
    {
      "a": rs.normal(size=200),
      "grade": pd.Categorical(
        rs.choice([1.5, 2.5, 3.5], 200),
        categories=[1.5, 2.5, 3.5],
        ordered=True,
      ),
    }
  )
  est = VineDensity(margins=Kde1d()).fit(df)
  assert len(est.distribution_.margins) == 2
