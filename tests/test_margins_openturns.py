"""Tests for the OpenTURNS margins.

Four contracts are pinned here. That the marshaling is right: OpenTURNS reads a
univariate argument as an ``(n, 1)`` ``Sample`` for a density and as a flat
``Point`` for a quantile, so the two conventions are opposites and only the
numbers can tell whether each was honored. That a discrete margin's left limit
is the *number* ``F(x^-)``, not merely a method that answers, including on a
support that is not the integer lattice. That the estimator and the selector
recover a family they were given data from, with every rejection recorded rather
than skipped. And that `pyvinecopulib.margins` imports, and `as_margin` keeps
working, without OpenTURNS installed at all.
"""

from __future__ import annotations


from typing import Any

import numpy as np
import pytest

from pyvinecopulib.core import Kde1d, MarginLike
from pyvinecopulib.margins import (
  OpenTURNSMargin,
  OpenTURNSSelector,
  as_margin,
)
from .helpers import run_without

openturns = pytest.importorskip("openturns")

#: The non-lattice discrete family, renamed in OpenTURNS 1.27.
FiniteDiscrete = (
  getattr(openturns, "FiniteDiscreteDistribution", None)
  or openturns.UserDefined
)


@pytest.fixture
def normal_sample() -> np.ndarray:
  """600 draws from a Normal(1, 2)."""
  return np.random.default_rng(0).normal(1.0, 2.0, size=600)


# --- marshaling ------------------------------------------------------------ #


def test_as_margin_matches_openturns_own_values() -> None:
  """`pdf` / `cdf` / `icdf` agree with OpenTURNS point by point."""
  dist = openturns.Normal(1.0, 2.0)
  margin: Any = as_margin(dist)
  x = np.array([-3.0, 0.0, 1.0, 2.5, 7.0])
  p = np.array([0.01, 0.1, 0.5, 0.9, 0.99])
  np.testing.assert_allclose(
    margin.pdf(x), [dist.computePDF(v) for v in x], rtol=0, atol=0
  )
  np.testing.assert_allclose(
    margin.cdf(x), [dist.computeCDF(v) for v in x], rtol=0, atol=0
  )
  np.testing.assert_allclose(
    margin.icdf(p),
    [dist.computeQuantile(q)[0] for q in p],
    rtol=0,
    atol=0,
  )


def test_as_margin_round_trips_cdf_and_icdf() -> None:
  """`cdf` inverts `icdf` on a continuous margin, and an empty batch survives."""
  margin: Any = as_margin(openturns.Gamma(2.5, 1.5, 0.0))
  p = np.array([0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999])
  np.testing.assert_allclose(margin.cdf(margin.icdf(p)), p, atol=1e-10)
  x = np.array([0.5, 1.0, 2.0, 6.0])
  np.testing.assert_allclose(margin.icdf(margin.cdf(x)), x, atol=1e-8)
  for empty in (margin.pdf(np.array([])), margin.icdf(np.array([]))):
    assert empty.shape == (0,)


def test_as_margin_reads_the_support_off_the_range() -> None:
  """Support bounds follow the range's finiteness flags, not its numbers.

  `support` is an optional capability rather than a `MarginLike` member, so it
  is read off a loosely typed binding -- the same idiom the other margin tests
  use.
  """
  bounded: Any = as_margin(openturns.Uniform(-1.0, 3.0))
  assert bounded.support == (-1.0, 3.0)
  positive: Any = as_margin(openturns.Gamma(2.0, 1.5, 0.0))
  assert positive.support == (0.0, float("inf"))
  real_line: Any = as_margin(openturns.Normal())
  assert real_line.support == (float("-inf"), float("inf"))
  counts: Any = as_margin(openturns.Binomial(10, 0.3))
  assert counts.support == (0.0, 10.0)


def test_as_margin_is_idempotent() -> None:
  """Coercing a margin returns it unchanged rather than re-wrapping it."""
  margin: Any = as_margin(openturns.Normal(1.0, 2.0))
  assert isinstance(margin, MarginLike)
  assert as_margin(margin) is margin


def test_as_margin_still_rejects_an_unrelated_object() -> None:
  """The OpenTURNS predicate does not claim objects from elsewhere."""
  with pytest.raises(TypeError, match="cannot use 'object' as a margin"):
    as_margin(object())


def test_as_margin_rejects_a_multivariate_distribution() -> None:
  """A margin is univariate, so a 2-d OpenTURNS distribution is refused."""
  with pytest.raises(ValueError, match="dimension 2"):
    as_margin(
      openturns.Normal([0.0, 0.0], [1.0, 1.0], openturns.CorrelationMatrix(2))
    )


# --- discrete left limits --------------------------------------------------- #


def test_discrete_left_limit_is_the_mass_below_the_atom() -> None:
  """`cdf_left` is `F(x^-)` numerically, at an atom and between atoms."""
  dist = openturns.Poisson(3.0)
  margin: Any = as_margin(dist)
  assert margin.var_type == "d"
  atoms = np.array([0.0, 1.0, 3.0, 9.0])
  np.testing.assert_allclose(
    margin.cdf_left(atoms),
    [0.0] + [dist.computeCDF(v - 1.0) for v in atoms[1:]],
    atol=1e-12,
  )
  # F(3^-) = P(X <= 2) for a Poisson(3).
  np.testing.assert_allclose(
    margin.cdf_left(np.array([3.0])), 0.4231900811, atol=1e-10
  )
  # Away from an atom there is no jump, so the left limit is the value.
  off = np.array([2.5, 7.25])
  np.testing.assert_allclose(margin.cdf_left(off), margin.cdf(off), atol=0)


def test_discrete_left_limit_does_not_assume_the_integer_lattice() -> None:
  """A support that is not the integers still gets the exact left limit."""
  dist = FiniteDiscrete(
    openturns.Sample([[0.5], [1.25], [3.0]]), [0.2, 0.3, 0.5]
  )
  margin: Any = as_margin(dist)
  assert margin.var_type == "d"
  np.testing.assert_allclose(
    margin.cdf_left(np.array([0.5, 1.25, 3.0])), [0.0, 0.2, 0.5], atol=1e-12
  )
  np.testing.assert_allclose(margin.cdf(np.array([1.25])), [0.5], atol=1e-12)


# --- OpenTURNSMargin -------------------------------------------------------- #


def test_estimator_recovers_the_generating_parameters(
  normal_sample: np.ndarray,
) -> None:
  """Fitting the true family returns its parameters and counts them all."""
  margin = OpenTURNSMargin("Normal").fit(normal_sample)
  assert margin.family_name == "Normal"
  assert margin.parameter_names == ("mu_0", "sigma_0")
  np.testing.assert_allclose(margin.parameters, (1.0, 2.0), atol=0.2)
  assert margin.n_parameters == 2.0
  np.testing.assert_allclose(
    margin.loglik(), np.sum(margin.logpdf(normal_sample)), atol=0
  )


def test_estimator_accepts_a_factory_object(count_sample: np.ndarray) -> None:
  """A `DistributionFactory` works wherever its name does."""
  by_name = OpenTURNSMargin("Poisson").fit(count_sample)
  by_object = OpenTURNSMargin(openturns.PoissonFactory()).fit(count_sample)
  np.testing.assert_allclose(by_object.parameters, by_name.parameters, atol=0)
  np.testing.assert_allclose(
    count_sample.mean(), by_name.parameters[0], atol=1e-8
  )


def test_estimator_knows_its_variable_type_before_it_is_fitted() -> None:
  """Discreteness is a property of the family, so it needs no data."""
  assert OpenTURNSMargin("Poisson").var_type == "d"
  assert OpenTURNSMargin("Gamma").var_type == "c"
  unfitted = OpenTURNSMargin("Gamma")
  assert not unfitted.is_fitted
  assert unfitted.support == (float("-inf"), float("inf"))
  with pytest.raises(RuntimeError, match="is not fitted"):
    unfitted.pdf(np.array([1.0]))


def test_estimator_from_distribution_estimated_nothing() -> None:
  """A margin given its parameters is fitted, and has none to its name."""
  margin = OpenTURNSMargin.from_distribution(openturns.Normal(1.0, 2.0))
  assert margin.is_fitted
  assert margin.n_parameters == 0.0
  with pytest.raises(TypeError, match="already carries its parameters"):
    margin.fit(np.zeros(10))


def test_estimator_samples_reproducibly(normal_sample: np.ndarray) -> None:
  """Seeded draws repeat, and land inside the fitted support."""
  margin = OpenTURNSMargin("Gamma").fit(np.abs(normal_sample) + 0.1)
  first = margin.sample(50, seeds=[7])
  np.testing.assert_allclose(margin.sample(50, seeds=[7]), first, atol=0)
  lo, hi = margin.support
  assert np.all((first >= lo) & (first <= hi))


def test_estimator_rejects_weights_and_unknown_families(
  normal_sample: np.ndarray,
) -> None:
  """Weights and unknown names fail loudly rather than quietly."""
  with pytest.raises(TypeError, match="cannot use observation weights"):
    OpenTURNSMargin("Normal").fit(normal_sample, weights=np.ones(600))
  with pytest.raises(ValueError, match="unknown OpenTURNS factory"):
    OpenTURNSMargin("NotAFamily")
  with pytest.raises(ValueError, match="not both and not neither"):
    OpenTURNSMargin()


# --- OpenTURNSSelector ------------------------------------------------------ #


def test_selector_recovers_the_generating_family(
  normal_sample: np.ndarray,
) -> None:
  """Among competing candidates, the family the data came from wins."""
  selector = OpenTURNSSelector(
    ["Normal", "Uniform", "Logistic", "Laplace"], criterion="bic"
  ).fit(normal_sample)
  assert selector.family_name == "Normal"
  assert selector.var_type == "c"
  np.testing.assert_allclose(
    selector.selected_.parameters, (1.0, 2.0), atol=0.2
  )
  np.testing.assert_allclose(
    selector.cdf(selector.icdf(np.array([0.25, 0.75]))),
    [0.25, 0.75],
    atol=1e-10,
  )
  assert len(selector.sample(20, seeds=[3])) == 20


def test_selector_recovers_a_discrete_family(count_sample: np.ndarray) -> None:
  """Counts are fitted with the discrete families, and Poisson data pick it."""
  selector = OpenTURNSSelector(criterion="bic", name="k").fit(count_sample)
  assert selector.family_name == "Poisson"
  assert selector.var_type == "d"
  np.testing.assert_allclose(
    selector.cdf_left(np.array([4.0])),
    selector.cdf(np.array([3.0])),
    atol=1e-12,
  )
  assert {row["column"] for row in selector.report_} == {"k"}


def test_selector_criteria_are_on_the_usual_scale(
  count_sample: np.ndarray,
) -> None:
  """OpenTURNS reports per observation; the rows carry the total."""
  selector = OpenTURNSSelector(criterion="bic").fit(count_sample)
  (winner,) = [row for row in selector.report_ if row["selected"]]
  n, k = count_sample.size, winner["n_parameters"]
  np.testing.assert_allclose(
    winner["bic"], -2.0 * winner["loglik"] + k * np.log(n), rtol=1e-9
  )
  np.testing.assert_allclose(
    winner["aic"], -2.0 * winner["loglik"] + 2.0 * k, rtol=1e-9
  )
  np.testing.assert_allclose(
    winner["aicc"],
    winner["aic"] + 2.0 * k * (k + 1.0) / (n - k - 1.0),
    rtol=1e-9,
  )
  assert winner["bic"] == min(row["bic"] for row in selector.report_)


def test_selector_records_every_rejection(count_sample: np.ndarray) -> None:
  """A candidate that cannot be built gets a row with its reason."""
  selector = OpenTURNSSelector(criterion="bic").fit(count_sample)
  status = {row["family"]: row["status"] for row in selector.report_}
  # Poisson(4) counts are neither 0/1 nor constant nor strictly positive.
  for family in ("Bernoulli", "Dirac", "Geometric"):
    assert family in status
    assert status[family] not in ("ok", "selected")
  assert status["Poisson"] == "selected"


def test_selector_never_compares_masses_against_densities(
  count_sample: np.ndarray,
) -> None:
  """A candidate on the wrong side of the discrete split cannot win."""
  with pytest.warns(UserWarning, match="no OpenTURNS family was admissible"):
    selector = OpenTURNSSelector(
      ["Normal"], criterion="bic", var_type="d", name="k"
    ).fit(count_sample)
  (row,) = [r for r in selector.report_ if r["family"] == "Normal"]
  assert "not comparable" in row["status"]
  assert isinstance(selector.selected_, Kde1d)
  assert selector.report_[-1]["status"] == "fallback"
  assert selector.var_type == "d"


def test_selector_rejects_bad_arguments(normal_sample: np.ndarray) -> None:
  """Argument checks happen before any fitting."""
  with pytest.raises(ValueError, match="unknown criterion"):
    OpenTURNSSelector(criterion="mdl")
  with pytest.raises(ValueError, match="unknown var_type"):
    OpenTURNSSelector(var_type="count")
  with pytest.raises(ValueError, match="unknown candidates"):
    OpenTURNSSelector("everything").fit(normal_sample)
  with pytest.raises(TypeError, match="cannot use observation weights"):
    OpenTURNSSelector().fit(normal_sample, weights=np.ones(600))


# --- absence of OpenTURNS --------------------------------------------------- #


def test_margins_import_and_coerce_without_openturns() -> None:
  """Without OpenTURNS, the package imports and `as_margin` still works."""
  run_without(
    "openturns",
    "import numpy as np\n"
    "from pyvinecopulib.core import Kde1d\n"
    "from pyvinecopulib.margins import as_margin\n"
    "m = as_margin(Kde1d().fit(np.arange(10.0)))\n"
    "assert m.cdf(np.array([5.0])).size == 1\n"
    "try:\n"
    "  as_margin(object())\n"
    "except TypeError:\n"
    "  sys.exit(0)\n"
    "sys.exit(3)\n",
  )


def test_openturns_margin_names_the_extra_without_openturns() -> None:
  """Without OpenTURNS, constructing one says which extra provides it."""
  run_without(
    "openturns",
    "from pyvinecopulib.margins import OpenTURNSMargin\n"
    "try:\n"
    "  OpenTURNSMargin('Normal')\n"
    "except ImportError as e:\n"
    "  sys.exit(0 if 'pyvinecopulib[openturns]' in str(e) else 2)\n"
    "sys.exit(3)\n",
  )
