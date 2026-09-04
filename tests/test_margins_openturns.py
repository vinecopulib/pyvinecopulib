"""Tests for the OpenTURNS margins.

Five contracts are pinned here. That the marshaling is right: OpenTURNS reads a
univariate argument as an ``(n, 1)`` ``Sample`` for a density and as a flat
``Point`` for a quantile, so the two conventions are opposites and only the
numbers can tell whether each was honored. That a discrete margin's left limit
is the *number* ``F(x^-)``, not merely a method that answers, including on a
support that is not the integer lattice. That `fit` recovers the parameters of
the family it was named and `select` recovers the family itself, refusing every
candidate on the wrong side of the discrete split and naming the ones it
refused rather than skipping them. That `margin_controls` reaches each
variable's own search when one of these margins is a vine distribution's. And
that `pyvinecopulib.margins` imports, and `as_margin` keeps working, without
OpenTURNS installed at all.
"""

from __future__ import annotations


from typing import Any, Callable

import numpy as np
import pytest

from pyvinecopulib.core import Kde1d, MarginLike, Vinedist
from pyvinecopulib.margins import (
  FitControlsMargin,
  OpenTURNSMargin,
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


# --- one named family ------------------------------------------------------- #


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
  with pytest.raises(ValueError, match="not both"):
    OpenTURNSMargin("Normal", distribution=openturns.Normal())


@pytest.mark.parametrize(
  "fitter",
  [
    lambda y: OpenTURNSMargin("Normal").fit(y),
    lambda y: OpenTURNSMargin().select(
      y, FitControlsMargin(family_set=["Normal"])
    ),
  ],
  ids=["fit", "select"],
)
@pytest.mark.parametrize("shape", [(4, 1), (2, 2)])
def test_openturns_fitters_require_a_univariate_shape(
  fitter: Callable[[np.ndarray], Any], shape: tuple[int, int]
) -> None:
  """Column matrices must not be flattened into a pooled sample."""
  with pytest.raises(ValueError, match=r"y must have shape \(n,\)"):
    fitter(np.arange(np.prod(shape), dtype=float).reshape(shape))


# --- family selection ------------------------------------------------------- #


def test_an_unnamed_margin_has_no_family_until_it_selects_one(
  normal_sample: np.ndarray,
) -> None:
  """Constructed with neither argument, the family is `select`'s to choose."""
  margin = OpenTURNSMargin()
  assert not margin.is_fitted
  with pytest.raises(RuntimeError, match="no family yet"):
    margin.family_name
  with pytest.raises(RuntimeError, match="no family yet"):
    margin.pdf(np.array([1.0]))
  chosen = margin.select(
    normal_sample, FitControlsMargin(family_set=["Normal", "Laplace"])
  )
  # `select` returns the margin it was called on, now carrying the winner:
  # the margin *is* the selected family rather than holding one.
  assert chosen is margin
  assert margin.family_name == "Normal"
  assert margin.is_fitted


def test_select_recovers_the_generating_family(
  normal_sample: np.ndarray,
) -> None:
  """Among competing candidates, the family the data came from wins."""
  families = ["Normal", "Uniform", "Logistic", "Laplace"]
  margin = OpenTURNSMargin().select(
    normal_sample,
    FitControlsMargin(family_set=families, selection_criterion="bic"),
  )
  assert margin.family_name == "Normal"
  assert margin.var_type == "c"
  np.testing.assert_allclose(margin.parameters, (1.0, 2.0), atol=0.2)
  np.testing.assert_allclose(
    margin.cdf(margin.icdf(np.array([0.25, 0.75]))),
    [0.25, 0.75],
    atol=1e-10,
  )
  assert len(margin.sample(20, seeds=[3])) == 20
  # The winner is the argmin on the criterion that was asked for, which is
  # checkable against the candidates fitted one at a time.
  scored = {
    family: OpenTURNSMargin(family).fit(normal_sample).bic(normal_sample)
    for family in families
  }
  assert min(scored, key=lambda family: scored[family]) == margin.family_name


def test_select_recovers_a_discrete_family(count_sample: np.ndarray) -> None:
  """Counts search the discrete registry, and Poisson data pick Poisson."""
  margin = OpenTURNSMargin().select(
    count_sample, FitControlsMargin(selection_criterion="bic")
  )
  assert margin.family_name == "Poisson"
  assert margin.var_type == "d"
  np.testing.assert_allclose(
    margin.cdf_left(np.array([4.0])), margin.cdf(np.array([3.0])), atol=1e-12
  )


def test_criteria_are_on_the_usual_scale(count_sample: np.ndarray) -> None:
  """OpenTURNS reports per observation; a margin reports the sample total."""
  margin = OpenTURNSMargin("Poisson").fit(count_sample)
  sample = openturns.Sample(count_sample.reshape(-1, 1))
  n, k = float(count_sample.size), int(margin.n_parameters)
  for name, method in (("aic", "AIC"), ("bic", "BIC"), ("aicc", "AICC")):
    test = getattr(openturns.FittingTest, method)
    per_observation = float(test(sample, margin.distribution, k))
    np.testing.assert_allclose(
      getattr(margin, name)(count_sample), n * per_observation, rtol=1e-9
    )
  loglik = float(margin.loglik())
  np.testing.assert_allclose(margin.aic(), -2.0 * loglik + 2.0 * k, rtol=1e-12)
  np.testing.assert_allclose(
    margin.bic(count_sample), -2.0 * loglik + k * np.log(n), rtol=1e-12
  )
  np.testing.assert_allclose(
    margin.aicc(count_sample),
    margin.aic() + 2.0 * k * (k + 1.0) / (n - k - 1.0),
    rtol=1e-12,
  )
  # The fit records its sample size, so a penalized criterion answers from it.
  assert margin.nobs == count_sample.size
  np.testing.assert_allclose(margin.bic(), margin.bic(count_sample), rtol=0)


def test_select_honors_a_declared_variable_type(
  count_sample: np.ndarray,
) -> None:
  """Caller schema wins over the integer-valued-data heuristic."""
  families = ["Normal", "Poisson"]
  as_counts = OpenTURNSMargin().select(
    count_sample, FitControlsMargin(family_set=families, var_type="d")
  )
  assert (as_counts.family_name, as_counts.var_type) == ("Poisson", "d")
  as_continuous = OpenTURNSMargin().select(
    count_sample, FitControlsMargin(family_set=families, var_type="c")
  )
  assert (as_continuous.family_name, as_continuous.var_type) == ("Normal", "c")


def test_select_never_compares_masses_against_densities(
  count_sample: np.ndarray,
) -> None:
  """A candidate on the wrong side of the discrete split cannot win.

  A probability mass and a Lebesgue density are not comparable on one
  criterion, so a continuous family on a variable declared discrete is refused
  by name instead of being ranked against the families that could have won.
  """
  with pytest.raises(ValueError, match="not comparable") as refusal:
    OpenTURNSMargin().select(
      count_sample, FitControlsMargin(family_set=["Normal"], var_type="d")
    )
  assert "Normal" in str(refusal.value)


def test_select_names_every_refused_candidate(
  count_sample: np.ndarray,
) -> None:
  """A candidate that cannot be built is reported with its reason."""
  with pytest.raises(
    ValueError, match="every candidate was refused"
  ) as refused:
    OpenTURNSMargin().select(
      count_sample,
      FitControlsMargin(
        family_set=["Bernoulli", "Dirac", "Geometric"], var_type="d"
      ),
    )
  message = str(refused.value)
  # Poisson(4) counts are neither 0/1 nor constant nor strictly positive, so
  # none of the three can be built and each says why.
  for family in ("Bernoulli", "Dirac", "Geometric"):
    assert family in message
  assert "narrow family_set" in message


def test_from_data_selects_a_family_where_fit_keeps_the_named_one(
  count_sample: np.ndarray,
) -> None:
  """`from_data` reaches `select`, so a wrong family does not survive it."""
  chosen = OpenTURNSMargin.from_data(count_sample)
  assert chosen.is_fitted
  assert chosen.family_name == "Poisson"
  # `fit` estimates the family it was given, however wrong; `select` replaces
  # it. That difference is the whole reason both verbs exist.
  families = FitControlsMargin(family_set=["Normal", "Poisson"])
  assert OpenTURNSMargin("Normal").fit(count_sample).family_name == "Normal"
  assert (
    OpenTURNSMargin("Normal").select(count_sample, families).family_name
    == "Poisson"
  )


def test_select_rejects_bad_arguments(normal_sample: np.ndarray) -> None:
  """Argument checks happen before any fitting."""
  with pytest.raises(ValueError, match="unknown OpenTURNS factory"):
    OpenTURNSMargin().select(
      normal_sample, FitControlsMargin(family_set=["everything"])
    )
  with pytest.raises(TypeError, match="cannot use observation weights"):
    OpenTURNSMargin().select(normal_sample, weights=np.ones(600))


# --- inside a vine distribution --------------------------------------------- #


def test_margin_controls_configure_each_variable_separately() -> None:
  """A mapping keyed by variable name configures one variable's search.

  The two columns are on different families, so a specification that was
  broadcast to both, or routed to the wrong one, would name the wrong family.
  """
  rng = np.random.default_rng(5)
  y = np.column_stack(
    [rng.gamma(2.5, 1.5, size=500), rng.normal(1.0, 2.0, size=500)]
  )
  dist = Vinedist.from_data(
    y,
    margins=OpenTURNSMargin(),
    margin_controls={
      "positive": FitControlsMargin(family_set=["Gamma"]),
      "real": FitControlsMargin(
        family_set=["Normal", "Laplace"], selection_criterion="bic"
      ),
    },
    names=("positive", "real"),
  )
  margins: tuple[Any, ...] = dist.margins
  assert [margin.family_name for margin in margins] == ["Gamma", "Normal"]
  # The fitted Gamma bounds the variable below, and the draws respect it.
  lo, _ = margins[0].support
  assert np.all(dist.sample(200, seeds=[2])[:, 0] >= lo)
  assert np.all(np.isfinite(dist.logpdf(y[:20])))


def test_fallback_substitutes_a_kernel_density_margin() -> None:
  """`on_failure="fallback"` answers an impossible search with a `Kde1d`.

  A search that refuses every candidate is decided where the margin for the
  variable is chosen, rather than inside a margin that would have to stop
  being an OpenTURNS one to decide it -- so the substitution is reachable
  through a fit, with one warning naming the cause.
  """
  rng = np.random.default_rng(1)
  y = np.column_stack(
    [rng.poisson(4.0, size=400).astype(float), rng.normal(1.0, 2.0, size=400)]
  )
  substituted = "kernel-density margin was substituted"
  with pytest.warns(UserWarning, match=substituted):
    dist = Vinedist.from_data(
      y,
      margins=OpenTURNSMargin(),
      margin_controls={
        0: FitControlsMargin(
          family_set=["Normal"], var_type="d", on_failure="fallback"
        ),
        1: FitControlsMargin(family_set=["Normal"]),
      },
    )
  margins: tuple[Any, ...] = dist.margins
  # The declared type survives the substitution: the counts keep their atoms.
  assert isinstance(margins[0], Kde1d)
  assert margins[0].var_type == "d"
  assert margins[1].family_name == "Normal"


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
