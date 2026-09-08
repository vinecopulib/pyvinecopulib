"""The margin contract, checked once against every ecosystem adapter.

`SciPyMargin` and `OpenTURNSMargin` wrap different family registries but answer
the same questions, so the answers are pinned here over the `ECOSYSTEMS` table
and a third adapter inherits the suite by adding a row.

What stays in `test_margins_scipy.py` / `test_margins_openturns.py` is what is
genuinely one ecosystem's: SciPy's curated candidate set, its traps and its
admissibility gate; OpenTURNS' marshaling conventions and its
`DistributionFactory`; and each lane's own refusal messages -- including the
discrete split, which SciPy pre-filters with a warning where OpenTURNS refuses
candidate by candidate.
"""

from __future__ import annotations

from typing import Any, NamedTuple

import numpy as np
import pytest

from pyvinecopulib.core import Kde1d, MarginLike, Vinedist
from pyvinecopulib.margins import (
  FitControlsMargin,
  OpenTURNSMargin,
  SciPyMargin,
)
from .helpers import widen


class Ecosystem(NamedTuple):
  """One adapter, and the family names its own registry spells its own way.

  ``real`` is a continuous family the positive sample did not come from,
  ``true`` the one it did, ``count`` what count data select, ``competitors``
  three families a criterion ranks; ``impossible`` is a ``family_set`` no
  column of the collapsed trap can satisfy, and ``reason`` the cause the
  fallback warning then carries.
  """

  cls: Any
  module: str
  real: str
  true: str
  count: str
  competitors: list[str]
  impossible: list[str]
  reason: str


ECOSYSTEMS = [
  Ecosystem(
    cls=SciPyMargin,
    module="scipy.stats",
    real="norm",
    true="gamma",
    count="poisson",
    competitors=["norm", "logistic", "laplace"],
    impossible=["t"],
    reason="t: degenerate parameter scale",
  ),
  Ecosystem(
    cls=OpenTURNSMargin,
    module="openturns",
    real="Normal",
    true="Gamma",
    count="Poisson",
    competitors=["Normal", "Logistic", "Laplace"],
    impossible=["Poisson"],
    reason="Poisson: TypeError",
  ),
]


@pytest.fixture(params=ECOSYSTEMS, ids=[eco.module for eco in ECOSYSTEMS])
def eco(request: pytest.FixtureRequest) -> Ecosystem:
  """One ecosystem's adapter, skipped when its package is not installed."""
  row: Ecosystem = request.param
  pytest.importorskip(row.module)
  return row


#: 500 draws from a gamma(2.5, 1.5), the running example of the design.
POSITIVE = np.random.default_rng(0).gamma(2.5, 1.5, size=500)


# --- a margin, fitted or not ------------------------------------------------ #


@pytest.mark.parametrize("verb", ["fit", "select"])
def test_a_fitted_margin_conforms_to_the_contract(
  eco: Ecosystem, verb: str
) -> None:
  """However a margin reached its family, it satisfies the margin contract.

  It inverts its own cdf, draws reproducibly inside its support, and reports a
  log-likelihood -- as a total over data handed to it, and as the fitted value
  when handed none. One method, two jobs, as on `Bicop`: a property would
  shadow the method and make the polymorphic `margin.loglik(sample)` call raise
  `TypeError`, which is why the fitted value is reached through the same name.
  """
  controls = FitControlsMargin(family_set=[eco.true])
  margin = (
    eco.cls(eco.true).fit(POSITIVE)
    if verb == "fit"
    else eco.cls().select(POSITIVE, controls)
  )
  assert isinstance(margin, MarginLike)
  assert margin.family_name == eco.true
  assert margin.var_type == "c"
  p = np.array([0.1, 0.5, 0.9])
  np.testing.assert_allclose(margin.cdf(margin.icdf(p)), p, atol=1e-10)

  total = margin.loglik(POSITIVE)
  assert np.ndim(total) == 0
  np.testing.assert_allclose(total, np.sum(margin.logpdf(POSITIVE)), atol=0)
  assert isinstance(margin.loglik(), float)

  draws = margin.sample(50, seeds=[7])
  np.testing.assert_array_equal(margin.sample(50, seeds=[7]), draws)
  lo, hi = margin.support
  assert np.all((draws >= lo) & (draws <= hi))


def test_an_unfitted_margin_answers_nothing(eco: Ecosystem) -> None:
  """A named margin has no numbers, and no support, until it is fitted."""
  margin = eco.cls(eco.real)
  assert not margin.is_fitted
  assert margin.support == (float("-inf"), float("inf"))
  with pytest.raises(RuntimeError, match="is not fitted"):
    margin.pdf(np.array([1.0]))
  with pytest.raises(RuntimeError, match="is not fitted"):
    margin.parameters
  with pytest.raises(RuntimeError, match="only defined after"):
    margin.loglik()


def test_an_unnamed_margin_has_no_family_until_select(eco: Ecosystem) -> None:
  """`cls()` is a request to choose a family, not a broken margin.

  It is what a "choose one" margin specification resolves to, so everything
  that reads a family has to say so rather than fail obscurely -- and the type
  it can represent is settled by the family it chooses.
  """
  margin = eco.cls()
  assert not margin.is_fitted
  with pytest.raises(RuntimeError, match="no family yet"):
    margin.family_name

  chosen = margin.select(
    POSITIVE, FitControlsMargin(family_set=[eco.real, eco.true])
  )
  # `select` returns the margin it was called on, now carrying the winner: the
  # margin *is* the selected family rather than holding one.
  assert chosen is margin
  assert margin.is_fitted
  assert margin.family_name == eco.true
  # Continuous data selected a continuous family; count data would have
  # selected a count one, which is what "either kind until chosen" means.
  assert margin.var_type == "c"


# --- fit, select and the criteria ------------------------------------------- #


def test_select_replaces_a_wrong_family_and_fit_keeps_it(
  eco: Ecosystem,
) -> None:
  """One method estimates the family it was given, the other also chooses it.

  The same margin object and the same controls: `fit` has nothing to search, so
  the `family_set` is inert there, while `select` adopts the winner and
  *becomes* it. That difference is the whole reason both verbs exist, and the
  log-likelihood is where it shows.
  """
  controls = FitControlsMargin(family_set=[eco.true])
  kept = eco.cls(eco.real).fit(POSITIVE, controls)
  assert kept.family_name == eco.real

  replaced = eco.cls(eco.real).select(POSITIVE, controls)
  assert replaced.family_name == eco.true
  assert replaced.is_fitted
  assert replaced.loglik() > kept.loglik()


def test_select_on_a_named_margin_keeps_the_family(eco: Ecosystem) -> None:
  """Naming a family *is* the choice, so `select` reduces to `fit`.

  It matters because `fit_margin` calls `select` by default, so without it a
  named margin comes back as whatever won the registry search -- answering a
  specification with a different model. `family_set` is how a caller asks for
  the search back on one.
  """
  from pyvinecopulib.core._resolve import fit_margin

  named = eco.cls(eco.real)
  assert not named.is_fitted
  assert named.select(POSITIVE).family_name == eco.real

  # And through the resolution path a vine distribution actually takes.
  resolved = widen(fit_margin(eco.cls(eco.real), POSITIVE))
  assert isinstance(resolved, eco.cls)
  assert resolved.family_name == eco.real

  reopened = eco.cls(eco.real).select(
    POSITIVE, FitControlsMargin(family_set=[eco.real, eco.true])
  )
  assert reopened.family_name == eco.true


def test_criteria_match_their_definitions(eco: Ecosystem) -> None:
  """Each criterion is its own formula in the fitted log-likelihood.

  The fit records its sample size, so a penalized criterion answers from it
  rather than needing the sample again.
  """
  m = eco.cls(eco.true).fit(POSITIVE)
  loglik, k, n = m.loglik(), m.n_parameters, float(m.nobs or 0)
  assert n == POSITIVE.size
  assert m.aic() == pytest.approx(-2.0 * loglik + 2.0 * k)
  assert m.bic() == pytest.approx(-2.0 * loglik + k * np.log(n))
  assert m.aicc() == pytest.approx(m.aic() + 2.0 * k * (k + 1.0) / (n - k - 1))
  assert m.bic() > m.aic()
  np.testing.assert_allclose(m.bic(), m.bic(POSITIVE), rtol=0)


@pytest.mark.parametrize("criterion", ["aic", "bic", "aicc"])
def test_the_selection_criterion_drives_the_winner(
  eco: Ecosystem, criterion: str
) -> None:
  """Whichever criterion is asked for is the one the search minimizes.

  The winner is the argmin over the candidates fitted one at a time, and
  reports that minimum as its own.
  """
  by_hand = {
    family: getattr(eco.cls(family).fit(POSITIVE), criterion)()
    for family in eco.competitors
  }
  winner = eco.cls().select(
    POSITIVE,
    FitControlsMargin(
      family_set=eco.competitors, selection_criterion=criterion
    ),
  )
  assert winner.family_name == min(by_hand, key=lambda f: by_hand[f])
  assert getattr(winner, criterion)() == pytest.approx(min(by_hand.values()))


def test_counts_select_a_count_family(
  eco: Ecosystem, count_sample: np.ndarray
) -> None:
  """Counts search the discrete candidates, and a declared type overrides them.

  `from_data` reaches `select`, so it chooses the family too, and the winner's
  `pdf` is a probability mass -- the jump `F(k) - F(k^-)`. A caller who
  declares the variable type beats the integer-valued-data heuristic in either
  direction, the curated set being partitioned by that type.
  """
  chosen = eco.cls.from_data(count_sample)
  assert chosen.family_name == eco.count
  assert chosen.var_type == "d"
  k = np.array([0.0, 1.0, 4.0])
  np.testing.assert_allclose(
    chosen.pdf(k), chosen.cdf(k) - chosen.cdf_left(k), atol=1e-12
  )

  forced_d = eco.cls().select(count_sample, FitControlsMargin(var_type="d"))
  assert forced_d.family_name == eco.count
  forced_c = eco.cls().select(count_sample, FitControlsMargin(var_type="c"))
  assert forced_c.var_type == "c"


def test_margin_controls_address_each_variable(eco: Ecosystem) -> None:
  """A mapping keyed by variable name configures one variable's own search.

  The two columns are on different families, so a specification that was
  broadcast to both, or routed to the wrong one, would name the wrong family --
  while a specification that *is* meant for both is broadcast by giving it
  once. The fitted positive family bounds its variable below, and the draws
  respect it.
  """
  rng = np.random.default_rng(5)
  y = np.column_stack(
    [rng.gamma(2.5, 1.5, size=500), rng.normal(1.0, 2.0, size=500)]
  )
  dist = Vinedist.from_data(
    y,
    margins=eco.cls(),
    margin_controls={
      "positive": FitControlsMargin(family_set=[eco.true]),
      "real": FitControlsMargin(
        family_set=eco.competitors, selection_criterion="bic"
      ),
    },
    names=("positive", "real"),
  )
  margins = [widen(m) for m in dist.margins]
  assert [m.family_name for m in margins] == [eco.true, eco.real]
  lo, _ = margins[0].support
  assert np.all(dist.sample(200, seeds=[2])[:, 0] >= lo)
  assert np.all(np.isfinite(dist.logpdf(y[:20])))

  shared = Vinedist.from_data(
    y,
    margins=eco.cls(),
    margin_controls=FitControlsMargin(
      family_set=[eco.real, eco.true], selection_criterion="bic"
    ),
  )
  assert [widen(m).family_name for m in shared.margins] == [eco.true, eco.real]


# --- refusals --------------------------------------------------------------- #


@pytest.mark.parametrize("verb", ["fit", "select"])
def test_margin_fitters_reject_weights(eco: Ecosystem, verb: str) -> None:
  """Neither registry's estimator weights observations, so asking must raise."""
  margin = eco.cls(eco.real) if verb == "fit" else eco.cls()
  with pytest.raises(TypeError, match="cannot use observation weights"):
    getattr(margin, verb)(POSITIVE, weights=np.ones_like(POSITIVE))


@pytest.mark.parametrize("verb", ["fit", "select"])
@pytest.mark.parametrize("shape", [(4, 1), (2, 2)])
def test_margin_fitters_require_a_univariate_shape(
  eco: Ecosystem, verb: str, shape: tuple[int, int]
) -> None:
  """Column matrices must not be flattened into a pooled sample."""
  margin = eco.cls(eco.real) if verb == "fit" else eco.cls()
  with pytest.raises(ValueError, match=r"y must have shape \(n,\)"):
    getattr(margin, verb)(np.arange(np.prod(shape), dtype=float).reshape(shape))


def test_on_failure_fallback_substitutes_a_kernel_density(
  eco: Ecosystem,
) -> None:
  """When nothing is admissible the fallback is nonparametric, never a family.

  Two columns, half exact zeros each, on which every candidate asked for is
  refused. Substituting another kind of margin is a decision about which margin
  the column gets, so it is made where the margin is chosen rather than inside
  one that would have to stop being parametric to make it -- which is what puts
  it on this path, with one warning per column carrying its own cause.
  """
  rng = np.random.default_rng
  trap = np.column_stack(
    [np.concatenate([np.zeros(150), rng(s).normal(size=150)]) for s in (0, 1)]
  )
  controls = FitControlsMargin(family_set=eco.impossible, on_failure="fallback")
  with pytest.warns(UserWarning, match="kernel-density margin was") as caught:
    dist = Vinedist.from_data(trap, margins=eco.cls(), margin_controls=controls)
  assert all(isinstance(m, Kde1d) for m in dist.margins)
  assert len(caught) == 2  # one per column, and no more
  assert eco.reason in str(caught[0].message)
  assert np.all(np.isfinite(dist.logpdf(trap)))
