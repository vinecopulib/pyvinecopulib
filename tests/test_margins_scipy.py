"""Tests for `SciPyMargin`, its family search, and `FitControlsMargin`.

Four contracts are pinned here. That a parametric margin counts only the
parameters it actually estimated, since SciPy's legacy fitter always returns
`loc` and `scale` and every miscounted parameter moves AIC by 2. That the
curated candidate set is *ours*: the documented traps -- `vonmises` winning a
blind sweep on clean gamma data, a `weibull_min` fit whose `loc` overshoots the
smallest observation, a Student `t` collapsing onto an atom, and counts ranked
against densities -- must not reach the caller, and every refusal names its
cause rather than being skipped. That the marginal configuration reaches each
variable's own search. And that `pyvinecopulib.margins` imports without SciPy
at all.

What every ecosystem adapter owes -- `fit` versus `select`, the criteria, the
weights and shape refusals -- is in `test_margins_contract.py`.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest

from pyvinecopulib.core import Kde1d, Vinedist
from pyvinecopulib.margins import (
  FitControlsMargin,
  SciPyMargin,
  resolve_margin_controls,
  resolve_margins,
)
from .helpers import run_without

pytest.importorskip("scipy.stats")


@pytest.fixture
def gamma_sample() -> np.ndarray:
  """500 draws from a gamma(2.5, 1.5), the running example of the design."""
  return np.random.default_rng(0).gamma(2.5, 1.5, size=500)


@pytest.fixture
def lognormal_pair() -> np.ndarray:
  """400 draws of two dependent lognormal variables."""
  z = np.random.default_rng(5).normal(size=(400, 2))
  return np.column_stack([np.exp(z[:, 0]), np.exp(z[:, 0] + 0.5 * z[:, 1])])


@pytest.fixture
def collapsed_trap() -> np.ndarray:
  """Two columns, half exact zeros, on which every `t` fit collapses."""

  def column(seed: int) -> np.ndarray:
    zeros = np.zeros(150)
    return np.concatenate([zeros, np.random.default_rng(seed).normal(size=150)])

  return np.column_stack([column(0), column(1)])


def parametric_margins(dist: Vinedist) -> list[SciPyMargin]:
  """The parametric margins of a fitted distribution, in column order."""
  return [m for m in dist.margins if isinstance(m, SciPyMargin)]


# --- SciPyMargin ------------------------------------------------------- #


def test_parametric_margin_counts_only_free_parameters(
  gamma_sample: np.ndarray,
) -> None:
  """A pinned `loc` is not a parameter, so it does not enter AIC."""
  free = SciPyMargin("gamma").fit(gamma_sample)
  pinned = SciPyMargin("gamma", floc=0.0).fit(gamma_sample)
  # SciPy returns three numbers either way; only one of the two estimated three.
  assert len(free.parameters) == len(pinned.parameters) == 3
  assert free.n_parameters == 3.0
  assert pinned.n_parameters == 2.0
  assert pinned.fixed_parameters == {"loc": 0.0}
  # And pinning it is what makes the support the family's own half-line.
  assert pinned.support == (0.0, float("inf"))


def test_parametric_margin_fully_specified_is_already_fitted() -> None:
  """Parameters given at construction need no data, and estimate nothing."""
  m = SciPyMargin("norm", (0.0, 1.0))
  assert m.is_fitted
  assert m.n_parameters == 0.0
  np.testing.assert_allclose(m.cdf(np.array([0.0])), [0.5], atol=1e-12)


def test_parametric_margin_logpdf_beats_log_of_pdf() -> None:
  """`logpdf` comes from SciPy, so the far tail does not underflow to `-inf`."""
  m = SciPyMargin("norm", (0.0, 1.0))
  far = np.array([40.0])
  assert np.isfinite(m.logpdf(far)).all()
  with np.errstate(divide="ignore"):
    assert not np.isfinite(np.log(m.pdf(far))).all()


def test_discrete_fit_is_reproducible(count_sample: np.ndarray) -> None:
  """`scipy.stats.fit` optimizes stochastically; the seed makes it repeatable.

  Without a pinned RNG the same counts give a different parameter vector on
  every call, which makes a search over the count families rank them
  differently run to run.
  """
  fits = [SciPyMargin("nbinom").fit(count_sample).parameters for _ in range(3)]
  for other in fits[1:]:
    np.testing.assert_array_equal(fits[0], other)


def test_parametric_margin_pins_the_discrete_lattice_offset(
  count_sample: np.ndarray,
) -> None:
  """A count family estimates its shape only; `loc` is not identified."""
  m = SciPyMargin("poisson").fit(count_sample)
  assert m.var_type == "d"
  assert m.fixed_parameters == {"loc": 0.0}
  assert m.n_parameters == 1.0
  # `pdf` is the probability mass, and `cdf_left` steps back one lattice point.
  k = np.array([0.0, 1.0, 4.0])
  np.testing.assert_allclose(m.pdf(k), m.cdf(k) - m.cdf_left(k), atol=1e-12)


def test_parametric_margin_fits_when_nothing_is_free() -> None:
  """Pinning every parameter skips the optimizer, which would refuse it."""
  x = np.random.default_rng(2).uniform(1.0, 3.0, size=100)
  m = SciPyMargin("uniform", floc=1.0, fscale=2.0).fit(x)
  assert m.parameters == (1.0, 2.0)
  assert m.n_parameters == 0.0
  assert m.loglik() == pytest.approx(100 * np.log(0.5))


def test_parametric_margin_rejects_an_unknown_family() -> None:
  """A typo in the family name fails at construction, not at fit."""
  with pytest.raises(ValueError, match="unknown scipy.stats family"):
    SciPyMargin("gaussian")


def test_parametric_margin_rejects_an_unknown_fixed_parameter() -> None:
  """A pinning keyword must name a parameter the family has."""
  with pytest.raises(ValueError, match="does not name a parameter"):
    SciPyMargin("gamma", fmu=0.0)


def test_parametric_margin_rejects_a_wrong_length_parameter_vector() -> None:
  """`params` must give every parameter of the family, in SciPy's order."""
  with pytest.raises(ValueError, match="takes 2 parameters"):
    SciPyMargin("norm", (0.0,))


@pytest.mark.parametrize("verb", ["fit", "select"])
def test_parametric_margin_rejects_an_empty_sample(verb: str) -> None:
  """A column of NaNs is not something to fit."""
  margin = SciPyMargin("norm") if verb == "fit" else SciPyMargin()
  with pytest.raises(ValueError, match="no usable observation"):
    getattr(margin, verb)(np.array([np.nan, np.nan]))


def test_parametric_margin_needs_bounds_for_an_uncurated_discrete_family() -> (
  None
):
  """`scipy.stats.fit` needs a bound per free parameter, and cannot guess one."""
  with pytest.raises(ValueError, match="needs search bounds for"):
    SciPyMargin("randint").fit(np.array([1.0, 2.0, 3.0]))


def test_an_unnamed_margin_refuses_what_needs_a_family(
  gamma_sample: np.ndarray,
) -> None:
  """`SciPyMargin()` is a request to choose a family, not a broken one.

  It is what `margins="parametric"` resolves to, so everything that reads a
  family has to say so rather than fail obscurely. There is nothing for `fit`
  to estimate either -- only `select` can supply the family, which is what
  `test_margins_contract.py` pins for both adapters.
  """
  m = SciPyMargin()
  assert "unfitted" in repr(m)
  with pytest.raises(RuntimeError, match="has no family yet"):
    m.fit(gamma_sample)
  with pytest.raises(ValueError, match="params= was given without a family"):
    SciPyMargin(None, (0.0, 1.0))


def test_parametric_margin_round_trips_through_pickle(
  gamma_sample: np.ndarray,
) -> None:
  """The margin pickles as a family name plus a parameter tuple."""
  import pickle

  m = SciPyMargin("gamma", floc=0.0).fit(gamma_sample)
  clone = pickle.loads(pickle.dumps(m))
  assert clone.parameters == m.parameters
  np.testing.assert_allclose(clone.pdf(gamma_sample), m.pdf(gamma_sample))


def test_a_selected_margin_round_trips_through_pickle(
  gamma_sample: np.ndarray,
) -> None:
  """A margin that chose its own family pickles like any other."""
  import pickle

  m = SciPyMargin().select(gamma_sample)
  clone = pickle.loads(pickle.dumps(m))
  assert clone.family_name == m.family_name
  assert clone.parameters == m.parameters
  assert clone.n_parameters == m.n_parameters
  assert clone.nobs == m.nobs


# --- fit versus select ------------------------------------------------------ #


def test_from_data_honors_a_named_family_and_chooses_an_unnamed_one(
  lognormal_pair: np.ndarray,
) -> None:
  """Naming a family *is* the choice, so a two-step fit keeps it.

  The columns are lognormal. A specification that names `norm` asks for a
  normal and gets one -- replacing it would answer the specification with a
  different model. Leaving the family out is what asks for the search.
  """
  named = Vinedist.from_data(lognormal_pair, margins=SciPyMargin("norm"))
  assert [m.family_name for m in parametric_margins(named)] == ["norm", "norm"]
  assert np.all(np.isfinite(named.logpdf(lognormal_pair)))

  chosen = Vinedist.from_data(lognormal_pair, margins="parametric")
  assert [m.family_name for m in parametric_margins(chosen)] == [
    "lognorm",
    "lognorm",
  ]

  # ... and `family_set` is how a caller asks to re-search a named family.
  research = Vinedist.from_data(
    lognormal_pair,
    margins=SciPyMargin("norm"),
    margin_controls=FitControlsMargin(family_set=["lognorm", "expon"]),
  )
  assert [m.family_name for m in parametric_margins(research)] == [
    "lognorm",
    "lognorm",
  ]


# --- the curated candidate set ---------------------------------------------- #


def test_select_picks_the_true_family(gamma_sample: np.ndarray) -> None:
  """Gamma data select a gamma, ahead of ten other admissible candidates."""
  m = SciPyMargin().select(gamma_sample)
  assert m.family_name == "gamma"
  assert m.n_parameters == 2.0
  assert m.fixed_parameters == {"loc": 0.0}
  assert m.nobs == gamma_sample.size


def test_the_curated_set_never_proposes_vonmises(
  gamma_sample: np.ndarray,
) -> None:
  """The `vonmises` trap: it wins a blind sweep, so it is not in the set.

  A sweep over all 110 SciPy continuous families ranks `vonmises` first on this
  very sample -- ahead of the true gamma by 306 AIC units here, and by 701 with
  its `loc` and `scale` estimated the way a sweep would -- because it advertises
  the whole real line as its support while its density integrates to 63.75
  there. The candidate set is curated precisely so that cannot happen.
  """
  would_have_won = SciPyMargin("vonmises").fit(gamma_sample)
  selected = SciPyMargin().select(gamma_sample)
  assert selected.family_name == "gamma"
  assert would_have_won.aic() < selected.aic()
  # It stays reachable on purpose, as an explicit opt-in.
  assert would_have_won.family_name == "vonmises"


def test_parametric_margin_documents_every_exclusion() -> None:
  """The excluded families are rendered into the docstring, not restated in it."""
  doc = SciPyMargin.__doc__ or ""
  assert "{excluded}" not in doc
  for family in ("vonmises", "levy_stable", "chi2", "binom"):
    assert f"``{family}``" in doc


def test_select_adds_the_unit_group_for_data_inside_the_unit_interval() -> None:
  """Data in `(0, 1)` are compatible with the real, positive and unit groups."""
  x = np.random.default_rng(2).beta(2.0, 5.0, size=400)
  m = SciPyMargin().select(x)
  assert m.family_name == "beta"
  # The unit group pins both endpoints, so only the two shapes are estimated.
  assert m.fixed_parameters == {"loc": 0.0, "scale": 1.0}
  assert m.n_parameters == 2.0


def test_declared_bounds_anchor_the_candidates() -> None:
  """Given `(a, b)`, the candidates live on `[a, b]` and are pinned there."""
  x = np.random.default_rng(3).uniform(1.0, 3.0, size=400)
  m = SciPyMargin().select(x, FitControlsMargin(support=(1.0, 3.0)))
  assert m.family_name == "uniform"
  assert m.support == (1.0, 3.0)
  assert m.n_parameters == 0.0


def test_a_half_bounded_support_is_not_a_bounded_one() -> None:
  """`support` selects the bounded families, which need a finite interval.

  A `(0, inf)` declaration would otherwise send a positive variable to
  `uniform` and a rescaled `beta`, neither of which can be pinned to an
  infinite endpoint.
  """
  x = np.random.default_rng(5).gamma(2.0, 1.5, size=400)
  m = SciPyMargin().select(x, FitControlsMargin(support=(0.0, None)))
  assert m.family_name not in ("uniform", "beta")


def test_the_variable_type_is_inferred_from_the_sample() -> None:
  """Non-negative integers read as counts; the same values signed do not.

  A declared `var_type` overrides the inference either way, which is the
  contract's and is pinned in `test_margins_contract.py`; where the inference
  lands with nothing declared is SciPy's own.
  """
  x = np.round(np.random.default_rng(2).normal(0.0, 20.0, size=400))
  assert SciPyMargin().select(x).var_type == "c"  # negative, so not counts
  assert SciPyMargin().select(np.abs(x)).var_type == "d"


def test_a_declared_var_type_supplies_what_the_sample_cannot_show() -> None:
  """A caller's declared type beats the margin's own inference.

  An ordered categorical with negative levels is the sharp case: the
  data-driven test needs non-negative integers, so the sample alone reads the
  column as continuous even though the input declared its levels. Declaring
  `"d"` makes the count families the candidates, and they are then refused by
  name rather than quietly replaced by a density.
  """
  x = np.asarray(
    np.random.default_rng(4).integers(-2, 3, size=400), dtype=float
  )
  assert SciPyMargin().select(x).var_type == "c"

  with pytest.raises(ValueError) as caught:
    SciPyMargin().select(x, FitControlsMargin(var_type="d"))
  message = str(caught.value)
  for family in ("poisson", "nbinom", "geom"):
    assert f"{family}: support" in message


def test_counts_and_densities_are_never_ranked_together(
  count_sample: np.ndarray,
) -> None:
  """A probability mass and a Lebesgue density are not one comparison.

  They are normalized against different reference measures, so an information
  criterion cannot rank them -- and the density wins nearly every time. The
  curated set partitions on the variable type, and an explicit `family_set`
  is partitioned the same way rather than taken at face value.
  """
  selected = SciPyMargin().select(count_sample)
  assert selected.family_name in ("poisson", "nbinom", "geom")
  assert selected.var_type == "d"

  with pytest.warns(UserWarning, match="not comparable"):
    mixed = SciPyMargin().select(
      count_sample,
      FitControlsMargin(family_set=["poisson", "norm", "t"], var_type="d"),
    )
  assert mixed.family_name == "poisson"

  with pytest.raises(ValueError, match="every family in family_set is of the"):
    SciPyMargin().select(
      count_sample, FitControlsMargin(family_set=["norm"], var_type="d")
    )


# --- the admissibility gate ------------------------------------------------- #


def test_the_gate_rejects_a_support_that_excludes_the_data() -> None:
  """The `weibull_min` trap: a free `loc` lands above the smallest observation.

  The fit succeeds and returns a plausible triple, but every log-density is then
  `-inf`. Reaching it needs an unanchored candidate, which no candidate set
  produces any more -- both the curated search and a named `family_set` pin
  `loc` on positive data -- so the gate is exercised directly.
  """
  from pyvinecopulib.margins.scipy import _reject

  x = np.random.default_rng(0).gamma(2.5, 1.5, size=200)
  escaped = SciPyMargin("weibull_min").fit(x)
  assert escaped.support[0] > np.min(x), "the fit did not escape after all"
  reason = _reject(escaped, x)
  assert reason is not None
  assert "does not contain the data range" in reason

  # Anchored the way every candidate set anchors it, the same family is fine.
  anchored = SciPyMargin().select(
    x, FitControlsMargin(family_set=["weibull_min"])
  )
  assert anchored.fixed_parameters == {"loc": 0.0}
  assert _reject(anchored, x) is None


def test_a_collapsed_scale_is_rejected() -> None:
  """A continuous family shrinking onto an atom must not win the criterion.

  Half exact zeros and half N(0, 1): a `t` fitted there comes back with
  `scale` around 1e-20 and a log-likelihood in the thousands, because a density
  concentrating on a repeated value diverges. It beats every honest candidate by
  more than 12,000 AIC units, so the gate -- not the criterion -- is what keeps
  it out.
  """
  y = np.concatenate([np.zeros(250), np.random.default_rng(0).normal(size=250)])
  with pytest.raises(ValueError, match="collapsing onto an atom") as caught:
    SciPyMargin().select(y, FitControlsMargin(family_set=["t"]))
  assert "degenerate parameter scale" in str(caught.value)

  honest = SciPyMargin().select(y)
  collapsed = SciPyMargin("t").fit(y)
  assert honest.family_name != "t"
  assert honest.aic() - collapsed.aic() > 12000
  # And the winner is scored on a finite, honest likelihood.
  assert np.isfinite(honest.loglik())


def test_select_rejects_a_degenerate_parameter() -> None:
  """A Student `t` fitted to normal data runs `df` to ~1e10; that is not a fit.

  It would also *win*: the limit is the normal, so it matches the true density
  and spends one more parameter to say so, which AIC forgives whenever the
  sample is large enough. Rejecting it is what leaves `norm` the winner.
  """
  x = np.random.default_rng(5).normal(0.0, 1.0, size=1000)
  assert SciPyMargin().select(x).family_name == "norm"
  with pytest.raises(ValueError, match="degenerate parameter df") as caught:
    SciPyMargin().select(x, FitControlsMargin(family_set=["t"]))
  assert "outside (0.0, 200.0)" in str(caught.value)


def test_select_reports_geom_as_inadmissible_on_zeros(
  count_sample: np.ndarray,
) -> None:
  """`geom` counts trials, so data containing 0 fall outside its support."""
  with pytest.raises(ValueError, match="does not contain the data range"):
    SciPyMargin().select(count_sample, FitControlsMargin(family_set=["geom"]))


def test_select_reports_a_candidate_that_cannot_be_fitted(
  count_sample: np.ndarray,
) -> None:
  """A candidate whose fitter raises is a refusal with its cause attached."""
  with pytest.raises(ValueError) as caught:
    SciPyMargin().select(
      count_sample, FitControlsMargin(family_set=["randint"])
    )
  assert "randint: ValueError: fitting 'randint' needs search bounds" in str(
    caught.value
  )


def test_the_gate_catches_what_a_family_set_cannot_express() -> None:
  """Two refusals need a candidate no `family_set` can name.

  `family_set` names families, so it cannot pin a parameter to a NaN or a
  density to zero mass on an observation -- yet both are ways a fitter comes
  back with something inadmissible, and neither is implied by the support check
  or by the parameter check alone. They are reached at the gate itself, which
  is the only place they are observable.
  """
  from pyvinecopulib.margins.scipy import _reject

  x = np.random.default_rng(0).normal(size=50)
  nan_pinned = SciPyMargin("norm", floc=float("nan"), fscale=1.0).fit(x)
  assert _reject(nan_pinned, x) == "non-finite parameter loc=nan"

  # A beta pinned to [0, 1] is admissible on data in [0, 1], and its density
  # still vanishes at an exact 0 -- which the support check cannot see.
  b = np.random.default_rng(2).beta(2.0, 5.0, size=400)
  b[0] = 0.0
  zero_mass = SciPyMargin("beta", fa=2.0, fb=5.0, floc=0.0, fscale=1.0).fit(b)
  assert _reject(zero_mass, b) == "non-finite log-likelihood (-inf)"

  assert _reject(SciPyMargin("norm").fit(x), x) is None


def test_candidates_that_would_tie_are_deduplicated(
  count_sample: np.ndarray,
) -> None:
  """The same family with the same pins twice would tie exactly.

  Two such candidates fit identical numbers, so the winner becomes an artifact
  of iteration order. The pins and the search bounds are part of the identity:
  the same family anchored differently is a different candidate, which is what
  `unit` and `bounded` do with `beta`. Deduplication happens where the
  candidate list is built, before anything is fitted, so it is checked there.
  """
  from pyvinecopulib.margins.scipy import _dedupe

  def families(candidates: list) -> list[str]:
    return [c.family_name for c in candidates]

  collapsed = _dedupe(
    [
      SciPyMargin("norm"),
      SciPyMargin("norm"),
      SciPyMargin("logistic"),
    ]
  )
  assert families(collapsed) == ["norm", "logistic"]

  # Same family, different pins: both are kept.
  kept = _dedupe([SciPyMargin("norm"), SciPyMargin("norm", floc=0.0)])
  assert families(kept) == ["norm", "norm"]

  # And distinct optimizer domains for one family are two candidates.
  bounded = _dedupe(
    [
      SciPyMargin("nbinom", param_bounds={"n": (1.0, 2.0), "p": (0.01, 0.99)}),
      SciPyMargin(
        "nbinom", param_bounds={"n": (10.0, 20.0), "p": (0.01, 0.99)}
      ),
    ]
  )
  assert families(bounded) == ["nbinom", "nbinom"]

  # Two already-fitted models with no exposed parameter vector have unreadable
  # identities. Neither may be discarded merely because the family names tie.
  fitted_kdes = [Kde1d().fit(count_sample), Kde1d().fit(count_sample + 3.0)]
  assert len(_dedupe(fitted_kdes)) == 2


# --- criteria --------------------------------------------------------------- #


def test_criteria_take_a_sample_or_read_the_fitted_value(
  gamma_sample: np.ndarray,
) -> None:
  """Given data they score it, and a fixed margin has only that option.

  A margin built from its parameters estimated nothing, so it costs no
  parameters and recorded no sample size: every criterion is then just the
  deviance, and asking for one without data says what is missing.
  """
  fixed = SciPyMargin("norm", (0.0, 1.0))
  deviance = -2.0 * float(np.sum(fixed.logpdf(gamma_sample)))
  assert fixed.n_parameters == 0.0
  assert fixed.aic(gamma_sample) == pytest.approx(deviance)
  assert fixed.bic(gamma_sample) == pytest.approx(deviance)
  assert fixed.aicc(gamma_sample) == pytest.approx(deviance)
  with pytest.raises(RuntimeError, match="only defined after"):
    fixed.bic()


# --- FitControlsMargin ------------------------------------------------------ #


@pytest.mark.parametrize(
  ("kwargs", "match"),
  [
    ({"selection_criterion": "cv"}, "selection_criterion='cv' is not one of"),
    ({"on_failure": "ignore"}, "on_failure='ignore' is not one of"),
    ({"var_type": "zi2"}, "var_type='zi2' is not one of"),
    ({"family_set": []}, "family_set is empty"),
    ({"family_set": [3]}, "must name families as strings"),
    ({"support": (1.0,)}, r"support must be a \(lo, hi\) pair"),
    ({"support": (3.0, 1.0)}, "is not an increasing interval"),
    ({"support": (1.0, 1.0)}, "is not an increasing interval"),
  ],
)
def test_fit_controls_margin_validates_its_arguments(
  kwargs: dict, match: str
) -> None:
  """Every setting is checked at construction, where the caller can see it."""
  with pytest.raises(ValueError, match=match):
    FitControlsMargin(**kwargs)


def test_fit_controls_margin_defaults_and_to_dict() -> None:
  """`to_dict` is what makes it a `ControlsLike`, and the defaults search."""
  assert FitControlsMargin().to_dict() == {
    "family_set": None,
    "selection_criterion": "aic",
    "var_type": None,
    "support": None,
    "on_failure": "raise",
  }
  # A one-sided bound is a legal declaration; it is simply not a bounded one.
  assert FitControlsMargin(support=(0.0, None)).support == (0.0, None)
  assert FitControlsMargin(family_set=("gamma",)).to_dict()["family_set"] == [
    "gamma"
  ]


def test_select_raises_rather_than_silently_falling_back(
  collapsed_trap: np.ndarray,
) -> None:
  """A parametric request answered nonparametrically is a silent downgrade.

  The fallback is still available, but it is opt-in: a warning is easy to lose
  in a loop over fifty columns, so the default names every family that was
  tried, why each failed, and how to ask for the fallback.
  """
  controls = FitControlsMargin(family_set=["t"])
  with pytest.raises(ValueError) as caught:
    SciPyMargin().select(collapsed_trap[:, 0], controls)
  message = str(caught.value)
  assert "no parametric family fits this variable" in message
  assert 'on_failure="fallback"' in message

  # Substituting another kind of margin is a decision about which margin the
  # column gets, so it is made where the margin is chosen, not inside one.
  with pytest.raises(ValueError, match="no parametric family fits"):
    Vinedist.from_data(
      collapsed_trap, margins="parametric", margin_controls=controls
    )


def test_family_set_is_refused_by_a_margin_that_cannot_search(
  lognormal_pair: np.ndarray,
) -> None:
  """A declaration is a default; an instruction to search is not.

  A kernel density has no family to choose, so honoring `family_set` is
  impossible and ignoring it would look like a choice.
  """
  with pytest.raises(TypeError, match="cannot select a family"):
    Vinedist.from_data(
      lognormal_pair,
      margins="kde",
      margin_controls=FitControlsMargin(family_set=["gamma"]),
    )


# --- per-variable resolution ------------------------------------------------ #


def test_resolve_margin_controls_addresses_each_variable() -> None:
  """The `margins=` rules, applied to the marginal configuration."""
  bounded = FitControlsMargin(support=(0.0, None))
  assert resolve_margin_controls(None, 3) == [None, None, None]
  # Broadcast is shared, not copied: controls are read, never written to.
  assert all(c is bounded for c in resolve_margin_controls(bounded, 3))
  assert resolve_margin_controls({"b": bounded}, 3, names=["a", "b", "c"]) == [
    None,
    bounded,
    None,
  ]
  assert resolve_margin_controls({1: bounded}, 3) == [None, bounded, None]
  assert resolve_margin_controls([bounded, None, bounded], 3) == [
    bounded,
    None,
    bounded,
  ]
  with pytest.raises(ValueError, match="has length 1, but there are 3"):
    resolve_margin_controls([bounded], 3)
  with pytest.raises(ValueError, match="names 'z', which is not a variable"):
    resolve_margin_controls({"z": bounded}, 3, names=["a", "b", "c"])


def test_declared_bounds_reach_the_default_margin(
  lognormal_pair: np.ndarray,
) -> None:
  """A declared support bounds the kernel density the library builds.

  The margin does not exist yet, so there is nothing for the declaration to
  override -- and a bound learned after the fit is a bound already padded past.
  The variable left unaddressed keeps the unbounded grid, whose draws go
  negative where the bounded one's cannot.
  """
  dist = Vinedist.from_data(
    lognormal_pair,
    names=["bounded", "free"],
    margin_controls={"bounded": FitControlsMargin(support=(0.0, None))},
  )
  draws = dist.sample(2000, seeds=[3])
  assert draws[:, 0].min() >= 0.0
  assert draws[:, 1].min() < 0.0


# --- wiring ----------------------------------------------------------------- #


def test_resolve_margins_knows_the_parametric_alias() -> None:
  """`margins="parametric"` means "choose a family", one margin per variable."""
  resolved = resolve_margins("parametric", 3)
  assert len(resolved) == 3
  assert all(isinstance(m, SciPyMargin) for m in resolved)
  for margin in resolved:
    with pytest.raises(RuntimeError, match="has no family yet"):
      margin.family_name
  # Copied, not shared: fitting one must not fit the others.
  assert len({id(m) for m in resolved}) == 3


def test_vinedist_from_data_selects_margins(
  lognormal_pair: np.ndarray,
) -> None:
  """The alias reaches `Vinedist.from_data`, which selects one per column."""
  dist = Vinedist.from_data(
    lognormal_pair, margins="parametric", names=["a", "b"]
  )
  chosen = parametric_margins(dist)
  assert len(chosen) == dist.dim
  assert all(m.is_fitted for m in chosen)
  assert np.all(np.isfinite(dist.logpdf(lognormal_pair)))


def test_margins_imports_without_scipy() -> None:
  """`pyvinecopulib.margins` must not need the `scipy` extra to import.

  The default margin is `pyvinecopulib.core.Kde1d`, so `core` has to stay
  SciPy-free as well, and so does `FitControlsMargin`.
  """
  code = (
    "import sys; sys.modules['scipy'] = None; "
    "import pyvinecopulib.core as c; import pyvinecopulib.margins as m; "
    "c.Kde1d(); m.FitControlsMargin(); "
    "sys.exit(0 if 'scipy.stats' not in sys.modules else 1)"
  )
  result = subprocess.run(  # noqa: S603
    [sys.executable, "-c", code], capture_output=True, text=True
  )
  assert result.returncode == 0, result.stderr


def test_parametric_margin_names_the_extra_without_scipy() -> None:
  """Without SciPy, constructing one says which extra provides it."""
  run_without(
    "scipy",
    "from pyvinecopulib.margins import SciPyMargin\n"
    "try:\n"
    "  SciPyMargin('norm')\n"
    "except ImportError as e:\n"
    "  sys.exit(0 if 'pyvinecopulib[scipy]' in str(e) else 2)\n"
    "sys.exit(3)\n",
  )


@pytest.mark.parametrize(
  ("data_kind", "expect_pinned"),
  [("real", False), ("positive", True), ("unit", True), ("count", True)],
)
def test_a_selected_margin_survives_a_json_round_trip(
  data_kind: str, expect_pinned: bool
) -> None:
  """A pinned parameter has two spellings, and the payload used the other one.

  `_fixed` holds normalized names (``loc``) while the constructor takes
  SciPy's pin spelling (``floc``), so writing one and reading it as the other
  raised for every margin with a pin -- which is every curated candidate
  outside the ``real`` group, and therefore the `margins="parametric"`
  default on positive, unit, bounded or count data. Only unpinned ``norm``
  round-tripped, which is what the one existing test used.
  """
  from pyvinecopulib.core import margin_from_json, margin_to_json

  rng = np.random.default_rng(0)
  data = {
    "real": rng.normal(size=300),
    "positive": np.abs(rng.normal(size=300)) + 0.1,
    "unit": rng.uniform(0.05, 0.95, size=300),
    "count": rng.poisson(3.0, size=300).astype(float),
  }[data_kind]

  margin = SciPyMargin().select(data)
  assert bool(margin.fixed_parameters) is expect_pinned, margin.family_name

  back = margin_from_json(margin_to_json(margin))
  assert back.family_name == margin.family_name
  assert back.fixed_parameters == margin.fixed_parameters
  assert back.parameters == margin.parameters
  np.testing.assert_allclose(back.pdf(data[:20]), margin.pdf(data[:20]))
  for name in ("loglik", "aic", "bic", "aicc"):
    assert getattr(back, name)() == pytest.approx(getattr(margin, name)())
