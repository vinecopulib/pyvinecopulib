"""Tests for `ParametricMargin` and `MarginSelector`.

Four contracts are pinned here. That a parametric margin counts only the
parameters it actually estimated, since SciPy's legacy fitter always returns
`loc` and `scale` and every miscounted parameter moves AIC by 2. That the
curated candidate set is *ours*: the two documented traps -- `vonmises` winning
a blind sweep on clean gamma data, and a `weibull_min` fit whose `loc` overshoots
the smallest observation -- must not reach the caller. That every rejection is
recorded rather than skipped, and that the fallback when everything fails is a
kernel-density margin, never `norm`. And that `pyvinecopulib.margins` imports
without SciPy at all.
"""

from __future__ import annotations

import subprocess
import sys
from typing import Any

import numpy as np
import pytest

from pyvinecopulib.core import MarginLike
from pyvinecopulib.core import Kde1d
from pyvinecopulib.margins import (
  MarginSelector,
  ParametricMargin,
  resolve_margins,
)

pytest.importorskip("scipy.stats")


@pytest.fixture
def gamma_sample() -> np.ndarray:
  """500 draws from a gamma(2.5, 1.5), the running example of the design."""
  return np.random.default_rng(0).gamma(2.5, 1.5, size=500)


@pytest.fixture
def count_sample() -> np.ndarray:
  """400 Poisson(4) counts, including zeros."""
  return np.random.default_rng(1).poisson(4.0, size=400).astype(float)


# --- ParametricMargin ------------------------------------------------------- #


def test_parametric_margin_conforms(gamma_sample: np.ndarray) -> None:
  """A fitted parametric margin satisfies the margin contract."""
  m = ParametricMargin("gamma", floc=0.0).fit(gamma_sample)
  assert isinstance(m, MarginLike)
  assert m.family_name == "gamma"
  assert m.var_type == "c"
  assert m.support == (0.0, float("inf"))
  np.testing.assert_allclose(
    m.cdf(m.icdf(np.array([0.1, 0.5, 0.9]))),
    [0.1, 0.5, 0.9],
    atol=1e-10,
  )


def test_parametric_margin_counts_only_free_parameters(
  gamma_sample: np.ndarray,
) -> None:
  """A pinned `loc` is not a parameter, so it does not enter AIC."""
  free = ParametricMargin("gamma").fit(gamma_sample)
  pinned = ParametricMargin("gamma", floc=0.0).fit(gamma_sample)
  # SciPy returns three numbers either way; only one of the two estimated three.
  assert len(free.parameters) == len(pinned.parameters) == 3
  assert free.n_parameters == 3.0
  assert pinned.n_parameters == 2.0
  assert pinned.fixed_parameters == {"loc": 0.0}


def test_parametric_margin_fully_specified_is_already_fitted() -> None:
  """Parameters given at construction need no data, and estimate nothing."""
  m = ParametricMargin("norm", (0.0, 1.0))
  assert m.is_fitted
  assert m.n_parameters == 0.0
  np.testing.assert_allclose(m.cdf(np.array([0.0])), [0.5], atol=1e-12)


def test_parametric_margin_logpdf_beats_log_of_pdf() -> None:
  """`logpdf` comes from SciPy, so the far tail does not underflow to `-inf`."""
  m = ParametricMargin("norm", (0.0, 1.0))
  far = np.array([40.0])
  assert np.isfinite(m.logpdf(far)).all()
  with np.errstate(divide="ignore"):
    assert not np.isfinite(np.log(m.pdf(far))).all()


def test_discrete_fit_is_reproducible(count_sample: np.ndarray) -> None:
  """`scipy.stats.fit` optimizes stochastically; the seed makes it repeatable.

  Without a pinned RNG the same counts give a different parameter vector on
  every call, which makes a `MarginSelector` count-group `report_` differ run
  to run.
  """
  fits = [
    ParametricMargin("nbinom").fit(count_sample).parameters for _ in range(3)
  ]
  for other in fits[1:]:
    np.testing.assert_array_equal(fits[0], other)


def test_parametric_margin_pins_the_discrete_lattice_offset(
  count_sample: np.ndarray,
) -> None:
  """A count family estimates its shape only; `loc` is not identified."""
  m = ParametricMargin("poisson").fit(count_sample)
  assert m.var_type == "d"
  assert m.fixed_parameters == {"loc": 0.0}
  assert m.n_parameters == 1.0
  # `pdf` is the probability mass, and `cdf_left` steps back one lattice point.
  k = np.array([0.0, 1.0, 4.0])
  np.testing.assert_allclose(m.pdf(k), m.cdf(k) - m.cdf_left(k), atol=1e-12)


def test_parametric_margin_fits_when_nothing_is_free() -> None:
  """Pinning every parameter skips the optimizer, which would refuse it."""
  x = np.random.default_rng(2).uniform(1.0, 3.0, size=100)
  m = ParametricMargin("uniform", floc=1.0, fscale=2.0).fit(x)
  assert m.parameters == (1.0, 2.0)
  assert m.n_parameters == 0.0
  assert m.loglik() == pytest.approx(100 * np.log(0.5))


def test_parametric_margin_rejects_weights(gamma_sample: np.ndarray) -> None:
  """SciPy cannot fit with weights, so asking must raise, not ignore them."""
  with pytest.raises(TypeError, match="cannot use observation weights"):
    ParametricMargin("gamma", floc=0.0).fit(
      gamma_sample, weights=np.ones_like(gamma_sample)
    )


def test_parametric_margin_rejects_an_unknown_family() -> None:
  """A typo in the family name fails at construction, not at fit."""
  with pytest.raises(ValueError, match="unknown scipy.stats family"):
    ParametricMargin("gaussian")


def test_parametric_margin_rejects_an_unknown_fixed_parameter() -> None:
  """A pinning keyword must name a parameter the family has."""
  with pytest.raises(ValueError, match="does not name a parameter"):
    ParametricMargin("gamma", fmu=0.0)


def test_parametric_margin_rejects_a_wrong_length_parameter_vector() -> None:
  """`params` must give every parameter of the family, in SciPy's order."""
  with pytest.raises(ValueError, match="takes 2 parameters"):
    ParametricMargin("norm", (0.0,))


def test_parametric_margin_rejects_an_empty_sample() -> None:
  """A column of NaNs is not something to fit."""
  with pytest.raises(ValueError, match="no usable observation"):
    ParametricMargin("norm").fit(np.array([np.nan, np.nan]))


def test_parametric_margin_needs_bounds_for_an_uncurated_discrete_family() -> (
  None
):
  """`scipy.stats.fit` needs a bound per free parameter, and cannot guess one."""
  with pytest.raises(ValueError, match="needs search bounds for"):
    ParametricMargin("randint").fit(np.array([1.0, 2.0, 3.0]))


def test_parametric_margin_raises_before_fit() -> None:
  """An unfitted margin has neither parameters nor a log-likelihood."""
  m = ParametricMargin("gamma", floc=0.0)
  assert not m.is_fitted
  assert m.support == (float("-inf"), float("inf"))
  with pytest.raises(RuntimeError, match="is not fitted"):
    m.parameters
  with pytest.raises(RuntimeError, match="only defined after"):
    m.loglik()


def test_parametric_margin_round_trips_through_pickle(
  gamma_sample: np.ndarray,
) -> None:
  """The margin pickles as a family name plus a parameter tuple."""
  import pickle

  m = ParametricMargin("gamma", floc=0.0).fit(gamma_sample)
  clone = pickle.loads(pickle.dumps(m))
  assert clone.parameters == m.parameters
  np.testing.assert_allclose(clone.pdf(gamma_sample), m.pdf(gamma_sample))


def test_parametric_margin_simulates_reproducibly() -> None:
  """Same seeds, same draws."""
  m = ParametricMargin("norm", (0.0, 1.0))
  np.testing.assert_array_equal(
    m.sample(16, seeds=[5]), m.sample(16, seeds=[5])
  )


# --- MarginSelector --------------------------------------------------------- #


def test_selector_picks_the_true_family(gamma_sample: np.ndarray) -> None:
  """Gamma data select a gamma, ahead of ten other admissible candidates."""
  sel = MarginSelector().fit(gamma_sample)
  assert sel.selected_.family_name == "gamma"
  assert sel.family_name == "gamma"
  assert sel.n_parameters == 2.0
  assert len(sel.report_) == 11


def test_selector_forwards_evaluation_to_the_winner(
  gamma_sample: np.ndarray,
) -> None:
  """A selector is a margin: it evaluates as whatever it selected."""
  sel = MarginSelector().fit(gamma_sample)
  grid = np.array([0.5, 1.0, 3.0])
  for name in ("pdf", "cdf", "logpdf", "cdf_left"):
    np.testing.assert_array_equal(
      getattr(sel, name)(grid), getattr(sel.selected_, name)(grid)
    )
  np.testing.assert_array_equal(
    sel.icdf(np.array([0.3])), sel.selected_.icdf(np.array([0.3]))
  )
  assert sel.support == sel.selected_.support
  assert sel.var_type == sel.selected_.var_type
  assert sel.loglik() == sel.selected_.loglik()


def test_loglik_evaluates_or_reports_the_fitted_value(
  gamma_sample: np.ndarray,
) -> None:
  """One method, two jobs, as on `Bicop`: data evaluates, no data reports.

  A property would shadow the method and make the polymorphic
  `margin.loglik(sample)` call raise `TypeError`, which is why the fitted value
  is reached through the same name rather than a second one.
  """
  for margin in (
    ParametricMargin("norm").fit(gamma_sample),
    MarginSelector(candidates=["norm", "gamma"]).fit(gamma_sample),
  ):
    total = margin.loglik(gamma_sample)
    assert np.ndim(total) == 0
    np.testing.assert_allclose(total, np.sum(margin.logpdf(gamma_sample)))
    assert isinstance(margin.loglik(), float)


def test_selector_never_returns_vonmises(gamma_sample: np.ndarray) -> None:
  """The `vonmises` trap: it wins a blind sweep, so it is not in the set.

  A sweep over all 110 SciPy continuous families ranks `vonmises` first on this
  very sample, ahead of the true gamma by 701 AIC units, because it advertises
  the whole real line as its support while its density integrates to 63.75
  there. The candidate set is curated precisely so that cannot happen.
  """
  sel = MarginSelector().fit(gamma_sample)
  tried = {row["family"] for row in sel.report_}
  assert "vonmises" not in tried
  assert "wrapcauchy" not in tried
  assert "levy_stable" not in tried
  # It stays reachable on purpose, as an explicit opt-in.
  assert ParametricMargin("vonmises").family_name == "vonmises"


def test_selector_documents_every_exclusion() -> None:
  """The excluded families are rendered into the docstring, not restated in it."""
  doc = MarginSelector.__doc__ or ""
  assert "{excluded}" not in doc
  for family in ("vonmises", "levy_stable", "chi2", "binom"):
    assert f"``{family}``" in doc


def test_selector_rejects_a_support_that_excludes_the_data() -> None:
  """The `weibull_min` trap: a free `loc` lands above the smallest observation.

  The fit succeeds and returns a plausible triple, but every log-density is then
  `-inf`. The rejection must be recorded, not silently skipped.
  """
  x = np.random.default_rng(0).gamma(2.5, 1.5, size=200)
  with pytest.warns(UserWarning, match="no parametric family was admissible"):
    sel = MarginSelector(candidates=["weibull_min"]).fit(x)
  row = sel.report_[0]
  assert row["family"] == "weibull_min"
  assert "does not contain the data range" in row["status"]
  assert row["loglik"] == float("-inf")
  # `loc` really did overshoot, which is what the gate saw.
  assert row["support"][0] > x.min()


def test_selector_falls_back_to_a_kernel_density_never_to_norm() -> None:
  """When nothing is admissible the fallback is nonparametric, and warns once."""
  x = np.random.default_rng(0).gamma(2.5, 1.5, size=200)
  with pytest.warns(UserWarning) as caught:
    sel = MarginSelector(candidates=["weibull_min"]).fit(x)
  assert len(caught) == 1
  assert isinstance(sel.selected_, Kde1d)
  assert sel.report_[-1]["status"] == "fallback"
  assert sel.report_[-1]["selected"] is True


def test_selector_rejects_a_degenerate_parameter() -> None:
  """A Student `t` fitted to normal data runs `df` to ~1e10; that is not a fit.

  It would also *win*: the limit is the normal, so it matches the true density
  and spends one more parameter to say so, which AIC forgives whenever the
  sample is large enough. Rejecting it is what leaves `norm` the winner.
  """
  x = np.random.default_rng(5).normal(0.0, 1.0, size=1000)
  sel = MarginSelector().fit(x)
  row = next(r for r in sel.report_ if r["family"] == "t")
  assert "degenerate parameter df" in row["status"]
  # The log-likelihood it reached is still reported, so the row explains itself.
  assert np.isfinite(row["loglik"])
  assert sel.selected_.family_name == "norm"


def test_selector_adds_the_unit_group_for_data_inside_the_unit_interval() -> (
  None
):
  """Data in `(0, 1)` are compatible with the real, positive and unit groups."""
  x = np.random.default_rng(2).beta(2.0, 5.0, size=400)
  sel = MarginSelector().fit(x)
  tried = {row["family"] for row in sel.report_}
  assert {"norm", "gamma", "beta"} <= tried
  assert sel.selected_.family_name == "beta"


def test_selector_rejects_a_non_finite_log_likelihood() -> None:
  """A fit whose support contains the data can still put zero mass on some of it.

  A beta pinned to `[0, 1]` is admissible on data in `[0, 1]`, and its density
  still vanishes at an exact 0 -- which the support check cannot see and the
  log-likelihood check can.
  """
  x = np.random.default_rng(2).beta(2.0, 5.0, size=400)
  x[0] = 0.0
  candidate = ParametricMargin("beta", fa=2.0, fb=5.0, floc=0.0, fscale=1.0)
  with pytest.warns(UserWarning, match="no parametric family was admissible"):
    sel = MarginSelector(candidates=[candidate]).fit(x)
  assert sel.report_[0]["status"] == "non-finite log-likelihood (-inf)"


def test_selector_rejects_a_non_finite_parameter() -> None:
  """A parameter that is not a number is rejected before anything is evaluated."""
  candidate = ParametricMargin("norm", floc=float("nan"), fscale=1.0)
  x = np.random.default_rng(0).normal(size=50)
  with pytest.warns(UserWarning, match="no parametric family was admissible"):
    sel = MarginSelector(candidates=[candidate]).fit(x)
  assert sel.report_[0]["status"] == "non-finite parameter loc=nan"


def test_selector_simulates_through_the_winner(
  gamma_sample: np.ndarray,
) -> None:
  """Sampling a selector samples whatever it selected."""
  sel = MarginSelector().fit(gamma_sample)
  draws = sel.sample(32, seeds=[7])
  assert draws.shape == (32,)
  np.testing.assert_array_equal(draws, sel.selected_.sample(32, seeds=[7]))


def test_selector_keeps_counts_out_of_the_continuous_comparison(
  count_sample: np.ndarray,
) -> None:
  """Count data are fitted with count families only."""
  sel = MarginSelector().fit(count_sample)
  assert {row["family"] for row in sel.report_} == {
    "poisson",
    "nbinom",
    "geom",
  }
  assert sel.selected_.family_name == "poisson"
  assert sel.var_type == "d"


def test_selector_reports_geom_as_inadmissible_on_zeros(
  count_sample: np.ndarray,
) -> None:
  """`geom` counts trials, so data containing 0 fall outside its support."""
  sel = MarginSelector().fit(count_sample)
  row = next(r for r in sel.report_ if r["family"] == "geom")
  assert "does not contain the data range" in row["status"]


def test_selector_var_type_overrides_the_inference() -> None:
  """Integer-valued data can still be declared continuous."""
  x = np.round(np.random.default_rng(2).normal(0.0, 20.0, size=400))
  assert MarginSelector().fit(x).var_type == "c"  # negative, so not counts
  forced = MarginSelector(var_type="c").fit(np.abs(x))
  assert forced.var_type == "c"
  assert MarginSelector(var_type="d").fit(np.abs(x)).var_type == "d"


def test_selector_bounded_support_anchors_the_candidates() -> None:
  """Given `(a, b)`, the candidates live on `[a, b]` and are pinned there."""
  x = np.random.default_rng(3).uniform(1.0, 3.0, size=400)
  sel = MarginSelector(bounds=(1.0, 3.0)).fit(x)
  assert {row["family"] for row in sel.report_} == {"uniform", "beta"}
  assert sel.support == (1.0, 3.0)
  assert sel.selected_.family_name == "uniform"


def test_declare_supplies_what_the_sample_cannot_show() -> None:
  """A caller's resolved type and support beat the selector's own inference.

  An ordered categorical with negative levels is the sharp case: the data-driven
  test needs non-negative integers, so the sample alone reads the column as
  continuous even though the input declared its levels.
  """
  x = np.asarray(
    np.random.default_rng(4).integers(-2, 3, size=400), dtype=float
  )
  assert MarginSelector().fit(x).var_type == "c"

  declared = MarginSelector().declare(var_type="d", support=(-2.0, 2.0))
  assert declared.bounds == (-2.0, 2.0)
  with pytest.warns(UserWarning, match="no parametric family"):
    assert declared.fit(x).var_type == "d"


def test_declare_leaves_a_pinned_argument_alone() -> None:
  """A constructor argument is the caller's instruction; a schema is a default."""
  sel = MarginSelector(var_type="c", bounds=(-5.0, 5.0)).declare(
    var_type="d", support=(-2.0, 2.0)
  )
  assert sel.bounds == (-5.0, 5.0)
  assert sel.fit(np.arange(-3.0, 4.0)).var_type == "c"


def test_declare_ignores_a_half_bounded_support() -> None:
  """`bounds` selects the bounded families, which need a finite interval.

  A `(0, inf)` support would otherwise send a positive variable to `uniform` and
  a rescaled `beta`, neither of which can be pinned to an infinite endpoint.
  """
  sel = MarginSelector().declare(support=(0.0, float("inf")))
  assert sel.bounds is None
  x = np.random.default_rng(5).gamma(2.0, 1.5, size=400)
  assert {row["family"] for row in sel.fit(x).report_} != {"uniform", "beta"}


@pytest.mark.parametrize("criterion", ["aic", "bic", "aicc"])
def test_selector_criteria_are_all_reported(
  criterion: str, gamma_sample: np.ndarray
) -> None:
  """Every criterion is in every row, and the chosen one drives the winner."""
  sel = MarginSelector(criterion=criterion).fit(gamma_sample)
  admissible = [r for r in sel.report_ if r["status"] in ("ok", "selected")]
  winner = next(r for r in sel.report_ if r["selected"])
  assert winner[criterion] == min(r[criterion] for r in admissible)
  for row in sel.report_:
    assert {"aic", "bic", "aicc"} <= set(row)


def test_selector_bic_penalizes_more_than_aic() -> None:
  """BIC's `log(n)` penalty is the only difference from AIC's 2."""
  from pyvinecopulib.margins.selection import _criteria

  values = _criteria(loglik=-100.0, k=3.0, n=1000)
  assert values["bic"] > values["aic"]
  assert values["aicc"] > values["aic"]
  undefined = _criteria(loglik=-1.0, k=10.0, n=5)
  assert undefined["aicc"] == float("inf")
  assert _criteria(loglik=float("-inf"), k=1.0, n=10)["aic"] == float("inf")


def test_selector_report_carries_the_column_name(
  gamma_sample: np.ndarray,
) -> None:
  """Rows from several variables concatenate into one tidy table."""
  sel = MarginSelector(name="claims").fit(gamma_sample)
  assert all(row["column"] == "claims" for row in sel.report_)


def test_selector_accepts_ready_made_candidates(
  gamma_sample: np.ndarray,
) -> None:
  """Candidates may be margins, not only family names."""
  sel = MarginSelector(candidates=[Kde1d(), "gamma"]).fit(gamma_sample)
  assert [row["family"] for row in sel.report_] == ["kde1d", "gamma"]


def test_selector_records_a_failed_fit_rather_than_raising() -> None:
  """A candidate that cannot be fitted at all is a row, not an exception."""

  class Broken(Kde1d):
    """A margin whose fit always fails."""

    def fit(self, x: Any, weights: Any = None) -> "Broken":
      """Fail.

      Parameters
      ----------
      x : array, shape (n,), dtype float
          Ignored.
      weights : array, shape (n,), or None, optional
          Ignored.

      Returns
      -------
      Broken
          Never returns.

      Raises
      ------
      RuntimeError
          Always.
      """
      raise RuntimeError("nope")

  x = np.random.default_rng(4).normal(size=50)
  with pytest.warns(UserWarning, match="no parametric family was admissible"):
    sel = MarginSelector(candidates=[Broken()]).fit(x)
  assert sel.report_[0]["status"] == "RuntimeError: nope"
  assert np.isnan(sel.report_[0]["loglik"])


def test_selector_rejects_weights(gamma_sample: np.ndarray) -> None:
  """Weights cannot reach the parametric candidates, so they are refused."""
  with pytest.raises(TypeError, match="cannot use observation weights"):
    MarginSelector().fit(gamma_sample, weights=np.ones_like(gamma_sample))


def test_selector_validates_its_arguments() -> None:
  """Bad `criterion` / `var_type` / `candidates` fail loudly."""
  with pytest.raises(ValueError, match="unknown criterion"):
    MarginSelector(criterion="cv")
  with pytest.raises(ValueError, match="unknown var_type"):
    MarginSelector(var_type="zi")
  with pytest.raises(ValueError, match="unknown candidates"):
    MarginSelector(candidates="everything").fit(np.array([1.0, 2.0, 3.0]))
  with pytest.raises(ValueError, match="no usable observation"):
    MarginSelector().fit(np.array([np.nan]))


def test_selector_raises_before_fit() -> None:
  """An unfitted selector has no winner and an empty report."""
  sel = MarginSelector()
  assert not sel.is_fitted
  assert sel.report_ == []
  assert "unfitted" in repr(sel)
  with pytest.raises(RuntimeError, match="is not fitted"):
    sel.selected_


def test_selector_round_trips_through_pickle(
  gamma_sample: np.ndarray,
) -> None:
  """A fitted selector pickles together with its winner and its report."""
  import pickle

  sel = MarginSelector().fit(gamma_sample)
  clone = pickle.loads(pickle.dumps(sel))
  assert clone.family_name == sel.family_name
  assert len(clone.report_) == len(sel.report_)
  np.testing.assert_allclose(clone.pdf(gamma_sample), sel.pdf(gamma_sample))


# --- wiring ----------------------------------------------------------------- #


def test_resolve_margins_knows_the_parametric_alias() -> None:
  """`margins="parametric"` means "select a family", one selector per variable."""
  resolved = resolve_margins("parametric", 3)
  assert len(resolved) == 3
  assert all(isinstance(m, MarginSelector) for m in resolved)
  # Copied, not shared: fitting one must not fit the others.
  assert len({id(m) for m in resolved}) == 3


def test_vinedist_from_data_selects_margins() -> None:
  """The alias reaches `Vinedist.from_data`, which fits one selector per column."""
  import pyvinecopulib as pv
  from pyvinecopulib.core import Vinedist

  rng = np.random.default_rng(5)
  z = rng.normal(size=(300, 2))
  x = np.column_stack([np.exp(z[:, 0]), z[:, 0] + 0.5 * z[:, 1]])
  dist = Vinedist.from_data(x, margins="parametric", names=["a", "b"])
  assert all(isinstance(m, MarginSelector) for m in dist.margins)
  assert np.all(np.isfinite(dist.logpdf(x)))
  del pv


def test_margins_imports_without_scipy() -> None:
  """`pyvinecopulib.margins` must not need the `scipy` extra to import.

  The default margin is `pyvinecopulib.core.Kde1d`, so `core` has to stay
  SciPy-free as well.
  """
  code = (
    "import sys; sys.modules['scipy'] = None; "
    "import pyvinecopulib.core as c; import pyvinecopulib.margins as m; "
    "c.Kde1d(); m.MarginSelector(); "
    "sys.exit(0 if 'scipy.stats' not in sys.modules else 1)"
  )
  result = subprocess.run(  # noqa: S603
    [sys.executable, "-c", code], capture_output=True, text=True
  )
  assert result.returncode == 0, result.stderr


def test_parametric_margin_names_the_extra_without_scipy() -> None:
  """Without SciPy, constructing one says which extra provides it."""
  code = (
    "import sys, builtins\n"
    "real = builtins.__import__\n"
    "def blocked(name, *a, **k):\n"
    "  if name.split('.')[0] == 'scipy':\n"
    "    raise ImportError('no scipy')\n"
    "  return real(name, *a, **k)\n"
    "builtins.__import__ = blocked\n"
    "from pyvinecopulib.margins import ParametricMargin\n"
    "try:\n"
    "  ParametricMargin('norm')\n"
    "except ImportError as e:\n"
    "  sys.exit(0 if 'pyvinecopulib[scipy]' in str(e) else 2)\n"
    "sys.exit(3)\n"
  )
  result = subprocess.run(  # noqa: S603
    [sys.executable, "-c", code], capture_output=True, text=True
  )
  assert result.returncode == 0, result.stdout + result.stderr


def test_a_collapsed_scale_is_rejected() -> None:
  """A continuous family shrinking onto an atom must not win the criterion.

  Half exact zeros and half N(0, 1): a `t` fitted there comes back with
  `scale` around 1e-20 and a log-likelihood in the thousands, because a density
  concentrating on a repeated value diverges. It beat every honest candidate by
  more than 12,000 AIC units.
  """
  y = np.concatenate([np.zeros(250), np.random.default_rng(0).normal(size=250)])
  sel = MarginSelector().fit(y)
  rows = {row["family"]: row for row in sel.report_}
  assert "below" in rows["t"]["status"]
  assert sel.selected_.family_name != "t"
  # And the winner is scored on a finite, honest likelihood.
  assert np.isfinite(sel.loglik()) and sel.loglik() < rows["t"]["loglik"]


def test_counts_and_densities_are_never_ranked_together() -> None:
  """An explicit candidate list is partitioned by kind, like `"auto"` is.

  A probability mass and a Lebesgue density are not comparable on one
  information criterion, and the density wins nearly every time.
  """
  y = np.random.default_rng(1).poisson(0.7, 500).astype(float)
  with pytest.warns(UserWarning, match="not comparable"):
    sel = MarginSelector(candidates=["poisson", "norm", "t"], var_type="d").fit(
      y
    )
  assert sel.selected_.family_name == "poisson"
  assert [row["family"] for row in sel.report_] == ["poisson"]

  with pytest.raises(ValueError, match="every candidate is of the other kind"):
    MarginSelector(candidates=["norm"], var_type="d").fit(y)
