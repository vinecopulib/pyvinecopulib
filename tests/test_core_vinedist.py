"""Tests for `pyvinecopulib.core.Vinedist`.

The claims worth pinning are the identities. That `logpdf` really is
`log c(F(x)) + sum_j log pdf_j(x_j)`, checked against a hand computation rather
than a tolerance. And that the `(n, d + k)` layout a discrete margin
needs is assembled for the user, since assembling it by hand is what
`examples/04_discrete_variables.ipynb` currently has to do.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import (
  FitControlsKde1d,
  Kde1d,
  MarginBase,
  Vinedist,
  VinedistBase,
  VinedistLike,
)
from pyvinecopulib.margins import FitControlsMargin, SciPyMargin

from .conftest import GaussianBicop, HostedVinecop
from .helpers import FlatMargin, ShiftedNormalMargin, widen

# The discrete cascade owns these; the end-to-end test at the bottom reuses them
# rather than duplicating the type patterns and the parity bound.
from .test_core_discrete_cascade import _D, _assert_parity, _both

stats = pytest.importorskip("scipy.stats")


class ParametricVinedist(Vinedist):
  """A `Vinedist` whose margins are parametric: the `margin_class` route."""

  margin_class = SciPyMargin


def fit_one(margin: Any, column: np.ndarray, x: Any = None) -> Any:
  """Fit one margin to its column, leaving an already-fitted one alone."""
  from pyvinecopulib.margins import as_margin

  m: Any = as_margin(margin)
  if getattr(m, "is_fitted", True):
    return m
  verb = getattr(m, "select", None) or m.fit
  kwargs = (
    {"x": x}
    if x is not None and getattr(m, "supports_covariates", False)
    else {}
  )
  verb(column, **kwargs)
  return m


def from_margins(
  y: Any, margins: list[Any], x: Any = None, **kwargs: Any
) -> Vinedist:
  """Fit a copula onto margins the caller supplies: the constructor route.

  The two-step estimator written out, for the tests that need the copula half
  configured directly. A caller who only wants their own margins fitted writes
  `Vinedist(copula, margins).select(y)` instead, which leaves a fixed margin
  alone and re-estimates the rest.
  """
  data = np.asarray(y, dtype=float)
  fitted = [fit_one(m, data[:, j], x) for j, m in enumerate(margins)]
  vinecop = pv.Vinecop.from_data(
    Vinedist.copula_data(fitted, data, x=x),
    var_types=Vinedist.copula_var_types(fitted),
    **kwargs,
  )
  return Vinedist(vinecop, fitted)


@pytest.fixture
def data() -> np.ndarray:
  """Three dependent columns: real-line, positive, and a count."""
  rng = np.random.default_rng(0)
  cov = [[1.0, 0.6, 0.3], [0.6, 1.0, 0.4], [0.3, 0.4, 1.0]]
  z = rng.multivariate_normal([0.0, 0.0, 0.0], cov, size=400)
  return np.column_stack(
    [z[:, 0], np.exp(z[:, 1]), rng.poisson(np.exp(0.5 * z[:, 2])) * 1.0]
  )


@pytest.fixture
def continuous(data: np.ndarray) -> np.ndarray:
  """The two continuous columns."""
  return data[:, :2]


# --- the Sklar identity ----------------------------------------------------- #


def _sklar_logpdf(dist: Any, y: np.ndarray) -> np.ndarray:
  """`log c(layout) + sum_j log pdf_j(y_j)`, by hand.

  The identity holds verbatim for atoms too, because `pdf` is the density with
  respect to each margin's *own* reference measure -- a mass where there is an
  atom. Pinning it rather than only finiteness is what makes a change to the
  discrete branch visible: `logpdf` stays finite through a wrong quotient.
  """
  manual = np.log(np.asarray(dist.vinecop.pdf(dist.copula_layout(y))))
  for j, m in enumerate(dist.margins):
    manual = manual + np.log(np.asarray(m.pdf(y[:, j])))
  return manual


def test_logpdf_is_the_sklar_factorization(continuous: np.ndarray) -> None:
  """`logpdf` equals the copula term plus the marginal log-densities."""
  dist = pv.Vinedist.from_data(continuous)
  np.testing.assert_allclose(
    dist.logpdf(continuous), _sklar_logpdf(dist, continuous), atol=0.0
  )


def test_pdf_is_the_exponential_of_logpdf(continuous: np.ndarray) -> None:
  """The two are consistent, and `logpdf` is the primitive."""
  dist = pv.Vinedist.from_data(continuous)
  np.testing.assert_allclose(
    dist.pdf(continuous), np.exp(dist.logpdf(continuous)), rtol=1e-12
  )


def test_loglik_sums_logpdf(continuous: np.ndarray) -> None:
  """`loglik` is the total, as a 0-d array."""
  dist = pv.Vinedist.from_data(continuous)
  total = dist.loglik(continuous)
  assert np.ndim(total) == 0
  np.testing.assert_allclose(total, dist.logpdf(continuous).sum(), rtol=1e-12)


def test_logpdf_preserves_an_extreme_tail_copula_density() -> None:
  """The copula term is logged without replacing valid tiny densities."""
  pair = pv.Bicop.from_family(
    pv.families.gaussian, parameters=np.array([[0.99]])
  )
  structure = pv.RVineStructure.from_order([1, 2])
  copula = HostedVinecop([[pair]], structure)
  dist = Vinedist(copula, [stats.uniform(), stats.uniform()])
  y = np.array([[1e-6, 1 - 1e-6], [1e-5, 1 - 1e-5]])
  copula_density = np.asarray(pair.pdf(y))
  assert np.all(copula_density < 1e-300)
  np.testing.assert_allclose(dist.logpdf(y), np.log(copula_density), rtol=1e-14)


def test_logpdf_survives_a_copula_density_that_underflows() -> None:
  """The copula half is read in log space, as the marginal half already was.

  Closes vinecopulib/pyvinecopulib#335: by the time a density arrives the
  product over the vine's edges has already underflowed to zero, and no
  logarithm applied afterwards can recover it. A one-truncated D-vine of
  strongly dependent Gaussian pairs, evaluated with the arguments alternating
  between the two tails, puts the log-density near ``-766`` -- past ``log`` of
  the smallest subnormal, ``-744.44``.
  """
  d = 10
  pairs = [
    [
      pv.Bicop.from_family(pv.families.gaussian, parameters=np.array([[0.9]]))
      for _ in range(d - 1)
    ]
  ]
  structure = pv.DVineStructure(list(range(1, d + 1)), trunc_lvl=1)
  copula = pv.Vinecop.from_structure(structure=structure, pair_copulas=pairs)
  dist = Vinedist(copula, [stats.uniform()] * d)
  y = np.tile(np.array([0.999, 0.001]), d // 2)[None, :]

  assert np.asarray(copula.pdf(y))[0] == 0.0  # noqa: RUF069 - an exact guarantee, not a computed approximation
  logpdf = np.asarray(dist.logpdf(y))[0]
  assert np.isfinite(logpdf)
  # Uniform margins contribute nothing, so the joint log-density is the
  # copula's own.
  np.testing.assert_allclose(logpdf, copula.logpdf(y)[0], rtol=1e-12)
  assert np.isfinite(np.asarray(dist.loglik(y)))


# --- conditional sampling --------------------------------------------------- #


def test_sample_conditional_matches_the_copula_scale(
  continuous: np.ndarray,
) -> None:
  """The data scale is the copula scale plus the marginal transforms, exactly.

  Same seeds, so the base uniforms are the same draw and the two agree to the
  bit rather than in distribution.
  """
  dist = pv.Vinedist.from_data(continuous)
  tail = int(dist.vinecop.structure.order[-1])
  y_cond = np.full((7, 1), float(np.median(continuous[:, tail - 1])))

  got = dist.sample_conditional(y_cond, seeds=[1, 2, 3])

  margin: Any = dist.margins[tail - 1]
  u_cond = np.asarray(margin.cdf(y_cond[:, 0])).reshape(-1, 1)
  reference = dist.marginal_icdf(
    np.asarray(widen(dist.vinecop).sample_conditional(u_cond, seeds=[1, 2, 3]))
  )
  np.testing.assert_array_equal(got, reference)


def test_sample_conditional_returns_the_conditioners_it_was_given(
  data: np.ndarray,
) -> None:
  """A conditioner comes back as itself, up to its own `cdf`/`icdf` round trip.

  Two of three variables are held, which is the widest set the copula accepts --
  conditioning on all of them would leave nothing to draw.
  """
  dist = from_margins(data, [Kde1d(), Kde1d(), stats.poisson(3.0)])
  y_cond = np.column_stack(
    [np.linspace(-1.0, 1.0, 6), np.linspace(0.5, 2.0, 6)]
  )
  out = dist.sample_conditional(y_cond, conditioning_set=[1, 2])
  assert out.shape == (6, 3)
  np.testing.assert_allclose(out[:, :2], y_cond, atol=1e-6)


def test_sample_conditional_derives_a_discrete_left_limit(
  data: np.ndarray,
) -> None:
  """A discrete conditioner needs no left-limit column from the caller.

  On the copula scale that column is mandatory; here it comes from the
  variable's own margin, so `y_cond` stays one column per conditioner.
  """
  dist = from_margins(data, [Kde1d(), Kde1d(), stats.poisson(3.0)])
  assert dist.var_types[2] == "d"
  out = dist.sample_conditional(
    np.full((8, 1), 3.0), conditioning_set=[3], seeds=[4, 5, 6]
  )
  # Reproduced up to its atom, since the margin's `icdf` lands on the lattice.
  np.testing.assert_array_equal(out[:, 2], np.full(8, 3.0))

  # The width the copula scale would demand is rejected here: two columns name
  # two conditioning variables, not one variable and its left limit.
  with pytest.raises(ValueError, match="names 1 variables but y_cond has 2"):
    dist.sample_conditional(np.full((8, 2), 3.0), conditioning_set=[3])


def test_sample_conditional_validates_its_arguments(
  continuous: np.ndarray,
) -> None:
  """A bad conditioning specification is refused, not guessed at."""
  dist = pv.Vinedist.from_data(continuous)
  with pytest.raises(ValueError, match="must be two-dimensional"):
    dist.sample_conditional(np.zeros(5))
  with pytest.raises(ValueError, match=r"must be in 1, \.\.\., 2"):
    dist.sample_conditional(np.zeros((5, 1)), conditioning_set=[3])
  with pytest.raises(ValueError, match="invalid number of columns"):
    dist.sample_conditional(np.zeros((5, 2)))


# --- reporting -------------------------------------------------------------- #


def test_margin_summary_has_a_row_per_variable(data: np.ndarray) -> None:
  """Every variable is described, whether its margin chose a family or not."""
  dist = from_margins(data, [Kde1d(), stats.norm(0.0, 1.0), stats.poisson(3.0)])
  rows = dist.margin_summary()
  assert [row["variable"] for row in rows] == [0, 1, 2]
  assert [row["var_type"] for row in rows] == dist.var_types
  assert [row["family"] for row in rows] == ["kde1d", "norm", "poisson"]
  # A fitted margin reports the log-likelihood it attained; a fixed one has no
  # fit to report, and says so with None rather than a number.
  assert isinstance(rows[0]["loglik"], float)
  assert rows[1]["loglik"] is None


# --- discrete margins ------------------------------------------------------- #


def test_discrete_margin_builds_the_compact_layout(data: np.ndarray) -> None:
  """One extra column per variable with atoms, appended after the first block."""
  dist = from_margins(data, [Kde1d(), Kde1d(), stats.poisson(3.0)])
  assert dist.var_types == ["c", "c", "d"]
  layout: Any = dist.copula_layout(data)
  assert layout.shape == (data.shape[0], 4)
  # The first block is the cdf values; the trailing column is the left limit.
  np.testing.assert_allclose(layout[:, :3], dist.marginal_cdf(data), atol=1e-12)
  assert np.all(layout[:, 3] <= layout[:, 2] + 1e-12)
  np.testing.assert_allclose(
    dist.logpdf(data), _sklar_logpdf(dist, data), atol=0.0
  )


def test_all_discrete_margins(data: np.ndarray) -> None:
  """A fully discrete model needs `2d` columns and still evaluates."""
  counts = np.round(np.abs(data)).astype(float)
  dist = from_margins(
    counts, [stats.poisson(1.0), stats.poisson(2.0), stats.poisson(1.0)]
  )
  assert dist.var_types == ["d", "d", "d"]
  assert dist.copula_layout(counts).shape == (counts.shape[0], 6)
  np.testing.assert_allclose(
    dist.logpdf(counts), _sklar_logpdf(dist, counts), atol=0.0
  )


def test_simulate_respects_a_discrete_margin(data: np.ndarray) -> None:
  """Draws through a count margin land on the lattice."""
  dist = from_margins(data, [Kde1d(), Kde1d(), stats.poisson(3.0)])
  drawn = dist.sample(200, seeds=[1, 2, 3])
  assert drawn.shape == (200, 3)
  np.testing.assert_array_equal(drawn[:, 2], np.round(drawn[:, 2]))


def test_cdf_accepts_a_margin_with_atoms(data: np.ndarray) -> None:
  """A copula with atoms validates the whole layout, so the layout is what it
  gets -- even though a distribution function reads only the first block."""
  dist = from_margins(data, [Kde1d(), Kde1d(), stats.poisson(3.0)])
  values = np.asarray(dist.cdf(data[:20], N=2000, seeds=[1, 2, 3]))
  assert values.shape == (20,)
  assert np.all((values >= 0.0) & (values <= 1.0))


def test_copula_data_needs_no_copula(data: np.ndarray) -> None:
  """The layout and the var_types are available before a copula exists."""
  margins: list[Any] = [
    Kde1d().fit(data[:, 0]),
    Kde1d().fit(data[:, 1]),
    stats.poisson(3.0),
  ]
  assert pv.Vinedist.copula_var_types(margins) == ["c", "c", "d"]
  layout = pv.Vinedist.copula_data(margins, data)
  assert layout.shape == (data.shape[0], 4)

  # Which is exactly the workflow of fitting your own copula and wrapping it.
  copula = pv.Vinecop.from_data(layout, var_types=["c", "c", "d"])
  dist = pv.Vinedist(copula, margins)
  np.testing.assert_array_equal(layout, dist.copula_layout(data))


def test_left_limit_above_the_cdf_is_refused_on_both_paths() -> None:
  """A margin that reports `F(x^-) > F(x)` is caught at the boundary.

  `copula_data` refused it while `sample_conditional` assembled the same block
  by hand without the check -- so a bad left limit was caught everywhere
  except the one path where it puts a conditioner outside its own atom. Both
  go through `copula_data` now, so one margin drives both.
  """

  class _Broken(FlatMargin):
    var_type = "d"

    def cdf(self, y: Any, /, *, x: Any | None = None) -> Any:
      return np.full_like(np.asarray(y, dtype=float), 0.3)

    def cdf_left(self, y: Any, /, *, x: Any | None = None) -> Any:
      return np.full_like(np.asarray(y, dtype=float), 0.9)

  # A discrete copula, so the var_types cross-check passes and the layout
  # builder is what rejects the margin.
  u = np.random.default_rng(0).uniform(size=(50, 2))
  layout = np.column_stack([u, u * 0.9])
  copula = pv.Vinecop.from_data(layout, var_types=["d", "d"])
  dist = pv.Vinedist(copula, [_Broken(), _Broken()])
  calls = (
    lambda: dist.copula_layout(np.ones((5, 2), dtype=float)),
    lambda: dist.sample_conditional(np.ones((5, 1), dtype=float)),
  )
  for call in calls:
    with pytest.raises(ValueError, match="cdf_left > cdf"):
      call()


# --- construction ----------------------------------------------------------- #


def test_margins_may_mix_fitted_and_unfitted(continuous: np.ndarray) -> None:
  """A fixed margin stays fixed; an unfitted one gets estimated."""
  fixed = stats.norm(0.0, 1.0)
  dist = from_margins(continuous, [fixed, Kde1d()])
  # The fixed margin is untouched, so its cdf is still the standard normal's.
  np.testing.assert_allclose(
    np.asarray(dist.margins[0].cdf(np.array([0.0]))), 0.5, atol=1e-12
  )
  fitted: Any = dist.margins[1]
  assert fitted.is_fitted


def test_each_variable_gets_its_own_margin(continuous: np.ndarray) -> None:
  """One `margin_class` must not carry a fit between variables."""
  dist = pv.Vinedist.from_data(continuous)
  first: Any = dist.margins[0]
  second: Any = dist.margins[1]
  assert first is not second
  # Distinct bandwidths prove they were fitted independently.
  assert first.bandwidth != second.bandwidth


def _assert_each_margin_matches_its_own_column(
  dist: Any, y: np.ndarray
) -> None:
  """Every margin must be the closest fit to its own column, not another's.

  Stated as "closest of the candidates" rather than as a tolerance on the
  median, because a kernel density on a skewed column is legitimately some way
  off its sample median while still being unmistakably that column's margin.
  """
  medians = [float(np.median(y[:, k])) for k in range(y.shape[1])]
  for j, margin in enumerate(dist.margins):
    center = float(margin.icdf(np.array([0.5]))[0])
    closest = min(range(len(medians)), key=lambda k: abs(center - medians[k]))
    assert closest == j, (
      f"margin {j} centered at {center}, closest to column {closest}"
    )


def test_the_constructor_copies_an_unfitted_broadcast_margin(
  continuous: np.ndarray,
) -> None:
  """An unfitted prototype cannot be shared: `fit` estimates it in place.

  Sharing one made every column carry the fit from the *last* column, since
  each estimate overwrote the previous one, and `logpdf` was then `-inf`
  wherever the columns were on different scales.
  """
  copula = pv.Vinecop.from_data(pv.to_pseudo_obs(continuous))
  dist = Vinedist(copula, Kde1d())
  assert len({id(m) for m in dist.margins}) == 2

  dist.fit(continuous)
  first: Any = dist.margins[0]
  second: Any = dist.margins[1]
  assert first is not second
  assert np.isfinite(dist.logpdf(continuous)).all()
  _assert_each_margin_matches_its_own_column(dist, continuous)


def test_the_constructor_shares_a_fitted_broadcast_margin(
  continuous: np.ndarray,
) -> None:
  """A fitted margin standing for every variable ties their parameters."""
  copula = pv.Vinecop.from_data(pv.to_pseudo_obs(continuous))
  shared = Kde1d().fit(continuous[:, 0])
  dist = Vinedist(copula, shared)
  assert all(m is shared for m in dist.margins)


@pytest.mark.parametrize("verb", ["fit", "select"])
def test_reestimating_separates_margins_tied_at_construction(
  continuous: np.ndarray, verb: str
) -> None:
  """Each column is estimated from its own data, so tied margins come apart."""
  copula = pv.Vinecop.from_data(pv.to_pseudo_obs(continuous))
  shared = Kde1d().fit(continuous[:, 0])
  dist = Vinedist(copula, shared)

  getattr(dist, verb)(continuous)
  assert len({id(m) for m in dist.margins}) == 2
  assert np.isfinite(dist.logpdf(continuous)).all()
  _assert_each_margin_matches_its_own_column(dist, continuous)


def test_an_explicitly_aliased_sequence_is_taken_as_given(
  continuous: np.ndarray,
) -> None:
  """Repeating one entry is how tied parameters are requested explicitly."""
  copula = pv.Vinecop.from_data(pv.to_pseudo_obs(continuous))
  shared = Kde1d().fit(continuous[:, 0])
  dist = Vinedist(copula, [shared, shared])
  assert dist.margins[0] is dist.margins[1] is shared


def test_each_half_is_weighted_by_its_own_controls(
  continuous: np.ndarray,
) -> None:
  """Weights ride in the controls, so each half carries its own -- or none.

  There is no propagation rule, which is the point: one controls object per
  part means a caller may weight the margins and the copula differently, or
  weight one and leave the other alone.
  """
  w = np.where(continuous[:, 0] > 0, 3.0, 1.0)
  plain = pv.Vinedist.from_data(continuous)
  copula_only = pv.Vinedist.from_data(
    continuous, pv.FitControlsVinecop(weights=w)
  )
  margins_only = pv.Vinedist.from_data(
    continuous, margin_controls=FitControlsKde1d(weights=w)
  )

  # Weighting the copula leaves the margins exactly where they were, and moves
  # the copula; weighting the margins moves the margins.
  np.testing.assert_array_equal(
    copula_only.marginal_cdf(continuous), plain.marginal_cdf(continuous)
  )
  assert not np.allclose(
    widen(copula_only.vinecop).get_pair_copula(0, 0).parameters,
    widen(plain.vinecop).get_pair_copula(0, 0).parameters,
  )
  assert not np.allclose(
    margins_only.marginal_cdf(continuous), plain.marginal_cdf(continuous)
  )


def test_weights_on_a_margin_that_cannot_use_them_raises(
  continuous: np.ndarray,
) -> None:
  """Silently dropping weights would fit a different model than requested."""

  class Unweighted(Vinedist):
    margin_class = _needs_fitting

  with pytest.raises(TypeError, match="honors no observation weights"):
    Unweighted.from_data(
      continuous,
      margin_controls=FitControlsMargin(weights=np.ones(continuous.shape[0])),
    )


def _needs_fitting() -> Any:
  """An unfitted margin that does not accept weights.

  Both declarations are the premise of the test above rather than scenery, so
  they are made here instead of on the shared `FlatMargin`.
  """

  class _Unweighted(FlatMargin):
    supports_weights = False

    @property
    def is_fitted(self) -> bool:
      return False

  return _Unweighted()


# --- transforms and consistency --------------------------------------------- #


def test_marginal_transforms_round_trip(continuous: np.ndarray) -> None:
  """`marginal_icdf` inverts `marginal_cdf` on continuous margins."""
  dist = pv.Vinedist.from_data(continuous)
  back = dist.marginal_icdf(dist.marginal_cdf(continuous))
  np.testing.assert_allclose(back, continuous, rtol=1e-4)


def test_rosenblatt_round_trip(continuous: np.ndarray) -> None:
  """The Rosenblatt transform inverts on the original scale."""
  dist = pv.Vinedist.from_data(continuous)
  head = continuous[:40]
  np.testing.assert_allclose(
    dist.inverse_rosenblatt(dist.rosenblatt(head)), head, rtol=1e-4
  )


def test_cdf_is_a_distribution_function(continuous: np.ndarray) -> None:
  """Monotone in each argument and within the unit interval."""
  dist = pv.Vinedist.from_data(continuous)
  grid = np.column_stack([np.linspace(-2.0, 2.0, 9), np.linspace(0.2, 6.0, 9)])
  values = np.asarray(dist.cdf(grid, N=20000, seeds=[1, 2, 3]))
  assert np.all((values >= 0.0) & (values <= 1.0))
  assert np.all(np.diff(values) >= -1e-3)


def test_dimension_mismatch_is_refused(continuous: np.ndarray) -> None:
  """The margin count must match the copula's dimension."""
  copula = pv.Vinecop.from_data(np.asarray(pv.to_pseudo_obs(continuous)))
  with pytest.raises(ValueError, match="2-dimensional copula"):
    pv.Vinedist(copula, [Kde1d().fit(continuous[:, 0])])


def test_var_type_mismatch_with_the_copula_is_refused(
  continuous: np.ndarray,
) -> None:
  """A continuous copula cannot host a margin with atoms."""
  copula = pv.Vinecop.from_data(np.asarray(pv.to_pseudo_obs(continuous)))
  with pytest.raises(ValueError, match="var_types"):
    pv.Vinedist(copula, [stats.poisson(3.0), stats.poisson(3.0)])


def test_wrong_column_count_is_refused(continuous: np.ndarray) -> None:
  """Evaluation validates the shape it was handed."""
  dist = pv.Vinedist.from_data(continuous)
  with pytest.raises(ValueError, match=r"shape \(n, 2\)"):
    dist.logpdf(np.ones((10, 3)))
  with pytest.raises(ValueError, match=r"shape \(n, 2\)"):
    pv.Vinedist.copula_data(dist.margins, np.ones((10, 3)))


def test_repr_names_the_margin_families(continuous: np.ndarray) -> None:
  """`repr` shows the dimension and what each margin is."""
  dist = pv.Vinedist.from_data(continuous)
  assert repr(dist) == "Vinedist(dim=2, margins=[kde1d, kde1d])"


# --- family selection ------------------------------------------------------- #


def _families(dist: Any) -> list[str]:
  """The family each margin settled on.

  Read off `Any`: `family_name` is an optional capability rather than part of
  the margin contract.
  """
  return [margin.family_name for margin in dist.margins]


def test_from_data_reads_dataframe_column_names(data: np.ndarray) -> None:
  """A DataFrame's own columns become the fitted distribution's names."""
  pd = pytest.importorskip("pandas")
  df = pd.DataFrame(data, columns=["real", "positive", "count"])
  dist = pv.Vinedist.from_data(df)
  assert dist.dim == 3
  assert dist.var_names == ["real", "positive", "count"]


def test_parametric_margins_choose_a_family_per_variable(
  data: np.ndarray,
) -> None:
  """Each variable's family is chosen from its own column, not from the first.

  The first column is normal and the second is the exponential of one, so the
  two land in different support groups from the same call -- and
  `margin_summary` is what makes that readable, one row per variable.
  """
  dist = ParametricVinedist.from_data(data[:, :2], names=["real", "positive"])
  assert all(isinstance(m, SciPyMargin) for m in dist.margins)
  rows = dist.margin_summary()
  assert [row["variable"] for row in rows] == [0, 1]
  assert [row["family"] for row in rows] == ["norm", "lognorm"]


def test_from_data_honors_a_named_family_and_chooses_an_unnamed_one(
  continuous: np.ndarray,
) -> None:
  """A named family is what the caller asked for; `from_data` keeps it.

  The second column is lognormal, so the two spellings separate here: naming
  `norm` for it gets a normal, and leaving the family out gets the search.
  """
  named = from_margins(continuous, [SciPyMargin("norm"), SciPyMargin("norm")])
  assert _families(named) == ["norm", "norm"]

  chosen = from_margins(continuous, [SciPyMargin(), SciPyMargin()])
  assert _families(chosen) == ["norm", "lognorm"]


def test_a_margin_that_only_selects_is_not_treated_as_fixed(
  continuous: np.ndarray,
) -> None:
  """Overriding `select` and leaving `fit` raising is an estimator too.

  A margin choosing between *kinds* of model -- a parametric family against a
  kernel density, say -- has nothing for `fit` to re-estimate, so it overrides
  `select` alone. Reading "no `fit`" as "fixed" left it unfitted and the
  distribution evaluated a margin that had chosen nothing.
  """

  class _Chooses(MarginBase[np.ndarray]):
    """Chooses a kernel density; `fit` stays the raising base one."""

    def __init__(self) -> None:
      self.chosen: Any = None

    @property
    def is_fitted(self) -> bool:
      return self.chosen is not None

    def pdf(self, y: np.ndarray, /, *, x: Any = None) -> np.ndarray:
      return np.asarray(self.chosen.pdf(y), dtype=float)

    def cdf(self, y: np.ndarray, /, *, x: Any = None) -> np.ndarray:
      return np.asarray(self.chosen.cdf(y), dtype=float)

    def select(
      self,
      y: np.ndarray,
      /,
      controls: Any = None,
      *,
      var_type: str | None = None,
      support: tuple[float | None, float | None] | None = None,
      x: Any = None,
      weights: Any = None,
    ) -> _Chooses:
      del controls, var_type, support, x, weights
      self.chosen = Kde1d().fit(np.asarray(y, dtype=float))
      return self

  class _Chooser(Vinedist):
    margin_class = _Chooses

  dist = _Chooser.from_data(continuous)
  assert all(cast("Any", m).is_fitted for m in dist.margins)
  # And `fit`, which has nothing to re-estimate here, says so rather than
  # quietly leaving the margins where they were.
  with pytest.raises(
    NotImplementedError, match=r"_Chooses\.fit is not defined"
  ):
    dist.fit(continuous)


def test_fit_keeps_the_family_where_select_replaces_it(
  continuous: np.ndarray,
) -> None:
  """The data-scale analog of `Vinecop.fit` versus `Vinecop.select`.

  Both verbs keep a family the caller named, since naming it is the choice.
  What separates them is `family_set`: it asks `select` to search again, and
  `fit` has nothing to search with -- so `fit` refuses it rather than dropping
  it and returning a model the controls do not describe.
  """
  margins = [SciPyMargin("norm").fit(continuous[:, j]) for j in range(2)]
  # Already fitted, so `from_data` leaves both alone and the wrong family on
  # the second column survives to be re-estimated below.
  dist = from_margins(continuous, margins)
  assert _families(dist) == ["norm", "norm"]

  wider = FitControlsMargin(family_set=["norm", "lognorm"])
  with pytest.raises(TypeError, match="family_set= would be ignored"):
    dist.fit(continuous, margin_controls=wider)
  assert _families(dist.select(continuous, margin_controls=wider)) == [
    "norm",
    "lognorm",
  ]

  # A controls object naming no family set asks for no search, so `fit` has
  # nothing to refuse and re-estimates the families the margins now hold.
  quiet = FitControlsMargin(selection_criterion="bic")
  assert _families(dist.fit(continuous, margin_controls=quiet)) == [
    "norm",
    "lognorm",
  ]


def test_a_margin_given_per_variable_is_fitted_in_place(
  continuous: np.ndarray,
) -> None:
  """One margin per variable is the caller's own object, estimated in place.

  A margin is both the specification and the fitted object, so the sequence
  form hands ownership over -- unlike the broadcast form, which has to copy.
  """
  spec = SciPyMargin()
  dist = from_margins(continuous, [spec, Kde1d()])
  assert dist.margins[0] is spec
  assert spec.is_fitted and spec.family_name == "norm"


@pytest.mark.parametrize(
  ("spec", "expected"),
  [
    (FitControlsMargin(family_set=["logistic"]), ["logistic", "logistic"]),
    (
      [
        FitControlsMargin(family_set=["norm"]),
        FitControlsMargin(family_set=["logistic"]),
      ],
      ["norm", "logistic"],
    ),
    (
      {"positive": FitControlsMargin(family_set=["logistic"])},
      ["norm", "logistic"],
    ),
    ({1: FitControlsMargin(family_set=["logistic"])}, ["norm", "logistic"]),
  ],
  ids=["broadcast", "sequence", "by-name", "by-index"],
)
def test_margin_controls_are_resolved_per_variable(
  continuous: np.ndarray, spec: Any, expected: list[str]
) -> None:
  """`margin_controls` resolves by the same four shapes, form for form.

  A mapping configures the variables it addresses and leaves the rest to the
  curated search, which is what lets one call constrain the one variable whose
  family is known.
  """
  dist = ParametricVinedist.from_data(
    continuous,
    margin_controls=spec,
    names=["real", "positive"],
  )
  assert _families(dist) == expected


def test_declared_supports_bound_the_default_kde_margin(
  continuous: np.ndarray,
) -> None:
  """A declared support reaches a margin the library itself constructs.

  Which is what makes a bounded default reachable without naming a class: the
  second column is positive, and an unbounded kernel density pads its grid past
  zero, so the draws go where nothing can occur. The declaration is per
  variable and keyword-only, exactly as `var_types` is on `Vinecop.from_data`.
  """
  plain = pv.Vinedist.from_data(continuous)
  bounded = pv.Vinedist.from_data(continuous, supports=[None, (0.0, None)])
  margin: Any = bounded.margins[1]
  assert isinstance(margin, Kde1d) and margin.xmin == 0.0  # noqa: RUF069 - the value this test set, read back
  assert plain.sample(500, seeds=[1, 2, 3])[:, 1].min() < 0.0
  assert bounded.sample(500, seeds=[1, 2, 3])[:, 1].min() >= 0.0


def test_margin_controls_resolution_is_refused_by_name(
  continuous: np.ndarray,
) -> None:
  """A misresolved marginal configuration says which argument it came from."""
  with pytest.raises(ValueError, match="margin_controls mapping names 'z'"):
    pv.Vinedist.from_data(
      continuous,
      margin_controls={"z": FitControlsMargin()},
      names=["a", "b"],
    )
  with pytest.raises(ValueError, match="margin_controls has length 1"):
    pv.Vinedist.from_data(continuous, margin_controls=[FitControlsMargin()])


def test_margin_controls_fallback_substitutes_a_kde_margin(
  continuous: np.ndarray,
) -> None:
  """An impossible family set fails loudly, or falls back once and says so.

  `beta` lives on the unit interval, so neither column admits it. Answering a
  parametric request nonparametrically is a downgrade the caller has to ask
  for, which is why the default is the refusal.
  """
  with pytest.raises(ValueError, match="no parametric family fits"):
    ParametricVinedist.from_data(
      continuous,
      margin_controls=FitControlsMargin(family_set=["beta"]),
    )

  with pytest.warns(UserWarning, match="kernel-density margin was substituted"):
    fallen_back = ParametricVinedist.from_data(
      continuous,
      margin_controls=FitControlsMargin(
        family_set=["beta"], on_failure="fallback"
      ),
    )
  assert all(isinstance(m, Kde1d) for m in fallen_back.margins)


# --- exogenous covariates ---------------------------------------------------- #


class _ConditionalVine(HostedVinecop):
  """A hosted vine whose pair copulas read external covariates."""

  supports_covariates = True


def _conditional_dist() -> tuple[Vinedist, GaussianBicop]:
  """Two conditional-normal margins and one externally conditional pair."""
  pair = GaussianBicop(scale=0.7, rho_max=0.75)
  structure = pv.RVineStructure.from_order([1, 2])
  copula = _ConditionalVine([[pair]], structure)
  return Vinedist(copula, [ShiftedNormalMargin(), ShiftedNormalMargin()]), pair


def test_full_y_given_x_matches_an_analytic_bivariate_normal() -> None:
  """The same row-aligned X reaches both margins and conditional dependence."""
  dist, _ = _conditional_dist()
  x = np.array([[-0.8], [-0.2], [0.3], [0.9]])
  y = np.array([[-1.1, 0.2], [0.4, -0.7], [1.2, 0.0], [1.6, 2.1]])
  z = y - x
  rho = 0.75 * np.tanh(0.7 * x[:, 0])
  one_minus = 1.0 - rho * rho
  expected = (
    -math.log(2.0 * math.pi)
    - 0.5 * np.log(one_minus)
    - (z[:, 0] ** 2 - 2.0 * rho * z[:, 0] * z[:, 1] + z[:, 1] ** 2)
    / (2.0 * one_minus)
  )
  np.testing.assert_allclose(dist.logpdf(y, x=x), expected, rtol=2e-12)

  # The inverse Rosenblatt transform is the sampling map. Check it against the
  # analytic conditional-normal representation, independently of the cascade.
  w = np.array([[0.2, 0.3], [0.4, 0.7], [0.6, 0.25], [0.8, 0.9]])
  z1, z2 = stats.norm.ppf(w[:, 0]), stats.norm.ppf(w[:, 1])
  expected_sample = np.column_stack(
    [x[:, 0] + rho * z2 + np.sqrt(one_minus) * z1, x[:, 0] + z2]
  )
  np.testing.assert_allclose(
    dist.inverse_rosenblatt(w, x=x), expected_sample, rtol=2e-9, atol=2e-9
  )


def test_full_y_given_x_sampler_matches_its_analytic_base_uniform_map() -> None:
  """`sample` forwards X through both the copula inverse and marginal quantiles."""
  dist, _ = _conditional_dist()
  n = 40
  x = np.linspace(-0.9, 0.9, n)[:, None]
  base = pv.utils.sample_uniform(n, 2, seeds=[17])
  z1, z2 = stats.norm.ppf(base[:, 0]), stats.norm.ppf(base[:, 1])
  rho = 0.75 * np.tanh(0.7 * x[:, 0])
  expected = np.column_stack(
    [
      x[:, 0] + rho * z2 + np.sqrt(1.0 - rho * rho) * z1,
      x[:, 0] + z2,
    ]
  )
  np.testing.assert_allclose(
    dist.sample(n, x=x, seeds=[17]), expected, rtol=2e-9, atol=2e-9
  )


def test_conditional_vinedist_cdf_surfaces_the_base_limitation() -> None:
  """A conditional copula does not imply a generic per-X Monte-Carlo CDF."""
  dist, _ = _conditional_dist()
  with pytest.raises(NotImplementedError, match="Conditional cdf"):
    dist.cdf(np.zeros((3, 2)), x=np.zeros((3, 1)))


def test_conditional_entry_points_reject_broadcasting_covariates() -> None:
  """Every composition path requires one two-dimensional X row per input row."""
  dist, _ = _conditional_dist()
  y = np.zeros((4, 2))
  for x in (np.zeros(4), np.zeros((1, 1))):
    calls = (
      lambda x=x: dist.marginal_cdf(y, x=x),
      lambda x=x: dist.marginal_icdf(np.full_like(y, 0.5), x=x),
      lambda x=x: dist.logpdf(y, x=x),
      lambda x=x: dist.cdf(y, x=x),
      lambda x=x: dist.rosenblatt(y, x=x),
      lambda x=x: dist.inverse_rosenblatt(np.full_like(y, 0.5), x=x),
      lambda x=x: dist.sample(4, x=x),
      lambda x=x: dist.sample_conditional(
        np.zeros((4, 1)), conditioning_set=[2], x=x
      ),
    )
    for call in calls:
      with pytest.raises(ValueError, match=r"one row per observation|shape"):
        call()


def test_from_data_rejects_misaligned_covariates_before_fitting() -> None:
  """The two-step fitter cannot broadcast one conditional design row."""
  y = np.zeros((4, 2))
  with pytest.raises(ValueError, match="one row per observation"):
    Vinedist.from_data(
      y,
      x=np.zeros((1, 1)),
    )


def test_covariates_reach_the_margins(continuous: np.ndarray) -> None:
  """Conditioning the margins moves the joint density, through every entry point."""
  copula = pv.Vinecop.from_data(np.asarray(pv.to_pseudo_obs(continuous)))
  dist = pv.Vinedist(copula, [ShiftedNormalMargin(), ShiftedNormalMargin()])
  y = continuous[:10]
  cov = np.full((10, 1), 0.5)

  # The Sklar sum, with every term conditioned on the covariates.
  manual = np.log(np.asarray(copula.pdf(dist.marginal_cdf(y, x=cov))))
  for j, m in enumerate(dist.margins):
    # `logpdf` is an optional capability, so it is read off `Any`.
    margin: Any = m
    manual = manual + margin.logpdf(y[:, j], x=cov)
  np.testing.assert_allclose(dist.logpdf(y, x=cov), manual, atol=1e-12)
  assert not np.allclose(dist.logpdf(y, x=cov), dist.logpdf(y))
  np.testing.assert_allclose(dist.pdf(y, x=cov), np.exp(dist.logpdf(y, x=cov)))
  np.testing.assert_allclose(
    dist.loglik(y, x=cov), np.sum(dist.logpdf(y, x=cov))
  )
  # Shifting the margins right can only lower F, up to the layout's clamp.
  shifted, plain = dist.marginal_cdf(y, x=cov), dist.marginal_cdf(y)
  assert np.all(shifted <= plain) and np.any(shifted < plain)
  assert np.all(dist.cdf(y, x=cov, N=2000, seeds=[1]) <= 1.0)
  # Both directions of the transform read the covariates.
  w = dist.rosenblatt(y, x=cov)
  assert not np.allclose(w, dist.rosenblatt(y))
  u_grid = np.column_stack([[0.2, 0.5, 0.8], [0.4, 0.5, 0.6]])
  np.testing.assert_allclose(
    dist.inverse_rosenblatt(w, x=cov),
    dist.marginal_icdf(np.asarray(copula.inverse_rosenblatt(w)), x=cov),
    atol=0,
  )
  # A pure location shift, so conditioning moves every quantile by 0.5.
  np.testing.assert_allclose(
    dist.marginal_icdf(u_grid, x=cov[:3]),
    dist.marginal_icdf(u_grid) + 0.5,
    atol=1e-6,
  )


def test_an_unconditional_copula_is_never_handed_covariates(
  continuous: np.ndarray,
) -> None:
  """`Vinecop` takes no conditioning matrix, so the check must omit it."""
  copula = pv.Vinecop.from_data(np.asarray(pv.to_pseudo_obs(continuous)))
  dist = pv.Vinedist(copula, [ShiftedNormalMargin(), ShiftedNormalMargin()])
  # A `TypeError` here would mean `x=` reached the compiled copula.
  assert dist.logpdf(continuous[:5], x=np.zeros((5, 1))).shape == (5,)


def test_from_data_fits_conditional_margins_on_the_covariates() -> None:
  """`from_data` forwards covariates to the margins that read them, only."""
  rng = np.random.default_rng(0)
  cov = rng.normal(size=(200, 1))
  y = np.column_stack([cov[:, 0] + rng.normal(size=200), rng.normal(size=200)])

  seen: list[Any | None] = []

  class _Recording(ShiftedNormalMargin):
    @property
    def is_fitted(self) -> bool:
      return False

    def fit(
      self,
      data: Any,
      /,
      controls: Any = None,
      *,
      var_type: str | None = None,
      support: tuple[float | None, float | None] | None = None,
      x: Any | None = None,
      weights: Any = None,
    ) -> Any:
      del var_type, support
      seen.append(x)
      return self

  dist = from_margins(y, [_Recording(), Kde1d()], x=cov)
  assert len(seen) == 1 and seen[0] is not None
  # The kde margin never declared covariates, so it was fitted plainly.
  kde: Any = dist.margins[1]
  assert kde.var_type == "c"


def test_covariates_nothing_reads_are_refused(continuous: np.ndarray) -> None:
  """Silently returning the unconditional answer is the failure to avoid."""
  dist = pv.Vinedist.from_data(continuous)  # Kde1d margins: unconditional
  cov = np.zeros((continuous.shape[0], 1))
  for call in (dist.logpdf, dist.pdf, dist.marginal_cdf):
    with pytest.raises(ValueError, match="supports_covariates"):
      call(continuous, x=cov)
  with pytest.raises(ValueError, match="supports_covariates"):
    dist.marginal_icdf(np.full_like(continuous, 0.5), x=cov)
  with pytest.raises(ValueError, match="the fit would ignore"):
    pv.Vinedist.from_data(continuous, x=cov)


def test_from_data_leaves_the_caller_s_controls_alone(
  continuous: np.ndarray,
) -> None:
  """A fit reads the controls it is handed; it never writes to them.

  Compared against what each object holds rather than against the input:
  `FitControlsVinecop` rescales weights to average one as it stores them,
  which is its own business and happens before any fit sees them.
  """
  weights = np.linspace(0.5, 1.5, continuous.shape[0])
  controls = pv.FitControlsVinecop(weights=weights)
  margin_controls = FitControlsKde1d(weights=weights)
  stored = np.array(controls.weights, copy=True)
  margin_stored = np.array(margin_controls.weights, copy=True)
  pv.Vinedist.from_data(continuous, controls, margin_controls=margin_controls)
  np.testing.assert_array_equal(controls.weights, stored)
  np.testing.assert_array_equal(margin_controls.weights, margin_stored)


def test_the_two_halves_may_be_weighted_differently() -> None:
  """The copula reads the weights on its own controls and nobody else's.

  Two subsamples with opposite dependence: whichever one the *copula's*
  weights favor is the one its parameter takes after, whatever the margins
  were weighted by.
  """
  rng = np.random.default_rng(23)
  positive = rng.multivariate_normal(
    [0.0, 0.0], [[1.0, 0.9], [0.9, 1.0]], size=250
  )
  negative = rng.multivariate_normal(
    [0.0, 0.0], [[1.0, -0.9], [-0.9, 1.0]], size=250
  )
  y = stats.norm.cdf(np.vstack([positive, negative]))
  favor_positive = np.r_[np.full(250, 10.0), np.ones(250)]
  favor_negative = np.r_[np.ones(250), np.full(250, 10.0)]

  def fitted(copula_weights: np.ndarray, margin_weights: np.ndarray) -> float:
    dist = Vinedist.from_data(
      y,
      pv.FitControlsVinecop(
        family_set=[pv.families.gaussian], weights=copula_weights
      ),
      margin_controls=FitControlsKde1d(weights=margin_weights),
    )
    return float(widen(dist.vinecop).get_pair_copula(0, 0).parameters[0, 0])

  # The margins' weighting is varied against a fixed copula weighting, and the
  # copula's against a fixed marginal one. Only the latter flips the sign.
  assert fitted(favor_positive, favor_positive) > 0.0
  assert fitted(favor_positive, favor_negative) > 0.0
  assert fitted(favor_negative, favor_positive) < 0.0


# ---------------------------------------------------------------------------
# The discrete cascade, end to end through Vinedist
# ---------------------------------------------------------------------------

#: Mixed continuous / discrete type patterns, as `test_core_discrete_cascade`
#: spells them; imported rather than duplicated where the cascade owns them.


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
  margins = [Kde1d(var_type="d", xmin=0.0).fit(counts)] + [
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


def test_json_round_trip_is_exact_for_every_shipped_margin(
  unique_json_path: Path,
) -> None:
  """`Vinedist` persists to JSON like the copula classes it composes.

  `Bicop` / `Vinecop` / `RVineStructure` have carried `to_json` / `to_file`
  since well before 1.0 and `docs/concepts.rst` presents that as the way to
  store a model; the objects this release adds had only `pickle`.
  """
  pytest.importorskip("scipy")

  rng = np.random.default_rng(0)
  cov = [[1.0, 0.7, 0.3], [0.7, 1.0, 0.5], [0.3, 0.5, 1.0]]
  x = rng.multivariate_normal([0.0, 0.0, 0.0], cov, size=400)
  q = x[:6]
  specs = {
    "default": [pv.core.Kde1d() for _ in range(3)],
    "parametric": [SciPyMargin("norm") for _ in range(3)],
    "selected": [SciPyMargin() for _ in range(3)],
    "mixed": [pv.core.Kde1d(), SciPyMargin("norm"), SciPyMargin()],
  }
  for label, margins in specs.items():
    dist = from_margins(x, margins)
    restored = pv.core.Vinedist.from_json(dist.to_json())
    # Exactly, not approximately: the stored grid and parameters are the model.
    np.testing.assert_array_equal(restored.logpdf(q), dist.logpdf(q), label)
    np.testing.assert_array_equal(
      restored.rosenblatt(q), dist.rosenblatt(q), label
    )
    assert [type(m).__name__ for m in restored.margins] == [
      type(m).__name__ for m in dist.margins
    ]


def test_to_file_selects_cbor_by_extension(unique_json_path: Path) -> None:
  """The extension rule is the one the compiled classes already follow."""
  rng = np.random.default_rng(1)
  x = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.6], [0.6, 1.0]], size=300)
  dist = pv.core.Vinedist.from_data(x)
  base = str(unique_json_path)
  for path in (base, base.replace(".json", ".cbor")):
    dist.to_file(path)
    np.testing.assert_array_equal(
      pv.core.Vinedist.from_file(path).logpdf(x[:5]), dist.logpdf(x[:5])
    )


def test_a_margin_without_to_json_is_refused_by_name() -> None:
  """A custom margin has to opt in, and is told how."""
  rng = np.random.default_rng(2)
  x = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.6], [0.6, 1.0]], size=300)
  copula = pv.core.Vinedist.from_data(x).vinecop

  dist = pv.core.Vinedist(copula, [FlatMargin(), pv.core.Kde1d().fit(x[:, 1])])
  with pytest.raises(TypeError, match="register_margin_json"):
    dist.to_json()


def test_a_registered_custom_margin_round_trips() -> None:
  """`register_margin_json` is the documented hook, so it must work."""
  rng = np.random.default_rng(3)
  x = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.6], [0.6, 1.0]], size=300)
  copula = pv.core.Vinedist.from_data(x).vinecop

  class Uniform(FlatMargin):
    def to_json(self) -> dict[str, Any]:
      return {"kind": "_TestUniform"}

  pv.core.register_margin_json("_TestUniform", lambda payload: Uniform())
  dist = pv.core.Vinedist(copula, [Uniform(), pv.core.Kde1d().fit(x[:, 1])])
  restored = pv.core.Vinedist.from_json(dist.to_json())
  assert isinstance(restored.margins[0], pv.core.MarginBase)


def test_an_unknown_margin_kind_and_a_bad_version_both_raise() -> None:
  """A format change must fail loudly rather than build a wrong model."""
  from pyvinecopulib.core import margin_from_json

  with pytest.raises(ValueError, match="no reader registered"):
    margin_from_json({"kind": "NotAMargin", "version": 1})
  with pytest.raises(ValueError, match="unsupported margin JSON version"):
    margin_from_json({"kind": "Kde1d", "version": 999})


# ---------------------------------------------------------------------------
# The extension-point triad: VinedistLike / VinedistBase / Vinedist
# ---------------------------------------------------------------------------


def test_both_shipped_distributions_satisfy_the_contract() -> None:
  # The contract is what downstream code types against, so both routes must
  # satisfy it -- and the sklearn estimators publish one as
  # `distribution_`. The name says *both*, so check both: the torch half
  # went untested here, which is the lane where a `ModuleList` of margins and
  # an `nn.Module` copula could plausibly diverge from the protocol.
  copula = pv.Vinecop.from_data(
    pv.utils.to_pseudo_obs(np.random.default_rng(0).normal(size=(200, 2)))
  )
  dist = Vinedist(copula, [Kde1d().fit(np.zeros(5)), Kde1d().fit(np.zeros(5))])
  assert isinstance(dist, VinedistLike)
  assert isinstance(dist, VinedistBase)

  torch = pytest.importorskip("torch")
  from pyvinecopulib.torch import TorchVinedist

  rng = np.random.default_rng(1)
  y = torch.as_tensor(rng.normal(size=(200, 2)) + rng.normal(size=(200, 1)))
  torch_dist = TorchVinedist.from_data(y)
  assert isinstance(torch_dist, VinedistLike)
  assert isinstance(torch_dist, VinedistBase)


def test_a_minimal_vinedist_base_subclass_needs_no_hook_to_evaluate(
  random_state: Any,
) -> None:
  # The answer to "what do I subclass to build a new kind of vine
  # distribution?". Evaluation needs no hook at all: a vine distribution is
  # determined by its two halves, so installing them is the whole job.
  class MyDist(VinedistBase[Any]):
    pass

  y = random_state.normal(size=(300, 3))
  u = pv.utils.to_pseudo_obs(y)
  margins = [Kde1d().fit(y[:, j]) for j in range(3)]
  dist = MyDist(pv.Vinecop.from_data(u), margins)

  assert isinstance(dist, VinedistLike)
  assert dist.dim == 3
  assert repr(dist).startswith("MyDist(dim=3")
  assert np.all(np.isfinite(dist.logpdf(y)))
  np.testing.assert_allclose(dist.pdf(y), np.exp(dist.logpdf(y)))
  assert dist.copula_layout(y).shape == (300, 3)
  # The Rosenblatt round trip holds on the data scale.
  np.testing.assert_allclose(
    dist.inverse_rosenblatt(dist.rosenblatt(y)), y, rtol=1e-6, atol=1e-6
  )


def test_a_subclass_without_fit_hooks_refuses_to_fit() -> None:
  # Fitting is the only namespace-specific half, so it is the only thing that
  # needs hooks -- and their absence is reported, not guessed around.
  class MyDist(VinedistBase[Any]):
    pass

  with pytest.raises(NotImplementedError, match="_coerce_fit_data"):
    MyDist.from_data(np.random.default_rng(0).normal(size=(50, 2)))


def test_base_from_json_names_the_missing_part_not_the_method() -> None:
  """Reading back is a declaration, so the refusal is about the declaration.

  `from_json` itself is concrete -- it decodes, checks the version and checks
  the `kind` -- so what a subclass can fail to supply is the copula class, and
  that is what the message has to name.
  """
  from pyvinecopulib.core._json import dumps

  class MyDist(VinedistBase[Any]):
    pass

  payload = dumps(
    {"kind": "MyDist", "version": 1, "copula": "{}", "margins": []}
  )
  with pytest.raises(NotImplementedError, match="names no `vinecop_class`"):
    MyDist.from_json(payload)

  # Declared but unable to read its own JSON is the other half, and the one
  # `TorchVinedist` is in -- named, rather than reported as a missing method.
  class Opaque:
    pass

  class Declared(VinedistBase[Any]):
    vinecop_class = Opaque

  with pytest.raises(NotImplementedError, match="Opaque has no `from_json`"):
    Declared.from_json(dumps({"kind": "Declared", "version": 1}))


def test_a_payload_is_refused_by_version_and_by_class() -> None:
  """Both checks live in the base now, so every subclass inherits them."""
  from pyvinecopulib.core._json import dumps

  with pytest.raises(ValueError, match="unsupported Vinedist JSON version"):
    pv.core.Vinedist.from_json(dumps({"kind": "Vinedist", "version": 999}))
  # A subclass's payload read as its base is quietly the wrong model.
  with pytest.raises(ValueError, match="written by 'Other', not 'Vinedist'"):
    pv.core.Vinedist.from_json(dumps({"kind": "Other", "version": 1}))


def test_a_non_finite_float_survives_the_round_trip() -> None:
  """Why the decode belongs to the base and the codec stays private.

  JSON has no literal for a non-finite float, so one travels as a string. An
  override using `json.loads` would read that string back where a `-inf`
  belongs, with no error to notice.
  """
  import json

  from pyvinecopulib.core._json import dumps, read_payload

  written = dumps({"kind": "K", "version": 1, "loglik": float("-inf")})
  assert read_payload(written, "K", kind="K")["loglik"] == float("-inf")
  # The trap, spelled out: the standard library gets the tagged object, not a
  # float, and nothing about that reads as an error.
  assert isinstance(json.loads(written)["loglik"], dict)

  # All three values, and both signs of infinity.
  for value in (float("-inf"), float("inf")):
    payload = dumps({"kind": "K", "version": 1, "v": value})
    assert read_payload(payload, "K", kind="K")["v"] == value
  nan = dumps({"kind": "K", "version": 1, "v": float("nan")})
  assert np.isnan(read_payload(nan, "K", kind="K")["v"])

  # The reason the tag is an object and not a marked string: payloads carry
  # arbitrary user text, and a marked string is a value user data can spell.
  for text in ("__nonfinite__:0", "__pyvinecopulib_nonfinite__", "-inf"):
    round_tripped = read_payload(
      dumps({"kind": "K", "version": 1, "name": text}), "K", kind="K"
    )["name"]
    assert round_tripped == text, round_tripped


def test_copula_var_types_dispatches_through_the_subclass() -> None:
  # `copula_data` used to reach `Vinedist.copula_var_types` by name, which
  # silently bypassed an override; it goes through `cls` now.
  seen: list[int] = []

  class Counting(Vinedist):
    @classmethod
    def copula_var_types(cls, margins: Any) -> list[str]:
      seen.append(1)
      return super().copula_var_types(margins)

  Counting.copula_data([Kde1d().fit(np.zeros(5))], np.zeros((3, 1)))
  assert seen, "copula_data must dispatch copula_var_types through `cls`"


def test_a_vinedist_base_subclass_fits_from_declared_parts(
  random_state: Any,
) -> None:
  # The payoff of naming the parts: a subclass declares which classes its two
  # halves are and inherits the whole two-step fit, with no hook but the array
  # coercion.
  class MyDist(VinedistBase[Any]):
    vinecop_class = pv.Vinecop
    margin_class = Kde1d

    @classmethod
    def _coerce_fit_data(cls, y: Any, controls: Any) -> Any:
      del controls
      return np.asarray(y, dtype=float)

  y = random_state.normal(size=(400, 3))
  dist = MyDist.from_data(y)
  assert isinstance(dist, VinedistLike)
  assert isinstance(dist.vinecop, pv.Vinecop)
  assert all(isinstance(m, Kde1d) for m in dist.margins)
  assert np.all(np.isfinite(dist.logpdf(y)))
  assert repr(dist).startswith("MyDist(dim=3")


def test_a_subclass_that_declares_only_its_parts_honors_weights() -> None:
  """Naming two parts that honor weights is all it takes to honor them.

  There is nothing for the subclass to declare or override: the controls go to
  the parts, and each part answers for itself. A distribution-level flag
  restating the question is what made this answer `False` for a class composed
  of `Vinecop` and `Kde1d`, both of which weight perfectly well.
  """

  class MyDist(VinedistBase[Any]):
    vinecop_class = pv.Vinecop
    margin_class = Kde1d

    @classmethod
    def _coerce_fit_data(cls, y: Any, controls: Any) -> Any:
      del controls
      return np.asarray(y, dtype=float)

  rng = np.random.default_rng(7)
  y = rng.normal(size=(200, 3))
  w = rng.uniform(0.5, 2.0, size=200)
  plain = MyDist.from_data(y)
  weighted = MyDist.from_data(
    y,
    pv.FitControlsVinecop(weights=w),
    margin_controls=FitControlsKde1d(weights=w),
  )
  assert np.all(np.isfinite(weighted.logpdf(y)))
  assert not np.allclose(plain.logpdf(y), weighted.logpdf(y))


def test_declaring_no_vinecop_class_reports_it() -> None:
  class MyDist(VinedistBase[Any]):
    # It names its margins, so the report it owes is about the other half.
    margin_class = Kde1d

    @classmethod
    def _coerce_fit_data(cls, y: Any, controls: Any) -> Any:
      del controls
      return np.asarray(y, dtype=float)

  with pytest.raises(NotImplementedError, match="vinecop_class"):
    MyDist.from_data(np.random.default_rng(0).normal(size=(60, 2)))


def test_vinedist_refuses_torch_parts() -> None:
  # The mirror of `TorchVinedist` refusing a NumPy copula: this class
  # evaluates on NumPy, so a torch part would be detached from its graph.
  pytest.importorskip("torch")
  import pyvinecopulib.torch as torch_mod

  u = pv.utils.to_pseudo_obs(np.random.default_rng(0).normal(size=(200, 2)))
  copula = pv.Vinecop.from_data(u)
  margins = [Kde1d().fit(np.zeros(5)), Kde1d().fit(np.zeros(5))]

  lifted = torch_mod.TorchVinecop.from_vinecop(copula)
  with pytest.raises(TypeError, match="TorchVinedist"):
    # A torch vine satisfies `VinecopLike[Tensor]`, not the `[ndarray]` this
    # class takes, so the cast is what a caller who skips type checking does.
    Vinedist(cast("Any", lifted), margins)

  # And a torch margin, on an otherwise fine NumPy copula.
  with pytest.raises(TypeError, match="TorchVinedist"):
    Vinedist(copula, [torch_mod.TorchKde1d(), torch_mod.TorchKde1d()])


def test_logpdf_reads_the_parts_namespace_not_the_inputs() -> None:
  """A torch copula hosting NumPy margins is legal, and must evaluate.

  ``marginal_cdf`` and ``copula_data`` both state the rule and take their
  namespace from the columns the parts returned; ``logpdf`` -- the primitive
  the whole data-scale surface routes through -- took it from the input, so
  ``numpy.log`` was applied to a tensor. That is silent on the CPU and raises
  the moment the grid tracks a gradient, which is the entire point of the
  torch lane.
  """
  torch = pytest.importorskip("torch")
  from pyvinecopulib.torch import TorchVinecop

  class NumpyNormal(MarginBase[Any]):
    """A NumPy margin: it answers in ndarray whatever it is handed."""

    def pdf(self, y: Any, /, *, x: Any | None = None) -> Any:
      a = np.asarray(y, dtype=float)
      return np.exp(-0.5 * a * a) / math.sqrt(2.0 * math.pi)

    def cdf(self, y: Any, /, *, x: Any | None = None) -> Any:
      a = np.asarray(y, dtype=float)
      return np.array([0.5 * (1.0 + math.erf(v / math.sqrt(2.0))) for v in a])

  class MixedDist(VinedistBase[Any]):
    """Torch copula, NumPy margins -- the configuration the rule allows."""

  y = np.random.default_rng(0).normal(size=(40, 2))
  copula = TorchVinecop.from_data(pv.to_pseudo_obs(y))
  dist = MixedDist(copula, [NumpyNormal(), NumpyNormal()])

  plain = dist.logpdf(y)
  assert bool(torch.isfinite(torch.as_tensor(plain)).all())

  # With the grid tracking grad there is no silent path: either the namespace
  # is right or NumPy reaches for `__array__` on a tensor and raises.
  # Typed against the evaluation-only contract, which carries no grid.
  pair = cast("Any", copula.get_pair_copula(0, 0))
  pair.interp_grid.values.requires_grad_(True)
  with_grad = dist.logpdf(y)
  assert with_grad.requires_grad


# --- fit re-estimates the parts it holds ------------------------------------- #


class _SelfFittingVine(HostedVinecop):
  """A hosted vine that can refit its own pairs, so `Vinedist.fit` can ask."""

  bicop_class = pv.Bicop


def _hosted_dist(y: np.ndarray) -> tuple[pv.Vinedist, Any]:
  """A `Vinedist` over a caller's own vine class and kernel-density margins."""
  structure = pv.RVineStructure.from_order(list(range(1, y.shape[1] + 1)))
  vine = _SelfFittingVine.from_data(pv.to_pseudo_obs(y), structure=structure)
  margins = [Kde1d.from_data(y[:, j]) for j in range(y.shape[1])]
  return pv.Vinedist(vine, margins), vine


def test_fit_re_estimates_the_copula_it_holds(continuous: np.ndarray) -> None:
  """`fit` must not swap a caller's own vine for the default one.

  The copula half is asked to re-estimate *itself*, so a hosted `VinecopLike`
  keeps its class and its identity across a refit -- the same promise `fit`
  keeps for a margin's family. Building a fresh `vinecop_class` instead would
  silently replace the part the caller composed the distribution from.
  """
  dist, vine = _hosted_dist(continuous)
  dist.fit(continuous)
  assert dist.vinecop is vine
  assert type(dist.vinecop) is _SelfFittingVine


def test_select_re_estimates_the_copula_it_holds(
  continuous: np.ndarray,
) -> None:
  """`select` re-selects the held copula's structure, still in place."""
  dist, vine = _hosted_dist(continuous)
  dist.select(continuous)
  assert dist.vinecop is vine


def test_fit_reports_a_copula_that_cannot_re_estimate_itself(
  continuous: np.ndarray,
) -> None:
  """A vine with no pair fitter says so, rather than being replaced."""
  structure = pv.RVineStructure.from_order(
    list(range(1, continuous.shape[1] + 1))
  )

  def fit_edge(
    tree: int,
    edge: int,
    u_e: Any,
    x_e: Any,
    var_types: Any = ("c", "c"),
  ) -> Any:
    del tree, edge, x_e, var_types
    return pv.Bicop.from_data(np.asarray(u_e))

  vine = HostedVinecop.from_data(
    pv.to_pseudo_obs(continuous), structure=structure, fit_edge=fit_edge
  )
  margins = [
    Kde1d.from_data(continuous[:, j]) for j in range(continuous.shape[1])
  ]
  dist = pv.Vinedist(vine, margins)
  with pytest.raises(ValueError, match="fit_edge` is required"):
    dist.fit(continuous)


def test_a_weighted_refit_still_weights_both_halves(
  continuous: np.ndarray,
) -> None:
  """A refit reads the controls the same way `from_data` does.

  `fit` asks the copula it *holds* rather than the class, so this is the one
  path where the part answering for the weights need not be an instance of
  anything the distribution names.
  """
  weights = np.linspace(0.5, 1.5, continuous.shape[0])
  flat = pv.Vinedist.from_data(continuous)
  flat.fit(continuous)
  weighted = pv.Vinedist.from_data(continuous)
  weighted.fit(
    continuous,
    pv.FitControlsVinecop(weights=weights),
    margin_controls=FitControlsKde1d(weights=weights),
  )
  grid = continuous[:20]
  assert not np.allclose(flat.pdf(grid), weighted.pdf(grid))


def test_covariates_reach_the_parts_that_declare_them_only() -> None:
  """A conditional fit is per part, and the copula half declares for itself.

  `Vinedist`'s copula is a `Vinecop` of compiled pair copulas, which models no
  covariates and takes no `x` argument at all -- so `x` must not be forwarded
  to it. Reaching it anyway would raise instead of fitting something, and
  dropping it silently is only correct because `from_data` refuses covariates
  nothing on the lane reads at all.
  """
  rng = np.random.default_rng(3)
  cov = rng.normal(size=(200, 1))
  y = np.column_stack([cov[:, 0] + rng.normal(size=200), rng.normal(size=200)])

  class _Refittable(ShiftedNormalMargin):
    """`ShiftedNormalMargin` plus the estimator `fit` re-runs (nothing to fit)."""

    def fit(
      self,
      data: Any,
      /,
      controls: Any = None,
      *,
      var_type: str | None = None,
      support: tuple[float | None, float | None] | None = None,
      x: Any | None = None,
    ) -> Any:
      del data, controls, var_type, support, x
      return self

  assert getattr(pv.Vinecop, "supports_covariates", False) is False
  dist = from_margins(y, [_Refittable(), _Refittable()], x=cov)
  # Both halves are there and the conditional margins were used, so the fit
  # completed rather than raising on the copula's missing `x`.
  assert isinstance(dist.vinecop, pv.Vinecop)
  assert all(getattr(m, "supports_covariates", False) for m in dist.margins)

  # And `fit` -- which re-estimates the held copula -- follows the same rule.
  dist.fit(y, x=cov)
  assert isinstance(dist.vinecop, pv.Vinecop)


def test_fit_holds_the_copula_families_where_select_re_chooses_them() -> None:
  """`fit` estimates the shape it holds; `select` is what may change it.

  The margin half has always refused a family search inside `fit` -- it raises
  on `family_set=` there. The copula half went through `from_data`, which
  re-searched, so one call meant two different things: a vine fitted to
  gaussian pairs came back independent. Asking the held copula to refit itself
  is what makes the two halves agree.
  """
  rng = np.random.default_rng(0)
  y = rng.normal(size=(400, 3)) + rng.normal(size=(400, 1))

  def families(dist: Any) -> list[str]:
    vine = dist.vinecop
    return [
      str(vine.get_pair_copula(t, e).family).split(".")[-1]
      for t in range(int(vine.structure.trunc_lvl))
      for e in range(vine.dim - 1 - t)
    ]

  controls = pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  dist = pv.Vinedist.from_data(y, controls=controls)
  assert set(families(dist)) == {"gaussian"}
  dist.fit(y)
  assert set(families(dist)) == {"gaussian"}

  # And `select` is the call that may change it.
  dist.select(y)
  assert families(dist)


def test_a_json_payload_is_read_back_by_the_class_that_wrote_it() -> None:
  """`to_json` records the class and `from_json` now checks it.

  A subclass is a different distribution, so loading its payload as this class
  is quietly the wrong model -- and the field was written from the start and
  read by nothing.
  """
  import json as _json

  rng = np.random.default_rng(0)
  dist = pv.Vinedist.from_data(rng.normal(size=(150, 2)))
  raw = dist.to_json()
  assert _json.loads(raw)["kind"] == "Vinedist"
  assert isinstance(pv.Vinedist.from_json(raw), pv.Vinedist)

  foreign = _json.dumps({**_json.loads(raw), "kind": "SomeOtherDist"})
  with pytest.raises(ValueError, match="written by 'SomeOtherDist'"):
    pv.Vinedist.from_json(foreign)


def test_a_margin_keeps_every_criterion_across_a_json_round_trip() -> None:
  """`bic` and `aicc` need the sample size, which the payload now carries.

  `loglik` and `aic` survived without it, so the loss showed up only on the
  two criteria that penalize by ``n`` -- and as a raise, not a wrong number.
  """
  from pyvinecopulib.core import margin_from_json, margin_to_json

  rng = np.random.default_rng(0)
  margin = SciPyMargin("norm").fit(rng.normal(size=300))
  back = margin_from_json(margin_to_json(margin))
  for name in ("loglik", "aic", "bic", "aicc"):
    assert getattr(back, name)() == pytest.approx(getattr(margin, name)())


def test_margin_summary_survives_a_margin_that_declines_a_field() -> None:
  """Every field is optional, and declining is a way of declaring.

  The docstring promises `None` for whatever a margin does not contribute, but
  `name` / `family_name` / `support` / `n_parameters` were read with a bare
  `getattr(..., None)`, which absorbs only `AttributeError` -- so a property
  that *raises* took the whole summary down, while `loglik()` beside it was
  already guarded. A margin wrapping a regressor with no well-defined free
  parameter count is the case that hits it.
  """

  class _Declines(ShiftedNormalMargin):
    @property
    def n_parameters(self) -> float:
      raise NotImplementedError("no well-defined free-parameter count")

    @property
    def support(self) -> tuple[float, float]:
      raise RuntimeError("support depends on covariates")

  dist = pv.Vinedist(
    pv.Vinecop.from_data(
      pv.to_pseudo_obs(np.random.default_rng(3).normal(size=(80, 2)))
    ),
    [_Declines(), _Declines()],
  )
  rows = dist.margin_summary()
  assert len(rows) == 2
  for row in rows:
    assert row["n_parameters"] is None
    assert row["support"] is None
    assert row["margin"] == "_Declines"
