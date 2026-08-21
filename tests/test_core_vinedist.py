"""Tests for `pyvinecopulib.core.Vinedist`.

The claims worth pinning are the identities. That `logpdf` really is
`log c(F(x)) + sum_j log pdf_j(x_j)`, checked against a hand computation rather
than a tolerance. And that the `(n, d + k)` layout a discrete margin
needs is assembled for the user, since assembling it by hand is what
`examples/04_discrete_variables.ipynb` currently has to do.
"""

from __future__ import annotations

import math
from typing import Any, Optional

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import Kde1d, MarginBase, Vinedist
from pyvinecopulib.margins import MarginSelector

# The discrete cascade owns these; the end-to-end test at the bottom reuses them
# rather than duplicating the type patterns and the parity bound.
from .test_core_discrete_cascade import _D, _assert_parity, _both

stats = pytest.importorskip("scipy.stats")


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


def test_logpdf_is_the_sklar_factorization(continuous: np.ndarray) -> None:
  """`logpdf` equals the copula term plus the marginal log-densities."""
  dist = pv.Vinedist.from_data(continuous, margins="kde")
  u = dist.marginal_cdf(continuous)
  manual = np.log(np.asarray(dist.copula.pdf(u)))
  for j, m in enumerate(dist.margins):
    manual = manual + np.log(np.asarray(m.pdf(continuous[:, j])))
  np.testing.assert_allclose(dist.logpdf(continuous), manual, atol=0.0)


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


# --- conditional sampling --------------------------------------------------- #


def test_sample_conditional_matches_the_copula_scale(
  continuous: np.ndarray,
) -> None:
  """The data scale is the copula scale plus the marginal transforms, exactly.

  Same seeds, so the base uniforms are the same draw and the two agree to the
  bit rather than in distribution.
  """
  dist = pv.Vinedist.from_data(continuous, margins="kde")
  tail = int(dist.copula.structure.order[-1])
  y_cond = np.full((7, 1), float(np.median(continuous[:, tail - 1])))

  got = dist.sample_conditional(y_cond, seeds=[1, 2, 3])

  margin: Any = dist.margins[tail - 1]
  u_cond = np.asarray(margin.cdf(y_cond[:, 0])).reshape(-1, 1)
  reference = dist.marginal_icdf(
    np.asarray(dist.copula.sample_conditional(u_cond, seeds=[1, 2, 3]))
  )
  np.testing.assert_array_equal(got, reference)


def test_sample_conditional_returns_the_conditioners_it_was_given(
  data: np.ndarray,
) -> None:
  """A conditioner comes back as itself, up to its own `cdf`/`icdf` round trip.

  Two of three variables are held, which is the widest set the copula accepts --
  conditioning on all of them would leave nothing to draw.
  """
  dist = pv.Vinedist.from_data(data, margins=["kde", "kde", stats.poisson(3.0)])
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
  dist = pv.Vinedist.from_data(data, margins=["kde", "kde", stats.poisson(3.0)])
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
  dist = pv.Vinedist.from_data(continuous, margins="kde")
  with pytest.raises(ValueError, match="must be two-dimensional"):
    dist.sample_conditional(np.zeros(5))
  with pytest.raises(ValueError, match="must be in 1, ..., 2"):
    dist.sample_conditional(np.zeros((5, 1)), conditioning_set=[3])
  with pytest.raises(ValueError, match="invalid number of columns"):
    dist.sample_conditional(np.zeros((5, 2)))


# --- reporting -------------------------------------------------------------- #


def test_margin_summary_has_a_row_per_variable(data: np.ndarray) -> None:
  """Every variable is described, whether its margin selected a family or not.

  `selection_report` covers only margins that *chose*, so an all-fixed
  distribution reports nothing there and everything here.
  """
  dist = pv.Vinedist.from_data(
    data, margins=["kde", stats.norm(0.0, 1.0), stats.poisson(3.0)]
  )
  rows = dist.margin_summary()
  assert [row["variable"] for row in rows] == [0, 1, 2]
  assert [row["var_type"] for row in rows] == dist.var_types
  assert [row["family"] for row in rows] == ["kde1d", "norm", "poisson"]
  # A fitted margin reports the log-likelihood it attained; a fixed one has no
  # fit to report, and says so with None rather than a number.
  assert isinstance(rows[0]["loglik"], float)
  assert rows[1]["loglik"] is None
  assert dist.selection_report() == []


# --- discrete margins ------------------------------------------------------- #


def test_discrete_margin_builds_the_compact_layout(data: np.ndarray) -> None:
  """One extra column per variable with atoms, appended after the first block."""
  dist = pv.Vinedist.from_data(data, margins=["kde", "kde", stats.poisson(3.0)])
  assert dist.var_types == ["c", "c", "d"]
  layout: Any = dist._u_layout(data)
  assert layout.shape == (data.shape[0], 4)
  # The first block is the cdf values; the trailing column is the left limit.
  np.testing.assert_allclose(layout[:, :3], dist.marginal_cdf(data), atol=1e-12)
  assert np.all(layout[:, 3] <= layout[:, 2] + 1e-12)
  assert np.all(np.isfinite(dist.logpdf(data)))


def test_all_discrete_margins(data: np.ndarray) -> None:
  """A fully discrete model needs `2d` columns and still evaluates."""
  counts = np.round(np.abs(data)).astype(float)
  dist = pv.Vinedist.from_data(
    counts, margins=[stats.poisson(1.0), stats.poisson(2.0), stats.poisson(1.0)]
  )
  assert dist.var_types == ["d", "d", "d"]
  assert dist._u_layout(counts).shape == (counts.shape[0], 6)
  assert np.all(np.isfinite(dist.logpdf(counts)))


def test_simulate_respects_a_discrete_margin(data: np.ndarray) -> None:
  """Draws through a count margin land on the lattice."""
  dist = pv.Vinedist.from_data(data, margins=["kde", "kde", stats.poisson(3.0)])
  drawn = dist.sample(200, seeds=[1, 2, 3])
  assert drawn.shape == (200, 3)
  np.testing.assert_array_equal(drawn[:, 2], np.round(drawn[:, 2]))


def test_cdf_accepts_a_margin_with_atoms(data: np.ndarray) -> None:
  """A copula with atoms validates the whole layout, so the layout is what it
  gets -- even though a distribution function reads only the first block."""
  dist = pv.Vinedist.from_data(data, margins=["kde", "kde", stats.poisson(3.0)])
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
  np.testing.assert_array_equal(layout, dist._u_layout(data))


def test_left_limit_above_the_cdf_is_refused() -> None:
  """A margin that reports `F(x^-) > F(x)` is caught at the boundary."""
  from pyvinecopulib.core import MarginBase

  class _Broken(MarginBase[np.ndarray]):
    @property
    def var_type(self) -> str:
      return "d"

    def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return np.full_like(np.asarray(y, dtype=float), 0.5)

    def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return np.full_like(np.asarray(y, dtype=float), 0.3)

    def cdf_left(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return np.full_like(np.asarray(y, dtype=float), 0.9)

  # A discrete copula, so the var_types cross-check passes and the layout
  # builder is what rejects the margin.
  u = np.random.default_rng(0).uniform(size=(50, 2))
  layout = np.column_stack([u, u * 0.9])
  copula = pv.Vinecop.from_data(layout, var_types=["d", "d"])
  with pytest.raises(ValueError, match="cdf_left > cdf"):
    pv.Vinedist(copula, [_Broken(), _Broken()])._u_layout(
      np.ones((5, 2), dtype=float)
    )


# --- construction ----------------------------------------------------------- #


def test_margins_may_mix_fitted_and_unfitted(continuous: np.ndarray) -> None:
  """A fixed margin stays fixed; an unfitted one gets estimated."""
  fixed = stats.norm(0.0, 1.0)
  dist = pv.Vinedist.from_data(continuous, margins=[fixed, Kde1d()])
  # The fixed margin is untouched, so its cdf is still the standard normal's.
  np.testing.assert_allclose(
    np.asarray(dist.margins[0].cdf(np.array([0.0]))), 0.5, atol=1e-12
  )
  fitted: Any = dist.margins[1]
  assert fitted.is_fitted


def test_margins_mapping_is_keyed_by_name(continuous: np.ndarray) -> None:
  """A mapping addresses variables by name over a default."""
  dist = pv.Vinedist.from_data(
    continuous, margins={"b": MarginSelector()}, names=["a", "b"]
  )
  assert isinstance(dist.margins[0], Kde1d)
  assert isinstance(dist.margins[1], MarginSelector)


def test_margins_mapping_rejects_an_unknown_name(
  continuous: np.ndarray,
) -> None:
  """Naming a variable that does not exist is an error, not a silent no-op."""
  with pytest.raises(ValueError, match="not a variable"):
    pv.Vinedist.from_data(continuous, margins={"z": "kde"}, names=["a", "b"])


def test_margins_sequence_length_is_checked(continuous: np.ndarray) -> None:
  """A short sequence is caught rather than broadcast."""
  with pytest.raises(ValueError, match="length 1, but there are 2"):
    pv.Vinedist.from_data(continuous, margins=["kde"])


def test_broadcast_margin_is_copied_not_shared(continuous: np.ndarray) -> None:
  """One prototype must not carry a fit between variables."""
  dist = pv.Vinedist.from_data(continuous, margins=Kde1d())
  first: Any = dist.margins[0]
  second: Any = dist.margins[1]
  assert first is not second
  # Distinct bandwidths prove they were fitted independently.
  assert first.bandwidth != second.bandwidth


def test_callable_margin_is_used_as_a_fitter(continuous: np.ndarray) -> None:
  """A plain callable receives the column and returns a margin."""
  dist = pv.Vinedist.from_data(
    continuous, margins=lambda col: stats.norm(col.mean(), col.std())
  )
  np.testing.assert_allclose(
    np.asarray(dist.margins[0].cdf(np.array([continuous[:, 0].mean()]))),
    0.5,
    atol=1e-12,
  )


def test_weights_reach_both_halves(continuous: np.ndarray) -> None:
  """Weighting changes the fit rather than being ignored."""
  w = np.where(continuous[:, 0] > 0, 3.0, 1.0)
  plain = pv.Vinedist.from_data(continuous, margins="kde")
  tilted = pv.Vinedist.from_data(continuous, margins="kde", weights=w)
  assert not np.allclose(
    plain.marginal_cdf(continuous), tilted.marginal_cdf(continuous)
  )


def test_weights_on_a_margin_that_cannot_use_them_raises(
  continuous: np.ndarray,
) -> None:
  """Silently dropping weights would fit a different model than requested."""
  with pytest.raises(TypeError, match="cannot use observation weights"):
    pv.Vinedist.from_data(
      continuous,
      margins=[stats.norm(0, 1), _needs_fitting()],
      weights=np.ones(continuous.shape[0]),
    )


def _needs_fitting() -> Any:
  """An unfitted margin that does not accept weights."""
  from pyvinecopulib.core import MarginBase

  class _Unweighted(MarginBase[np.ndarray]):
    supports_weights = False

    @property
    def is_fitted(self) -> bool:
      return False

    def fit(
      self, y: Any, *, x: Optional[Any] = None, weights: Any = None
    ) -> Any:
      return self

    def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return np.ones_like(np.asarray(y, dtype=float))

    def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return np.clip(np.asarray(y, dtype=float), 0.0, 1.0)

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
  dist = pv.Vinedist.from_data(continuous, margins="kde")
  assert repr(dist) == "Vinedist(dim=2, margins=[kde1d, kde1d])"


# --- selection reporting ---------------------------------------------------- #


def test_from_data_reads_dataframe_column_names(data: np.ndarray) -> None:
  """A mapping keyed by column name works on a DataFrame without `names=`."""
  pd = pytest.importorskip("pandas")
  df = pd.DataFrame(data, columns=["real", "positive", "count"])
  dist = pv.Vinedist.from_data(df, margins={"real": Kde1d()})
  assert dist.dim == 3


def test_selection_report_names_each_variable(data: np.ndarray) -> None:
  """Report rows identify their variable, so a multi-column table is readable."""
  dist = pv.Vinedist.from_data(
    data[:, :2], margins="parametric", names=["real", "positive"]
  )
  report = dist.selection_report()
  assert report
  assert {row["column"] for row in report} == {"real", "positive"}
  selected = {row["column"] for row in report if row["selected"]}
  assert selected == {"real", "positive"}


def test_selection_report_is_empty_without_a_selector(data: np.ndarray) -> None:
  """Margins that were given rather than selected contribute no rows."""
  assert pv.Vinedist.from_data(data).selection_report() == []


def test_from_data_leaves_the_caller_s_specification_alone(
  data: np.ndarray,
) -> None:
  """Naming a selector must not mutate an object the caller still holds."""
  selector = MarginSelector(candidates=["norm", "logistic"])
  dist = pv.Vinedist.from_data(data[:, :1], margins=[selector], names=["real"])
  assert selector.name is None
  assert list(selector.report_) == []
  assert dist.selection_report()[0]["column"] == "real"


# --- exogenous covariates ---------------------------------------------------- #


class _ShiftedNormal(MarginBase[np.ndarray]):
  """A conditional margin: a standard normal shifted by the first covariate."""

  supports_covariates = True

  def _shift(self, x: Optional[Any]) -> Any:
    return 0.0 if x is None else np.asarray(x, dtype=float)[:, 0]

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    z = np.asarray(y, dtype=float) - self._shift(x)
    return np.exp(-0.5 * z**2) / np.sqrt(2.0 * np.pi)

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    z = np.asarray(y, dtype=float) - self._shift(x)
    return 0.5 * (1.0 + np.vectorize(math.erf)(z / np.sqrt(2.0)))


def test_covariates_reach_the_margins(continuous: np.ndarray) -> None:
  """Conditioning the margins moves the joint density, through every entry point."""
  copula = pv.Vinecop.from_data(np.asarray(pv.to_pseudo_obs(continuous)))
  dist = pv.Vinedist(copula, [_ShiftedNormal(), _ShiftedNormal()])
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
  """`Vinecop` takes no conditioning matrix, so the gate must omit it."""
  copula = pv.Vinecop.from_data(np.asarray(pv.to_pseudo_obs(continuous)))
  dist = pv.Vinedist(copula, [_ShiftedNormal(), _ShiftedNormal()])
  # A `TypeError` here would mean `x=` reached the compiled copula.
  assert dist.logpdf(continuous[:5], x=np.zeros((5, 1))).shape == (5,)


def test_from_data_fits_conditional_margins_on_the_covariates() -> None:
  """`from_data` forwards covariates to the margins that read them, only."""
  rng = np.random.default_rng(0)
  cov = rng.normal(size=(200, 1))
  y = np.column_stack([cov[:, 0] + rng.normal(size=200), rng.normal(size=200)])

  seen: list[Optional[Any]] = []

  class _Recording(_ShiftedNormal):
    @property
    def is_fitted(self) -> bool:
      return False

    def fit(
      self, data: Any, *, x: Optional[Any] = None, weights: Any = None
    ) -> Any:
      seen.append(x)
      return self

  dist = pv.Vinedist.from_data(y, x=cov, margins=[_Recording(), "kde"])
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
  with pytest.raises(ValueError, match="no margin reads them"):
    pv.Vinedist.from_data(continuous, x=cov)


def test_from_data_leaves_the_caller_s_controls_alone(
  continuous: np.ndarray,
) -> None:
  """Weights must not be written into the controls object the caller owns."""
  controls = pv.FitControlsVinecop()
  weights = np.linspace(0.5, 1.5, continuous.shape[0])
  pv.Vinedist.from_data(continuous, weights=weights, controls=controls)
  assert len(controls.weights) == 0


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
