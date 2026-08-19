"""Tests for `pyvinecopulib.core.Vinedist`.

The claims worth pinning are the identities. That `logpdf` really is
`log c(F(x)) + sum_j log pdf_j(x_j)`, checked against a hand computation rather
than a tolerance. That a `Vinedist` with empirical margins is the same model as
a `Vinecop` fitted to `to_pseudo_obs`, which is what makes the long-standing
workflow a special case. And that the `(n, d + k)` layout a discrete margin
needs is assembled for the user, since assembling it by hand is what
`examples/04_discrete_variables.ipynb` currently has to do.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import Kde1d
from pyvinecopulib.margins import (
  EmpiricalMargin,
  MarginSelector,
)

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


# --- the empirical-margin identity ------------------------------------------ #


def test_empirical_margins_reproduce_a_copula_on_pseudo_obs(
  continuous: np.ndarray,
) -> None:
  """`Vinedist` with empirical margins is `Vinecop` on `to_pseudo_obs`.

  Same copula data, so the same structure and the same pair copulas. The joint
  density is not defined — an empirical margin has none — so `logpdf` raises
  rather than reporting a number on the original scale.
  """
  dist = pv.Vinedist.from_data(continuous, margins="empirical")
  u = np.asarray(pv.to_pseudo_obs(continuous))
  reference: Any = pv.Vinecop.from_data(u)

  assert dist.var_types == ["c", "c"]
  np.testing.assert_array_equal(
    np.asarray(dist.copula.structure.matrix),
    np.asarray(reference.structure.matrix),
  )
  np.testing.assert_array_equal(
    np.asarray(dist.copula.pdf(u)), np.asarray(reference.pdf(u))
  )

  # The identity at its source: the margins hand the copula exactly the
  # pseudo-observations. Ties would separate the two (`to_pseudo_obs` averages
  # ranks where the margin counts upward); continuous data has none.
  np.testing.assert_allclose(
    pv.Vinedist.copula_data(dist.margins, continuous), u, atol=1e-12
  )

  with pytest.raises(NotImplementedError, match="Kde1d"):
    dist.logpdf(continuous)


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

    def pdf(self, x: Any) -> Any:
      return np.full_like(np.asarray(x, dtype=float), 0.5)

    def cdf(self, x: Any) -> Any:
      return np.full_like(np.asarray(x, dtype=float), 0.3)

    def cdf_left(self, x: Any) -> Any:
      return np.full_like(np.asarray(x, dtype=float), 0.9)

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
    continuous, margins={"b": "empirical"}, names=["a", "b"]
  )
  assert isinstance(dist.margins[0], Kde1d)
  assert isinstance(dist.margins[1], EmpiricalMargin)


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

    def fit(self, x: Any, *, weights: Any = None) -> Any:
      return self

    def pdf(self, x: Any) -> Any:
      return np.ones_like(np.asarray(x, dtype=float))

    def cdf(self, x: Any) -> Any:
      return np.clip(np.asarray(x, dtype=float), 0.0, 1.0)

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
