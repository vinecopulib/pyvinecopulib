"""Tests for ``VineRegressor``'s quadrature over the response.

What this file pins is where the nodes come from. ``use_grid=True`` integrates
on the probability scale --- nodes at ``margin.icdf(p)``, weights the copula
density alone --- so the rule asks the response margin for nothing but its
inverse CDF: every continuous margin can be the response, the nodes cannot
leave its support, and a bounded margin is no less accurate than an unbounded
one. It also pins the arithmetic: the weights carry the node spacing and *not*
the marginal density, which is what makes the substitution correct.

The reference throughout is ``use_grid=False``, the exact weighted sum over the
training rows, which approximates the same integral by Monte Carlo and shares no
code with the quadrature.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

pytest.importorskip("sklearn")

import pyvinecopulib as pv  # noqa: E402
from pyvinecopulib.core import Kde1d  # noqa: E402
from pyvinecopulib.sklearn import VineRegressor  # noqa: E402

from .helpers import AtomicMargin  # noqa: E402

BETA = np.array([1.0, -0.7])


def _data(n: int, seed: int = 20260818, heavy: bool = False):
  """Covariates and a linear response, optionally heavy-tailed."""
  rng = np.random.default_rng(seed)
  X = rng.standard_normal((n, 2))
  noise = rng.standard_t(3, size=n) if heavy else rng.standard_normal(n)
  return X, X @ BETA + 0.5 * noise


def _bounded_kde(column: np.ndarray) -> Kde1d:
  """A margin bounded by the observed range, so its own grid is non-uniform."""
  column = np.asarray(column, dtype=float)
  return Kde1d(xmin=float(column.min()), xmax=float(column.max())).fit(column)


def _relative_rmse(a: np.ndarray, b: np.ndarray, scale: float) -> float:
  return float(np.sqrt(np.mean((a - b) ** 2)) / scale)


def test_the_weights_are_the_copula_density_times_the_node_spacing() -> None:
  """No marginal density factor: the substitution to ``p`` cancels it."""
  X, y = _data(300)
  est = VineRegressor(n_nodes=41).fit(X, y)
  dist = est.distribution_

  # The documented rule: p = Phi(z) on a uniform z grid spanning five standard
  # deviations either side of the median. Retuning that half-width moves this
  # expectation with it.
  z = np.linspace(-5.0, 5.0, 41)
  p = np.array([0.5 * (1.0 + math.erf(zi / math.sqrt(2.0))) for zi in z])
  nodes = np.asarray(dist.margins[0].icdf(p), dtype=float)
  spacing = np.exp(-0.5 * z**2)

  expected = []
  for row in X[:5]:
    ux = np.asarray(
      pv.Vinedist.copula_data(dist.margins[1:], row.reshape(1, -1))
    )
    u = np.column_stack([p, np.repeat(ux, p.size, axis=0)])
    w = np.asarray(dist.vinecop.pdf(u)) * spacing
    expected.append((w / w.sum()) @ nodes)

  np.testing.assert_allclose(
    est.predict(X[:5]), np.asarray(expected), rtol=1e-12, atol=1e-12
  )


def test_the_node_count_is_n_nodes() -> None:
  """``n_nodes`` sets the width of the weight matrix, whatever the sample."""
  X, y = _data(300)
  est = VineRegressor(n_nodes=17).fit(X, y)
  assert est._weights_for_batch(X[:4]).shape == (4, 17)


@pytest.mark.parametrize("heavy", [False, True])
def test_the_quadrature_converges(heavy: bool) -> None:
  """Refining the grid stops moving the prediction well before the default."""
  X, y = _data(400, heavy=heavy)
  scale = float(np.std(y))
  fine = VineRegressor(n_nodes=4001).fit(X, y).predict(X[:25])
  default = VineRegressor().fit(X, y).predict(X[:25])
  coarse = VineRegressor(n_nodes=51).fit(X, y).predict(X[:25])

  assert _relative_rmse(default, fine, scale) < 1e-3
  # A coarse grid is visibly worse, so the default is not accidentally exact.
  assert _relative_rmse(coarse, fine, scale) > _relative_rmse(
    default, fine, scale
  )


def test_it_integrates_the_inverse_cdf_against_the_copula() -> None:
  """The rule computes the integral it claims to, on an independent grid.

  The expectation here is a midpoint rule uniform in ``p`` -- a different node
  set and no probit substitution -- written out from the fitted distribution
  alone, so agreeing with it says the estimator integrates
  :math:`F_Y^{-1}(p)\\, c(p, u_x)` and not something proportional to it.
  """
  X, y = _data(400)
  scale = float(np.std(y))
  est = VineRegressor(n_nodes=4001).fit(X, y)
  dist = est.distribution_

  n_reference = 8001
  p = (np.arange(n_reference) + 0.5) / n_reference
  nodes = np.asarray(dist.margins[0].icdf(p), dtype=float)
  expected = []
  for row in X[:8]:
    ux = np.asarray(
      pv.Vinedist.copula_data(dist.margins[1:], row.reshape(1, -1))
    )
    u = np.column_stack([p, np.repeat(ux, n_reference, axis=0)])
    w = np.asarray(dist.vinecop.pdf(u))
    expected.append((w / w.sum()) @ nodes)

  assert _relative_rmse(est.predict(X[:8]), np.asarray(expected), scale) < 1e-4


def test_the_quadrature_agrees_with_the_exact_weighted_sum() -> None:
  """Both paths estimate the same integral, so they must nearly agree."""
  X, y = _data(2000)
  scale = float(np.std(y))
  grid = VineRegressor().fit(X, y).predict(X[:25])
  exact = VineRegressor(use_grid=False).fit(X, y).predict(X[:25])
  # The bound is the Monte-Carlo error of the weighted sum, whose effective
  # sample size is a fraction of `n`; the quadrature's own error is three
  # orders smaller (`test_it_integrates_the_inverse_cdf_against_the_copula`).
  assert _relative_rmse(grid, exact, scale) < 0.05


def test_the_grid_reaches_past_the_observed_responses() -> None:
  """A heavy-tailed response makes the two paths differ, by design.

  The quadrature integrates over the whole probability scale, so it sees the
  response margin wherever the margin extrapolates; the weighted sum over
  training rows cannot go past the largest observation. The two therefore
  disagree by more than Monte-Carlo error on heavy tails -- and it is not
  quadrature error, which is what refining the grid here shows.
  """
  X, y = _data(2000, heavy=True)
  scale = float(np.std(y))
  default = VineRegressor().fit(X, y).predict(X[:25])
  fine = VineRegressor(n_nodes=4001).fit(X, y).predict(X[:25])
  exact = VineRegressor(use_grid=False).fit(X, y).predict(X[:25])

  assert _relative_rmse(default, fine, scale) < 1e-3
  assert _relative_rmse(default, exact, scale) > 0.01


def test_a_bounded_response_margin_is_as_accurate_as_an_unbounded_one() -> None:
  """Bounding the response leaves the quadrature where it was.

  A bounded ``Kde1d`` grids its own support non-uniformly, which is what
  made reading nodes and a density off that grid a biased rule.
  """
  X, y = _data(1500)
  scale = float(np.std(y))
  exact = (
    VineRegressor(margins=_bounded_kde, use_grid=False)
    .fit(X, y)
    .predict(X[:25])
  )
  bounded = VineRegressor(margins=_bounded_kde).fit(X, y).predict(X[:25])
  assert _relative_rmse(bounded, exact, scale) < 0.05


def test_the_nodes_stay_inside_the_response_support() -> None:
  """Predictions are convex combinations of points the margin can produce."""
  X, y = _data(400)
  est = VineRegressor(quantiles=[0.05, 0.95], margins=_bounded_kde).fit(X, y)
  pred = est.predict(X[:50])
  assert pred.min() >= y.min()
  assert pred.max() <= y.max()


@pytest.mark.parametrize("spec", ["kde", "parametric"])
def test_any_continuous_response_margin_works_on_the_grid(spec: str) -> None:
  """The rule needs an inverse CDF and nothing else, so every margin serves."""
  if spec == "parametric":
    pytest.importorskip("scipy")
  X, y = _data(1200)
  scale = float(np.std(y))
  grid = VineRegressor(margins=spec).fit(X, y).predict(X[:25])
  exact = VineRegressor(margins=spec, use_grid=False).fit(X, y).predict(X[:25])
  assert np.all(np.isfinite(grid))
  # Correlation rather than a distance: a margin that extrapolates past the
  # data legitimately moves the conditional mean away from the weighted sum
  # (`test_the_grid_reaches_past_the_observed_responses`).
  assert np.corrcoef(grid, exact)[0, 1] > 0.99
  assert _relative_rmse(grid, exact, scale) < 0.2


def test_a_step_response_margin_degenerates_to_its_atoms() -> None:
  """A step inverse CDF returns order statistics, so the quadrature becomes a
  weighted sum over the observed responses."""
  X, y = _data(300)
  est = VineRegressor(margins=AtomicMargin(), n_nodes=201).fit(X, y)
  nodes = np.asarray(est.distribution_.margins[0].icdf(np.linspace(0.01, 0.99)))
  assert np.all(np.isin(np.round(nodes, 12), np.round(y, 12)))
  assert np.all(np.isfinite(est.predict(X[:10])))
