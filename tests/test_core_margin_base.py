"""Tests for the backend-neutral ``pyvinecopulib.core`` margin base.

Exercises :class:`pyvinecopulib.core.MarginBase` — the canonical,
array-backend-agnostic partial implementation of the ``MarginLike`` contract —
purely on NumPy, so it also confirms that the neutral ``core`` layer runs
without PyTorch. A conformance test pins that :class:`pyvinecopulib.utils.Kde1d`
satisfies ``MarginLike`` directly: it is the library's default margin, and the
contract was named after its surface so that it needs no wrapper.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import MarginBase, MarginLike


class _ShiftedExp(MarginBase[np.ndarray]):
  """Shifted exponential, implementing only the abstract surface.

  ``icdf`` / ``logpdf`` / ``cdf_left`` / ``loglik`` therefore come from
  :class:`MarginBase` and are what these tests exercise. The support is
  half-bounded, so ``icdf`` has to widen the bracket outward.
  """

  def __init__(self, rate: float = 1.0, shift: float = 0.0) -> None:
    self.rate = rate
    self.shift = shift

  @property
  def support(self) -> tuple[float, float]:
    return (self.shift, float("inf"))

  def pdf(self, x: Any) -> Any:
    inside = x >= self.shift
    return np.where(
      inside, self.rate * np.exp(-self.rate * (x - self.shift)), 0.0
    )

  def cdf(self, x: Any) -> Any:
    inside = x >= self.shift
    return np.where(inside, 1.0 - np.exp(-self.rate * (x - self.shift)), 0.0)

  def exact_icdf(self, p: Any) -> Any:
    """The closed form, for comparison against the inherited bisection."""
    return self.shift - np.log1p(-p) / self.rate


class _Uniform01(MarginBase[np.ndarray]):
  """Bounded support, so ``icdf`` needs no bracket search at all."""

  @property
  def support(self) -> tuple[float, float]:
    return (0.0, 1.0)

  def pdf(self, x: Any) -> Any:
    return np.where((x >= 0.0) & (x <= 1.0), 1.0, 0.0)

  def cdf(self, x: Any) -> Any:
    return np.clip(x, 0.0, 1.0)


class _Geometricish(MarginBase[np.ndarray]):
  """A discrete margin on the non-negative integers (``p = 1/2``)."""

  @property
  def var_type(self) -> str:
    return "d"

  @property
  def support(self) -> tuple[float, float]:
    return (0.0, float("inf"))

  def pdf(self, x: Any) -> Any:
    return np.where(x >= 0.0, 0.5 ** (x + 1.0), 0.0)

  def cdf(self, x: Any) -> Any:
    return np.where(x >= 0.0, 1.0 - 0.5 ** (x + 1.0), 0.0)


class _ZeroInflated(MarginBase[np.ndarray]):
  """An atom at zero plus an exponential body — the mixed case."""

  def __init__(self, prob0: float = 0.3) -> None:
    self.prob0 = prob0

  @property
  def var_type(self) -> str:
    return "zi"

  @property
  def support(self) -> tuple[float, float]:
    return (0.0, float("inf"))

  def pdf(self, x: Any) -> Any:
    body = (1.0 - self.prob0) * np.exp(-x)
    return np.where(x == 0.0, self.prob0, np.where(x > 0.0, body, 0.0))

  def cdf(self, x: Any) -> Any:
    body = (1.0 - self.prob0) * (1.0 - np.exp(-x))
    return np.where(x >= 0.0, self.prob0 + body, 0.0)


def test_kde1d_satisfies_the_margin_contract() -> None:
  """`Kde1d` is a `MarginLike` with no adapter, fitted or not."""
  kde = pv.utils.Kde1d()
  assert isinstance(kde, MarginLike)
  kde.fit(np.random.default_rng(0).normal(size=300))
  assert isinstance(kde, MarginLike)


def test_subclass_satisfies_the_contract() -> None:
  """Two methods are enough to conform."""
  assert isinstance(_ShiftedExp(), MarginLike)


@pytest.mark.parametrize("rate,shift", [(1.0, 0.0), (2.0, 1.0), (0.5, -3.0)])
def test_icdf_matches_the_closed_form(rate: float, shift: float) -> None:
  """The inherited bisection inverts a half-bounded cdf to machine precision."""
  m = _ShiftedExp(rate=rate, shift=shift)
  p = np.array([1e-6, 0.01, 0.25, 0.5, 0.9, 0.999])
  np.testing.assert_allclose(m.icdf(p), m.exact_icdf(p), atol=1e-9)


def test_icdf_on_a_bounded_support() -> None:
  """A finite bracket needs no widening and stays inside the support."""
  m = _Uniform01()
  p = np.linspace(0.0, 1.0, 11)
  got = m.icdf(p)
  np.testing.assert_allclose(got, p, atol=1e-12)
  assert np.all(got >= 0.0) and np.all(got <= 1.0)


def test_icdf_roundtrips_with_cdf() -> None:
  """`cdf(icdf(p)) == p` on the continuous margin."""
  m = _ShiftedExp(rate=1.5, shift=2.0)
  p = np.linspace(0.01, 0.99, 25)
  np.testing.assert_allclose(m.cdf(m.icdf(p)), p, atol=1e-9)


def test_logpdf_is_minus_inf_outside_the_support() -> None:
  """A vanishing density gives `-inf`, without a numpy warning."""
  m = _ShiftedExp(shift=1.0)
  with np.errstate(divide="raise", invalid="raise"):
    got = m.logpdf(np.array([0.0, 0.5, 2.0]))
  assert got[0] == -np.inf and got[1] == -np.inf
  assert np.isfinite(got[2])


def test_logpdf_matches_log_of_pdf_where_positive() -> None:
  """No floor is applied where the density is positive."""
  m = _ShiftedExp(rate=2.0)
  x = np.array([0.1, 1.0, 5.0])
  np.testing.assert_allclose(m.logpdf(x), np.log(m.pdf(x)), atol=0.0)


def test_loglik_is_a_zero_d_array_and_weightable() -> None:
  """`loglik` sums log-densities, optionally weighted."""
  m = _ShiftedExp()
  x = np.array([0.5, 1.0, 2.0])
  total = m.loglik(x)
  assert np.ndim(total) == 0
  np.testing.assert_allclose(total, np.sum(m.logpdf(x)))
  w = np.array([1.0, 2.0, 0.0])
  np.testing.assert_allclose(m.loglik(x, weights=w), np.sum(m.logpdf(x) * w))
  # Weighting by 2 equals duplicating the observation.
  dup = np.array([0.5, 1.0, 1.0])
  np.testing.assert_allclose(
    m.loglik(dup), m.loglik(np.array([0.5, 1.0]), weights=np.array([1.0, 2.0]))
  )


def test_cdf_left_defaults_to_cdf_when_continuous() -> None:
  """A continuous margin has no atoms, so the left limit coincides."""
  m = _ShiftedExp()
  x = np.array([0.5, 1.0, 2.0])
  assert np.array_equal(m.cdf_left(x), m.cdf(x))
  assert m.var_type == "c"


def test_cdf_left_steps_back_a_lattice_point_when_discrete() -> None:
  """`F(k^-)` equals both `F(k-1)` and `F(k) - pmf(k)`."""
  m = _Geometricish()
  k = np.arange(0.0, 8.0)
  left = m.cdf_left(k)
  np.testing.assert_allclose(left, m.cdf(k - 1.0), atol=0.0)
  np.testing.assert_allclose(left, m.cdf(k) - m.pdf(k), atol=1e-15)
  assert left[0] == 0.0
  assert np.all(left <= m.cdf(k))


def test_cdf_left_rejects_non_integer_discrete_input() -> None:
  """The default lattice rule is guarded rather than silently wrong."""
  with pytest.raises(ValueError, match="integer-valued"):
    _Geometricish().cdf_left(np.array([1.5]))


def test_cdf_left_removes_the_atom_when_zero_inflated() -> None:
  """`F(0) - F(0^-)` recovers exactly the point mass."""
  m = _ZeroInflated(prob0=0.3)
  x = np.array([0.0, 1.0, 2.0])
  left = m.cdf_left(x)
  np.testing.assert_allclose(m.cdf(x)[0] - left[0], 0.3, atol=1e-12)
  # Away from the atom the margin is continuous, so nothing is removed.
  np.testing.assert_allclose(left[1:], m.cdf(x)[1:], atol=0.0)


def test_fit_raises_by_default_and_is_fitted_is_true() -> None:
  """A fully-specified margin is already fitted and needs no estimator."""
  m = _ShiftedExp()
  assert m.is_fitted is True
  assert m.supports_weights is False
  with pytest.raises(NotImplementedError, match="fit is not defined"):
    m.fit(np.array([1.0, 2.0]))


def test_simulate_requires_the_draw_hook() -> None:
  """Sampling needs an array namespace, so the hook raises until overridden."""
  with pytest.raises(NotImplementedError, match="_simulate_uniform"):
    _ShiftedExp().simulate(5)


def test_simulate_inverts_uniforms_once_the_hook_exists() -> None:
  """With the hook supplied, `simulate` is `icdf` of uniforms."""

  class _Seeded(_ShiftedExp):
    def _simulate_uniform(self, n: int, seeds: list[int]) -> Any:
      rng = np.random.default_rng(seeds[0] if seeds else 0)
      return rng.uniform(size=n)

  m = _Seeded(rate=2.0, shift=1.0)
  draws = m.simulate(500, seeds=[42])
  assert draws.shape == (500,)
  assert np.all(draws >= 1.0)
  # The sample cdf of the draws should look uniform.
  u = m.cdf(draws)
  assert 0.4 < float(np.mean(u)) < 0.6


def test_repr_is_structural() -> None:
  """`repr` names the class and its variable type."""
  assert repr(_ShiftedExp()) == "_ShiftedExp(var_type='c')"
  assert repr(_Geometricish()) == "_Geometricish(var_type='d')"
