"""Tests for the backend-neutral ``pyvinecopulib.core`` margin base.

Exercises :class:`pyvinecopulib.core.MarginBase` — the canonical,
array-backend-agnostic partial implementation of the ``MarginLike`` contract —
purely on NumPy, so it also confirms that the neutral ``core`` layer runs
without PyTorch. A conformance test pins that :class:`pyvinecopulib.core.Kde1d`
satisfies ``MarginLike`` directly: it is the library's default margin, and the
contract was named after its surface so that it needs no wrapper.
"""

from __future__ import annotations

import math
from typing import Any, Optional

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import MarginBase, MarginLike
from pyvinecopulib.core.margin_base import _reject_covariates


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

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    inside = y >= self.shift
    return np.where(
      inside, self.rate * np.exp(-self.rate * (y - self.shift)), 0.0
    )

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    inside = y >= self.shift
    return np.where(inside, 1.0 - np.exp(-self.rate * (y - self.shift)), 0.0)

  def exact_icdf(self, p: Any) -> Any:
    """The closed form, for comparison against the inherited bisection."""
    return self.shift - np.log1p(-p) / self.rate


class _Uniform01(MarginBase[np.ndarray]):
  """Bounded support, so ``icdf`` needs no bracket search at all."""

  @property
  def support(self) -> tuple[float, float]:
    return (0.0, 1.0)

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    return np.where((y >= 0.0) & (y <= 1.0), 1.0, 0.0)

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    return np.clip(y, 0.0, 1.0)


class _Geometricish(MarginBase[np.ndarray]):
  """A discrete margin on the non-negative integers (``p = 1/2``)."""

  @property
  def var_type(self) -> str:
    return "d"

  @property
  def support(self) -> tuple[float, float]:
    return (0.0, float("inf"))

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    return np.where(y >= 0.0, 0.5 ** (y + 1.0), 0.0)

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    return np.where(y >= 0.0, 1.0 - 0.5 ** (y + 1.0), 0.0)


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

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    body = (1.0 - self.prob0) * np.exp(-y)
    return np.where(y == 0.0, self.prob0, np.where(y > 0.0, body, 0.0))

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    body = (1.0 - self.prob0) * (1.0 - np.exp(-y))
    return np.where(y >= 0.0, self.prob0 + body, 0.0)


def test_kde1d_satisfies_the_margin_contract() -> None:
  """`Kde1d` is a `MarginLike` with no adapter, fitted or not."""
  kde = pv.core.Kde1d()
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
  with pytest.raises(NotImplementedError, match="_sample_uniform"):
    _ShiftedExp().sample(5)


def test_simulate_inverts_uniforms_once_the_hook_exists() -> None:
  """With the hook supplied, `sample` is `icdf` of uniforms."""

  class _Seeded(_ShiftedExp):
    def _sample_uniform(self, n: int, seeds: list[int]) -> Any:
      rng = np.random.default_rng(seeds[0] if seeds else 0)
      return rng.uniform(size=n)

  m = _Seeded(rate=2.0, shift=1.0)
  draws = m.sample(500, seeds=[42])
  assert draws.shape == (500,)
  assert np.all(draws >= 1.0)
  # The sample cdf of the draws should look uniform.
  u = m.cdf(draws)
  assert 0.4 < float(np.mean(u)) < 0.6


def test_repr_is_structural() -> None:
  """`repr` names the class and its variable type."""
  assert repr(_ShiftedExp()) == "_ShiftedExp(var_type='c')"
  assert repr(_Geometricish()) == "_Geometricish(var_type='d')"


# --- exogenous covariates ---------------------------------------------------- #


class _Recording(MarginBase[np.ndarray]):
  """Records whether covariates reached each primitive.

  A location-shift model, so `pdf` / `cdf` genuinely change with `x` and a
  forwarding bug shows up in the numbers as well as in the log.
  """

  supports_covariates = True

  def __init__(self) -> None:
    self.seen: list[str] = []

  def _shift(self, x: Optional[Any]) -> Any:
    return 0.0 if x is None else np.asarray(x, dtype=float)[:, 0]

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    self.seen.append("pdf" if x is not None else "pdf-bare")
    z = np.asarray(y, dtype=float) - self._shift(x)
    return np.exp(-0.5 * z**2) / np.sqrt(2.0 * np.pi)

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    self.seen.append("cdf" if x is not None else "cdf-bare")
    z = np.asarray(y, dtype=float) - self._shift(x)
    return 0.5 * (1.0 + np.vectorize(math.erf)(z / np.sqrt(2.0)))


def test_covariates_reach_every_derived_member() -> None:
  """A margin that declares them sees them through the inherited members."""
  y = np.array([0.0, 1.0, 2.0])
  cov = np.array([[0.5], [0.5], [0.5]])
  m = _Recording()

  np.testing.assert_allclose(m.logpdf(y, x=cov), np.log(m.pdf(y, x=cov)))
  np.testing.assert_allclose(m.cdf_left(y, x=cov), m.cdf(y, x=cov))
  np.testing.assert_allclose(m.loglik(y, x=cov), np.sum(m.logpdf(y, x=cov)))
  # `icdf` inverts the shifted cdf, so the covariate has to reach the closure.
  np.testing.assert_allclose(
    m.icdf(np.array([0.5, 0.5, 0.5]), x=cov), cov[:, 0], atol=1e-6
  )
  assert "pdf-bare" not in m.seen
  assert "cdf-bare" not in m.seen


def test_covariates_are_not_forwarded_to_an_unconditional_margin() -> None:
  """The gate omits `x` entirely, which is what keeps `pdf(self, y)` valid."""

  class _Unconditional(_Recording):
    supports_covariates = False

  m = _Unconditional()
  cov = np.array([[0.5], [0.5]])
  y = np.array([0.0, 1.0])
  # Same answer as without covariates, and the primitives were called bare.
  np.testing.assert_array_equal(m.logpdf(y, x=cov), m.logpdf(y))
  assert set(m.seen) == {"pdf-bare"}


def test_fit_refuses_covariates_it_cannot_read() -> None:
  """Silently fitting `f(y)` when `f(y | x)` was asked for is the bad outcome."""

  class _Fittable(MarginBase[np.ndarray]):
    def fit(
      self, y: Any, *, x: Optional[Any] = None, weights: Any = None
    ) -> Any:
      _reject_covariates(self, x)
      return self

    def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return np.ones_like(np.asarray(y, dtype=float))

    def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return np.clip(np.asarray(y, dtype=float), 0.0, 1.0)

  m = _Fittable()
  assert m.fit(np.array([0.5])) is m
  with pytest.raises(ValueError, match="supports_covariates"):
    m.fit(np.array([0.5]), x=np.array([[1.0]]))


def test_loglik_without_data_needs_a_fitted_value() -> None:
  """The no-argument form reports what `fit` attained, or says it has none."""
  m = _ShiftedExp()
  with pytest.raises(NotImplementedError, match="no fitted log-likelihood"):
    m.loglik()
  with pytest.raises(ValueError, match="only meaningful with data"):
    m.loglik(weights=np.array([1.0]))

  class _Stored(_ShiftedExp):
    @property
    def _fitted_loglik(self) -> float:
      return -12.5

  assert _Stored().loglik() == -12.5


def test_a_discrete_icdf_lands_on_its_own_lattice() -> None:
  """`icdf` must return values its own `cdf_left` accepts.

  Bisection converges *to* a jump rather than onto it, so an answer a few ulp
  below an integer used to come back — and `cdf_left`, which steps back one
  lattice point, then rejected it as non-integral.
  """
  m = _Geometricish()
  drawn = m.icdf(np.array([0.05, 0.3, 0.5, 0.9, 0.999]))
  np.testing.assert_array_equal(drawn, np.round(drawn))
  assert np.all(np.isfinite(m.cdf_left(drawn)))


def test_icdf_resolves_a_quantile_far_below_the_bracket() -> None:
  """Bisection halves an absolute interval, so the step count has to allow for
  an answer many orders of magnitude below the support's scale.

  With the h-inverses' 50 steps, distinct probabilities here all returned the
  same 4.44e-16 — a relative error that grows without bound as the target
  shrinks. The floor stays absolute, so this widens the resolvable range rather
  than guaranteeing exactness; `icdf` is documented as overridable for that.
  """

  class _Tiny(MarginBase[np.ndarray]):
    """`cdf(y) = y ** 0.05` on the unit interval: a very steep left tail."""

    @property
    def support(self) -> tuple[float, float]:
      return (0.0, 1.0)

    def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return 0.05 * np.asarray(y, dtype=float) ** -0.95

    def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
      return np.asarray(y, dtype=float) ** 0.05

  p = np.array([0.2, 0.1, 0.05])
  got = _Tiny().icdf(p)
  exact = p**20.0  # 1.0e-14 down to 9.5e-27, above the resolvable floor
  np.testing.assert_allclose(got, exact, rtol=1e-6)
  # The defect was a collapse, not mere imprecision: distinct probabilities
  # all came back as the same value.
  assert np.all(np.diff(got) < 0.0)


def test_declare_is_a_no_op_that_chains() -> None:
  """A margin fixed by construction ignores the caller's schema.

  `declare` exists so a caller with schema knowledge can hand it over before
  fitting; a margin whose type and support are settled needs nothing from it, so
  the base implementation must accept the call and change nothing.
  """
  m = _ShiftedExp(rate=2.0, shift=1.0)
  assert m.declare(var_type="d", support=(0.0, 10.0)) is m
  assert m.var_type == "c"
  assert m.support == (1.0, math.inf)
