"""Tests for `pyvinecopulib.margins`.

Three things are pinned here. That `to_pseudo_obs` is scale-invariant in its
weights, which is what keeps a weighted fit on the same scale as an
unweighted one. That `Kde1d` re-selects its bandwidth on every fit, which a
bare `Kde1d` does not. And that `as_margin`
coerces each supported ecosystem to the *right numbers* -- especially for
discrete distributions, where satisfying the contract's member names is not the
same as satisfying its semantics.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import MarginLike
from pyvinecopulib.core import Kde1d
from pyvinecopulib.margins import (
  ParametricMargin,
  as_margin,
  register_margin_adapter,
  resolve_margins,
)

scipy_stats = pytest.importorskip("scipy.stats")

# `scipy.stats.Binomial` is the new-API *discrete* class, added in SciPy 1.17.
# SciPy 1.17 requires Python >= 3.11, so a 3.10 environment resolves 1.15 or
# 1.16, which ship `Normal` but no discrete counterpart. Skips are applied per
# parameter, so the legacy half of each parametrized test stays live.
needs_scipy_discrete = pytest.mark.skipif(
  not hasattr(scipy_stats, "Binomial"),
  reason="scipy.stats.Binomial requires SciPy >= 1.17",
)


@pytest.fixture
def sample() -> np.ndarray:
  """A tie-free positive sample."""
  return np.random.default_rng(0).gamma(2.0, 1.5, size=400)


# --- to_pseudo_obs ---------------------------------------------------------- #


def test_to_pseudo_obs_is_scale_invariant_in_the_weights() -> None:
  """Only weight *ratios* matter, so rescaling them changes nothing.

  Weights summing to 1 must not shrink every value toward `1 / (n + 1)`. The
  transform feeds the copula data of every weighted fit, so this property is
  what keeps a weighted fit on the same scale as an unweighted one.
  """
  y = np.array([1.0, 2.0, 3.0]).reshape(-1, 1)
  w = np.array([1.0, 2.0, 1.0])
  reference = np.asarray(pv.to_pseudo_obs(y, weights=w))

  for scale in (0.1, 1.0, 1000.0):
    np.testing.assert_allclose(
      np.asarray(pv.to_pseudo_obs(y, weights=w * scale)), reference, atol=1e-12
    )

  # Uniform weights of any scale reproduce the unweighted transform.
  np.testing.assert_allclose(
    np.asarray(pv.to_pseudo_obs(y, weights=np.full(3, 0.25))),
    np.asarray(pv.to_pseudo_obs(y)),
    atol=1e-12,
  )


# --- Kde1d ------------------------------------------------------------ #


def test_kde1d_margin_conforms_and_round_trips(sample: np.ndarray) -> None:
  """A fitted KDE margin inverts its own cdf."""
  m = Kde1d(xmin=0.0).fit(sample)
  assert isinstance(m, MarginLike)
  assert m.is_fitted and m.var_type == "c" and m.support == (0.0, float("inf"))
  p = np.linspace(0.05, 0.95, 19)
  np.testing.assert_allclose(m.cdf(m.icdf(p)), p, atol=1e-6)


def test_kde1d_margin_refit_reselects_the_bandwidth(
  sample: np.ndarray,
) -> None:
  """Refitting re-selects the bandwidth instead of reusing the previous fit's.

  Before kde1d-cpp#28 an estimator fed its last selected bandwidth back into the
  selector, so a prototype copied across columns kept the first column's scale.
  """
  m = Kde1d()
  m.fit(sample)
  narrow = m.bandwidth
  m.fit(sample * 10.0)
  wide = m.bandwidth
  assert wide / narrow == pytest.approx(10.0, rel=0.05)


def test_kde1d_margin_supports_weights(sample: np.ndarray) -> None:
  """Weights reach the underlying estimator and change the fit."""
  assert Kde1d.supports_weights is True
  w = np.where(sample > np.median(sample), 3.0, 1.0)
  plain = Kde1d().fit(sample)
  tilted = Kde1d().fit(sample, weights=w)
  assert tilted.icdf(np.array([0.5]))[0] > plain.icdf(np.array([0.5]))[0]


def test_kde1d_margin_discrete_left_limit() -> None:
  """The inherited `cdf_left` steps back a lattice point for counts."""
  counts = np.random.default_rng(1).poisson(3.0, size=500).astype(float)
  m = Kde1d(type="discrete", xmin=0.0, xmax=15.0).fit(counts)
  assert m.var_type == "d"
  k = np.arange(0.0, 6.0)
  np.testing.assert_allclose(m.cdf_left(k), m.cdf(k - 1.0), atol=0.0)
  assert np.all(m.cdf(k) - m.cdf_left(k) >= 0.0)


def test_kde1d_margin_zero_inflated_left_limit() -> None:
  """For a zero-inflated margin the jump at 0 is the point mass."""
  rng = np.random.default_rng(2)
  data = np.where(rng.uniform(size=600) < 0.3, 0.0, rng.exponential(size=600))
  m = Kde1d(type="zero-inflated", xmin=0.0).fit(data)
  assert m.var_type == "zi"
  jump = m.cdf(np.array([0.0])) - m.cdf_left(np.array([0.0]))
  np.testing.assert_allclose(jump, m.prob0, atol=1e-12)


def test_kde1d_margin_rejects_unknown_type() -> None:
  """The type is validated at construction, where the mistake is."""
  with pytest.raises(ValueError, match="variable type .* unknown"):
    Kde1d(type="ordinal")


def test_kde1d_margin_raises_before_fit() -> None:
  """`kde1d` is unavailable until fitted."""
  m = Kde1d()
  assert m.is_fitted is False
  with pytest.raises(RuntimeError, match="must first fit"):
    m.pdf(np.array([0.0]))


def test_kde1d_is_passed_through_by_as_margin(sample: np.ndarray) -> None:
  """A `Kde1d` needs no adapter: it satisfies the contract already."""
  kde = Kde1d(xmin=0.0).fit(sample)
  assert as_margin(kde) is kde
  assert as_margin(as_margin(kde)) is kde


# --- resolve_margins -------------------------------------------------------- #


def test_resolve_margins_falls_back_to_the_given_default() -> None:
  """An unaddressed variable takes the caller's default, not the library's."""
  default = [Kde1d(type="discrete"), Kde1d(type="zero-inflated")]
  resolved = resolve_margins({0: Kde1d()}, 2, default=default)
  assert resolved[0].type == "continuous"
  assert resolved[1].type == "zero-inflated"
  assert resolve_margins(None, 2, default=default)[0].type == "discrete"


def test_resolve_margins_checks_the_default_length() -> None:
  """A default is per variable, so its length is checked like a sequence's."""
  with pytest.raises(ValueError, match="default has length 1"):
    resolve_margins(None, 2, default=[Kde1d()])


# --- as_margin -------------------------------------------------------------- #


def test_as_margin_is_idempotent(sample: np.ndarray) -> None:
  """Anything this library made is returned unchanged."""
  m = Kde1d().fit(sample)
  assert as_margin(m) is m


def test_as_margin_adopts_a_raw_kde1d(sample: np.ndarray) -> None:
  """A bare fitted `Kde1d` becomes a `Kde1d`."""
  kde = pv.core.Kde1d()
  kde.fit(sample)
  assert isinstance(as_margin(kde), Kde1d)


@pytest.mark.parametrize(
  "factory,var_type",
  [
    (lambda: scipy_stats.Normal(mu=0.0, sigma=1.0), "c"),
    (lambda: scipy_stats.gamma(2.0, scale=1.5), "c"),
    (lambda: scipy_stats.poisson(3.0), "d"),
    pytest.param(
      lambda: scipy_stats.Binomial(n=10, p=0.3),
      "d",
      marks=needs_scipy_discrete,
    ),
  ],
)
def test_as_margin_coerces_scipy(factory, var_type: str) -> None:
  """Both SciPy generations coerce, with the right variable type."""
  m: Any = as_margin(factory())
  assert isinstance(m, MarginLike)
  assert m.var_type == var_type


@pytest.mark.parametrize(
  "factory,ref",
  [
    (
      lambda: scipy_stats.poisson(3.0),
      lambda k: scipy_stats.poisson.pmf(k, 3.0),
    ),
    pytest.param(
      lambda: scipy_stats.Binomial(n=10, p=0.3),
      lambda k: scipy_stats.binom.pmf(k, 10, 0.3),
      marks=needs_scipy_discrete,
    ),
  ],
)
def test_as_margin_discrete_pdf_is_the_mass_not_the_lebesgue_density(
  factory, ref
) -> None:
  """The trap: a modern SciPy discrete `pdf` is `+inf` at every atom.

  Coercion must route to `pmf`, and the jump `F(k) - F(k^-)` must equal it.
  """
  m: Any = as_margin(factory())
  k = np.arange(0.0, 6.0)
  np.testing.assert_allclose(m.pdf(k), ref(k), rtol=1e-10)
  assert np.all(np.isfinite(m.pdf(k)))
  np.testing.assert_allclose(m.cdf(k) - m.cdf_left(k), m.pdf(k), atol=1e-12)


@needs_scipy_discrete
def test_as_margin_scipy_new_discrete_pdf_would_have_been_inf() -> None:
  """Pin the hazard itself, so the reason for the wrapper stays visible."""
  raw = scipy_stats.Binomial(n=10, p=0.3)
  assert np.all(np.isinf(raw.pdf(np.arange(0.0, 4.0))))
  assert np.all(np.isfinite(as_margin(raw).pdf(np.arange(0.0, 4.0))))


def test_as_margin_scipy_icdf_matches_the_source() -> None:
  """`ppf` / `icdf` are wired to the contract's `icdf`."""
  p = np.array([0.1, 0.5, 0.9])
  legacy = as_margin(scipy_stats.gamma(2.0, scale=1.5))
  np.testing.assert_allclose(
    legacy.icdf(p), scipy_stats.gamma(2.0, scale=1.5).ppf(p), rtol=1e-12
  )
  modern = as_margin(scipy_stats.Normal(mu=0.0, sigma=1.0))
  np.testing.assert_allclose(
    modern.icdf(p), scipy_stats.norm.ppf(p), atol=1e-10
  )


def test_as_margin_rejects_the_unknown() -> None:
  """An unrecognized object names both escape hatches."""
  with pytest.raises(TypeError, match="MarginBase"):
    as_margin(object())


def test_register_margin_adapter_takes_precedence() -> None:
  """A later registration wins, so a user can override a built-in."""

  class _Sentinel:
    pass

  sentinel = _Sentinel()
  target = ParametricMargin("norm", (0.0, 1.0))
  register_margin_adapter(lambda o: isinstance(o, _Sentinel), lambda o: target)
  assert as_margin(sentinel) is target


# --- torch ------------------------------------------------------------------ #


def test_as_margin_coerces_torch_distributions() -> None:
  """`log_prob` becomes `pdf`; `support` becomes a plain pair."""
  torch = pytest.importorskip("torch")
  d = torch.distributions.Normal(0.0, 1.0)
  m: Any = as_margin(d)
  x = torch.tensor([-1.0, 0.0, 1.0], dtype=torch.float64)
  np.testing.assert_allclose(
    m.pdf(x).numpy(), scipy_stats.norm.pdf(x.numpy()), rtol=1e-6
  )
  np.testing.assert_allclose(
    m.cdf(x).numpy(), scipy_stats.norm.cdf(x.numpy()), rtol=1e-6
  )
  assert m.support == (float("-inf"), float("inf"))


def test_as_margin_torch_icdf_falls_back_to_bisection() -> None:
  """`Gamma` implements `cdf` but not `icdf`; the wrapper inverts numerically."""
  torch = pytest.importorskip("torch")
  raw = torch.distributions.Gamma(
    torch.tensor(2.0, dtype=torch.float64),
    torch.tensor(1.0, dtype=torch.float64),
  )
  with pytest.raises(NotImplementedError):
    raw.icdf(torch.tensor([0.5], dtype=torch.float64))
  m = as_margin(raw)
  p = torch.tensor([0.25, 0.5, 0.75], dtype=torch.float64)
  np.testing.assert_allclose(
    m.icdf(p).numpy(), scipy_stats.gamma(2.0).ppf(p.numpy()), atol=1e-6
  )


def test_as_margin_torch_autograd_survives() -> None:
  """Gradients flow through the wrapped density, which is the point of torch."""
  torch = pytest.importorskip("torch")
  mu = torch.tensor(0.5, dtype=torch.float64, requires_grad=True)
  m = as_margin(torch.distributions.Normal(mu, 1.0))
  m.pdf(torch.tensor([0.2, 0.7], dtype=torch.float64)).sum().backward()
  assert mu.grad is not None and torch.isfinite(mu.grad)
