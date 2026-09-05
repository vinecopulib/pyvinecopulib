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
  SciPyMargin,
  as_margin,
  register_margin_adapter,
  resolve_margins,
)

scipy_stats = pytest.importorskip("scipy.stats")

# `scipy.stats.Binomial` is the new-API *discrete* class, added in SciPy 1.17,
# one release above the `[scipy]` extra's floor -- 1.16 ships `Normal` but no
# discrete counterpart. Skips are applied per parameter, so the legacy half of
# each parametrized test stays live.
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


def test_resolve_margins_defers_a_default_no_variable_needs() -> None:
  """A specification naming every variable must not build the default at all.

  The `default` parameter documents this, and it matters because building one
  can legitimately raise -- the sklearn estimators pass a callable that reads
  the variable types off the data, and `Kde1d` refuses a categorical whose
  levels are not integers. Every branch has to honor it, mapping included.
  """
  calls = {"n": 0}

  def default() -> Any:
    calls["n"] += 1
    raise AssertionError("built a default no variable needed")

  assert len(resolve_margins({0: Kde1d(), 1: Kde1d()}, 2, default=default)) == 2
  assert len(resolve_margins([Kde1d(), Kde1d()], 2, default=default)) == 2
  assert len(resolve_margins(Kde1d(), 2, default=default)) == 2
  assert len(resolve_margins("kde", 2, default=default)) == 2
  assert calls["n"] == 0

  # And it *is* built when a variable is genuinely left over.
  resolved = resolve_margins(
    {0: Kde1d()}, 2, default=lambda: [Kde1d(), Kde1d(type="discrete")]
  )
  assert resolved[1].type == "discrete"


def test_resolve_margins_checks_the_default_length() -> None:
  """A default is per variable, so its length is checked like a sequence's."""
  with pytest.raises(ValueError, match="default has length 1"):
    resolve_margins(None, 2, default=[Kde1d()])


@pytest.mark.parametrize("key", [0.9, 1.0, np.float64(1.0)])
def test_resolve_margins_rejects_noninteger_mapping_keys(key: Any) -> None:
  """A numeric key must be an integer position, never silently truncated."""
  with pytest.raises(ValueError, match="integer position"):
    resolve_margins({key: Kde1d()}, 2)

  resolved = resolve_margins({np.int64(1): Kde1d(type="discrete")}, 2)
  assert resolved[1].type == "discrete"


def test_callable_margin_specifications_receive_weights() -> None:
  """A callable owns its fitting, so observation weights must reach it."""
  from pyvinecopulib.core import Vinedist

  seen: list[dict[str, Any]] = []

  def make_margin(y: Any, **kwargs: Any) -> Any:
    seen.append(kwargs)
    return SciPyMargin("norm", (float(np.mean(y)), 1.0))

  rng = np.random.default_rng(8)
  data = rng.normal(size=(40, 2))
  weights = np.linspace(1.0, 2.0, data.shape[0])
  Vinedist.from_data(data, margins=make_margin, weights=weights)

  assert len(seen) == data.shape[1]
  for kwargs in seen:
    np.testing.assert_array_equal(kwargs["weights"], weights)


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


def test_as_margin_accepts_a_structural_margin() -> None:
  """The documented `MarginLike` seam does not require a library base class."""

  class StructuralMargin:
    def pdf(self, y: Any, /, *, x: Any = None) -> Any:
      return np.exp(-(np.asarray(y) ** 2) / 2) / np.sqrt(2 * np.pi)

    def cdf(self, y: Any, /, *, x: Any = None) -> Any:
      return scipy_stats.norm.cdf(y)

    def icdf(self, p: Any, /, *, x: Any = None) -> Any:
      return scipy_stats.norm.ppf(p)

  margin = StructuralMargin()
  assert isinstance(margin, MarginLike)
  assert as_margin(margin) is margin


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


def test_as_margin_scipy_discrete_accepts_a_shifted_lattice() -> None:
  """The legacy adapter derives the left limit without assuming integers."""
  raw = scipy_stats.poisson(3.0, loc=0.5)
  margin: Any = as_margin(raw)
  points = np.array([0.5, 1.0, 1.5, 2.25])
  expected = raw.cdf(points) - raw.pmf(points)
  np.testing.assert_allclose(margin.cdf_left(points), expected, atol=0)
  assert margin.cdf_left(np.array([0.5]))[0] == pytest.approx(0.0, abs=1e-15)


@pytest.mark.parametrize(
  "raw",
  [
    pytest.param(scipy_stats.norm(), id="legacy"),
    pytest.param(scipy_stats.Normal(), id="modern"),
  ],
)
def test_as_margin_scipy_forwards_native_log_density(raw: Any) -> None:
  """A finite native tail log-density must not underflow through `log(pdf)`."""
  far = np.array([40.0])
  margin: Any = as_margin(raw)
  np.testing.assert_allclose(margin.logpdf(far), raw.logpdf(far), atol=0)
  assert np.isfinite(margin.logpdf(far)).all()
  assert margin.pdf(far)[0] == 0.0


def test_as_margin_rejects_the_unknown() -> None:
  """An unrecognized object names both escape hatches."""
  with pytest.raises(TypeError, match="MarginBase"):
    as_margin(object())


def test_register_margin_adapter_takes_precedence() -> None:
  """A later registration wins, so a user can override a built-in."""

  class _Sentinel:
    pass

  sentinel = _Sentinel()
  target = SciPyMargin("norm", (0.0, 1.0))
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


def test_as_margin_torch_forwards_native_log_density() -> None:
  """Torch's stable `log_prob` remains available in far tails."""
  torch = pytest.importorskip("torch")
  raw = torch.distributions.Normal(
    torch.tensor(0.0, dtype=torch.float64),
    torch.tensor(1.0, dtype=torch.float64),
  )
  far = torch.tensor([40.0], dtype=torch.float64)
  margin = as_margin(raw)
  margin_with_logpdf: Any = margin
  torch.testing.assert_close(margin_with_logpdf.logpdf(far), raw.log_prob(far))
  assert torch.isfinite(margin_with_logpdf.logpdf(far)).all()
  assert margin.pdf(far).item() == 0.0


@pytest.mark.parametrize("family", ["Poisson", "Bernoulli"])
def test_as_margin_rejects_discrete_torch_distributions(family: str) -> None:
  """Torch margins with atoms lack the cdf-left contract vines require."""
  torch = pytest.importorskip("torch")
  raw = getattr(torch.distributions, family)(torch.tensor(0.4))
  with pytest.raises(TypeError, match="discrete torch distribution.*Kde1d"):
    as_margin(raw)


@pytest.mark.parametrize("family", ["Beta", "StudentT"])
def test_as_margin_rejects_torch_families_without_a_cdf(family: str) -> None:
  """A recognized ecosystem object is refused until it can meet the contract."""
  torch = pytest.importorskip("torch")
  args = (2.0, 3.0) if family == "Beta" else (5.0,)
  raw = getattr(torch.distributions, family)(*args)
  with pytest.raises(TypeError, match="cdf is not implemented.*MarginBase"):
    as_margin(raw)


@pytest.mark.parametrize("family", ["Normal", "Gamma", "LogNormal"])
def test_as_margin_supports_torch_families_with_a_cdf(family: str) -> None:
  """Documented continuous Torch families round-trip through their cdf."""
  torch = pytest.importorskip("torch")
  args = (0.0, 1.0) if family != "Gamma" else (2.0, 1.0)
  raw = getattr(torch.distributions, family)(
    *(torch.tensor(v, dtype=torch.float64) for v in args)
  )
  margin = as_margin(raw)
  p = torch.tensor([0.25, 0.5, 0.75], dtype=torch.float64)
  torch.testing.assert_close(margin.cdf(margin.icdf(p)), p, atol=1e-6, rtol=0)
