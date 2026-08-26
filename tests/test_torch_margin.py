"""Tests for `pyvinecopulib.torch.TorchMargin`.

Skipped without PyTorch. The claims worth pinning are the ones that motivated
the class rather than the arithmetic, which is torch's:

- the parameters are *registered*, so `state_dict` captures them, `.to()` moves
  them and an optimizer sees them, while the distribution object — which has
  none of those — is rebuilt per call and never stored;
- a Python float reaches the margin at full precision, rather than through
  torch's default `float32`;
- coverage is uneven and the class says so out loud: `Gamma` / `Chi2` are
  rescued by inverting `cdf`, `Beta` / `StudentT` cannot be rescued at all, and
  discrete families are out of scope on this lane;
- `Gamma.cdf` has no derivative in its shape parameter, so a Gamma margin is
  fittable by its own likelihood but not through a copula.

Every family is checked against the equivalent SciPy distribution, which is the
independent implementation available.
"""

from __future__ import annotations

import pickle
from typing import Any

import numpy as np
import pytest

torch = pytest.importorskip("torch")
stats = pytest.importorskip("scipy.stats")

from pyvinecopulib.torch import TorchMargin  # noqa: E402
from .helpers import widen  # noqa: E402

_D = torch.distributions
_F64 = torch.float64

#: Quantile levels every family is probed at; safely inside every support.
_LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])

#: Family -> (factory, parameters, equivalent frozen SciPy distribution). Covers
#: every `torch.distributions` family with a usable `cdf`: the nine that also
#: have `icdf`, plus `Gamma` / `Chi2`, which reach it by inversion.
_FAMILIES: dict[str, tuple[Any, dict[str, float], Any]] = {
  "Normal": (_D.Normal, {"loc": 0.3, "scale": 1.4}, stats.norm(0.3, 1.4)),
  "Uniform": (
    _D.Uniform,
    {"low": -1.0, "high": 2.0},
    stats.uniform(-1.0, 3.0),
  ),
  "Exponential": (_D.Exponential, {"rate": 1.7}, stats.expon(scale=1 / 1.7)),
  "Cauchy": (_D.Cauchy, {"loc": 0.2, "scale": 0.9}, stats.cauchy(0.2, 0.9)),
  "Laplace": (
    _D.Laplace,
    {"loc": -0.4, "scale": 1.1},
    stats.laplace(-0.4, 1.1),
  ),
  "LogNormal": (
    _D.LogNormal,
    {"loc": 0.1, "scale": 0.7},
    stats.lognorm(0.7, scale=float(np.exp(0.1))),
  ),
  "Weibull": (
    _D.Weibull,
    {"scale": 1.3, "concentration": 2.1},
    stats.weibull_min(2.1, scale=1.3),
  ),
  "Gumbel": (_D.Gumbel, {"loc": 0.5, "scale": 1.2}, stats.gumbel_r(0.5, 1.2)),
  "HalfNormal": (_D.HalfNormal, {"scale": 1.6}, stats.halfnorm(0.0, 1.6)),
  "Gamma": (
    _D.Gamma,
    {"concentration": 2.4, "rate": 1.5},
    stats.gamma(2.4, scale=1 / 1.5),
  ),
  "Chi2": (_D.Chi2, {"df": 3.0}, stats.chi2(3.0)),
}

#: Families with no `cdf` at all, so no margin and no possible fallback.
_NO_CDF = {
  "Beta": _D.Beta(2.0, 3.0),
  "StudentT": _D.StudentT(5.0),
}

#: Families whose support is discrete; out of scope on the torch lane.
_DISCRETE = {
  "Poisson": _D.Poisson(3.0),
  "Binomial": _D.Binomial(10, 0.3),
  "Categorical": _D.Categorical(torch.tensor([0.2, 0.3, 0.5])),
}


def _normal(loc: float = 0.0, scale: float = 1.0, **kwargs: Any) -> TorchMargin:
  """A normal margin, the workhorse of these tests."""
  return TorchMargin(_D.Normal, {"loc": loc, "scale": scale}, **kwargs)


# --- agreement with SciPy --------------------------------------------------- #


@pytest.mark.parametrize("family", sorted(_FAMILIES))
def test_evaluation_matches_scipy(family: str) -> None:
  """`pdf` / `logpdf` / `cdf` / `icdf` agree with the SciPy equivalent."""
  factory, parameters, ref = _FAMILIES[family]
  margin = TorchMargin(factory, parameters)
  x = ref.ppf(_LEVELS)
  x_t = torch.as_tensor(x, dtype=_F64)

  np.testing.assert_allclose(
    margin.pdf(x_t).detach().numpy(), ref.pdf(x), rtol=1e-12, atol=0.0
  )
  np.testing.assert_allclose(
    margin.logpdf(x_t).detach().numpy(), ref.logpdf(x), rtol=1e-12, atol=0.0
  )
  np.testing.assert_allclose(
    margin.cdf(x_t).detach().numpy(), _LEVELS, rtol=1e-12, atol=0.0
  )
  got = margin.icdf(torch.as_tensor(_LEVELS, dtype=_F64)).detach().numpy()
  np.testing.assert_allclose(got, x, rtol=1e-11, atol=1e-11)


@pytest.mark.parametrize("family", ["Gamma", "Chi2"])
def test_icdf_falls_back_to_bisection(family: str) -> None:
  """These two implement `cdf` but not `icdf`, and are rescued by inverting it.

  The bisection is the generic one on `MarginBase`, and it reaches the same
  accuracy as a closed form on the families that have one — which is why the
  `cdf`-only column of `torch.distributions` is usable at all.
  """
  factory, parameters, ref = _FAMILIES[family]
  margin = TorchMargin(factory, parameters)
  with pytest.raises(NotImplementedError):
    margin.distribution.icdf(torch.tensor([0.5], dtype=_F64))
  got = margin.icdf(torch.as_tensor(_LEVELS, dtype=_F64)).detach().numpy()
  np.testing.assert_allclose(got, ref.ppf(_LEVELS), rtol=1e-11, atol=1e-11)


@pytest.mark.parametrize("family", sorted(_NO_CDF))
def test_rejects_families_without_a_cdf(family: str) -> None:
  """No `cdf` means no copula scale, and nothing to invert; say so up front."""
  with pytest.raises(NotImplementedError, match="does not implement cdf"):
    TorchMargin.from_distribution(_NO_CDF[family])


@pytest.mark.parametrize("family", sorted(_DISCRETE))
def test_rejects_discrete_families(family: str) -> None:
  """Discrete margins need a left-limit cdf the torch lane does not carry."""
  with pytest.raises(NotImplementedError, match="continuous-only"):
    TorchMargin.from_distribution(_DISCRETE[family])


# --- the nn.Module half ----------------------------------------------------- #


def test_parameters_are_registered_and_trainable() -> None:
  """`state_dict` sees them and autograd reaches them; the point of the class.

  A `torch.distributions.Distribution` held as a plain attribute contributes
  nothing to `state_dict` and has no `.to()`, which is why the parameters are
  registered here and the distribution is rebuilt from them.
  """
  margin = _normal(0.3, 1.4)
  assert sorted(margin.state_dict()) == ["loc", "scale"]
  assert {name for name, _ in margin.named_parameters()} == {"loc", "scale"}
  assert all(p.requires_grad for p in margin.parameters())
  assert margin.parameter_names == ("loc", "scale")


def test_trainable_false_registers_buffers_instead() -> None:
  """A frozen margin still travels in `state_dict`, but carries no gradient."""
  margin = _normal(0.3, 1.4, trainable=False)
  assert sorted(margin.state_dict()) == ["loc", "scale"]
  assert list(margin.parameters()) == []
  assert {name for name, _ in margin.named_buffers()} == {"loc", "scale"}


def test_python_floats_keep_full_precision() -> None:
  """A Python float must not be rounded through torch's default `float32`.

  `torch.as_tensor(-0.2)` is `float32`, and widening that back to `float64`
  leaves the parameter wrong in its eighth digit — enough to move a fitted
  log-density by 1e-7, which is far above every parity tolerance in this suite.
  """
  margin = _normal(-0.2, 0.8)
  assert widen(margin.loc).dtype == _F64
  assert widen(margin.loc).item() == -0.2
  assert widen(margin.scale).item() == 0.8


def test_the_distribution_is_rebuilt_from_the_current_parameters() -> None:
  """A fresh object per call, so an in-place update is picked up at once."""
  margin = _normal()
  assert margin.distribution is not margin.distribution

  x = torch.tensor([0.5], dtype=_F64)
  before = margin.cdf(x).item()
  with torch.no_grad():
    widen(margin.loc).fill_(0.5)
  after = margin.cdf(x).item()
  assert before == pytest.approx(stats.norm(0.0, 1.0).cdf(0.5))
  assert after == pytest.approx(0.5)


def test_state_dict_round_trip() -> None:
  """A fresh margin loads another's parameters and reproduces its output."""
  fitted = _normal(0.3, 1.4)
  fresh = _normal()
  fresh.load_state_dict(fitted.state_dict(), strict=True)
  x = torch.tensor([-1.0, 0.0, 1.0], dtype=_F64)
  torch.testing.assert_close(fresh.cdf(x), fitted.cdf(x))


def test_to_device_round_trip(device: str) -> None:
  """`.to()` moves the registered tensors, and the rebuilt object follows.

  ``y`` is built on ``device`` because ``TorchMargin`` hands it straight to
  ``torch.distributions``, unlike ``TorchVinecop`` / ``TorchVinedist``, which
  coerce their inputs in ``_prep``.
  """
  want = torch.device(device).type
  margin = _normal(0.3, 1.4)
  moved = margin.to(device)
  assert moved is margin
  assert widen(margin.loc).device.type == want
  assert widen(margin.distribution).loc.device.type == want
  x = torch.tensor([0.0], dtype=_F64, device=device)
  assert torch.isfinite(margin.cdf(x)).all()


def test_pickle_round_trip() -> None:
  """Margins round-trip through `pickle`, as every other fitted object does."""
  margin = _normal(0.3, 1.4)
  restored = pickle.loads(pickle.dumps(margin))
  x = torch.tensor([-1.0, 0.5], dtype=_F64)
  torch.testing.assert_close(restored.logpdf(x), margin.logpdf(x))


def test_support_follows_the_parameters() -> None:
  """`Uniform` derives its support from its parameters, so it is not cached."""
  margin = TorchMargin(_D.Uniform, {"low": -1.0, "high": 2.0})
  assert margin.support == (-1.0, 2.0)
  with torch.no_grad():
    widen(margin.high).fill_(5.0)
  assert margin.support == (-1.0, 5.0)


def test_a_factory_registering_no_parameters_still_evaluates() -> None:
  """A factory closing over its own tensors is placed by construction alone."""
  loc = torch.tensor(0.0, dtype=_F64)
  margin = TorchMargin(lambda: _D.Normal(loc, torch.tensor(1.0, dtype=_F64)))
  assert margin.parameter_names == ()
  assert margin.state_dict() == {}
  draws = margin.sample(5, seeds=[3])
  assert draws.shape == (5,) and draws.dtype == _F64


# --- autograd --------------------------------------------------------------- #


def test_gradients_reach_the_parameters() -> None:
  """Both `logpdf` and `cdf` are differentiable in a normal's parameters."""
  margin = _normal(0.3, 1.4)
  x = torch.tensor([-0.5, 0.2, 1.1], dtype=_F64)
  (margin.logpdf(x).sum() + margin.cdf(x).sum()).backward()
  for parameter in margin.parameters():
    assert parameter.grad is not None
    assert torch.isfinite(parameter.grad).all()


def test_gamma_marginal_likelihood_is_differentiable() -> None:
  """A Gamma margin can be fitted by its own likelihood: `log_prob` has grads."""
  margin = TorchMargin(_D.Gamma, {"concentration": 2.0, "rate": 1.0})
  margin.logpdf(torch.tensor([0.5, 1.5], dtype=_F64)).sum().backward()
  assert widen(margin.concentration).grad is not None
  assert torch.isfinite(widen(margin.concentration).grad).all()


@pytest.mark.xfail(
  raises=NotImplementedError,
  strict=True,
  reason=(
    "torch has no derivative for igamma's first argument, so Gamma.cdf carries "
    "no dF/dconcentration. A Gamma margin is therefore fittable by its own "
    "marginal likelihood but not end-to-end through a copula, which is a "
    "concrete argument for keeping two-step estimation the default."
  ),
)
def test_gamma_shape_gradient_through_the_copula_scale() -> None:
  """Optimizing a Gamma margin through `cdf` — the copula scale — cannot work."""
  margin = TorchMargin(_D.Gamma, {"concentration": 2.0, "rate": 1.0})
  margin.cdf(torch.tensor([0.5, 1.5], dtype=_F64)).sum().backward()


def test_gamma_cdf_is_differentiable_in_its_argument() -> None:
  """What does work: the derivative in the data, which is the density."""
  margin = TorchMargin(
    _D.Gamma, {"concentration": 2.0, "rate": 1.0}, trainable=False
  )
  x = torch.tensor([0.5, 1.5], dtype=_F64, requires_grad=True)
  margin.cdf(x).sum().backward()
  np.testing.assert_allclose(
    x.grad.numpy(),
    stats.gamma(2.0).pdf([0.5, 1.5]),
    rtol=1e-12,
  )


# --- the rest of the MarginBase surface ------------------------------------- #


def test_log_prob_is_an_alias_for_logpdf() -> None:
  """`log_prob` is the torch spelling; `logpdf` is the library's."""
  margin = _normal(0.3, 1.4)
  x = torch.tensor([-1.0, 0.4], dtype=_F64)
  torch.testing.assert_close(margin.log_prob(x), margin.logpdf(x))


def test_simulate_is_seeded_and_inverts_the_cdf() -> None:
  """Draws are reproducible and land back on the levels they came from."""
  margin = _normal(0.3, 1.4)
  first = margin.sample(64, seeds=[7])
  torch.testing.assert_close(margin.sample(64, seeds=[7]), first)
  assert not torch.allclose(margin.sample(64, seeds=[8]), first)
  torch.testing.assert_close(
    margin.icdf(margin.cdf(first)), first, rtol=1e-10, atol=1e-10
  )


def test_var_type_and_is_fitted_are_the_continuous_defaults() -> None:
  """Nothing here defers parameters to `fit`, and nothing has atoms."""
  margin = _normal()
  assert margin.var_type == "c"
  assert margin.is_fitted


def test_loglik_sums_the_log_density() -> None:
  """`loglik` stays a 0-d tensor, so it remains differentiable."""
  margin = _normal(0.3, 1.4)
  x = torch.tensor([-0.5, 0.2, 1.1], dtype=_F64)
  total = margin.loglik(x)
  assert total.ndim == 0
  torch.testing.assert_close(total, margin.logpdf(x).sum())


def test_repr_names_the_family_and_its_parameters() -> None:
  """The repr has to survive a grad-tracking parameter without detaching it."""
  assert repr(_normal(0.3, 1.4)) == "TorchMargin(Normal(loc=0.3, scale=1.4))"


# --- from_distribution ------------------------------------------------------ #


def test_from_distribution_copies_the_parameters() -> None:
  """Lifting reads `arg_constraints`, and copies rather than shares."""
  source = _D.Normal(
    torch.tensor(0.3, dtype=_F64), torch.tensor(1.4, dtype=_F64)
  )
  margin = TorchMargin.from_distribution(source)
  assert margin.parameter_names == ("loc", "scale")
  with torch.no_grad():
    widen(source.loc).fill_(99.0)
  assert widen(margin.loc).item() == pytest.approx(0.3)


def test_from_distribution_rejects_unreadable_parameters() -> None:
  """A family declaring a parameter it does not expose cannot be lifted."""

  class _Opaque(_D.Distribution):
    arg_constraints = {"hyperparameter": _D.constraints.real}

  with pytest.raises(ValueError, match="arg_constraints"):
    TorchMargin.from_distribution(_Opaque(validate_args=False))


def test_validate_args_is_forwarded_to_the_family() -> None:
  """Opting into validation turns an out-of-support evaluation into an error.

  Left to the family's own default it is silent, which is the wrong answer for
  an optimizer that has stepped a parameter out of its domain.
  """
  strict = TorchMargin(
    _D.Uniform, {"low": 0.0, "high": 1.0}, validate_args=True
  )
  with pytest.raises(ValueError):
    strict.logpdf(torch.tensor([2.0], dtype=_F64))

  lax = TorchMargin(_D.Uniform, {"low": 0.0, "high": 1.0}, validate_args=False)
  assert lax.logpdf(torch.tensor([2.0], dtype=_F64)).item() == float("-inf")
