"""Tests for `pyvinecopulib.torch.TorchKde1d`.

The load-bearing claim is **parity**: a lifted grid evaluates to the same
numbers as the compiled `Kde1d` it came from, across all three variable types,
bounded and unbounded, at every local-polynomial degree. Three of those numbers
come from C++ behavior that looks like a bug and is not, so they are pinned
separately: the interpolant drops by `exp(-0.5)` exactly at the grid's right
end, the unnormalized integral carries no Gaussian-tail mass, and the discrete
masses are normalized by the *raw* interpolation rather than the clamped
density. `icdf` is compared with `assert_array_equal` rather than a tolerance,
because reproducing the 35-step bisection makes it exact.

Beyond parity: the grid is differentiable when a caller asks for it, and the
module round-trips through `state_dict`, `pickle` and `.to()`.
"""

from __future__ import annotations

import pickle
from typing import Any

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from pyvinecopulib.core import Kde1d, MarginLike  # noqa: E402
from pyvinecopulib.torch import TorchKde1d  # noqa: E402

_PROBS = np.array([1e-6, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1 - 1e-6])


def _sample(kind: str, n: int = 600, seed: int = 0) -> np.ndarray:
  rng = np.random.default_rng(seed)
  if kind == "continuous":
    return rng.gamma(2.0, 1.5, n)
  if kind == "unit":
    return rng.beta(2.0, 5.0, n)
  if kind == "discrete":
    return rng.poisson(3.0, n).astype(float)
  return np.where(rng.random(n) < 0.3, 0.0, rng.gamma(2.0, 1.5, n))


def _fitted(kind: str, **kwargs: Any) -> tuple[Kde1d, TorchKde1d, np.ndarray]:
  """One fit, lifted, plus the data it was fitted to."""
  y = _sample(kind)
  kde = Kde1d(**kwargs)
  kde.fit(y)
  return kde, TorchKde1d.from_kde1d(kde), y


def _t(values: np.ndarray) -> Any:
  return torch.as_tensor(values, dtype=torch.float64)


# --- parity ----------------------------------------------------------------- #


@pytest.mark.parametrize("degree", [0, 1, 2])
@pytest.mark.parametrize(
  ("kind", "kwargs"),
  [
    ("continuous", {}),
    ("continuous", {"xmin": 0.0}),
    ("unit", {"xmin": 0.0, "xmax": 1.0}),
    ("discrete", {"type": "discrete", "xmin": 0.0}),
    ("discrete", {"type": "discrete"}),
    ("zi", {"type": "zero-inflated", "xmin": 0.0}),
  ],
)
def test_pdf_and_cdf_match_the_compiled_estimator(
  kind: str, kwargs: dict, degree: int
) -> None:
  """Every type, bounded and unbounded, at every degree."""
  kde, lifted, y = _fitted(kind, degree=degree, **kwargs)
  grid = np.asarray(kde.grid_points, dtype=float)
  # Inside, on the knots, and outside on both sides.
  q = np.concatenate(
    [
      np.linspace(grid[0], grid[-1], 97),
      grid[:4],
      grid[-4:],
      np.unique(np.round(y)),
      [grid[0] - 1.0, grid[-1] + 1.0],
    ]
  )
  np.testing.assert_allclose(
    lifted.pdf(_t(q)).numpy(), np.asarray(kde.pdf(q)), rtol=1e-12, atol=1e-12
  )
  np.testing.assert_allclose(
    lifted.cdf(_t(q)).numpy(), np.asarray(kde.cdf(q)), rtol=1e-12, atol=1e-12
  )


@pytest.mark.parametrize(
  ("kind", "kwargs"),
  [
    ("continuous", {}),
    ("discrete", {"type": "discrete", "xmin": 0.0}),
    ("zi", {"type": "zero-inflated", "xmin": 0.0}),
  ],
)
def test_icdf_matches_the_compiled_estimator_exactly(
  kind: str, kwargs: dict
) -> None:
  """Not a tolerance: the 35-step bisection is reproduced step for step."""
  kde, lifted, _ = _fitted(kind, **kwargs)
  np.testing.assert_array_equal(
    lifted.icdf(_t(_PROBS)).numpy(), np.asarray(kde.icdf(_PROBS))
  )


def test_both_tails_leave_the_grid_continuously() -> None:
  """Each tail decays from the end it leaves, so neither jumps there.

  The Gaussian tail is written in the end cell's normalized coordinate, and
  upstream measures it from the grid's own end -- `t` on the left, `t - 1` on
  the right. Get the right one wrong and `pdf(grid[-1])` drops by a factor
  `exp(-0.5)` while the left end looks fine, because `exp(0)` is 1.

  The fit has to be **bounded** for this to bite: on an unbounded one the grid
  runs out to where the density has vanished, so `values[-1]` is exactly zero
  and both sides of the comparison are trivially `0`.
  """
  y = np.random.default_rng(3).uniform(0.0, 1.0, 600)
  kde = Kde1d(xmin=0.0, xmax=1.0)
  kde.fit(y)
  lifted = TorchKde1d.from_kde1d(kde)
  grid = np.asarray(kde.grid_points, dtype=float)
  values = np.asarray(kde.values, dtype=float)
  assert values[-1] > 1e-3, "the case would be vacuous with a vanished tail"

  for end, expected in ((grid[-1:], values[-1]), (grid[:1], values[0])):
    got = float(lifted.pdf(_t(end)).numpy()[0])
    np.testing.assert_allclose(got, expected, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(
      got, float(np.asarray(kde.pdf(end))[0]), rtol=1e-12, atol=1e-12
    )

  # And the decay away from each end matches the compiled estimator too, which
  # is what distinguishes a re-centred tail from a merely continuous one.
  outside = np.array([grid[0] - 0.3, grid[-1] + 0.3])
  np.testing.assert_allclose(
    lifted.pdf(_t(outside)).numpy(),
    np.asarray(kde.pdf(outside)),
    rtol=1e-12,
    atol=1e-12,
  )


def test_the_cdf_saturates_at_the_grid_mass() -> None:
  """`integrate` adds no tail contribution, so the normalized cdf hits 1 there."""
  kde, lifted, _ = _fitted("continuous")
  grid = np.asarray(kde.grid_points, dtype=float)
  far = _t(np.array([grid[-1], grid[-1] + 10.0]))
  np.testing.assert_allclose(
    lifted.cdf(far).numpy(), np.asarray(kde.cdf(far.numpy())), atol=1e-12
  )
  np.testing.assert_allclose(lifted.cdf(far).numpy(), [1.0, 1.0], atol=1e-12)


def test_nan_in_gives_nan_out() -> None:
  """As `unaryExpr_or_nan` does, rather than propagating a cell index."""
  _, lifted, _ = _fitted("continuous")
  q = _t(np.array([1.0, np.nan, 2.0]))
  assert bool(torch.isnan(lifted.pdf(q)[1]))
  assert bool(torch.isnan(lifted.cdf(q)[1]))
  assert bool(torch.isnan(lifted.icdf(_t(np.array([0.5, np.nan])))[1]))


def test_discrete_masses_sum_to_one_over_the_lattice() -> None:
  """The normalization is the whole point of the discrete branch."""
  _, lifted, _ = _fitted("discrete", type="discrete", xmin=0.0)
  levels = _t(np.arange(-2.0, 40.0))
  assert float(lifted.pdf(levels).sum()) == pytest.approx(1.0, abs=1e-10)
  # Off-lattice points carry no mass at all.
  assert float(lifted.pdf(_t(np.array([1.5, 2.25]))).abs().sum()) == 0.0


def test_cdf_left_is_derived_for_a_discrete_margin() -> None:
  """Inherited from `MarginBase`, which is what lets the copula difference it."""
  kde, lifted, _ = _fitted("discrete", type="discrete", xmin=0.0)
  y = _t(np.arange(0.0, 8.0))
  np.testing.assert_allclose(
    lifted.cdf_left(y).numpy(),
    np.asarray(kde.cdf_left(y.numpy())),
    rtol=1e-12,
    atol=1e-12,
  )
  assert bool(torch.all(lifted.cdf_left(y) <= lifted.cdf(y)))


# --- the contract ----------------------------------------------------------- #


def test_it_satisfies_the_margin_contract() -> None:
  """Structurally, so it drops into a `Vinedist` with no adapter."""
  _, lifted, _ = _fitted("continuous")
  assert isinstance(lifted, MarginLike)
  assert isinstance(lifted, torch.nn.Module)
  assert lifted.supports_weights is True
  assert lifted.supported_var_types == ("c", "d", "zi")


@pytest.mark.parametrize(
  ("kwargs", "expected"),
  [
    ({}, "c"),
    ({"type": "discrete"}, "d"),
    ({"type": "zero-inflated"}, "zi"),
  ],
)
def test_var_type_maps_the_compiled_spelling(
  kwargs: dict, expected: str
) -> None:
  """`kde_type` keeps the hyphenated compiled name; `var_type` is the contract's."""
  margin = TorchKde1d(**kwargs)
  assert margin.var_type == expected
  assert margin.kde_type == kwargs.get("type", "continuous")


def test_type_does_not_shadow_the_module_dtype_cast() -> None:
  """`nn.Module.type` has to keep working, which is why `kde_type` exists."""
  _, lifted, _ = _fitted("continuous")
  assert lifted.type(torch.float32).grid_points.dtype is torch.float32


def test_support_reports_the_declared_bounds() -> None:
  _, bounded, _ = _fitted("unit", xmin=0.0, xmax=1.0)
  assert bounded.support == (0.0, 1.0)
  _, unbounded, _ = _fitted("continuous")
  assert unbounded.support == (float("-inf"), float("inf"))


def test_an_unfitted_margin_refuses_to_evaluate() -> None:
  margin = TorchKde1d()
  assert margin.is_fitted is False
  with pytest.raises(RuntimeError, match="not fitted"):
    margin.pdf(_t(np.array([0.0])))


def test_an_unknown_type_is_refused_at_construction() -> None:
  with pytest.raises(ValueError, match="unknown type"):
    TorchKde1d(type="zero_inflated")  # the underscore spelling never existed


def test_covariates_are_refused_at_fit_time() -> None:
  """A kernel density reads none, so accepting them would fit the wrong model."""
  with pytest.raises(ValueError, match="supports_covariates"):
    TorchKde1d().fit(_t(_sample("continuous")), x=_t(np.zeros((600, 2))))


def test_weights_change_the_fit() -> None:
  y = _sample("continuous")
  w = np.linspace(0.1, 2.0, y.size)
  plain = TorchKde1d().fit(_t(y))
  weighted = TorchKde1d().fit(_t(y), weights=_t(w))
  assert not np.allclose(
    plain.values.numpy(), weighted.values.numpy(), atol=1e-8
  )
  # And they agree with the compiled fit on the same weights.
  reference = Kde1d()
  reference.fit(y, w)
  np.testing.assert_allclose(
    weighted.values.numpy(),
    np.asarray(reference.values),
    rtol=1e-12,
    atol=1e-12,
  )


def test_probabilities_outside_the_unit_interval_are_refused() -> None:
  _, lifted, _ = _fitted("continuous")
  with pytest.raises(ValueError, match="must lie in"):
    lifted.icdf(_t(np.array([-0.1])))


def test_loglik_and_n_parameters_come_from_the_fit() -> None:
  kde, lifted, y = _fitted("continuous")
  assert lifted.loglik() == pytest.approx(float(kde.loglik()), abs=1e-12)
  assert lifted.n_parameters == pytest.approx(float(kde.edf), abs=1e-12)
  # With data it evaluates instead, as `MarginBase.loglik` promises.
  assert float(lifted.loglik(_t(y))) == pytest.approx(
    float(kde.loglik(y)), rel=1e-10
  )


def test_a_grid_supplied_directly_reports_no_fitted_loglik() -> None:
  kde, _, _ = _fitted("continuous")
  built = TorchKde1d.from_grid(
    _t(np.asarray(kde.grid_points, dtype=float)),
    _t(np.asarray(kde.values, dtype=float)),
  )
  assert built.is_fitted
  assert np.isnan(built.n_parameters)
  with pytest.raises(RuntimeError, match="supplied directly"):
    built.loglik()


def test_from_grid_validates_its_input() -> None:
  with pytest.raises(ValueError, match="same length"):
    TorchKde1d.from_grid(_t(np.arange(4.0)), _t(np.arange(3.0)))
  with pytest.raises(ValueError, match="strictly ascending"):
    TorchKde1d.from_grid(_t(np.array([0.0, 2.0, 1.0])), _t(np.ones(3)))


def test_from_kde1d_refuses_an_unfitted_estimator() -> None:
  with pytest.raises(ValueError, match="not fitted"):
    TorchKde1d.from_kde1d(Kde1d())


# --- autograd and placement ------------------------------------------------- #


def test_gradients_reach_the_grid_values() -> None:
  """Opt-in, as `TorchBicop` does with its grid: the density is fitted, not learned."""
  _, lifted, y = _fitted("continuous")
  assert lifted.values.requires_grad is False
  lifted.values.requires_grad_(True)

  loss = -lifted.logpdf(_t(y[:64])).mean()
  loss.backward()
  grad = lifted.values.grad
  assert grad is not None
  assert bool(torch.isfinite(grad).all())
  assert float(grad.abs().sum()) > 0.0


def test_the_quantile_carries_an_exact_gradient() -> None:
  """The forward value stays the bisection's; the gradient comes from Newton.

  `dq/dtheta = -(dF/dtheta) / f(q)` by the implicit function theorem, which one
  Newton step expresses exactly -- so differentiating the correction while
  returning the bisection gives a usable gradient without moving the value.
  """
  kde, lifted, _ = _fitted("continuous")
  reference = np.asarray(kde.icdf(_PROBS))
  lifted.values.requires_grad_(True)

  q = lifted.icdf(_t(_PROBS))
  np.testing.assert_array_equal(q.detach().numpy(), reference)

  q.sum().backward()
  grad = lifted.values.grad
  assert grad is not None
  assert bool(torch.isfinite(grad).all())
  assert float(grad.abs().sum()) > 0.0


def test_the_quantile_is_differentiable_in_the_probability() -> None:
  """`d icdf / d p = 1 / f(q)`, and it must not depend on the grid being learned.

  The Newton correction that supplies the gradient was gated on
  `values.requires_grad`, so `d icdf/d p` was dead for a fitted, fixed grid --
  the common case, since the density is fitted rather than learned. The
  correction is exact, so this is an equality rather than a tolerance.
  """
  kde, lifted, _ = _fitted("continuous")
  assert lifted.values.requires_grad is False

  p = _t(_PROBS).requires_grad_(True)
  q = lifted.icdf(p)
  assert q.requires_grad
  # The value is still the bisection's, bit for bit.
  np.testing.assert_array_equal(
    q.detach().numpy(), np.asarray(kde.icdf(_PROBS))
  )

  (grad,) = torch.autograd.grad(q.sum(), p)
  np.testing.assert_allclose(
    grad.numpy(), 1.0 / lifted.pdf(q.detach()).numpy(), rtol=1e-12, atol=0.0
  )

  # Under `no_grad` the value is unchanged, so the correction cannot move it.
  with torch.no_grad():
    np.testing.assert_array_equal(
      lifted.icdf(_t(_PROBS)).numpy(), np.asarray(kde.icdf(_PROBS))
    )


def test_state_dict_round_trip() -> None:
  kde, lifted, _ = _fitted("discrete", type="discrete", xmin=0.0)
  restored = TorchKde1d(type="discrete", xmin=0.0)
  restored.load_state_dict(lifted.state_dict())
  q = _t(np.arange(0.0, 10.0))
  np.testing.assert_array_equal(restored.pdf(q).numpy(), lifted.pdf(q).numpy())
  assert set(lifted.state_dict()) == {"grid_points", "values", "prob0"}


def test_pickle_round_trip() -> None:
  kde, lifted, _ = _fitted("zi", type="zero-inflated", xmin=0.0)
  restored = pickle.loads(pickle.dumps(lifted))
  q = _t(np.array([0.0, 0.5, 1.0, 4.0]))
  np.testing.assert_array_equal(restored.pdf(q).numpy(), lifted.pdf(q).numpy())
  np.testing.assert_array_equal(restored.cdf(q).numpy(), lifted.cdf(q).numpy())
  assert restored.kde_type == "zero-inflated"
  assert restored.support == lifted.support == (0.0, float("inf"))


def test_to_moves_every_buffer() -> None:
  _, lifted, _ = _fitted("continuous")
  moved = lifted.to(torch.float32)
  assert moved.grid_points.dtype is torch.float32
  assert moved.values.dtype is torch.float32
  assert moved.prob0.dtype is torch.float32


def test_sample_lands_in_the_support() -> None:
  _, lifted, _ = _fitted("unit", xmin=0.0, xmax=1.0)
  draws = lifted.sample(256, seeds=[1])
  assert draws.shape == (256,)
  assert bool(torch.all((draws >= 0.0) & (draws <= 1.0)))


def test_a_discrete_sample_lands_on_the_lattice() -> None:
  _, lifted, _ = _fitted("discrete", type="discrete", xmin=0.0)
  draws = lifted.sample(256, seeds=[2])
  np.testing.assert_array_equal(draws.numpy(), np.round(draws.numpy()))
