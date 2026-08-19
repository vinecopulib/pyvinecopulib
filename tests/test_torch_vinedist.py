"""Tests for `pyvinecopulib.torch.TorchVinedist`.

Skipped without PyTorch. `Vinedist` is already array-agnostic, so the cascade
needs no separate check here — what this file pins is the `nn.Module` half and
the boundaries:

- parity of `logpdf` with the NumPy `Vinedist` on the *same* fitted model, which
  is what makes the torch lane a port rather than a second model;
- the copula and every margin are registered children, so `state_dict` captures
  all of them and nothing derived leaks in;
- the joint log-density is differentiable back to the marginal parameters, which
  is the reason to assemble the thing in torch at all;
- the parts have to be torch parts: a SciPy margin or a compiled `Vinecop`
  would silently detach every gradient, so both are refused up front.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

import pyvinecopulib as pv

torch = pytest.importorskip("torch")
stats = pytest.importorskip("scipy.stats")

from pyvinecopulib.core import MarginBase  # noqa: E402
from pyvinecopulib.torch import (  # noqa: E402
  TorchMargin,
  TorchVinecop,
  TorchVinedist,
)

_D = torch.distributions
_F64 = torch.float64

#: Marginal parameters used on both sides of the parity check.
_PARAMS = [(0.0, 1.0), (0.5, 1.5), (-0.2, 0.8)]


@pytest.fixture
def data() -> np.ndarray:
  """Three dependent columns on the real line."""
  rng = np.random.default_rng(0)
  base = rng.standard_normal((600, 1))
  return 0.6 * base + 0.4 * rng.standard_normal((600, 3))


@pytest.fixture
def copula(data: np.ndarray) -> pv.Vinecop:
  """A TLL vine fitted to the pseudo-observations of `data`."""
  return pv.Vinecop.from_data(
    pv.to_pseudo_obs(data),
    controls=pv.FitControlsVinecop(family_set=[pv.families.tll], num_threads=1),
  )


def _wide(value: object) -> Any:
  """Widen a value to ``Any``.

  ``nn.Module.__getattr__`` hands back a ``Tensor | Parameter | Module`` union
  that no static checker can narrow, and ``Vinedist.margins`` is typed by the
  margin protocol rather than by the concrete class. Both are read through this.
  """
  return value


def _margins() -> list[TorchMargin]:
  """One normal margin per column, matching `_PARAMS`."""
  return [
    TorchMargin(_D.Normal, {"loc": loc, "scale": scale})
    for loc, scale in _PARAMS
  ]


@pytest.fixture
def dist(copula: pv.Vinecop) -> TorchVinedist:
  """The distribution under test: a lifted copula plus normal margins."""
  return TorchVinedist(
    TorchVinecop.from_vinecop(copula, cache_integrals=False), _margins()
  )


# --- parity with the NumPy Vinedist ----------------------------------------- #


def test_logpdf_matches_the_numpy_vinedist(
  data: np.ndarray, copula: pv.Vinecop, dist: TorchVinedist
) -> None:
  """The same model evaluated on either lane gives the same log-density.

  The tolerance is the one `test_torch_vinecop.py` pins for the copula term:
  what separates the two sides is `TorchBicop`'s bilinear grid against the C++
  on-the-fly cascade, since the marginal terms are closed forms that agree to
  machine precision.
  """
  reference = pv.Vinedist(
    copula, [stats.norm(loc, scale) for loc, scale in _PARAMS]
  )
  got = dist.logpdf(torch.as_tensor(data, dtype=_F64)).detach().numpy()
  np.testing.assert_allclose(
    got, reference.logpdf(data), atol=1e-10, rtol=1e-10
  )


def test_marginal_cdf_matches_the_numpy_vinedist(
  data: np.ndarray, copula: pv.Vinecop, dist: TorchVinedist
) -> None:
  """The copula-scale transform agrees too, which is where a dtype slip shows."""
  reference = pv.Vinedist(
    copula, [stats.norm(loc, scale) for loc, scale in _PARAMS]
  )
  got = dist.marginal_cdf(torch.as_tensor(data, dtype=_F64)).detach().numpy()
  np.testing.assert_allclose(
    got, reference.marginal_cdf(data), atol=1e-13, rtol=1e-13
  )


# --- the nn.Module half ----------------------------------------------------- #


def test_the_parts_are_registered_children(dist: TorchVinedist) -> None:
  """The copula and every margin travel in `state_dict`.

  A plain tuple of modules is invisible to `nn.Module`, so the margins go into a
  `ModuleList`; without it no optimizer, checkpoint or `.to()` would reach their
  parameters.
  """
  keys = set(dist.state_dict())
  assert {f"_margins.{j}.{p}" for j in range(3) for p in ("loc", "scale")} <= (
    keys
  )
  assert any(key.startswith("_copula.") for key in keys)
  assert {name for name, _ in dist.named_parameters()} == {
    f"_margins.{j}.{p}" for j in range(3) for p in ("loc", "scale")
  }


def test_margins_property_reads_the_module_list(dist: TorchVinedist) -> None:
  """The public accessor still hands back a tuple of the very same objects."""
  assert isinstance(dist.margins, tuple)
  assert len(dist.margins) == dist.dim == 3
  assert all(isinstance(m, TorchMargin) for m in dist.margins)
  # The very same objects the registered `ModuleList` holds, not copies.
  registered = [m for m in dist.modules() if isinstance(m, TorchMargin)]
  assert [id(m) for m in dist.margins] == [id(m) for m in registered]


def test_state_dict_round_trip_and_no_derived_cache_leak(
  data: np.ndarray, copula: pv.Vinecop, dist: TorchVinedist
) -> None:
  """Evaluation adds no keys, and a fresh distribution loads them all.

  The copula's batched state is a memo derived from the pair copulas; if it
  entered `state_dict` every checkpoint taken after a batched call would be
  rejected by a fresh model as carrying unexpected keys.
  """
  x = torch.as_tensor(data, dtype=_F64)
  keys_before = set(dist.state_dict())
  dist.pdf(x)
  dist.rosenblatt(x)
  dist.copula.pdf(dist.marginal_cdf(x), batched=True)
  assert set(dist.state_dict()) == keys_before

  fresh = TorchVinedist(
    TorchVinecop.from_vinecop(copula, cache_integrals=False),
    [TorchMargin(_D.Normal, {"loc": 0.0, "scale": 1.0}) for _ in range(3)],
  )
  fresh.load_state_dict(dist.state_dict(), strict=True)
  torch.testing.assert_close(fresh.logpdf(x), dist.logpdf(x))


def test_to_device_round_trip(data: np.ndarray, dist: TorchVinedist) -> None:
  """`.to()` walks the registered children and the evaluation still runs."""
  x = torch.as_tensor(data, dtype=_F64)
  moved = dist.to("cpu")
  assert moved is dist
  assert all(_wide(m).loc.device.type == "cpu" for m in dist.margins)
  assert torch.isfinite(dist.logpdf(x)).all()


def test_backward_reaches_the_margin_parameters(
  data: np.ndarray, dist: TorchVinedist
) -> None:
  """The joint negative log-likelihood is differentiable in the margins.

  Gradients arrive through both routes at once: the marginal log-densities, and
  the copula evaluated at `F_j(x_j)`.
  """
  loss = -dist.logpdf(torch.as_tensor(data, dtype=_F64)).mean()
  loss.backward()
  for margin in dist.margins:
    for parameter in _wide(margin).parameters():
      assert parameter.grad is not None
      assert torch.isfinite(parameter.grad).all()
      assert not torch.allclose(parameter.grad, torch.zeros_like(parameter))


def test_an_optimizer_step_lowers_the_loss(
  data: np.ndarray, dist: TorchVinedist
) -> None:
  """End to end: the margins really are optimizable through the copula."""
  x = torch.as_tensor(data, dtype=_F64)
  optimizer = torch.optim.Adam(dist.parameters(), lr=1e-2)
  before = -dist.logpdf(x).mean().item()
  for _ in range(5):
    optimizer.zero_grad()
    loss = -dist.logpdf(x).mean()
    loss.backward()
    optimizer.step()
  assert -dist.logpdf(x).mean().item() < before


# --- the inherited surface -------------------------------------------------- #


def test_log_prob_is_an_alias_for_logpdf(
  data: np.ndarray, dist: TorchVinedist
) -> None:
  """`log_prob` is the torch spelling, and only exists on this subclass."""
  x = torch.as_tensor(data[:20], dtype=_F64)
  torch.testing.assert_close(dist.log_prob(x), dist.logpdf(x))
  assert not hasattr(pv.Vinedist, "log_prob")


def test_pdf_is_the_exponential_of_logpdf(
  data: np.ndarray, dist: TorchVinedist
) -> None:
  """The two stay consistent, and `logpdf` is the primitive."""
  x = torch.as_tensor(data[:50], dtype=_F64)
  torch.testing.assert_close(dist.pdf(x), dist.logpdf(x).exp())


def test_marginal_transforms_round_trip(
  data: np.ndarray, dist: TorchVinedist
) -> None:
  """`marginal_icdf` inverts `marginal_cdf` column by column."""
  x = torch.as_tensor(data, dtype=_F64)
  torch.testing.assert_close(
    dist.marginal_icdf(dist.marginal_cdf(x)), x, rtol=1e-10, atol=1e-10
  )


def test_rosenblatt_round_trips_through_the_data_scale(
  data: np.ndarray, dist: TorchVinedist
) -> None:
  """`inverse_rosenblatt` lands back on the observations it started from."""
  x = torch.as_tensor(data[:100], dtype=_F64)
  w = dist.rosenblatt(x)
  assert w.shape == x.shape
  torch.testing.assert_close(
    dist.inverse_rosenblatt(w), x, rtol=1e-6, atol=1e-6
  )


def test_simulate_returns_data_scale_samples(dist: TorchVinedist) -> None:
  """Samples are on the original scale and reproducible from a seed."""
  first = dist.sample(200, seeds=[7])
  assert first.shape == (200, 3)
  torch.testing.assert_close(dist.sample(200, seeds=[7]), first)
  # Each column should sit near its own margin rather than on [0, 1].
  for j, (loc, scale) in enumerate(_PARAMS):
    assert abs(first[:, j].mean().item() - loc) < 0.5 * scale


def test_loglik_sums_the_log_density(
  data: np.ndarray, dist: TorchVinedist
) -> None:
  """`loglik` stays 0-d, so it remains differentiable."""
  x = torch.as_tensor(data, dtype=_F64)
  total = dist.loglik(x)
  assert total.ndim == 0
  torch.testing.assert_close(total, dist.logpdf(x).sum())


def test_a_single_margin_is_broadcast_across_the_variables(
  copula: pv.Vinecop, data: np.ndarray
) -> None:
  """One margin standing for every variable ties their parameters together."""
  shared = TorchMargin(_D.Normal, {"loc": 0.0, "scale": 1.0})
  dist = TorchVinedist(
    TorchVinecop.from_vinecop(copula, cache_integrals=False), shared
  )
  assert dist.dim == 3
  assert all(m is shared for m in dist.margins)
  assert {name for name, _ in dist.named_parameters()} == {
    "_margins.0.loc",
    "_margins.0.scale",
  }
  assert torch.isfinite(dist.logpdf(torch.as_tensor(data, dtype=_F64))).all()


# --- boundaries ------------------------------------------------------------- #


def test_rejects_a_non_module_margin(copula: pv.Vinecop) -> None:
  """A SciPy margin would detach every gradient, so it is refused here."""
  torch_copula = TorchVinecop.from_vinecop(copula)
  with pytest.raises(TypeError, match="torch.nn.Module"):
    TorchVinedist(torch_copula, [stats.norm(0, 1) for _ in range(3)])


def test_rejects_a_compiled_vinecop(copula: pv.Vinecop) -> None:
  """The compiled copula evaluates on NumPy; point at the lift instead."""
  with pytest.raises(TypeError, match="from_vinecop"):
    TorchVinedist(copula, _margins())


def test_from_data_stays_on_the_numpy_lane(
  data: np.ndarray,
) -> None:
  """There is no torch marginal estimator, so `from_data` cannot land here.

  It would fit `Kde1d` margins and a compiled copula, neither of which
  this class can hold, so the refusal names the two routes that do work.
  """
  with pytest.raises(NotImplementedError, match="assembled, not fitted"):
    TorchVinedist.from_data(data)


def test_rejects_a_discrete_margin(copula: pv.Vinecop) -> None:
  """Atoms need a left-limit cdf the torch cascade does not carry yet."""

  class _CountMargin(MarginBase[Any], torch.nn.Module):
    """A margin that is a module, but declares atoms."""

    def __init__(self) -> None:
      torch.nn.Module.__init__(self)

    @property
    def var_type(self) -> str:
      return "d"

    def pdf(self, y: Any, *, x: Any = None) -> Any:
      return torch.full_like(y, 0.5)

    def cdf(self, y: Any, *, x: Any = None) -> Any:
      return torch.clamp(y, 0.0, 1.0)

  torch_copula = TorchVinecop.from_vinecop(copula)
  with pytest.raises(NotImplementedError, match="continuous-only"):
    TorchVinedist(torch_copula, [_CountMargin() for _ in range(3)])
