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

from typing import Any, cast

import numpy as np
import pytest

import pyvinecopulib as pv

from .helpers import widen

torch = pytest.importorskip("torch")
stats = pytest.importorskip("scipy.stats")

from pyvinecopulib.torch import (  # noqa: E402
  FitControlsTorchVinecop,
  TorchKde1d,
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
  assert any(key.startswith("_vinecop.") for key in keys)
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
  widen(dist.vinecop).pdf(dist.marginal_cdf(x), batched=True)
  assert set(dist.state_dict()) == keys_before

  fresh = TorchVinedist(
    TorchVinecop.from_vinecop(copula, cache_integrals=False),
    [TorchMargin(_D.Normal, {"loc": 0.0, "scale": 1.0}) for _ in range(3)],
  )
  fresh.load_state_dict(dist.state_dict(), strict=True)
  torch.testing.assert_close(fresh.logpdf(x), dist.logpdf(x))


def test_to_device_round_trip(
  device: str, data: np.ndarray, dist: TorchVinedist
) -> None:
  """`.to()` walks the registered children and the evaluation still runs.

  ``x`` deliberately stays on the host: coercing it is ``_prep``'s job, and
  this is the test that would notice if it stopped doing it.
  """
  want = torch.device(device).type
  x = torch.as_tensor(data, dtype=_F64)
  moved = dist.to(device)
  assert moved is dist
  assert all(widen(m).loc.device.type == want for m in dist.margins)
  out = dist.logpdf(x)
  assert torch.isfinite(out).all()
  assert out.device.type == want


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
    for parameter in widen(margin).parameters():
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


def test_from_data_fits_end_to_end_in_torch(data: np.ndarray) -> None:
  """Margins and copula both fitted on tensors, on one device, in one dtype.

  This is what `TorchKde1d` unlocked: before it there was no torch marginal
  estimator, so the two-step fit could only produce `Kde1d` margins and a
  compiled copula -- neither of which this class can hold.
  """
  y = torch.as_tensor(data, dtype=torch.float64)
  dist = TorchVinedist.from_data(y)

  assert isinstance(dist.vinecop, TorchVinecop)
  assert all(isinstance(m, TorchKde1d) for m in dist.margins)
  assert dist.var_types == ["c"] * y.shape[1]

  logpdf = dist.logpdf(y[:32])
  assert logpdf.shape == (32,)
  assert bool(torch.isfinite(logpdf).all())
  # The Sklar identity, on the object's own terms.
  manual = torch.log(dist.vinecop.pdf(dist.marginal_cdf(y[:32])))
  for j, margin in enumerate(dist.margins):
    # `logpdf` is an optional capability, not a protocol member.
    lifted: Any = margin
    manual = manual + lifted.logpdf(y[:32, j])
  torch.testing.assert_close(logpdf, manual, rtol=1e-10, atol=1e-10)


def test_from_data_refuses_covariates(data: np.ndarray) -> None:
  """No torch margin reads them, so an unconditional fit would be a lie."""
  y = torch.as_tensor(data, dtype=torch.float64)
  with pytest.raises(NotImplementedError, match="takes no covariates"):
    TorchVinedist.from_data(y, x=torch.zeros(y.shape[0], 2))


def test_holds_margins_with_atoms() -> None:
  """A `TorchKde1d` with atoms is accepted, and matches the NumPy distribution.

  What a margin with atoms must supply is the left limit the copula's discrete
  cascade differences; `TorchKde1d` inherits `cdf_left` from `MarginBase`. The
  reference spells Sklar's factorization out by hand over the *same* margins, so
  what is checked is the assembly and the copula half, not the marginal fit.
  """
  rng = np.random.default_rng(11)
  n = 600
  base = rng.standard_normal((n, 1))
  y = np.column_stack(
    [
      rng.poisson(np.exp(0.4 * base.ravel() + 1.0)).astype(float),
      0.7 * base.ravel() + 0.7 * rng.standard_normal(n),
      0.7 * base.ravel() + 0.7 * rng.standard_normal(n),
    ]
  )
  # `from_data` does not guess a variable's type any more than the NumPy
  # `Vinedist.from_data` does -- the caller declares it on the margin.
  dist = TorchVinedist.from_data(
    torch.as_tensor(y, dtype=_F64),
    margins=[
      TorchKde1d(type="discrete", xmin=0.0),
      TorchKde1d(),
      TorchKde1d(),
    ],
  )
  assert dist.var_types == ["d", "c", "c"]

  y_t = torch.as_tensor(y, dtype=_F64)
  # The copula-scale layout the discrete cascade consumes: three value columns
  # plus one left limit for the count variable.
  u = dist.copula_data(dist.margins, y_t)
  assert u.shape == (n, 4)
  u_np = u.detach().numpy()

  # Sklar's factorization, spelled out: the compiled vine supplies the copula
  # term and the same torch margins supply theirs, so a wrong assembly -- a
  # missing left limit, a mis-ordered column -- shows up as a wrong number.
  # `structure=None`, as the torch fit used. Both selectors reuse the pairs they
  # fitted while selecting, so a pair may arrive at its slot flipped -- and a
  # flipped TLL fit is not the fit of the swapped arguments: the grid's margin
  # renormalization runs a fixed three alternating sweeps, which does not commute
  # with transposition (2.7e-4), and on a discrete edge the latent sample is
  # drawn per column position on top of that (4.2e-3 on the density). Both are
  # upstream's and reproduced exactly, so compare select to select.
  cop = pv.Vinecop.from_data(
    u_np,
    var_types=["d", "c", "c"],
    controls=pv.FitControlsVinecop(family_set=[pv.families.tll], num_threads=1),
  )
  assert np.array_equal(
    np.asarray(cop.structure.matrix),
    np.asarray(dist.vinecop.structure.matrix),
  )
  expected = np.log(np.asarray(cop.pdf(u_np)))
  for j, margin in enumerate(dist.margins):
    # `logpdf` is an optional capability on `MarginLike`, declared on the base.
    kde = cast(TorchKde1d, margin)
    expected = expected + kde.logpdf(y_t[:, j]).detach().numpy()

  np.testing.assert_allclose(
    dist.logpdf(y_t).detach().numpy(), expected, rtol=1e-11, atol=1e-11
  )


def test_rejects_a_margin_with_atoms_and_no_left_limit(
  copula: pv.Vinecop,
) -> None:
  """A margin declaring atoms must supply the left limit the cascade needs.

  `TorchMargin` cannot: `torch.distributions`' discrete families implement
  neither `cdf` nor `icdf`, so there is nothing to take a left limit of.
  """

  class _Atomic(TorchKde1d):
    """Declares atoms and hides the inherited left limit."""

    cdf_left = None  # type: ignore[assignment]

  discrete = _Atomic(type="discrete", xmin=0.0)
  discrete.fit(
    torch.as_tensor(
      np.random.default_rng(6).poisson(3.0, 400).astype(float), dtype=_F64
    )
  )
  assert discrete.var_type == "d"

  torch_copula = TorchVinecop.from_vinecop(copula)
  with pytest.raises(NotImplementedError, match="no `cdf_left`"):
    TorchVinedist(torch_copula, [discrete for _ in range(3)])


@pytest.mark.parametrize(
  ("device", "dtype"),
  [(None, None), (None, torch.float32), ("cpu", torch.float64)],
)
def test_from_data_puts_everything_on_one_device_and_dtype(
  device, dtype
) -> None:
  """`from_data` documents one device and one dtype for the whole object.

  The margins took theirs from `y` while the copula took `controls.device`, so
  `from_data(..., controls=FitControlsTorchVinecop(device="cuda"))` left every
  margin on the CPU: `state_dict` spanned two devices and `logpdf` raised.
  """
  rng = np.random.default_rng(0)
  y = torch.as_tensor(rng.normal(size=(200, 3)))
  controls = FitControlsTorchVinecop(device=device, dtype=dtype)
  dist = TorchVinedist.from_data(y, controls=controls)
  tensors = [v for v in dist.state_dict().values() if hasattr(v, "device")]
  assert len({t.device for t in tensors}) == 1
  assert len({t.dtype for t in tensors}) == 1
  if dtype is not None:
    assert tensors[0].dtype == dtype
  # And the object it produced can evaluate its own data.
  assert torch.isfinite(dist.log_prob(dist.sample(4))).all()


def test_from_data_does_not_take_its_dtype_from_integer_data() -> None:
  """An integer `y` must not give the margins an integer grid."""
  rng = np.random.default_rng(1)
  y = torch.as_tensor(rng.poisson(5.0, size=(200, 2)))
  assert not y.dtype.is_floating_point
  dist = TorchVinedist.from_data(y)
  tensors = [v for v in dist.state_dict().values() if hasattr(v, "dtype")]
  assert len({t.dtype for t in tensors}) == 1
  assert tensors[0].dtype.is_floating_point
