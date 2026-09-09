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

from pyvinecopulib.margins import FitControlsMargin  # noqa: E402
from pyvinecopulib.torch import (  # noqa: E402
  FitControlsTorchVinecop,
  TorchKde1d,
  TorchDistributionMargin,
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


def _margins() -> list[TorchDistributionMargin]:
  """One normal margin per column, matching `_PARAMS`."""
  return [
    TorchDistributionMargin(_D.Normal, {"loc": loc, "scale": scale})
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
  what separates the two sides is `TorchTllBicop`'s bilinear grid against the C++
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
  assert all(isinstance(m, TorchDistributionMargin) for m in dist.margins)
  # The very same objects the registered `ModuleList` holds, not copies.
  registered = [
    m for m in dist.modules() if isinstance(m, TorchDistributionMargin)
  ]
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
    [
      TorchDistributionMargin(_D.Normal, {"loc": 0.0, "scale": 1.0})
      for _ in range(3)
    ],
  )
  fresh.load_state_dict(dist.state_dict(), strict=True)
  torch.testing.assert_close(fresh.logpdf(x), dist.logpdf(x))


def test_to_device_round_trip(
  device: str, data: np.ndarray, dist: TorchVinedist
) -> None:
  """`.to()` walks the registered children and the evaluation still runs.

  ``x`` stays on the host: coercing it is ``_prep``'s job, and
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
  shared = TorchDistributionMargin(_D.Normal, {"loc": 0.0, "scale": 1.0})
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


def test_an_unfitted_broadcast_margin_is_copied_on_this_lane(
  copula: pv.Vinecop, data: np.ndarray
) -> None:
  """A torch margin is an `nn.Module`, so `callable(margin)` is `True`.

  The copy is guarded on there being no `cdf`, which is what separates a
  *fitter* -- a plain callable handed the column -- from a margin. Testing
  `callable` alone treated every torch margin as a fitter, so the copy never
  happened and all three columns shared one estimate.
  """
  torch_copula = TorchVinecop.from_vinecop(copula, cache_integrals=False)
  dist = TorchVinedist(torch_copula, TorchKde1d())
  assert len({id(m) for m in dist.margins}) == 3

  y = torch.as_tensor(data, dtype=_F64)
  dist.fit(y)
  assert torch.isfinite(dist.logpdf(y)).all()
  medians = [float(np.median(data[:, k])) for k in range(3)]
  for j, margin in enumerate(dist.margins):
    center = float(margin.icdf(torch.tensor([0.5], dtype=_F64))[0])
    closest = min(range(3), key=lambda k: abs(center - medians[k]))
    assert closest == j


def test_resolve_margins_copies_an_unfitted_torch_margin() -> None:
  """The same guard, at the resolver every `from_data` goes through."""
  from pyvinecopulib.margins import resolve_margins

  resolved = resolve_margins(TorchKde1d(), 3)
  assert len({id(m) for m in resolved}) == 3
  assert all(isinstance(m, TorchKde1d) for m in resolved)


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


def test_from_data_refuses_a_family_set_it_cannot_search(
  data: np.ndarray,
) -> None:
  """`TorchKde1d` reads no controls, so a `family_set` must be a refusal.

  Answering a parametric request with a kernel density is the silent downgrade
  the weights contract already refuses. The margin declares that it cannot
  search, which is what turns the request into an error -- introspection cannot
  answer it, since the fit accepts a `controls` argument either way.
  """
  from pyvinecopulib.margins import FitControlsMargin

  assert not TorchKde1d.supports_controls
  with pytest.raises(TypeError, match="cannot select a family"):
    TorchVinedist.from_data(
      torch.as_tensor(data, dtype=_F64),
      margin_controls=FitControlsMargin(family_set=["gamma"]),
    )
  # A declared type or support is a *default*, so it is still honored.
  fitted = TorchVinedist.from_data(
    torch.as_tensor(data, dtype=_F64),
    margin_controls=FitControlsMargin(support=(-10.0, 10.0)),
  )
  assert all(isinstance(m, TorchKde1d) for m in fitted.margins)


@pytest.mark.parametrize(
  ("declared", "expected_kde_type", "expected_var_type"),
  [
    ("c", "continuous", "c"),
    ("d", "discrete", "d"),
    ("zi", "zero-inflated", "d"),
  ],
)
def test_margin_controls_declare_the_variable_type(
  declared: str, expected_kde_type: str, expected_var_type: str
) -> None:
  """Every declared type reaches the torch margin's constructor.

  Parametrized over all three because the two lanes translate the declaration
  separately: the core `Kde1d` accepts either spelling of the zero-inflated
  type and `TorchKde1d` accepts only the hyphenated one, so a second copy of
  the mapping diverged silently on exactly that value.
  """
  rng = np.random.default_rng(0)
  y = torch.as_tensor(
    np.column_stack([rng.normal(size=300), rng.poisson(3.0, 300).astype(float)])
  )
  dist = TorchVinedist.from_data(
    y, margin_controls={1: FitControlsMargin(var_type=declared)}
  )
  margin = cast("Any", dist.margins[1])
  assert margin.kde_type == expected_kde_type
  assert dist.var_types[1] == expected_var_type
  assert torch.isfinite(dist.logpdf(y)).all()


def test_margin_controls_declare_a_bound() -> None:
  """A declared support bounds the margin the library builds."""
  rng = np.random.default_rng(1)
  y = torch.as_tensor(rng.gamma(2.0, 1.0, size=(400, 2)))
  bounded = TorchVinedist.from_data(
    y, margin_controls=FitControlsMargin(support=(0.0, None))
  )
  assert all(float(cast("Any", m).xmin) == 0.0 for m in bounded.margins)


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
  # term and the same torch margins theirs, so a wrong assembly -- a missing
  # left limit, a mis-ordered column -- is a wrong number. `structure=None`, as
  # the torch fit used: both selectors reuse the pairs they fitted, so a pair
  # can arrive at its slot flipped, and a flipped TLL fit is not the fit of the
  # swapped arguments (2.7e-4 from the renormalization sweeps, 4.2e-3 more on a
  # discrete edge). Both behaviors are upstream's, so compare select to
  # select.
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
    kde = cast("TorchKde1d", margin)
    expected = expected + kde.logpdf(y_t[:, j]).detach().numpy()

  np.testing.assert_allclose(
    dist.logpdf(y_t).detach().numpy(), expected, rtol=1e-11, atol=1e-11
  )


def test_rejects_a_margin_with_atoms_and_no_left_limit(
  copula: pv.Vinecop,
) -> None:
  """A margin declaring atoms must supply the left limit the cascade needs.

  `TorchDistributionMargin` cannot: `torch.distributions`' discrete families implement
  neither `cdf` nor `icdf`, so there is nothing to take a left limit of.
  """

  class _Atomic(TorchKde1d):
    """Declares atoms and hides the inherited left limit."""

    cdf_left = None

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
  device: str | None, dtype: torch.dtype | None
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


def test_fit_and_select_work_on_every_torch_vine_distribution() -> None:
  """``fit`` and ``select`` raised ``TypeError`` for every torch distribution.

  ``_bind_dist`` runs again on each refit, and the base stored the margins with
  a plain assignment. That is fine the first time -- ``_margins`` is not yet a
  registered child -- and fatal the second, when ``nn.Module`` refuses a tuple
  over the ``ModuleList`` it is tracking. The storage is a hook now, so the
  torch lane installs a ``ModuleList`` on every bind rather than repairing one
  the base already wrote.
  """
  y = torch.as_tensor(
    np.random.default_rng(4).normal(size=(300, 3)), dtype=torch.float64
  )
  dist = TorchVinedist.from_data(y)
  for verb in ("fit", "select"):
    returned = getattr(dist, verb)(y)
    assert returned is dist
    # Still registered children, or the parameters would be invisible to
    # `state_dict`, `.to()` and every optimizer.
    assert isinstance(dist._margins, torch.nn.ModuleList)
    assert any("_margins" in key for key in dist.state_dict())
    assert bool(torch.isfinite(dist.logpdf(y)).all())
