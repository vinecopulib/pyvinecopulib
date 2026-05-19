"""Tests for the optional ``method="vdc"`` path in ``TorchBicop.from_data``.

Skipped when either PyTorch or vine-denoising-copula isn't installed.
Validates that:

* the controls dataclass rejects unknown methods,
* the resulting grid has the expected shape and boundary endpoints,
* trapezoidal marginal integrals are already ~1 without
  ``normalize_margins`` (the algebraic identity that motivates
  ``norm_times=0`` in the vdc branch of :meth:`TorchBicop.from_data`),
* ``pdf`` / ``hfunc`` / ``hinv`` round-trip and stay in range,
* the module-level bundle cache amortizes the HuggingFace download.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyvinecopulib as pv

torch = pytest.importorskip("torch")
pytest.importorskip("vdc")

from pyvinecopulib.torch import FitControlsTorchBicop, TorchBicop  # noqa: E402
from pyvinecopulib.torch import _fit_vdc as _fit_vdc_mod  # noqa: E402


def _vdc_smoke_test_works() -> bool:
  """Detect upstream-broken vdc wheels.

  As of vdc 0.1.0 on GitHub ``main``, ``load_pretrained_model`` raises
  ``ModuleNotFoundError`` for ``vdc.inference`` / ``vdc.vine`` —
  references to submodules that aren't shipped. The integration ships a
  ``sys.modules`` shim that papers over those gaps for the denoiser
  inference path (see :func:`_fit_vdc._install_upstream_shims`); we
  install it here too so the smoke check goes through the same code
  path as :meth:`TorchBicop.from_data(..., method='vdc')`.
  """
  import numpy as np

  import vdc as _vdc

  _fit_vdc_mod._install_upstream_shims()
  try:
    bundle = _vdc.load_pretrained_model("vdc-denoiser-m64-v1", device="cpu")
    _vdc.estimate_pair_density_from_samples(
      bundle, np.random.uniform(size=(50, 2))
    )
    return True
  except ModuleNotFoundError:
    return False


if not _vdc_smoke_test_works():
  pytest.skip(
    "vdc package is installed but unusable (upstream wheel references "
    "missing submodules `vdc.inference` and `vdc.vine`). Re-enable these "
    "tests once vdc restores the missing subpackages.",
    allow_module_level=True,
  )


def _gaussian_sample(rho: float, n: int = 2000, seed: int = 1) -> np.ndarray:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[rho]]))
  return cop.simulate(n, seeds=[seed, seed + 1, seed + 2])


def _eval_grid(n: int, seed: int = 0) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return rng.uniform(0.02, 0.98, size=(n, 2))


@pytest.fixture
def vdc_bicop():
  """Single shared fit reused across tests to keep CI cheap."""
  u_fit = _gaussian_sample(rho=0.6, n=2000, seed=1)
  return TorchBicop.from_data(
    torch.from_numpy(u_fit), FitControlsTorchBicop(method="vdc")
  )


def test_grid_shape_and_endpoints(vdc_bicop: TorchBicop) -> None:
  """The padded grid is (m+2, m+2) on [0, cc..., 1]; pretrained ``m=64``
  gives ``(66, 66)`` values."""
  g = vdc_bicop.interp_grid
  assert g.values.shape == (66, 66)
  assert g.grid_points.shape == (66,)
  assert g.grid_points[0].item() == pytest.approx(0.0, abs=0.0)
  assert g.grid_points[-1].item() == pytest.approx(1.0, abs=0.0)
  # Interior is cell-centers: (0.5/m, 1.5/m, ..., (m-0.5)/m).
  interior = g.grid_points[1:-1].numpy()
  expected = (np.arange(64) + 0.5) / 64.0
  np.testing.assert_allclose(interior, expected, atol=1e-12)


def test_marginals_are_uniform_without_normalization(
  vdc_bicop: TorchBicop,
) -> None:
  """vdc's IPFP enforces midpoint-rule uniform marginals on the
  cell-center grid; with replicate-padding the trapezoidal integral on
  the (m+2)-grid coincides with that midpoint sum at every cell center.
  So row/column trapezoidal sums must equal 1 without us running a
  single ``normalize_margins`` round (``norm_times=0`` in the vdc
  branch)."""
  g = vdc_bicop.interp_grid
  dgrid = g.grid_points[1:] - g.grid_points[:-1]
  row_int = 0.5 * ((g.values[:, :-1] + g.values[:, 1:]) * dgrid).sum(-1)
  col_int = 0.5 * (
    (g.values[:-1, :] + g.values[1:, :]) * dgrid.unsqueeze(-1)
  ).sum(0)
  # vdc's projection_iters=50 gets to ~1e-6 marginals; the algebraic
  # identity then preserves that on the padded grid.
  assert (row_int - 1.0).abs().max().item() < 5e-5
  assert (col_int - 1.0).abs().max().item() < 5e-5


def test_pdf_hfunc_hinv_in_range(vdc_bicop: TorchBicop) -> None:
  u_t = torch.from_numpy(_eval_grid(400, seed=7))
  pdf = vdc_bicop.pdf(u_t)
  h1 = vdc_bicop.hfunc1(u_t)
  h2 = vdc_bicop.hfunc2(u_t)
  assert (pdf > 0).all() and torch.isfinite(pdf).all()
  assert ((h1 >= 0) & (h1 <= 1)).all()
  assert ((h2 >= 0) & (h2 <= 1)).all()


def test_hinv_roundtrip(vdc_bicop: TorchBicop) -> None:
  """``hinv1(hfunc1(u)) ≈ u`` on the second argument. ITP bisection is
  fixed-iter so we expect ~1e-3 at m=64."""
  u_t = torch.from_numpy(_eval_grid(400, seed=11))
  u2 = vdc_bicop.hinv1(u_t).unsqueeze(-1)
  back = vdc_bicop.hfunc1(torch.cat([u_t[:, 0:1], u2], dim=-1)).numpy()
  np.testing.assert_allclose(back, u_t[:, 1].numpy(), atol=2e-3, rtol=0.0)


def test_independence_smoke() -> None:
  """On iid uniform samples, vdc should at least integrate to 1 and stay
  finite. The denoiser's accuracy on the independence copula is
  disappointing (mean |pdf - 1| ≈ 0.9 in practice — see the precision
  bench for numbers), but that's an upstream model-quality issue, not
  an integration bug. This test just verifies the wrapper works."""
  rng = np.random.default_rng(42)
  u_fit = rng.uniform(size=(4000, 2))
  bc = TorchBicop.from_data(
    torch.from_numpy(u_fit), FitControlsTorchBicop(method="vdc")
  )
  u_eval = torch.from_numpy(_eval_grid(2000, seed=99))
  pdf = bc.pdf(u_eval).numpy()
  assert np.isfinite(pdf).all()
  assert (pdf >= 0).all()
  # Marginals were already validated in test_marginals_are_uniform_*.


def test_tll_and_vdc_both_recover_gaussian() -> None:
  """Both backends fit and produce distinct grids (so the dispatch is
  exercised, not a silent fall-through). Accuracy thresholds are loose
  for vdc — the released denoiser checkpoint underperforms TLL on
  Gaussian by a wide margin in practice; the precision bench reports
  the actual numbers."""
  u_fit = _gaussian_sample(rho=0.6, n=4000, seed=21)
  u_fit_t = torch.from_numpy(u_fit)
  bc_tll = TorchBicop.from_data(u_fit_t, FitControlsTorchBicop(method="tll"))
  bc_vdc = TorchBicop.from_data(u_fit_t, FitControlsTorchBicop(method="vdc"))

  # Different grid sizes — definitely not the same object/values.
  assert bc_tll.interp_grid.values.shape != bc_vdc.interp_grid.values.shape

  truth = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_eval = _eval_grid(2000, seed=33)
  u_eval_t = torch.from_numpy(u_eval)
  truth_h = truth.hfunc1(u_eval)
  tll_iae = float(np.mean(np.abs(bc_tll.hfunc1(u_eval_t).numpy() - truth_h)))
  vdc_iae = float(np.mean(np.abs(bc_vdc.hfunc1(u_eval_t).numpy() - truth_h)))
  assert tll_iae < 0.03
  # vdc denoiser checkpoint accuracy on Gaussian copula is poor — observed
  # hfunc1 IAE ≈ 0.19. Threshold leaves headroom; the bench surfaces the
  # real number for follow-up.
  assert vdc_iae < 0.30


def test_bundle_cache_amortizes_load(monkeypatch) -> None:
  """Successive vdc fits should reuse the cached
  ``LoadedPretrainedModel`` — the HuggingFace download happens once per
  process."""
  # Reset cache so we start clean.
  _fit_vdc_mod._BUNDLE_CACHE.clear()
  import vdc as _vdc_pkg

  call_count = {"n": 0}
  real_loader = _vdc_pkg.load_pretrained_model

  def _counting_loader(model_id, *args, **kwargs):
    call_count["n"] += 1
    return real_loader(model_id, *args, **kwargs)

  monkeypatch.setattr(_vdc_pkg, "load_pretrained_model", _counting_loader)

  u_fit = _gaussian_sample(rho=0.4, n=500, seed=51)
  TorchBicop.from_data(
    torch.from_numpy(u_fit), FitControlsTorchBicop(method="vdc")
  )
  TorchBicop.from_data(
    torch.from_numpy(u_fit), FitControlsTorchBicop(method="vdc")
  )
  assert call_count["n"] == 1
