"""Tests for the backend system under ``pyvinecopulib.sklearn.backends``."""

from __future__ import annotations

import copy
import subprocess
import sys

import numpy as np
import pytest

pytest.importorskip("sklearn")

import pyvinecopulib as pv  # noqa: E402
from pyvinecopulib.core import VinecopLike  # noqa: E402
from pyvinecopulib.sklearn import VineDensity, VineRegressor  # noqa: E402
from pyvinecopulib.sklearn.backends import (  # noqa: E402
  TorchVinecopBackend,
  VinecopBackend,
  resolve_backend,
)


# ---------------------------------------------------------------------------
# Protocol conformance
# ---------------------------------------------------------------------------


class TestVinecopLikeProtocol:
  """Both ``pv.Vinecop`` and ``pv.torch.TorchVinecop`` satisfy the canonical
  :class:`pyvinecopulib.core.VinecopLike` protocol structurally, so either
  backend's fitted vine is usable through the neutral contract."""

  def test_cpp_vinecop_satisfies_protocol(self):
    rng = np.random.default_rng(0)
    U = rng.uniform(0.001, 0.999, (200, 3))
    cop = pv.Vinecop.from_data(
      U,
      controls=pv.FitControlsVinecop(
        family_set=[pv.families.tll], num_threads=1
      ),
    )
    assert isinstance(cop, VinecopLike)

  def test_torch_vinecop_satisfies_protocol(self):
    torch = pytest.importorskip("torch")
    del torch
    from pyvinecopulib.torch import TorchVinecop

    rng = np.random.default_rng(0)
    U = rng.uniform(0.001, 0.999, (200, 3))
    cop = pv.Vinecop.from_data(
      U,
      controls=pv.FitControlsVinecop(
        family_set=[pv.families.tll], num_threads=1
      ),
    )
    tv = TorchVinecop.from_vinecop(cop)
    assert isinstance(tv, VinecopLike)


# ---------------------------------------------------------------------------
# resolve_backend
# ---------------------------------------------------------------------------


class TestResolveBackend:
  def test_none_defaults_to_vinecop_backend(self):
    b = resolve_backend(None)
    assert isinstance(b, VinecopBackend)

  def test_instance_passthrough(self):
    b = VinecopBackend()
    assert resolve_backend(b) is b


# ---------------------------------------------------------------------------
# with_* immutability — forest copy-on-write semantics
# ---------------------------------------------------------------------------


class TestVinecopBackendWith:
  def test_with_num_threads_returns_new_instance(self):
    # ``FitControlsVinecop`` clamps ``num_threads`` to the number of
    # available CPU cores, so use a value that's safe on any CI
    # runner (≥ 2 cores).
    b = VinecopBackend()
    b2 = b.with_num_threads(2)
    assert b is not b2
    assert b2.controls is not None
    assert b2.controls.num_threads == 2

  def test_with_random_structure(self):
    b = VinecopBackend()
    b2 = b.with_random_structure(3, seeds=[1, 2, 3, 4, 5])
    assert b is not b2
    assert b2.structure is not None
    assert b.structure is None

  def test_with_local_random_does_not_mutate_parent(self):
    parent_controls = pv.FitControlsVinecop(
      family_set=[pv.families.tll], num_threads=1
    )
    parent_algo = parent_controls.tree_algorithm
    b = VinecopBackend(controls=parent_controls)
    b2 = b.with_local_random([7, 8, 9])
    # `with_local_random` always materialises the controls; narrow
    # the type for the assertions below.
    new_controls = b2.controls
    assert new_controls is not None
    assert new_controls.tree_algorithm == "random_weighted"
    assert list(new_controls.seeds) == [7, 8, 9]
    # Parent controls instance must remain pristine.
    assert parent_controls.tree_algorithm == parent_algo


class TestTorchBackendWith:
  def test_with_num_threads_is_noop(self):
    pytest.importorskip("torch")
    b = TorchVinecopBackend()
    assert b.with_num_threads(8) is b

  def test_with_random_structure(self):
    pytest.importorskip("torch")
    b = TorchVinecopBackend()
    b2 = b.with_random_structure(3, seeds=[1, 2, 3, 4, 5])
    assert b2 is not b
    assert b2.structure is not None
    assert b.structure is None

  def test_with_local_random_threads_seeds_into_controls(self):
    pytest.importorskip("torch")
    b = TorchVinecopBackend()
    b2 = b.with_local_random([7, 8, 9])
    assert b2.structure is None
    # Seeds / tree algorithm land on the torch controls' native selection
    # fields that `TorchVinecop.from_data(structure=None)` consumes directly.
    assert b2.controls is not None
    assert b2.controls.tree_algorithm == "random_weighted"
    assert list(b2.controls.seeds) == [7, 8, 9]
    # Parent backend's controls remain untouched (copy-on-write).
    assert b.controls is None


# ---------------------------------------------------------------------------
# Lazy torch import
# ---------------------------------------------------------------------------


def test_lazy_torch_import():
  """Importing ``pyvinecopulib.sklearn`` (and constructing a default
  estimator + a :class:`VinecopBackend`) must not pull torch in."""
  code = (
    "import sys\n"
    "import pyvinecopulib.sklearn\n"
    "from pyvinecopulib.sklearn import VineDensity\n"
    "from pyvinecopulib.sklearn.backends import VinecopBackend\n"
    "assert 'torch' not in sys.modules, "
    "'torch imported on default-only path'\n"
    "VineDensity()\n"
    "VineDensity(backend=VinecopBackend())\n"
    "assert 'torch' not in sys.modules, "
    "'torch imported by cpp-only construction'\n"
  )
  subprocess.check_call([sys.executable, "-c", code])


# ---------------------------------------------------------------------------
# Estimator wiring — sklearn dev-guide compliance
# ---------------------------------------------------------------------------


@pytest.fixture
def small_data():
  rng = np.random.default_rng(0)
  return rng.standard_normal((200, 3))


class TestEstimatorWiring:
  def test_init_stores_params_verbatim(self):
    b = VinecopBackend(controls=pv.FitControlsVinecop(num_threads=4))
    est = VineDensity(backend=b, batch_size=50, random_state=7)
    assert est.backend is b
    assert est.batch_size == 50
    assert est.random_state == 7

  def test_init_does_not_validate(self):
    # Per the sklearn dev guide, ``__init__`` stores parameters
    # as-is; out-of-range values are caught later by
    # ``_validate_params()`` at ``fit`` time. Here we confirm
    # construction with an invalid (but correctly-typed) value
    # succeeds; the matching fit-time rejection is exercised in
    # ``tests/test_sklearn_density.py::test_constructor_validation``
    # and the regressor's analogue.
    est = VineDensity(batch_size=0)
    assert est.batch_size == 0
    reg = VineRegressor(quantiles=[1.5])  # 1.5 not in (0, 1)
    assert reg.quantiles == [1.5]

  def test_fit_resolves_backend_and_random_state(self, small_data):
    est = VineDensity(random_state=42).fit(small_data)
    assert isinstance(est.backend_, VinecopBackend)
    assert isinstance(est.random_state_, np.random.RandomState)

  def test_fit_sets_feature_names_in(self, small_data):
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame(small_data, columns=["a", "b", "c"])
    est = VineDensity().fit(df)
    assert list(est.feature_names_in_) == ["a", "b", "c"]

  def test_fit_sets_n_features_in(self, small_data):
    est = VineDensity().fit(small_data)
    assert est.n_features_in_ == 3

  def test_fit_sets_schema_underscore(self, small_data):
    est = VineDensity().fit(small_data)
    assert est.schema_ == {"kde1d_types": ["continuous"] * 3}

  def test_fit_sets_structure_underscore(self, small_data):
    est = VineDensity().fit(small_data)
    assert est.structure_ is not None
    assert est.structure_.dim == 3

  def test_clone_roundtrip_preserves_state(self):
    from sklearn.base import clone

    b = VinecopBackend(controls=pv.FitControlsVinecop(num_threads=2))
    est = VineDensity(backend=b, batch_size=25, random_state=3)
    est2 = clone(est)
    # backend stored as-is — clone() copies it but the underlying
    # FitControlsVinecop should carry the same num_threads.
    assert isinstance(est2.backend, VinecopBackend)
    assert est2.backend.controls.num_threads == 2
    assert est2.batch_size == 25
    assert est2.random_state == 3

  def test_check_is_fitted_raises_pre_fit(self, small_data):
    from sklearn.exceptions import NotFittedError

    est = VineDensity()
    with pytest.raises(NotFittedError):
      est.pdf(small_data[:3])

  def test_random_state_reproducible(self, small_data):
    est1 = VineDensity(random_state=42).fit(small_data)
    est2 = VineDensity(random_state=42).fit(small_data)
    np.testing.assert_allclose(
      est1.sample(50, random_state=11),
      est2.sample(50, random_state=11),
    )


# ---------------------------------------------------------------------------
# Cross-backend parity (skip if torch missing)
# ---------------------------------------------------------------------------


class TestCrossBackend:
  def test_density_pdf_parity(self, small_data):
    pytest.importorskip("torch")
    from pyvinecopulib.torch import FitControlsTorchVinecop

    est_cpp = VineDensity().fit(small_data)
    # Both backends select their structure independently, but the torch
    # selection is an exact port of Vinecop's (same structure, same reused
    # pairs), so the densities agree to TLL-fit precision. Pin
    # cache_integrals=False: the default cached evaluation trades ~1e-3 IAE
    # for speed.
    est_torch = VineDensity(
      backend=TorchVinecopBackend(
        controls=FitControlsTorchVinecop(cache_integrals=False)
      )
    ).fit(small_data)
    p_cpp = est_cpp.pdf(small_data[:10])
    p_torch = est_torch.pdf(small_data[:10])
    np.testing.assert_allclose(p_cpp, p_torch, rtol=1e-6)

  def test_cdf_works_on_both_backends(self, small_data):
    pytest.importorskip("torch")
    est_cpp = VineDensity().fit(small_data)
    est_torch = VineDensity(backend=TorchVinecopBackend()).fit(small_data)
    c_cpp = est_cpp.cdf(small_data[:5], N=5000, random_state=1)
    c_torch = est_torch.cdf(small_data[:5], N=5000, random_state=1)
    # Both are MC estimates with N=5000; agreement to ~5%.
    np.testing.assert_allclose(c_cpp, c_torch, atol=5e-2)

  def test_torch_discrete_raises(self):
    pytest.importorskip("torch")
    pd = pytest.importorskip("pandas")
    rng = np.random.default_rng(0)
    df = pd.DataFrame(
      {
        "a": pd.Categorical([0, 1] * 50, ordered=True),
        "b": rng.standard_normal(100),
      }
    )
    with pytest.raises(NotImplementedError, match="continuous-only"):
      VineDensity(backend=TorchVinecopBackend()).fit(df)


# ---------------------------------------------------------------------------
# Forest plumbing — backends propagate via base_params
# ---------------------------------------------------------------------------


class TestForestPlumbing:
  def test_forest_local_does_not_mutate_parent_backend(self, small_data):
    from pyvinecopulib.sklearn import VineForestDensity

    parent_controls = pv.FitControlsVinecop(
      family_set=[pv.families.tll], num_threads=1
    )
    parent_algo = parent_controls.tree_algorithm
    parent_backend = VinecopBackend(controls=parent_controls)
    VineForestDensity(
      base_params={"backend": parent_backend},
      n_vines=3,
      vines_sampling="local",
      n_jobs=1,
      val_fraction=0.2,
      random_state=42,
    ).fit(small_data)
    # The forest must not mutate the parent backend's controls in place.
    assert parent_controls.tree_algorithm == parent_algo

  def test_forest_with_torch_backend(self, small_data):
    pytest.importorskip("torch")
    from pyvinecopulib.sklearn import VineForestDensity

    f = VineForestDensity(
      base_params={"backend": TorchVinecopBackend()},
      n_vines=3,
      n_jobs=1,
      val_fraction=0.2,
      random_state=42,
    ).fit(small_data)
    assert f.pdf(small_data[:5]).shape == (5,)


# ---------------------------------------------------------------------------
# FitControlsVinecop copy.copy semantics (load-bearing for with_* methods)
# ---------------------------------------------------------------------------


def test_fitcontrols_copy_independent():
  c = pv.FitControlsVinecop(family_set=[pv.families.tll], num_threads=2)
  c2 = copy.copy(c)
  c2.num_threads = 8
  c2.tree_algorithm = "random_weighted"
  c2.seeds = [1, 2, 3]
  assert c.num_threads != c2.num_threads
  assert c.tree_algorithm != c2.tree_algorithm
  assert list(c.seeds) != [1, 2, 3]
