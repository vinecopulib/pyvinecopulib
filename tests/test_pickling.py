import pickle
from typing import Any

import numpy as np
import pytest

import pyvinecopulib as pv

from .helpers import (
  compare_bicop,
  compare_kde1d,
  compare_properties,
  compare_rvinestructure,
  compare_vinecop,
  random_data,
)


def _custom_criterion(data: np.ndarray, weights: np.ndarray) -> float:
  # Module-level (hence picklable) custom tree criterion.
  return float(pv.utils.wdm(data[:, 0], data[:, 1], "tau"))


def test_fitcontrolsbicop() -> None:
  original_controls = pv.FitControlsBicop()
  original_controls.family_set = pv.families.itau
  original_controls.parametric_method = "itau"

  # Serialize the object
  serialized = pickle.dumps(original_controls)

  # Deserialize the object
  deserialized_controls = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  attrs = [
    "family_set",
    "parametric_method",
    "nonparametric_method",
    "nonparametric_mult",
    "nonparametric_grid_size",
    "selection_criterion",
    "weights",
    "psi0",
    "preselect_families",
    "num_threads",
  ]
  compare_properties(original_controls, deserialized_controls, attrs)


def test_fitcontrolsvinecop() -> None:
  # Create an instance of FitControlsVinecop with some configuration
  original_controls = pv.FitControlsVinecop()
  original_controls.family_set = pv.families.itau
  original_controls.parametric_method = "itau"

  # Serialize the object
  serialized = pickle.dumps(original_controls)

  # Deserialize the object
  deserialized_controls = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  attrs = [
    "family_set",
    "parametric_method",
    "nonparametric_method",
    "weights",
    "nonparametric_mult",
    "nonparametric_grid_size",
    "trunc_lvl",
    "tree_criterion",
    "tree_criterion_function",
    "threshold",
    "selection_criterion",
    "psi0",
    "preselect_families",
    "select_trunc_lvl",
    "select_threshold",
    "select_families",
    "show_trace",
    "num_threads",
    "tree_algorithm",
    "allow_rotations",
    "seeds",
  ]
  compare_properties(original_controls, deserialized_controls, attrs)


def test_fitcontrolsvinecop_custom_criterion() -> None:
  # A module-level custom criterion round-trips by reference; the getter
  # returns the same object after unpickling.
  original_controls = pv.FitControlsVinecop(tree_criterion="custom")
  original_controls.tree_criterion_function = _custom_criterion

  deserialized_controls = pickle.loads(pickle.dumps(original_controls))
  assert deserialized_controls.tree_criterion == "custom"
  assert deserialized_controls.tree_criterion_function is _custom_criterion

  # A non-picklable callable (lambda) raises the standard pickling error.
  original_controls.tree_criterion_function = lambda data, weights: 0.0
  with pytest.raises((pickle.PicklingError, AttributeError)):
    pickle.dumps(original_controls)


def test_bicop() -> None:
  original_bicop = pv.Bicop(pv.families.gaussian)
  original_bicop.parameters = np.array([[0.5]])

  # Serialize the object
  serialized = pickle.dumps(original_bicop)

  # Deserialize the object
  deserialized_bicop = pickle.loads(serialized)

  # Assert that the deserialized object's properties match the original
  compare_bicop(original_bicop, deserialized_bicop)


def test_fitted_bicop_preserves_diagnostics() -> None:
  """A fitted pair round-trips its attained public estimator state."""
  u = pv.to_pseudo_obs(random_data(2, 300))
  original = pv.Bicop.from_data(u)
  restored = pickle.loads(pickle.dumps(original))
  compare_bicop(original, restored)
  assert restored.nobs == original.nobs
  for name in ("loglik", "aic", "bic"):
    assert getattr(restored, name)() == pytest.approx(getattr(original, name)())


def test_rvinestructure() -> None:
  # Create an instance of RVineStructure with some configuration
  original_structure = pv.RVineStructure.sample(5)

  # Serialize the object
  serialized = pickle.dumps(original_structure)

  # Deserialize the object
  deserialized_structure = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  compare_rvinestructure(original_structure, deserialized_structure)


def test_kde1d() -> None:
  # Test with unfitted Kde1d object first
  original_kde = pv.core.Kde1d(
    xmin=-5.0,
    xmax=5.0,
    type="continuous",
    multiplier=1.5,
    degree=1,
    bandwidth=0.1,
    grid_size=100,
  )

  # Serialize the unfitted object
  serialized = pickle.dumps(original_kde)

  # Deserialize the object
  deserialized_kde = pickle.loads(serialized)

  # Assert that the deserialized object's properties match the original (unfitted)
  compare_kde1d(original_kde, deserialized_kde)

  # Now test with fitted model
  np.random.seed(1234)
  x = np.random.normal(0, 1, 100)
  original_kde.fit(x)

  # Serialize the fitted object
  serialized_fitted = pickle.dumps(original_kde)

  # Deserialize the fitted object
  deserialized_fitted = pickle.loads(serialized_fitted)

  # Assert that the deserialized fitted object's properties match the original
  compare_kde1d(original_kde, deserialized_fitted)

  # Saved fitting controls govern subsequent fits as well as inspection.
  refit_data = np.random.default_rng(9).normal(1.0, 0.7, 120)
  original_kde.fit(refit_data)
  deserialized_fitted.fit(refit_data)
  compare_kde1d(original_kde, deserialized_fitted)


def test_vinecop() -> None:
  d = 5
  n = 1000
  u = pv.to_pseudo_obs(random_data(d, n))

  controls = pv.FitControlsVinecop(family_set=[pv.families.gaussian])
  assert controls.family_set == [pv.families.gaussian]
  original_cop = pv.Vinecop.from_data(u, controls=controls)

  # Serialize the object
  serialized = pickle.dumps(original_cop)

  # Deserialize the object
  deserialized_cop = pickle.loads(serialized)

  # Ensure the deserialized object has the same attributes as the original
  compare_vinecop(original_cop, deserialized_cop)
  assert deserialized_cop.nobs == original_cop.nobs
  for name in ("loglik", "aic", "bic", "mbicv"):
    assert getattr(deserialized_cop, name)() == pytest.approx(
      getattr(original_cop, name)()
    )


def _fitted_estimator(backend: object) -> tuple[Any, np.ndarray]:
  from pyvinecopulib.sklearn import VineDensity

  rng = np.random.default_rng(0)
  cov = np.full((3, 3), 0.5) + 0.5 * np.eye(3)
  x = rng.multivariate_normal(np.zeros(3), cov, size=300)
  return VineDensity(backend=backend, random_state=0).fit(x), x


@pytest.mark.parametrize("torch_backend", [False, True])
def test_vinedensity(torch_backend: bool) -> None:
  """A fitted estimator round-trips, including its ``distribution_``.

  On the torch backend that attribute is a ``TorchVinedist`` -- an
  ``nn.Module`` holding the copula and every margin as registered children --
  so this is the check that publishing it did not cost the pickling guarantee.
  """
  pytest.importorskip("sklearn")
  backend = None
  if torch_backend:
    pytest.importorskip("torch")
    from pyvinecopulib.sklearn.backends import TorchVinecopBackend

    backend = TorchVinecopBackend()

  original, x = _fitted_estimator(backend)
  restored = pickle.loads(pickle.dumps(original))

  np.testing.assert_allclose(restored.pdf(x[:50]), original.pdf(x[:50]))
  assert type(restored.distribution_) is type(original.distribution_)
  assert [type(m) for m in restored.distribution_.margins] == [
    type(m) for m in original.distribution_.margins
  ]


def test_vinedist() -> None:
  """A whole distribution, both halves, through `pickle`.

  `Vinedist` is new in this release and was covered here only indirectly, via
  the sklearn estimator that holds one.
  """
  rng = np.random.default_rng(0)
  y = rng.normal(size=(300, 3)) + rng.normal(size=(300, 1))
  dist = pv.Vinedist.from_data(y)
  back = pickle.loads(pickle.dumps(dist))

  grid = y[:20]
  np.testing.assert_allclose(back.logpdf(grid), dist.logpdf(grid))
  np.testing.assert_allclose(back.pdf(grid), dist.pdf(grid))
  assert [type(m).__name__ for m in back.margins] == [
    type(m).__name__ for m in dist.margins
  ]
  assert np.array_equal(
    np.asarray(back.vinecop.structure.matrix),
    np.asarray(dist.vinecop.structure.matrix),
  )


def test_scipy_margin() -> None:
  """A parametric margin, including the criteria its payload has to carry."""
  scipy_margin = pytest.importorskip("pyvinecopulib.margins").SciPyMargin
  rng = np.random.default_rng(1)
  y = rng.normal(size=250)
  margin = scipy_margin("norm").fit(y)
  back = pickle.loads(pickle.dumps(margin))

  q = np.array([-1.0, 0.0, 1.0])
  np.testing.assert_allclose(back.pdf(q), margin.pdf(q))
  np.testing.assert_allclose(back.cdf(q), margin.cdf(q))
  assert back.family_name == margin.family_name
  for name in ("loglik", "aic", "bic", "aicc"):
    assert getattr(back, name)() == pytest.approx(getattr(margin, name)())


def test_torch_tll_bicop() -> None:
  """The renamed torch pair copula, with its cache mode and revision."""
  torch = pytest.importorskip("torch")
  from pyvinecopulib.torch import TorchTllBicop

  rng = np.random.default_rng(2)
  u = pv.to_pseudo_obs(rng.normal(size=(400, 2)) + rng.normal(size=(400, 1)))
  pair = TorchTllBicop.from_data(torch.from_numpy(u))
  back = pickle.loads(pickle.dumps(pair))

  ut = torch.from_numpy(u[:30])
  for name in ("pdf", "cdf", "hfunc1", "hfunc2"):
    torch.testing.assert_close(getattr(back, name)(ut), getattr(pair, name)(ut))


def test_torch_distribution_margin() -> None:
  """The renamed torch margin: parameters are registered, so they travel."""
  torch = pytest.importorskip("torch")
  from pyvinecopulib.torch import TorchDistributionMargin

  margin = TorchDistributionMargin.from_distribution(
    torch.distributions.Normal(loc=torch.tensor(0.5), scale=torch.tensor(2.0))
  )
  back = pickle.loads(pickle.dumps(margin))
  q = torch.tensor([-1.0, 0.0, 1.0], dtype=torch.float32)
  torch.testing.assert_close(back.pdf(q), margin.pdf(q))
  torch.testing.assert_close(back.cdf(q), margin.cdf(q))
