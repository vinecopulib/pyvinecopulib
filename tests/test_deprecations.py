"""Tests for the `simulate` -> `sample` deprecation shims.

Three things are pinned. That each name released in 0.7.6 still resolves and still
warns, since a hard break there would strand existing callers. That the warning
message names both the old and the new path in full -- `pyproject.toml` promotes
warnings whose message starts with `` `pyvinecopulib. `` to errors, so the wording
is what keeps an internal caller from silently riding a shim. And that renaming
the RNG hook fails loudly for a third-party subclass, which is the one rename that
cannot fail visibly on its own.

Every test here runs inside `pytest.warns`, which re-arms `always` for its block
and so escapes that `error` filter; the redundant `filterwarnings` marks make the
intent explicit rather than relying on that.
"""

from typing import Any

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import BicopBase, RVineStructure, Vinecop, VinecopBase
from pyvinecopulib.utils import Kde1d


@pytest.fixture
def gaussian() -> pv.Bicop:
  """A fitted Gaussian pair copula."""
  return pv.Bicop.from_family(
    pv.BicopFamily.gaussian, parameters=np.array([[0.5]])
  )


# --- the five names that shipped in 0.7.6 ----------------------------------- #


@pytest.mark.filterwarnings("default:`pyvinecopulib.")
def test_bicop_simulate_warns_and_still_works(gaussian: pv.Bicop) -> None:
  """`Bicop.simulate` forwards to `sample`."""
  with pytest.warns(DeprecationWarning, match=r"Bicop\.simulate\(\)"):
    old = gaussian.simulate(8, seeds=[1])
  new = gaussian.sample(8, seeds=[1])
  np.testing.assert_array_equal(old, new)


@pytest.mark.filterwarnings("default:`pyvinecopulib.")
def test_vinecop_simulate_warns_and_still_works() -> None:
  """`Vinecop.simulate` forwards to `sample`."""
  u = np.random.default_rng(0).uniform(size=(80, 3))
  vine = pv.Vinecop.from_data(u)
  with pytest.warns(DeprecationWarning, match=r"Vinecop\.simulate\(\)"):
    old = vine.simulate(6, seeds=[2])
  np.testing.assert_array_equal(old, vine.sample(6, seeds=[2]))


@pytest.mark.filterwarnings("default:`pyvinecopulib.")
def test_rvinestructure_simulate_stays_static() -> None:
  """The alias for a `staticmethod` is reachable without an instance."""
  with pytest.warns(DeprecationWarning, match=r"RVineStructure\.simulate\(\)"):
    old = RVineStructure.simulate(4, seeds=[3])
  assert old.dim == RVineStructure.sample(4, seeds=[3]).dim


@pytest.mark.filterwarnings("default:`pyvinecopulib.")
def test_kde1d_simulate_warns_and_still_works() -> None:
  """`Kde1d.simulate` forwards to `sample`."""
  kde = Kde1d()
  kde.fit(np.random.default_rng(0).normal(size=200))
  with pytest.warns(DeprecationWarning, match=r"Kde1d\.simulate\(\)"):
    old = kde.simulate(5, seeds=[4])
  np.testing.assert_array_equal(old, kde.sample(5, seeds=[4]))


@pytest.mark.filterwarnings("default:`pyvinecopulib.")
def test_utils_simulate_uniform_warns_and_still_works() -> None:
  """The free function is served from a module `__getattr__`."""
  with pytest.warns(DeprecationWarning, match=r"utils\.simulate_uniform"):
    old = pv.utils.simulate_uniform(6, 2, seeds=[5])
  np.testing.assert_array_equal(old, pv.utils.sample_uniform(6, 2, seeds=[5]))


def test_deprecated_function_is_absent_from_the_public_surface() -> None:
  """The old free-function name is reachable but not advertised."""
  assert "sample_uniform" in pv.utils.__all__
  assert "simulate_uniform" not in pv.utils.__all__
  # `__dir__` still lists it, so tab completion finds the migration path.
  assert "simulate_uniform" in dir(pv.utils)


# --- what carries no alias -------------------------------------------------- #


def test_sample_conditional_has_no_alias() -> None:
  """It never shipped, so it is renamed outright rather than aliased."""
  assert hasattr(Vinecop, "sample_conditional")
  assert not hasattr(Vinecop, "simulate_conditional")


# --- the renamed RNG hook --------------------------------------------------- #


@pytest.mark.parametrize("base", [BicopBase, VinecopBase])
def test_overriding_the_old_rng_hook_name_is_rejected(base: Any) -> None:
  """A silent rename would leave the override ignored and the base raising."""
  with pytest.raises(TypeError, match=r"_simulate_uniform.*_sample_uniform"):

    class _Stale(base):  # type: ignore[misc, valid-type]
      def _simulate_uniform(self, n: int, qrng: bool, seeds: list[int]) -> Any:
        return np.zeros((n, 2))


@pytest.mark.parametrize("base", [BicopBase, VinecopBase])
def test_overriding_the_new_rng_hook_name_is_accepted(base: Any) -> None:
  """Defining only the current name is the supported case."""

  class _Fresh(base):  # type: ignore[misc, valid-type]
    def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> Any:
      return np.zeros((n, 2))

  assert hasattr(_Fresh, "_sample_uniform")


# --- the two coexisting `sample` conventions -------------------------------- #


def test_sklearn_sample_signature_is_pinned() -> None:
  """`check_estimator` requires `sample(n_samples, random_state)`; core takes `n`.

  The two conventions coexist deliberately. Pinned so a later "unification"
  cannot quietly break scikit-learn compatibility.
  """
  import inspect

  # Guard on the third-party package, as every other test here does:
  # `pyvinecopulib.sklearn` re-raises as a plain `ImportError`, which
  # `importorskip` does not treat as "absent" -- it only auto-skips
  # `ModuleNotFoundError`.
  pytest.importorskip("sklearn")
  from pyvinecopulib.sklearn import VineDensity

  params = list(inspect.signature(VineDensity.sample).parameters)
  assert params[:2] == ["self", "n_samples"]
  assert "random_state" in params
