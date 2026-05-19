"""Scikit-learn-compatible vine-copula estimators.

This subpackage exposes four non-parametric estimators that act as
thin scikit-learn wrappers around the core pyvinecopulib machinery:

- :class:`VineDensity` — single-vine joint-density estimator built on
  top of :class:`pyvinecopulib.core.Vinecop` and
  :class:`pyvinecopulib.utils.Kde1d`.
- :class:`VineRegressor` — single-vine conditional mean / quantile
  regressor on the same primitives.
- :class:`VineForestDensity` — ensemble of :class:`VineDensity`
  learners with random-search structure generation and MCS-based
  survivor selection.
- :class:`VineForestRegressor` — ensemble of :class:`VineRegressor`
  learners, same construction.

All four follow the standard pipeline of marginal kernel-density
estimation, transformation to pseudo-observations, and a vine-copula
fit on those pseudo-observations; the forest variants additionally
run a random search over vine structures and prune the candidate
pool via a model confidence set. The wrappers add the standard
``fit`` / ``predict``-style sklearn interface, ``DataFrame`` input
handling with auto-inferred ``continuous`` / ``discrete`` column
types, and batched evaluation. All low-level fitting knobs (pair
family, threading, structure-selection algorithm, …) are passed
through to :class:`pyvinecopulib.core.FitControlsVinecop` and
:class:`pyvinecopulib.core.RVineStructure` — reach for those
directly whenever you need control beyond the sklearn convenience
layer. See the class docstrings for the full methodology and
references.

Requires scikit-learn, pandas, joblib, and scipy. Install with
``pip install pyvinecopulib[sklearn]``.
"""

try:
  import sklearn  # noqa: F401
except ImportError as e:
  raise ImportError(
    "pyvinecopulib.sklearn requires scikit-learn. "
    "Install it with `pip install pyvinecopulib[sklearn]`."
  ) from e

from .density import VineDensity
from .forest_density import VineForestDensity
from .forest_regressor import VineForestRegressor
from .regressor import VineRegressor

__all__ = [
  "VineDensity",
  "VineForestDensity",
  "VineForestRegressor",
  "VineRegressor",
]
