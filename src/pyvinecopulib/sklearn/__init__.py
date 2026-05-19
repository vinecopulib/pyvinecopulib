"""Scikit-learn-compatible vine-copula estimators.

This subpackage exposes two non-parametric estimators that act as
thin scikit-learn wrappers around the core pyvinecopulib machinery:

- :class:`VineDensity` — joint-density estimator built on top of
  :class:`pyvinecopulib.core.Vinecop` and
  :class:`pyvinecopulib.utils.Kde1d`.
- :class:`VineRegressor` — conditional mean / quantile regressor on
  the same primitives.

The wrappers add the standard ``fit`` / ``predict``-style sklearn
interface, ``DataFrame`` input handling with auto-inferred
``continuous`` / ``discrete`` column types, and batched evaluation
for large prediction sets. All low-level fitting knobs (pair family,
threading, structure-selection algorithm, …) are passed through to
:class:`pyvinecopulib.core.FitControlsVinecop` and
:class:`pyvinecopulib.core.RVineStructure` — reach for those
directly whenever you need control beyond the sklearn convenience
layer.

Both estimators follow the standard pipeline of marginal
kernel-density estimation, transformation to pseudo-observations,
and a vine-copula fit on those pseudo-observations; see the class
docstrings for the full methodology and references.

Requires scikit-learn and pandas. Install with
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
from .regressor import VineRegressor

__all__ = ["VineDensity", "VineRegressor"]
