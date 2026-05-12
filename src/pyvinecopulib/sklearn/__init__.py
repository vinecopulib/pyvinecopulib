"""Scikit-learn-compatible vine-copula estimators.

This subpackage exposes two non-parametric estimators built on the
core :mod:`pyvinecopulib.core` and :mod:`pyvinecopulib.utils` machinery:

- :class:`VineDensity` — joint-density estimator.
- :class:`VineRegressor` — conditional mean / quantile regressor.

Both follow the standard pipeline of marginal kernel-density
estimation, transformation to pseudo-observations, and a vine-copula
fit on those pseudo-observations; see the class docstrings for the
full methodology and references.

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
