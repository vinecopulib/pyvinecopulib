"""Scikit-learn compatible vine estimators.

Requires scikit-learn. Install with `pip install pyvinecopulib[sklearn]`.
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
