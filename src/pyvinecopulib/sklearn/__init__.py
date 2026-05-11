"""Scikit-learn compatible vine estimators.

Requires scikit-learn. Install with `pip install pyvinecopulib[sklearn]`.

This subpackage is a placeholder; estimator classes ship in a follow-up PR.
"""

try:
  import sklearn  # noqa: F401
except ImportError as e:
  raise ImportError(
    "pyvinecopulib.sklearn requires scikit-learn. "
    "Install it with `pip install pyvinecopulib[sklearn]`."
  ) from e

__all__: list[str] = []
