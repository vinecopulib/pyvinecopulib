"""pyvinecopulib — Python interface to vinecopulib.

The public API is organized into four subpackages:

- :mod:`pyvinecopulib.core` — `Bicop`, `Vinecop`, fit controls, vine structures.
- :mod:`pyvinecopulib.families` — `BicopFamily` and family constants/groups.
- :mod:`pyvinecopulib.utils` — `Kde1d`, `wdm`, `to_pseudo_obs`, `pairs_copula_data`, etc.
- :mod:`pyvinecopulib.sklearn` — sklearn-compatible vine estimators (optional extra).

Core types (``Bicop``, ``Vinecop``, structures, controls, ``BicopFamily``,
``to_pseudo_obs``) are also re-exported at the top level for convenience. The
noisier aliases (family constants like ``pyvinecopulib.gaussian`` and utilities
like ``pyvinecopulib.Kde1d``) still resolve at the top level for backward
compatibility but emit a ``DeprecationWarning`` on access.
"""

from . import core, families, pyvinecopulib_ext, utils
from ._deprecations import _DEPRECATED_TOP_LEVEL, _resolve_deprecated
from .core import (
  Bicop,
  BicopFamily,
  CVineStructure,
  DVineStructure,
  FitControlsBicop,
  FitControlsVinecop,
  RVineStructure,
  Vinecop,
)
from .utils import to_pseudo_obs

__version__ = pyvinecopulib_ext.__version__

__all__ = [
  # Core types kept at the top level indefinitely
  "Bicop",
  "BicopFamily",
  "CVineStructure",
  "DVineStructure",
  "FitControlsBicop",
  "FitControlsVinecop",
  "RVineStructure",
  "Vinecop",
  "to_pseudo_obs",
  # Subpackages
  "core",
  "families",
  "utils",
  "sklearn",
  # Version
  "__version__",
]


def __getattr__(name: str):
  if name in _DEPRECATED_TOP_LEVEL:
    return _resolve_deprecated(name)
  if name == "sklearn":
    # Lazy import so that vanilla `import pyvinecopulib` doesn't require
    # scikit-learn — only `import pyvinecopulib.sklearn` does.
    import importlib

    return importlib.import_module("pyvinecopulib.sklearn")
  raise AttributeError(f"module 'pyvinecopulib' has no attribute {name!r}")


def __dir__():
  return sorted(set(__all__) | set(_DEPRECATED_TOP_LEVEL))
