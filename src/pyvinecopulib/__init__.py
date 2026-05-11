"""pyvinecopulib — Python interface to vinecopulib.

Public API is organized into four subpackages: :mod:`~pyvinecopulib.core`,
:mod:`~pyvinecopulib.families`, :mod:`~pyvinecopulib.utils`, and
:mod:`~pyvinecopulib.sklearn` (optional extra). Core types are re-exported
at the top level; noisier aliases (family constants, utilities) still
resolve there but emit a ``DeprecationWarning`` on access.
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
    # Lazy: only `import pyvinecopulib.sklearn` should trigger the extra.
    import importlib

    return importlib.import_module("pyvinecopulib.sklearn")
  raise AttributeError(f"module 'pyvinecopulib' has no attribute {name!r}")


def __dir__():
  return sorted(set(__all__) | set(_DEPRECATED_TOP_LEVEL))
