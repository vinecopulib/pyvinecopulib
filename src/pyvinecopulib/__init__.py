"""Python interface to vinecopulib.

Public API is organized into five subpackages: :mod:`~pyvinecopulib.core`,
:mod:`~pyvinecopulib.families`, :mod:`~pyvinecopulib.utils`,
:mod:`~pyvinecopulib.margins`, and
:mod:`~pyvinecopulib.sklearn` (optional extra). Core types are re-exported
at the top level; noisier aliases (family constants, utilities) still
resolve there but emit a ``DeprecationWarning`` on access.
"""

from . import core, families, margins, pyvinecopulib_ext, utils
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
  Vinedist,
)
from .utils import to_pseudo_obs

__version__ = pyvinecopulib_ext.__version__

__all__ = [
  "Bicop",
  "BicopFamily",
  "CVineStructure",
  "DVineStructure",
  "FitControlsBicop",
  "FitControlsVinecop",
  "RVineStructure",
  "Vinecop",
  "Vinedist",
  "to_pseudo_obs",
  "core",
  "families",
  "margins",
  "utils",
  "sklearn",
  "__version__",
]


def __getattr__(name: str):
  if name in _DEPRECATED_TOP_LEVEL:
    return _resolve_deprecated(name)
  if name == "sklearn":
    # Lazy: only `import pyvinecopulib.sklearn` triggers the extra.
    import importlib

    return importlib.import_module("pyvinecopulib.sklearn")
  raise AttributeError(f"module 'pyvinecopulib' has no attribute {name!r}")


def __dir__():
  return sorted(set(__all__) | set(_DEPRECATED_TOP_LEVEL))
