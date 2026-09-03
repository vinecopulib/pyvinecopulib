"""Python interface to vinecopulib.

Public API is organized into six subpackages: :mod:`~pyvinecopulib.core`,
:mod:`~pyvinecopulib.families`, :mod:`~pyvinecopulib.utils`,
:mod:`~pyvinecopulib.margins`, and
:mod:`~pyvinecopulib.sklearn` and :mod:`~pyvinecopulib.torch` (optional extras).
Core types are re-exported
at the top level; noisier aliases (family constants, utilities) still
resolve there but emit a ``DeprecationWarning`` on access.
"""

from ._cpu import require_x86_64_v3

require_x86_64_v3()

from . import core, families, margins, pyvinecopulib_ext, utils  # noqa: E402

# The CPU gate must run before importing the compiled extension.
from ._deprecations import _DEPRECATED_TOP_LEVEL, _resolve_deprecated  # noqa: E402
from .core import (  # noqa: E402
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
from .utils import to_pseudo_obs  # noqa: E402

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
