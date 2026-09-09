"""Python interface to vinecopulib.

Public API is organized into six subpackages: :mod:`~pyvinecopulib.core`,
:mod:`~pyvinecopulib.families`, :mod:`~pyvinecopulib.utils`,
:mod:`~pyvinecopulib.margins`, and
:mod:`~pyvinecopulib.sklearn` and :mod:`~pyvinecopulib.torch` (optional extras).
Core types are re-exported
at the top level; noisier aliases (family constants, utilities) still
resolve there but emit a ``DeprecationWarning`` on access.
"""

from typing import Any

from ._cpu import require_x86_64_v3

require_x86_64_v3()

from . import core, families, margins, pyvinecopulib_ext, utils  # noqa: E402

# The CPU check must run before importing the compiled extension.
from ._deprecations import _DEPRECATED_TOP_LEVEL, _resolve_deprecated  # noqa: E402
from .core import (  # noqa: E402
  Bicop,
  BicopFamily,
  CVineStructure,
  DVineStructure,
  FitControlsBicop,
  FitControlsMargin,
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
  "FitControlsMargin",
  "FitControlsVinecop",
  "RVineStructure",
  "Vinecop",
  "Vinedist",
  "to_pseudo_obs",
  "core",
  "families",
  "margins",
  "utils",
  "__version__",
]

#: Subpackages that need an optional dependency, reachable by attribute access
#: and by `import pyvinecopulib.<name>` -- but **out** of
#: `__all__`, because `from pyvinecopulib import *` resolves every name in it
#: and would then require every extra. `margins` is in `__all__` instead: it
#: imports with no extra, deferring SciPy to the margin class that needs it.
_LAZY_SUBPACKAGES = ("sklearn", "torch")


def __getattr__(name: str) -> Any:  # noqa: ANN401 - a name resolves to any
  # class, function or subpackage, which is what `Any` means here.
  if name in _DEPRECATED_TOP_LEVEL:
    return _resolve_deprecated(name)
  if name in _LAZY_SUBPACKAGES:
    # Lazy: attribute access is what triggers the extra, so neither importing
    # `pyvinecopulib` nor star-importing from it requires either one.
    import importlib

    return importlib.import_module(f"pyvinecopulib.{name}")
  raise AttributeError(f"module 'pyvinecopulib' has no attribute {name!r}")


def __dir__() -> list[str]:
  # The lazy subpackages belong here even though they are not in `__all__`:
  # `dir()` is discovery, which should name them, while `__all__` is what a
  # star-import binds, which must not require an extra.
  return sorted(
    set(__all__) | set(_LAZY_SUBPACKAGES) | set(_DEPRECATED_TOP_LEVEL)
  )
