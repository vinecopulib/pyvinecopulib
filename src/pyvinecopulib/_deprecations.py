"""Warn-on-access table for the top-level deprecation shim.

`pyvinecopulib.__getattr__` calls `_resolve_deprecated(name)`, which emits
a `DeprecationWarning` and returns the symbol from its canonical subpackage.
Slated for removal in the next major release.
"""

import importlib
import warnings
from typing import Any

# name -> (subpackage, attr_name)
_DEPRECATED_TOP_LEVEL: dict[str, tuple[str, str]] = {
  "indep": ("families", "indep"),
  "gaussian": ("families", "gaussian"),
  "student": ("families", "student"),
  "clayton": ("families", "clayton"),
  "gumbel": ("families", "gumbel"),
  "frank": ("families", "frank"),
  "joe": ("families", "joe"),
  "bb1": ("families", "bb1"),
  "bb6": ("families", "bb6"),
  "bb7": ("families", "bb7"),
  "bb8": ("families", "bb8"),
  "tawn": ("families", "tawn"),
  "tll": ("families", "tll"),
  "all": ("families", "all"),
  "parametric": ("families", "parametric"),
  "nonparametric": ("families", "nonparametric"),
  "one_par": ("families", "one_par"),
  "two_par": ("families", "two_par"),
  "three_par": ("families", "three_par"),
  "elliptical": ("families", "elliptical"),
  "archimedean": ("families", "archimedean"),
  "extreme_value": ("families", "extreme_value"),
  "bb": ("families", "bb"),
  "rotationless": ("families", "rotationless"),
  "lt": ("families", "lt"),
  "ut": ("families", "ut"),
  "itau": ("families", "itau"),
  "Kde1d": ("utils", "Kde1d"),
  "wdm": ("utils", "wdm"),
  "sobol": ("utils", "sobol"),
  "ghalton": ("utils", "ghalton"),
  "simulate_uniform": ("utils", "simulate_uniform"),
  "benchmark": ("utils", "benchmark"),
  "pairs_copula_data": ("utils", "pairs_copula_data"),
}


def _resolve_deprecated(name: str) -> Any:
  subpkg, attr = _DEPRECATED_TOP_LEVEL[name]
  warnings.warn(
    f"`pyvinecopulib.{name}` is deprecated; use "
    f"`pyvinecopulib.{subpkg}.{attr}` instead.",
    DeprecationWarning,
    stacklevel=3,
  )
  return getattr(importlib.import_module(f"pyvinecopulib.{subpkg}"), attr)
