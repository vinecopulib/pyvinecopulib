"""Deprecation table for module-level __getattr__ in pyvinecopulib/__init__.py.

When a user accesses a name in `_DEPRECATED_TOP_LEVEL` (e.g. `pyvinecopulib.gaussian`),
`_resolve_deprecated` emits a DeprecationWarning pointing at the new home
(e.g. `pyvinecopulib.families.gaussian`) and returns the resolved attribute.
Access via the canonical new path does NOT warn.

The intent is to phase out these top-level aliases in the next major release.
"""

import importlib
import warnings
from typing import Any

# name -> (subpackage, attr_name in subpackage)
_DEPRECATED_TOP_LEVEL: dict[str, tuple[str, str]] = {
  # Family constants -> families
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
  # Family groups -> families
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
  # Utilities -> utils. Kde1d is safe to deprecate: the C++ binding sets
  # Kde1d.__module__ = "pyvinecopulib.utils" so new pickles already use
  # the canonical path; only pre-1.0 pickles fall through this shim.
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
    stacklevel=3,  # skip __getattr__ and _resolve_deprecated frames
  )
  return getattr(importlib.import_module(f"pyvinecopulib.{subpkg}"), attr)
