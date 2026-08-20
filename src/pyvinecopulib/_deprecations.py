"""Warn-on-access table for the top-level deprecation shim.

Slated for removal in the next major release.
"""

import functools
import importlib
import warnings
from typing import Any

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
  "simulate_uniform": ("utils", "sample_uniform"),
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


def _method_alias(new: Any, old_name: str, qualname: str) -> Any:
  """Build a deprecated alias for a renamed method.

  Parameters
  ----------
  new : callable
      The method under its current name.
  old_name : str
      The deprecated spelling.
  qualname : str
      The owner's path below ``pyvinecopulib``, e.g. ``"core.Vinecop"``.

  Returns
  -------
  callable
      A wrapper that warns and forwards. Its docstring's first line carries the
      *old* name with the new signature, because the stub generator recovers
      signatures -- and ``@staticmethod`` -- by parsing that line.
  """
  new_name = getattr(new, "__name__", "sample")

  # `functools.wraps` is not optional here: without it Sphinx refuses the alias
  # with "Cannot handle as a local function". It also sets `__wrapped__`, so
  # `inspect.signature` reports the forwarded method's real signature.
  @functools.wraps(new)
  def alias(*args: Any, **kwargs: Any) -> Any:
    warnings.warn(
      f"`pyvinecopulib.{qualname}.{old_name}()` is deprecated; use "
      f"`pyvinecopulib.{qualname}.{new_name}()` instead.",
      DeprecationWarning,
      stacklevel=2,
    )
    return new(*args, **kwargs)

  head, _, rest = (new.__doc__ or f"{new_name}(*args, **kwargs)").partition(
    "\n"
  )
  alias.__name__ = old_name
  alias.__doc__ = "\n".join(
    [
      head.replace(f"{new_name}(", f"{old_name}(", 1),
      "",
      f"Deprecated alias for ``{new_name}``.",
      rest,
    ]
  )
  # `functools.wraps` also set `__wrapped__`, which makes `inspect.signature`
  # follow through to the nanobind method -- and `sphinx_autodoc_typehints`
  # cannot introspect one ("no signature found for builtin"). The docstring head
  # already carries the signature, so drop the link. Bound through `Any`
  # because a function's dunder attributes are not expressible in its type.
  tagged: Any = alias
  del tagged.__wrapped__
  return alias


def _reject_renamed_hook(cls: type, old: str, new: str) -> None:
  """Fail loudly when a subclass overrides a hook under its former name.

  A renamed hook is the one rename that cannot fail visibly on its own: the base
  class simply stops calling the old name, so the override is ignored and the
  inherited default raises as though nothing were overridden.

  Parameters
  ----------
  cls : type
      The subclass being defined.
  old : str
      The former hook name.
  new : str
      The current hook name.

  Raises
  ------
  TypeError
      If ``cls`` defines ``old`` and not ``new``.
  """
  if old in vars(cls) and new not in vars(cls):
    raise TypeError(
      f"{cls.__name__} defines `{old}`, which is now `{new}`. Rename the "
      f"override: the base class no longer calls `{old}`, so sampling would "
      "raise NotImplementedError rather than use it."
    )
