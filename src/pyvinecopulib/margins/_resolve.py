"""Turn a user's ``margins=`` specification into one margin per variable."""

from __future__ import annotations

import copy
from typing import Any, Optional, Sequence

from ..core import MarginLike
from ._adapters import as_margin
from .empirical import EmpiricalMargin
from .kde import Kde1dMargin

__all__ = ["resolve_margins", "fit_margin"]

#: String aliases for the built-in margin families.
_ALIASES = {
  "kde": Kde1dMargin,
  "empirical": EmpiricalMargin,
}


def _prototype(alias: str) -> Any:
  """Build a margin from a string alias.

  Parameters
  ----------
  alias : str
      One of the keys of the alias table.

  Returns
  -------
  MarginBase
      A fresh, unfitted margin.

  Raises
  ------
  ValueError
      If the alias is unknown.
  """
  if alias not in _ALIASES:
    raise ValueError(
      f"unknown margins={alias!r}; expected one of {sorted(_ALIASES)}, "
      "a margin, a sequence of margins, a mapping, or a callable"
    )
  return _ALIASES[alias]()


def _resolve_one(entry: Any) -> Any:
  """Resolve a single column's specification.

  Parameters
  ----------
  entry : object
      A string alias, a margin, a foreign distribution, or a callable.

  Returns
  -------
  object
      Either something to fit, or a callable to call with the column.
  """
  if isinstance(entry, str):
    return _prototype(entry)
  if callable(entry) and not hasattr(entry, "cdf"):
    # A fitter: it receives the column and returns a margin.
    return entry
  return entry


def resolve_margins(
  spec: Any,
  d: int,
  *,
  names: Optional[Sequence[str]] = None,
) -> list[Any]:
  """Expand ``spec`` into one specification per variable.

  Accepted forms, in the order they are recognized: ``None`` (the default
  kernel-density margin), a string alias (``"kde"``, ``"empirical"``), a mapping
  keyed by variable name or position over a default, a sequence of length ``d``
  mixing any of the above, a single margin broadcast to every variable, or a
  callable receiving a column and returning a margin.

  A broadcast margin is **copied** per variable rather than shared: margins are
  mutable, and one estimator fitted repeatedly would carry state between
  variables.

  Parameters
  ----------
  spec : object
      The user's ``margins=`` argument.
  d : int
      Number of variables.
  names : sequence of str, or None, optional
      Variable names, needed only to resolve a mapping keyed by name.

  Returns
  -------
  list
      ``d`` specifications, each either fittable or callable.

  Raises
  ------
  ValueError
      If a sequence has the wrong length, or a mapping names an unknown
      variable.
  """
  if spec is None:
    spec = "kde"

  if isinstance(spec, dict):
    resolved: list[Any] = [_prototype("kde") for _ in range(d)]
    lookup = {name: j for j, name in enumerate(names or [])}
    for key, value in spec.items():
      if isinstance(key, str):
        if key not in lookup:
          raise ValueError(
            f"margins mapping names {key!r}, which is not a variable"
            + (f"; known names are {sorted(lookup)}" if lookup else "")
          )
        index = lookup[key]
      else:
        index = int(key)
        if not 0 <= index < d:
          raise ValueError(f"margins mapping has out-of-range index {index}")
      resolved[index] = _resolve_one(value)
    return resolved

  if isinstance(spec, str):
    return [_prototype(spec) for _ in range(d)]

  # A sequence, but not a margin (which may itself be iterable in principle).
  if isinstance(spec, (list, tuple)) and not hasattr(spec, "cdf"):
    if len(spec) != d:
      raise ValueError(
        f"margins has length {len(spec)}, but there are {d} variables"
      )
    return [_resolve_one(entry) for entry in spec]

  one = _resolve_one(spec)
  # Copy per variable: a shared mutable estimator would leak state.
  return [one if callable(one) else copy.deepcopy(one) for _ in range(d)]


def fit_margin(
  entry: Any, x: Any, *, weights: Optional[Any] = None
) -> MarginLike[Any]:
  """Obtain a fitted margin for one column.

  Parameters
  ----------
  entry : object
      A specification produced by :func:`resolve_margins`.
  x : array, shape (n,), dtype float
      The column.
  weights : array, shape (n,), or None, optional
      Observation weights.

  Returns
  -------
  MarginLike
      A fitted margin.

  Raises
  ------
  TypeError
      If weights were supplied for a margin that cannot use them, since
      ignoring them would silently fit a different model than was asked for.
  """
  if callable(entry) and not hasattr(entry, "cdf"):
    return as_margin(entry(x))

  # Capability-based dispatch: `fit` / `is_fitted` / `supports_weights` are
  # optional members, so this is deliberately not narrowed to `MarginLike`.
  margin: Any = as_margin(entry)
  if getattr(margin, "is_fitted", True):
    return margin

  if weights is not None and not getattr(margin, "supports_weights", False):
    raise TypeError(
      f"{type(margin).__name__} cannot use observation weights; pass "
      "margins='kde' or margins='empirical' for a weighted fit, or drop "
      "weights="
    )
  if weights is None:
    margin.fit(x)
  else:
    margin.fit(x, weights=weights)
  return margin
