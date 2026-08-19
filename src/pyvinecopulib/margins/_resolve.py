"""Turn a user's ``margins=`` specification into one margin per variable."""

from __future__ import annotations

import copy
from typing import Any, Optional, Sequence

from ..core import MarginLike
from ._adapters import as_margin
from .empirical import EmpiricalMargin
from ..core import Kde1d
from .selection import MarginSelector

__all__ = ["resolve_margins", "fit_margin"]

#: String aliases for the built-in margin families. ``"parametric"`` is the
#: selector rather than one family: a parametric margin is only meaningful once
#: a family has been chosen, and choosing it from the data is the whole job.
_ALIASES = {
  "kde": Kde1d,
  "empirical": EmpiricalMargin,
  "parametric": MarginSelector,
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
  default: Optional[Sequence[Any]] = None,
) -> list[Any]:
  """Expand ``spec`` into one specification per variable.

  Accepted forms, in the order they are recognized: ``None`` (the default
  margin), a string alias (``"kde"``, ``"empirical"``, ``"parametric"``), a
  mapping keyed by variable name or position over the default, a sequence of
  length ``d`` mixing any of the above, a single margin broadcast to every
  variable, or a callable receiving a column and returning a margin.

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
  default : sequence, or None, optional
      What an unaddressed variable gets, one entry per variable. ``None`` means
      a kernel-density margin throughout. Worth setting whenever the caller
      knows something per variable that the library cannot: the sklearn
      estimators pass the variable types and bounds they inferred from the
      data, so a mapping that addresses one column does not silently retype the
      others.

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
  if default is None:
    base = [_prototype("kde") for _ in range(d)]
  elif len(default) != d:
    raise ValueError(
      f"default has length {len(default)}, but there are {d} variables"
    )
  else:
    base = [_resolve_one(entry) for entry in default]

  if spec is None:
    return base

  if isinstance(spec, dict):
    resolved: list[Any] = list(base)
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
  entry: Any, y: Any, *, x: Optional[Any] = None, weights: Optional[Any] = None
) -> MarginLike[Any]:
  """Obtain a fitted margin for one column.

  Parameters
  ----------
  entry : object
      A specification produced by :func:`resolve_margins`.
  y : array, shape (n,), dtype float
      The column.
  x : array, shape (n, k), or None, optional
      Exogenous covariates, forwarded only to a margin that declares
      ``supports_covariates``; a callable specification is handed them by
      keyword, so one that models no covariates fails loudly rather than
      fitting the wrong model.
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
    return as_margin(entry(y) if x is None else entry(y, x=x))

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
  kwargs: dict[str, Any] = {}
  if weights is not None:
    kwargs["weights"] = weights
  if x is not None and getattr(margin, "supports_covariates", False):
    kwargs["x"] = x
  margin.fit(y, **kwargs)
  return margin
