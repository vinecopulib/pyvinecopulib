"""Turn a user's ``margins=`` specification into one margin per variable."""

from __future__ import annotations

import copy
import operator
import warnings
from typing import Any, Callable, Optional, Sequence, Union

import numpy as _np

from ..core import MarginLike
from ._adapters import as_margin
from ..core import Kde1d

__all__ = [
  "resolve_margins",
  "resolve_margin_controls",
  "kde_from_controls",
  "fit_margin",
]


def _parametric_margin() -> Any:
  """Build a parametric margin with its family left to be selected.

  Deferred so that :mod:`pyvinecopulib.margins` imports without SciPy, which
  ``SciPyMargin`` needs and this module does not.

  Returns
  -------
  SciPyMargin
      An unfitted margin with no family yet.
  """
  from .scipy import SciPyMargin

  return SciPyMargin()


#: String aliases for the built-in margins. ``"parametric"`` names the class
#: without naming a family, which is what makes it an alias at all: choosing
#: the family from the data is `SciPyMargin.select`'s job.
_ALIASES = {
  "kde": Kde1d,
  "parametric": _parametric_margin,
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


#: What an unaddressed variable gets: one entry per variable, or a callable
#: producing them, invoked only if some variable is in fact unaddressed.
_Default = Union[
  Sequence[Any],
  Callable[[], Optional[Sequence[Any]]],
]


def _per_variable(
  spec: Any,
  d: int,
  *,
  names: Optional[Sequence[str]],
  default: Any,
  label: str,
  is_atom: Callable[[Any], bool],
) -> list[Any]:
  """Expand a per-variable specification, the one rule every layer follows.

  Recognized in this order: ``None`` (every variable gets ``default``), a
  mapping keyed by variable name or position over that default, a sequence of
  length ``d``, or a single value broadcast to every variable.

  Parameters
  ----------
  spec : object
      The user's specification.
  d : int
      Number of variables.
  names : sequence of str, or None
      Variable names, needed only to resolve a mapping keyed by name.
  default : object
      What an unaddressed variable gets.
  label : str
      The argument's name, quoted in the error messages.
  is_atom : callable
      ``value -> bool``, true for a single value that must not be read as a
      sequence of per-variable values.

  Returns
  -------
  list
      ``d`` entries.

  Raises
  ------
  ValueError
      If a sequence has the wrong length, or a mapping names an unknown
      variable or an out-of-range position.
  """
  if spec is None:
    return [default] * d

  if isinstance(spec, dict) and not is_atom(spec):
    resolved: list[Any] = [default] * d
    lookup = {name: j for j, name in enumerate(names or [])}
    for key, value in spec.items():
      if isinstance(key, str):
        if key not in lookup:
          raise ValueError(
            f"{label} mapping names {key!r}, which is not a variable"
            + (
              f"; known names are {sorted(lookup)}"
              if lookup
              else "; no names are known here, so key the mapping by position"
            )
          )
        index = lookup[key]
      else:
        try:
          index = operator.index(key)
        except TypeError as e:
          raise ValueError(
            f"{label} mapping key {key!r} is neither a variable name nor an "
            "integer position"
          ) from e
        if not 0 <= index < d:
          raise ValueError(f"{label} mapping has out-of-range index {index}")
      resolved[index] = value
    return resolved

  if isinstance(spec, (list, tuple)) and not is_atom(spec):
    if len(spec) != d:
      raise ValueError(
        f"{label} has length {len(spec)}, but there are {d} variables"
      )
    return list(spec)

  return [spec] * d


def resolve_margin_controls(
  spec: Any,
  d: int,
  *,
  names: Optional[Sequence[str]] = None,
) -> list[Any]:
  """Expand ``margin_controls=`` into one controls object per variable.

  Follows the same rules as ``margins=`` (see :func:`resolve_margins`): one
  controls object broadcast to every variable, a length-``d`` sequence, or a
  mapping keyed by variable name or position with the unaddressed variables
  left unconfigured. That is what lets one call give a bound to the two
  variables that need one and leave the rest alone.

  A broadcast object is **shared**, not copied: controls are read during a fit,
  never written to, so there is no state to leak between variables.

  Parameters
  ----------
  spec : object
      The user's ``margin_controls=`` argument.
  d : int
      Number of variables.
  names : sequence of str, or None, optional
      Variable names, needed only to resolve a mapping keyed by name.

  Returns
  -------
  list
      ``d`` entries, each a controls object or ``None``.

  Raises
  ------
  ValueError
      If a sequence has the wrong length, or a mapping names an unknown
      variable.

  See Also
  --------
  resolve_margins : The same rules, for the margins themselves.
  """
  return _per_variable(
    spec,
    d,
    names=names,
    default=None,
    label="margin_controls",
    # A controls object is a single value even though it may be a mapping or a
    # sequence in principle: `to_dict` is what makes it one.
    is_atom=lambda value: hasattr(value, "to_dict"),
  )


def resolve_margins(
  spec: Any,
  d: int,
  *,
  names: Optional[Sequence[str]] = None,
  default: Optional[_Default] = None,
) -> list[Any]:
  """Expand ``spec`` into one specification per variable.

  Accepted forms, in the order they are recognized: ``None`` (the default
  margin), a string alias (``"kde"``, ``"parametric"``), a
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
  default : sequence, callable, or None, optional
      What an unaddressed variable gets, one entry per variable. ``None`` means
      a kernel-density margin throughout. Worth setting whenever the caller
      knows something per variable that the library cannot: the sklearn
      estimators pass the variable types and bounds they inferred from the
      data, so a mapping that addresses one column does not silently retype the
      others. A callable is invoked only if some variable is actually
      unaddressed, which matters when building the default is expensive or can
      fail -- a specification that names every variable should not be held
      hostage to a default it never uses.

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

  def base_specs() -> list[Any]:
    if default is None:
      return [_prototype("kde") for _ in range(d)]
    computed = default() if isinstance(default, Callable) else default
    if computed is None:
      # A callable that answers `None` says what `default=None` says: this
      # module's own default. Deferred precisely so a specification naming
      # every variable is not held hostage to a default it never uses -- one
      # that can legitimately raise.
      return [_prototype("kde") for _ in range(d)]
    entries: Sequence[Any] = computed
    if len(entries) != d:
      raise ValueError(
        f"default has length {len(entries)}, but there are {d} variables"
      )
    return [_resolve_one(entry) for entry in entries]

  if spec is None:
    return base_specs()

  if isinstance(spec, dict):
    lookup = {name: j for j, name in enumerate(names or [])}
    addressed: dict[int, Any] = {}
    for key, value in spec.items():
      if isinstance(key, str):
        if key not in lookup:
          raise ValueError(
            f"margins mapping names {key!r}, which is not a variable"
            + (
              f"; known names are {sorted(lookup)}"
              if lookup
              else "; no names are known here, so key the mapping by position"
            )
          )
        index = lookup[key]
      else:
        try:
          index = operator.index(key)
        except TypeError as e:
          raise ValueError(
            f"margins mapping key {key!r} is neither a variable name nor an "
            "integer position"
          ) from e
        if not 0 <= index < d:
          raise ValueError(f"margins mapping has out-of-range index {index}")
      addressed[index] = _resolve_one(value)
    # Resolve the keys first, then build the default only if a variable is
    # actually left over -- the deferral the `default` parameter documents, and
    # the reason it matters is that building one can raise.
    if len(addressed) == d:
      return [addressed[j] for j in range(d)]
    resolved: list[Any] = list(base_specs())
    for index, value in addressed.items():
      resolved[index] = value
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


def kde_from_controls(controls: Optional[Any]) -> Kde1d:
  """Build a kernel-density margin honoring what the controls declare.

  The declaration is authoritative here, unlike on a margin the caller
  constructed: this margin does not exist yet, so there is nothing for it to
  override. It is what makes a bounded default reachable without naming a
  class -- ``margin_controls`` alone can say that one variable is positive and
  another is a proportion.

  Parameters
  ----------
  controls : FitControlsMargin, or None
      Read for ``var_type`` and ``support``; both optional.

  Returns
  -------
  Kde1d
      An unfitted margin, bounded and typed as declared.
  """
  kwargs: dict[str, Any] = {}
  var_type = getattr(controls, "var_type", None)
  if var_type == "d":
    kwargs["type"] = "discrete"
  elif var_type == "zi":
    kwargs["type"] = "zero_inflated"
  support = getattr(controls, "support", None)
  if support is not None:
    lo, hi = support
    if lo is not None and _np.isfinite(lo):
      kwargs["xmin"] = float(lo)
    if hi is not None and _np.isfinite(hi):
      kwargs["xmax"] = float(hi)
  return Kde1d(**kwargs)


def _fallback_kde(
  controls: Optional[Any], y: Any, weights: Optional[Any]
) -> Kde1d:
  """Fit the kernel-density margin substituted for a failed candidate set.

  Parameters
  ----------
  controls : FitControlsMargin, or None
      Read for the declared type and bounds.
  y : array, shape (n,), dtype float
      The column.
  weights : array, shape (n,), or None
      Observation weights.

  Returns
  -------
  Kde1d
      The fitted margin.
  """
  margin = kde_from_controls(controls)
  return margin.fit(y) if weights is None else margin.fit(y, weights)


def fit_margin(
  entry: Any,
  y: Any,
  *,
  x: Optional[Any] = None,
  weights: Optional[Any] = None,
  var_type: Optional[str] = None,
  support: Optional[tuple[float, float]] = None,
  controls: Optional[Any] = None,
  verb: str = "select",
  refit: bool = False,
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
      Observation weights. A callable specification receives them by keyword.
  var_type : str or None, optional
      The variable type the caller resolved, handed to a margin that
      implements ``declare``. Without it such a margin re-infers the type
      from the sample, which knows less than the caller does.
  support : tuple of float, or None, optional
      The declared bounds, handed over the same way.
  controls : object, or None, optional
      Fit configuration, forwarded to the margin's estimator.
  verb : str, optional
      Which estimator to call, ``"select"`` (the default, so a margin that
      searches a family set does) or ``"fit"`` (the current family only).
      ``MarginBase.select`` reduces to ``fit`` where there is nothing to
      choose, so the default is the weaker requirement.
  refit : bool, optional
      Re-estimate a margin that reports itself already fitted. Off by default,
      so a specification may mix fixed margins with ones to estimate.

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
    kwargs: dict[str, Any] = {}
    if x is not None:
      kwargs["x"] = x
    if weights is not None:
      kwargs["weights"] = weights
    return as_margin(entry(y, **kwargs))

  # Capability-based dispatch: `fit` / `is_fitted` / `supports_weights` are
  # optional members, so this is deliberately not narrowed to `MarginLike`.
  margin: Any = as_margin(entry)
  if not refit and getattr(margin, "is_fitted", True):
    return margin

  if weights is not None and not getattr(margin, "supports_weights", False):
    raise TypeError(
      f"{type(margin).__name__} cannot use observation weights; pass "
      "margins='kde' for a weighted fit, or drop weights="
    )
  kwargs: dict[str, Any] = {}
  if weights is not None:
    kwargs["weights"] = weights
  if x is not None and getattr(margin, "supports_covariates", False):
    kwargs["x"] = x
  declare = getattr(margin, "declare", None)
  if declare is not None and (var_type is not None or support is not None):
    declare(var_type=var_type, support=support)
  estimator = getattr(margin, verb, None) or margin.fit
  if controls is not None:
    if getattr(margin, "supports_controls", False):
      kwargs["controls"] = controls
    elif getattr(controls, "family_set", None) is not None:
      # A declared type or support is a *default*, so an estimator that cannot
      # read one has usually already been built with it. A family_set is not a
      # default -- it is an instruction to search -- so a margin that cannot
      # search must say so rather than fit one family and look like it chose.
      raise TypeError(
        f"{type(margin).__name__} cannot select a family, so family_set= "
        "would be ignored. Pass margins='parametric' (or a margin with a "
        "select method), or drop family_set"
      )
  try:
    estimator(y, **kwargs)
  except ValueError as e:
    # Substituting a different kind of margin is a decision about which margin
    # this column gets, so it is made here rather than inside a margin that
    # would have to stop being itself to make it.
    if getattr(controls, "on_failure", "raise") != "fallback":
      raise
    warnings.warn(
      f"{type(margin).__name__} could not fit this variable, so a "
      f"kernel-density margin was substituted. The cause was:\n{e}",
      UserWarning,
      stacklevel=2,
    )
    return as_margin(_fallback_kde(controls, y, weights))
  return margin
