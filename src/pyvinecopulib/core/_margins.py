"""Everything about a margin except its contract: coercion, resolution and
JSON persistence.

Three questions about a margin whose type is not known in advance, answered
the same way -- one table, resolved when it is needed:

* **Coercion.** :func:`as_margin` presents a distribution object from another
  ecosystem as a :class:`~pyvinecopulib.core.MarginLike`. Every margin the
  library receives passes through it, so a discrete SciPy object cannot slip
  past on a bare ``pdf`` (in SciPy's modern API that is ``+inf`` at an atom;
  the mass lives on ``pmf``).
* **Persistence.** :func:`margin_to_json` and :func:`margin_from_json` round a
  margin through a mapping carrying a ``kind`` naming its type and a
  ``version`` a reader checks, so a format change fails loudly rather than
  producing a wrong model.

* **Resolution.** :func:`resolve_margin_controls` expands ``margin_controls=``
  -- one controls object broadcast per variable, a length-``d`` sequence, or a
  mapping keyed by position or name -- into one entry per variable.

Each registry has a table of the margins this package ships and a registration
function for everything else (:func:`register_margin_adapter`,
:func:`register_margin_json`). A built-in entry imports lazily, because two of
the classes live behind an extra. Nothing registers itself from its own module,
so whether an object is recognized does not depend on some other module's
import list.

Everything here is about a margin whose type the caller chose and this library
did not. The contract itself is :mod:`~pyvinecopulib.core.margin_base`, and how
a payload is encoded is :mod:`~pyvinecopulib.core._json`.
"""

from __future__ import annotations

import copy
import operator
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

from ..pyvinecopulib_ext import Kde1d
from ._json import MODEL_JSON_VERSION, read_payload
from .margin_base import MarginBase, support_of
from .protocols import ArrayT, MarginLike

__all__ = [
  "as_margin",
  "margin_from_json",
  "margin_to_json",
  "register_margin_adapter",
  "register_margin_json",
  "resolve_margin_controls",
]


#: Adapters registered by a caller, newest first so a later registration can
#: take precedence over an earlier one -- and over ``_BUILTIN_ADAPTERS``, which
#: is consulted after this.
_ADAPTERS: list[
  tuple[Callable[[Any], bool], Callable[[Any], MarginLike[Any]]]
] = []


def register_margin_adapter(
  predicate: Callable[[Any], bool],
  adapter: Callable[[Any], MarginLike[Any]],
) -> None:
  """Teach :func:`as_margin` about another kind of distribution object.

  Parameters
  ----------
  predicate : callable
      Returns ``True`` for objects this adapter handles. It must not raise for
      unrelated objects, and it should avoid importing heavy modules until it
      has cheap evidence that the object belongs to them.
  adapter : callable
      Wraps a matching object into something satisfying
      :class:`~pyvinecopulib.core.MarginLike`.

  Returns
  -------
  None
  """
  _ADAPTERS.insert(0, (predicate, adapter))


class _WrappedMargin(MarginBase[ArrayT]):
  """A foreign distribution presented as a :class:`MarginBase`.

  Holds the wrapped object and the handful of callables that differ between
  ecosystems, so each adapter below is a declaration rather than code.

  Parameters
  ----------
  obj : object
      The wrapped distribution, kept reachable as ``wrapped``.
  pdf, cdf, icdf : callable
      The three primitives, already bound to ``obj`` and named as the
      :class:`~pyvinecopulib.core.MarginLike` contract expects.
  logpdf : callable, or None, optional
      Native log-density; ``None`` derives it from ``pdf``.
  var_type : str, default='c'
      Variable type of the wrapped distribution.
  cdf_left : callable, or None, optional
      Left-limit cdf; ``None`` derives it from ``var_type``.
  support : tuple of float, or None, optional
      Support bounds; ``None`` reads them off ``obj``.
  family_name : str, or None, optional
      Name to report in selection output; ``None`` uses the wrapped type's.
  """

  def __init__(
    self,
    obj: object,
    *,
    pdf: Callable[[ArrayT], ArrayT],
    cdf: Callable[[ArrayT], ArrayT],
    icdf: Callable[[ArrayT], ArrayT],
    logpdf: Callable[[ArrayT], ArrayT] | None = None,
    var_type: str = "c",
    cdf_left: Callable[[ArrayT], ArrayT] | None = None,
    support: tuple[float, float] | None = None,
    family_name: str | None = None,
  ) -> None:
    self.wrapped = obj
    self._pdf = pdf
    self._cdf = cdf
    self._icdf = icdf
    self._logpdf = logpdf
    self._var_type = var_type
    self._cdf_left = cdf_left
    self._support = support if support is not None else support_of(obj)
    self.family_name = family_name or type(obj).__name__

  @property
  def var_type(self) -> str:
    return self._var_type

  @property
  def support(self) -> tuple[float, float]:
    return self._support

  def pdf(self, y: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    return self._pdf(y)

  def logpdf(self, y: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    if self._logpdf is None:
      return super().logpdf(y)
    return self._logpdf(y)

  def cdf(self, y: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    return self._cdf(y)

  def icdf(self, p: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    return self._icdf(p)

  def cdf_left(self, y: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    if self._cdf_left is None:
      return super().cdf_left(y)
    return self._cdf_left(y)

  def __repr__(self) -> str:
    return f"as_margin({self.family_name}, var_type={self._var_type!r})"


def _is_scipy_new(obj: object) -> bool:
  """Whether ``obj`` is one of SciPy's modern distribution objects."""
  return any(
    base.__name__
    in (
      "ContinuousDistribution",
      "DiscreteDistribution",
      "UnivariateDistribution",
    )
    for base in type(obj).__mro__
  )


def _is_scipy_new_discrete(obj: object) -> bool:
  """Whether ``obj`` is a modern SciPy *discrete* distribution."""
  return any(
    base.__name__ == "DiscreteDistribution" for base in type(obj).__mro__
  )


def _adapt_scipy_new(obj: Any) -> MarginLike[np.ndarray]:  # noqa: ANN401
  """Adapt a modern SciPy distribution.

  A continuous one already matches the contract; a discrete one does not, and
  the mismatch is silent: ``pdf`` is the Lebesgue density, which is ``+inf`` at
  an atom, while the mass lives on ``pmf``.
  """
  if not _is_scipy_new_discrete(obj):
    return _WrappedMargin(
      obj,
      pdf=obj.pdf,
      logpdf=obj.logpdf,
      cdf=obj.cdf,
      icdf=obj.icdf,
      var_type="c",
    )
  return _WrappedMargin(
    obj,
    pdf=obj.pmf,
    logpdf=obj.logpmf,
    cdf=obj.cdf,
    icdf=obj.icdf,
    var_type="d",
    # `cdf` interpolates between atoms here, so stepping back a lattice point
    # would read a value strictly inside the jump. Subtract the mass instead.
    cdf_left=lambda x: obj.cdf(x) - obj.pmf(x),
  )


def _is_scipy_legacy(obj: object) -> bool:
  """Whether ``obj`` is a frozen legacy SciPy distribution."""
  dist = getattr(obj, "dist", None)
  return dist is not None and hasattr(obj, "ppf") and hasattr(dist, "name")


def _adapt_scipy_legacy(obj: Any) -> MarginLike[np.ndarray]:  # noqa: ANN401
  """Adapt a frozen legacy SciPy distribution (``ppf`` -> ``icdf``)."""
  discrete = hasattr(obj, "pmf") and not hasattr(obj, "pdf")
  name = getattr(obj.dist, "name", type(obj).__name__)
  if not discrete:
    return _WrappedMargin(
      obj,
      pdf=obj.pdf,
      logpdf=obj.logpdf,
      cdf=obj.cdf,
      icdf=obj.ppf,
      var_type="c",
      family_name=name,
    )
  return _WrappedMargin(
    obj,
    pdf=obj.pmf,
    logpdf=obj.logpmf,
    cdf=obj.cdf,
    icdf=obj.ppf,
    var_type="d",
    # This works on shifted and irregular numeric lattices as well as counts:
    # away from an atom `pmf` is zero, so the left limit equals `cdf`.
    cdf_left=lambda x: obj.cdf(x) - obj.pmf(x),
    family_name=name,
  )


def _is_torch_distribution(obj: object) -> bool:
  """Whether ``obj`` is a ``torch.distributions.Distribution``."""
  return any(
    base.__module__.startswith("torch.distributions")
    and base.__name__ == "Distribution"
    for base in type(obj).__mro__
  )


def _adapt_torch(obj: Any) -> MarginLike[Any]:  # noqa: ANN401
  """Adapt a ``torch.distributions`` object.

  ``log_prob`` is the only density it offers, and ``cdf`` / ``icdf`` are
  declared on the base class and raise unless the concrete family implements
  them — so conformance cannot be inferred from member names here. Continuous
  families without a cdf and every discrete family are rejected immediately;
  a missing ``icdf`` falls back to inverting an implemented ``cdf`` numerically.
  """
  support = getattr(obj, "support", None)
  if bool(getattr(support, "is_discrete", False)):
    raise TypeError(
      f"cannot adapt discrete torch distribution {type(obj).__name__!r}: "
      "torch.distributions does not provide the cdf and left-limit cdf a "
      "discrete vine margin needs. Use Kde1d, a SciPy or OpenTURNS margin, "
      "or implement MarginBase.cdf_left explicitly."
    )

  cdf = getattr(type(obj), "cdf", None)
  if cdf is None or (
    getattr(cdf, "__module__", "") == "torch.distributions.distribution"
    and getattr(cdf, "__qualname__", "") == "Distribution.cdf"
  ):
    raise TypeError(
      f"cannot adapt torch distribution {type(obj).__name__!r}: its cdf is "
      "not implemented. Use a continuous torch distribution with a cdf, "
      "provide a MarginBase implementation, or use a SciPy/OpenTURNS margin."
    )

  lo, hi = support_of(obj)

  def _icdf(p: Any) -> Any:  # noqa: ANN401 - a tensor, which `core` cannot name
    try:
      return obj.icdf(p)
    except NotImplementedError:
      from ..core._rootfind import solve_increasing

      return solve_increasing(obj.cdf, p, lo=lo, hi=hi)

  return _WrappedMargin(
    obj,
    pdf=lambda x: obj.log_prob(x).exp(),
    logpdf=obj.log_prob,
    cdf=obj.cdf,
    icdf=_icdf,
    var_type="c",
    support=(lo, hi),
  )


def as_margin(obj: object) -> MarginLike[Any]:
  """Present ``obj`` as a :class:`~pyvinecopulib.core.MarginLike`.

  Idempotent: anything this library produced, or a structural implementation of
  ``MarginLike``, is returned unchanged. Recognized foreign objects are checked
  and wrapped first, because satisfying the contract's member *names* is not the
  same as satisfying its semantics — a modern SciPy discrete distribution has
  ``pdf`` / ``cdf`` / ``icdf`` and would pass an ``isinstance`` check while
  reporting ``+inf`` for every mass.

  Parameters
  ----------
  obj : object
      A margin, or a distribution object from a supported ecosystem.

  Returns
  -------
  MarginLike
      ``obj`` itself, or a thin wrapper around it.

  Raises
  ------
  TypeError
      If nothing recognizes ``obj``. Subclass
      :class:`~pyvinecopulib.core.MarginBase` or call
      :func:`register_margin_adapter`.
  """
  if isinstance(obj, MarginBase):
    return obj
  for predicate, adapter in (*_ADAPTERS, *_BUILTIN_ADAPTERS):
    if predicate(obj):
      return adapter(obj)
  # Known foreign APIs are considered above, before this structural check.
  # Their matching names do not necessarily carry the contract's semantics
  # (notably SciPy's modern discrete `pdf`). A user-defined structural margin,
  # however, is the extension point promised by `MarginLike` and needs no base
  # class or registry entry.
  if isinstance(obj, MarginLike):
    return obj
  raise TypeError(
    f"cannot use {type(obj).__name__!r} as a margin. Subclass "
    "pyvinecopulib.core.MarginBase, or teach as_margin about it with "
    "pyvinecopulib.margins.register_margin_adapter."
  )


#: The ecosystems this package adapts, tried in order after anything
#: registered. Every predicate reads ``type(obj).__mro__``, so none of them
#: imports the ecosystem it recognizes. The two SciPy predicates are mutually
#: exclusive. ``Kde1d`` needs no entry: nothing here matches it, and it reaches
#: the structural check in :func:`as_margin` unchanged.
_BUILTIN_ADAPTERS: tuple[
  tuple[Callable[[Any], bool], Callable[[Any], MarginLike[Any]]], ...
] = (
  (_is_scipy_new, _adapt_scipy_new),
  (_is_scipy_legacy, _adapt_scipy_legacy),
  (_is_torch_distribution, _adapt_torch),
)


#: Bumped when the payload's shape changes incompatibly.

#: Readers registered by a caller. Consulted before ``_BUILTIN_READERS``, so a
#: caller may replace a shipped one.
_READERS: dict[str, Callable[[dict[str, Any]], Any]] = {}


def register_margin_json(
  kind: str, reader: Callable[[dict[str, Any]], Any]
) -> None:
  """Teach :func:`margin_from_json` how to rebuild one margin type.

  Parameters
  ----------
  kind : str
      The value the margin's ``to_json`` writes under ``"kind"``. Conventionally
      the class name.
  reader : callable
      ``reader(payload) -> margin``, receiving the mapping that
      :func:`margin_to_json` produced.

  Raises
  ------
  ValueError
      If ``kind`` is already registered to a different reader.
  """
  existing = _READERS.get(kind)
  if existing is not None and existing is not reader:
    raise ValueError(f"a reader for margin kind {kind!r} is already registered")
  _READERS[kind] = reader


def margin_to_json(margin: object) -> dict[str, Any]:
  """Return one margin's JSON payload.

  Parameters
  ----------
  margin : object
      The margin to serialize. It must provide ``to_json``, which every margin
      this package ships does.

  Returns
  -------
  dict
      A JSON-serializable mapping carrying ``kind`` and ``version``.

  Raises
  ------
  TypeError
      If the margin has no ``to_json``.
  """
  to_json = getattr(margin, "to_json", None)
  if to_json is None:
    raise TypeError(
      f"{type(margin).__name__} cannot be serialized: it has no `to_json`. "
      "Implement `to_json` returning a JSON-serializable mapping, and call "
      "`pyvinecopulib.core.register_margin_json` so it can be read back."
    )
  payload = to_json()
  if isinstance(payload, str):
    # A compiled margin (`Kde1d`) serializes itself to a JSON string.
    payload = {"kind": type(margin).__name__, "json": payload}
  else:
    payload = dict(payload)
    payload.setdefault("kind", type(margin).__name__)
  payload.setdefault("version", MODEL_JSON_VERSION)
  return payload


def margin_from_json(payload: dict[str, Any]) -> MarginLike[Any]:
  """Rebuild one margin from the payload :func:`margin_to_json` produced.

  Parameters
  ----------
  payload : dict
      A mapping carrying ``kind`` and ``version``.

  Returns
  -------
  MarginLike
      The reconstructed margin.

  Raises
  ------
  ValueError
      If ``kind`` is unknown or unregistered, or the version is unrecognized.
  """
  # No `kind=` here: this reader dispatches on `kind` itself and says far
  # more about an unknown one than a bare mismatch could.
  payload = read_payload(payload, "margin")
  kind = payload.get("kind")
  reader = _READERS.get(kind) or _BUILTIN_READERS.get(kind)
  if reader is None:
    known = ", ".join(sorted({*_READERS, *_BUILTIN_READERS})) or "(none)"
    raise ValueError(
      f"no reader registered for margin kind {kind!r}; known kinds: {known}. "
      "A margin from an optional extra registers its reader when its module "
      "is imported, so `import pyvinecopulib.margins` or "
      "`import pyvinecopulib.torch` first if the payload names one of those; "
      "otherwise call `pyvinecopulib.core.register_margin_json`."
    )
  return reader(payload)


def _read_kde1d(payload: dict[str, Any]) -> Kde1d:
  """Rebuild a ``Kde1d``, whose own JSON is a string rather than a mapping."""
  return Kde1d.from_json(payload["json"])


#: ``kind`` -> the reader for it. Only ``Kde1d`` is here: it is a ``core``
#: class, so ``core`` can name it. Every other margin registers its own
#: reader from its own module through :func:`register_margin_json`, which is
#: what keeps ``core`` from naming a class it must not import -- and makes the
#: first-party margins use the same hook a third party does.
_BUILTIN_READERS: dict[str, Callable[[dict[str, Any]], Any]] = {
  "Kde1d": _read_kde1d,
}


def _index_for(
  key: str | int, lookup: dict[str, int], d: int, label: str
) -> int:
  """Resolve one mapping key to a variable position.

  Shared by ``margins=`` and ``margin_controls=``, which accept the same four
  shapes and must therefore refuse the same keys the same way.

  Parameters
  ----------
  key : object
      A variable name or an integer position.
  lookup : dict
      Variable name to position; empty when the data carry no names.
  d : int
      Number of variables.
  label : str
      The argument's name, quoted in the error messages.

  Returns
  -------
  int
      The variable's position.

  Raises
  ------
  ValueError
      If the key names no variable, is neither a name nor an integer, or is
      out of range.
  """
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
    return lookup[key]
  try:
    index = operator.index(key)
  except TypeError as e:
    raise ValueError(
      f"{label} mapping key {key!r} is neither a variable name nor an "
      "integer position"
    ) from e
  if not 0 <= index < d:
    raise ValueError(f"{label} mapping has out-of-range index {index}")
  return index


def _per_variable(
  spec: object,
  d: int,
  *,
  names: Sequence[str] | None,
  default: object,
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
  names : sequence of str, or None, optional
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
      index = _index_for(key, lookup, d, label)
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
  spec: object,
  d: int,
  *,
  names: Sequence[str] | None = None,
) -> list[Any]:
  """Expand ``margin_controls=`` into one controls object per variable.

  Four shapes: ``None``, one controls object broadcast to every variable, a
  length-``d`` sequence, or a mapping keyed by variable name or position with
  the unaddressed variables left unconfigured. That is what lets one call bound
  the family search on the variables that need it and leave the rest alone.

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


def unshare(margins: Sequence[Any]) -> list[Any]:
  """Give every position its own margin, copying only what is aliased.

  Parameters
  ----------
  margins : sequence of MarginLike
      The margins to de-alias, one per variable.

  Returns
  -------
  list
      The same margins, with every repeated reference replaced by a copy.

  Notes
  -----
  Estimating a margin mutates it, so one margin standing at several positions
  would be estimated once per column -- each fit overwriting the last -- and
  every one of those columns would end up on the fit from the final column.
  Each *distinct* margin is left alone, so a caller holding one still sees it
  re-estimated where a lane promises to do that in place.
  """
  seen: set[int] = set()
  out: list[Any] = []
  for margin in margins:
    if id(margin) in seen:
      out.append(copy.deepcopy(margin))
    else:
      seen.add(id(margin))
      out.append(margin)
  return out
