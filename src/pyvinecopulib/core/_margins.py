"""The two margin registries: coercion, and JSON persistence.

Both answer a question about a margin whose type is not known in advance, and
both answer it the same way -- one table, resolved when it is needed:

* **Coercion.** :func:`as_margin` presents a distribution object from another
  ecosystem as a :class:`~pyvinecopulib.core.MarginLike`. Every margin the
  library receives passes through it, so a discrete SciPy object cannot slip
  past on a bare ``pdf`` (in SciPy's modern API that is ``+inf`` at an atom;
  the mass lives on ``pmf``).
* **Persistence.** :func:`margin_to_json` and :func:`margin_from_json` round a
  margin through a mapping carrying a ``kind`` naming its type and a
  ``version`` a reader checks, so a format change fails loudly rather than
  producing a wrong model.

Each registry has a table of the margins this package ships and a registration
function for everything else (:func:`register_margin_adapter`,
:func:`register_margin_json`). A built-in entry imports lazily, because two of
the classes live behind an extra. Nothing registers itself from its own module,
so whether an object is recognized does not depend on some other module's
import list.
"""

from __future__ import annotations

import json
import math
from typing import TYPE_CHECKING, Any, Callable, Optional, cast

import numpy as np

from .margin_base import MarginBase, support_of
from .protocols import ArrayT, MarginLike

if TYPE_CHECKING:
  # `core` imports without PyTorch; these names are read by the type checker
  # only. torch ships `py.typed`, so they are its own declarations rather than
  # a restatement of them.
  from torch import Tensor
  from torch.distributions import Distribution

  from ..pyvinecopulib_ext import Kde1d

__all__ = [
  "as_margin",
  "margin_from_json",
  "margin_to_json",
  "register_margin_adapter",
  "register_margin_json",
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
    logpdf: Optional[Callable[[ArrayT], ArrayT]] = None,
    var_type: str = "c",
    cdf_left: Optional[Callable[[ArrayT], ArrayT]] = None,
    support: Optional[tuple[float, float]] = None,
    family_name: Optional[str] = None,
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

  def pdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    return self._pdf(y)

  def logpdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    if self._logpdf is None:
      return super().logpdf(y)
    return self._logpdf(y)

  def cdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    return self._cdf(y)

  def icdf(self, p: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    return self._icdf(p)

  def cdf_left(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
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


def _adapt_torch(obj: Distribution) -> MarginLike[Tensor]:
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
  if (
    cdf is None
    or getattr(cdf, "__module__", "") == "torch.distributions.distribution"
    and getattr(cdf, "__qualname__", "") == "Distribution.cdf"
  ):
    raise TypeError(
      f"cannot adapt torch distribution {type(obj).__name__!r}: its cdf is "
      "not implemented. Use a continuous torch distribution with a cdf, "
      "provide a MarginBase implementation, or use a SciPy/OpenTURNS margin."
    )

  lo, hi = support_of(obj)

  def _icdf(p: Tensor) -> Tensor:
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


def _adapt_openturns(obj: Any) -> MarginLike[Any]:  # noqa: ANN401
  """Adapt an ``openturns`` distribution.

  Parameters
  ----------
  obj : openturns.Distribution
      The distribution to wrap.

  Returns
  -------
  MarginLike
      An ``OpenTURNSMargin`` around it, already fitted.
  """
  from ..margins.openturns import OpenTURNSMargin

  return OpenTURNSMargin.from_distribution(obj)


def _is_openturns_distribution(obj: object) -> bool:
  """Whether ``obj`` is an ``openturns`` distribution.

  Parameters
  ----------
  obj : object
      Any object.

  Returns
  -------
  bool
      ``True`` for a concrete OpenTURNS distribution and for the
      ``Distribution`` interface object a factory returns, which share no base
      class beyond ``Object``.
  """
  return any(
    str(getattr(base, "__module__", "")).startswith("openturns")
    and base.__name__ in ("Distribution", "DistributionImplementation")
    for base in type(obj).__mro__
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
  (_is_openturns_distribution, _adapt_openturns),
)


#: Bumped when the payload's shape changes incompatibly.
MARGIN_JSON_VERSION = 1

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
  payload.setdefault("version", MARGIN_JSON_VERSION)
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
  kind = payload.get("kind")
  version = payload.get("version")
  if version != MARGIN_JSON_VERSION:
    raise ValueError(
      f"unsupported margin JSON version {version!r}; this build reads "
      f"version {MARGIN_JSON_VERSION}"
    )
  reader = _READERS.get(kind) or _BUILTIN_READERS.get(kind)
  if reader is None:
    known = ", ".join(sorted({*_READERS, *_BUILTIN_READERS})) or "(none)"
    raise ValueError(
      f"no reader registered for margin kind {kind!r}; known kinds: {known}. "
      "Call `pyvinecopulib.core.register_margin_json` first."
    )
  return reader(payload)


def _encode_nonfinite(value: object) -> object:
  """Replace non-finite floats with strings, recursively.

  ``json.dumps`` writes ``Infinity`` / ``NaN``, which strict JSON has no
  literal for -- and which the reader behind :func:`write_file` rejects.
  They travel as strings and are restored on read, so a ``-inf`` log-likelihood
  in a selection report survives exactly rather than becoming ``null``.

  Parameters
  ----------
  value : object
      Any JSON-serializable structure.

  Returns
  -------
  object
      The same structure with non-finite floats replaced.
  """
  if isinstance(value, float) and not math.isfinite(value):
    return f"__nonfinite__:{value!r}"
  if isinstance(value, dict):
    return {k: _encode_nonfinite(v) for k, v in value.items()}
  if isinstance(value, (list, tuple)):
    return [_encode_nonfinite(v) for v in value]
  return value


def _decode_nonfinite(value: object) -> object:
  """Invert :func:`_encode_nonfinite`.

  Parameters
  ----------
  value : object
      A structure parsed from JSON.

  Returns
  -------
  object
      The same structure with the encoded floats restored.
  """
  if isinstance(value, str) and value.startswith("__nonfinite__:"):
    return float(value.split(":", 1)[1])
  if isinstance(value, dict):
    return {k: _decode_nonfinite(v) for k, v in value.items()}
  if isinstance(value, list):
    return [_decode_nonfinite(v) for v in value]
  return value


def dumps(payload: dict[str, Any]) -> str:
  """Serialize a payload to a JSON string.

  Parameters
  ----------
  payload : dict
      The mapping to serialize.

  Returns
  -------
  str
      Its JSON representation, with non-finite floats encoded as strings so the
      result is strict JSON.
  """
  return json.dumps(_encode_nonfinite(payload), allow_nan=False)


def loads(text: str) -> dict[str, Any]:
  """Parse a JSON string into a payload.

  Parameters
  ----------
  text : str
      A JSON object produced by :func:`dumps`.

  Returns
  -------
  dict
      The parsed mapping.

  Raises
  ------
  ValueError
      If the text is not a JSON object.
  """
  payload = _decode_nonfinite(json.loads(text))
  if not isinstance(payload, dict):
    raise ValueError("expected a JSON object")
  return payload


def _read_kde1d(payload: dict[str, Any]) -> Kde1d:
  """Rebuild a ``Kde1d``, whose own JSON is a string rather than a mapping."""
  from . import Kde1d

  return Kde1d.from_json(payload["json"])


def _read_scipy_margin(payload: dict[str, Any]) -> MarginLike[Any]:
  """Rebuild a ``SciPyMargin``, which lives behind the SciPy extra."""
  from ..margins.scipy import SciPyMargin

  return SciPyMargin.from_json_payload(payload)


#: ``kind`` -> the reader for it, for the margins this package ships. Consulted
#: after :func:`register_margin_json`'s table, so a caller may override one.
_BUILTIN_READERS: dict[str, Callable[[dict[str, Any]], Any]] = {
  "Kde1d": _read_kde1d,
  "SciPyMargin": _read_scipy_margin,
}


def write_file(filename: str, text: str) -> None:
  """Write a JSON payload, as CBOR when the name ends in ``.cbor``.

  The extension rule is the one ``Bicop.to_file`` / ``Vinecop.to_file`` follow
  -- the same helper, so the whole model surface reads and writes the same
  formats.

  Parameters
  ----------
  filename : str
      Path to write.
  text : str
      A JSON string.
  """
  from ..pyvinecopulib_ext import _json_to_file

  _json_to_file(filename, text)


def read_file(filename: str) -> str:
  """Read a JSON payload written by :func:`write_file`.

  Parameters
  ----------
  filename : str
      Path to read.

  Returns
  -------
  str
      The payload as a JSON string.
  """
  from ..pyvinecopulib_ext import _file_to_json

  return cast("str", _file_to_json(filename))
