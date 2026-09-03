"""Coercion of foreign distribution objects into the margin contract.

Every ecosystem spells the univariate surface differently — SciPy's modern
classes say ``icdf`` and ``sample``, its legacy frozen distributions say ``ppf``
and ``rvs``, PyTorch says ``log_prob`` and has no density at all — and no single
set of member names admits more than one of them. Rather than contort the
contract, :func:`as_margin` wraps what it recognizes; the adapters are a few
lines each.

Third-party ecosystems register their own with
:func:`register_margin_adapter`, so nothing outside this module needs to know
about them.
"""

from __future__ import annotations

from typing import Any, Callable, Optional, cast

from ..core import MarginBase, MarginLike
from ..core.margin_base import support_of

__all__ = ["as_margin", "register_margin_adapter"]

#: Registered ``(predicate, adapter)`` pairs, newest first so a later
#: registration can take precedence over an earlier one.
_ADAPTERS: list[tuple[Callable[[Any], bool], Callable[[Any], MarginLike]]] = []


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


class _WrappedMargin(MarginBase[Any]):
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
  logpdf : callable or None, optional
      Native log-density; ``None`` derives it from ``pdf``.
  var_type : str, optional
      Variable type of the wrapped distribution.
  cdf_left : callable or None, optional
      Left-limit cdf; ``None`` derives it from ``var_type``.
  support : tuple of float, or None, optional
      Support bounds; ``None`` reads them off ``obj``.
  family_name : str or None, optional
      Name to report in selection output; ``None`` uses the wrapped type's.
  """

  def __init__(
    self,
    obj: Any,
    *,
    pdf: Callable[[Any], Any],
    cdf: Callable[[Any], Any],
    icdf: Callable[[Any], Any],
    logpdf: Optional[Callable[[Any], Any]] = None,
    var_type: str = "c",
    cdf_left: Optional[Callable[[Any], Any]] = None,
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

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    return self._pdf(y)

  def logpdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    if self._logpdf is None:
      return super().logpdf(y)
    return self._logpdf(y)

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> Any:
    return self._cdf(y)

  def icdf(self, p: Any, *, x: Optional[Any] = None) -> Any:
    return self._icdf(p)

  def cdf_left(self, y: Any, *, x: Optional[Any] = None) -> Any:
    if self._cdf_left is None:
      return super().cdf_left(y)
    return self._cdf_left(y)

  def __repr__(self) -> str:
    return f"as_margin({self.family_name}, var_type={self._var_type!r})"


def _is_scipy_new(obj: Any) -> bool:
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


def _is_scipy_new_discrete(obj: Any) -> bool:
  """Whether ``obj`` is a modern SciPy *discrete* distribution."""
  return any(
    base.__name__ == "DiscreteDistribution" for base in type(obj).__mro__
  )


def _adapt_scipy_new(obj: Any) -> MarginLike[Any]:
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


def _is_scipy_legacy(obj: Any) -> bool:
  """Whether ``obj`` is a frozen legacy SciPy distribution."""
  dist = getattr(obj, "dist", None)
  return dist is not None and hasattr(obj, "ppf") and hasattr(dist, "name")


def _adapt_scipy_legacy(obj: Any) -> MarginLike[Any]:
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


def _is_torch_distribution(obj: Any) -> bool:
  """Whether ``obj`` is a ``torch.distributions.Distribution``."""
  return any(
    base.__module__.startswith("torch.distributions")
    and base.__name__ == "Distribution"
    for base in type(obj).__mro__
  )


def _adapt_torch(obj: Any) -> MarginLike[Any]:
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

  def _icdf(p: Any) -> Any:
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


def as_margin(obj: Any) -> MarginLike[Any]:
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
  for predicate, adapter in _ADAPTERS:
    if predicate(obj):
      return adapter(obj)
  # Known foreign APIs are considered above, before this structural check.
  # Their matching names do not necessarily carry the contract's semantics
  # (notably SciPy's modern discrete `pdf`). A user-defined structural margin,
  # however, is the extension seam promised by `MarginLike` and needs no base
  # class or registry entry.
  if isinstance(obj, MarginLike):
    return obj
  raise TypeError(
    f"cannot use {type(obj).__name__!r} as a margin. Subclass "
    "pyvinecopulib.core.MarginBase, or teach as_margin about it with "
    "pyvinecopulib.margins.register_margin_adapter."
  )


def _is_kde1d(obj: Any) -> bool:
  """Whether ``obj`` is a :class:`pyvinecopulib.core.Kde1d`."""
  from ..core import Kde1d

  return isinstance(obj, Kde1d)


def _adapt_kde1d(obj: Any) -> MarginLike[Any]:
  """Pass a ``Kde1d`` through: it satisfies the contract already.

  The pass-through is explicit rather than a structural check on
  ``MarginLike``, because a name-only test would also admit objects whose
  ``pdf`` means something else -- a modern SciPy discrete distribution
  reports ``+inf`` for every mass and would sail through.
  """
  return cast(MarginLike[Any], obj)


# Registered newest-first, so the order here is reversed on lookup: the two
# SciPy predicates are mutually exclusive, and the torch one is cheap.
register_margin_adapter(_is_kde1d, _adapt_kde1d)
register_margin_adapter(_is_torch_distribution, _adapt_torch)
register_margin_adapter(_is_scipy_legacy, _adapt_scipy_legacy)
register_margin_adapter(_is_scipy_new, _adapt_scipy_new)
