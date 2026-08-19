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


def _support_of(obj: Any) -> tuple[float, float]:
  """Read a support from an object that may spell it three different ways.

  SciPy exposes ``support()`` as a method, PyTorch a ``support`` property
  holding a ``Constraint``, and most other objects nothing at all.

  Parameters
  ----------
  obj : object
      The foreign distribution.

  Returns
  -------
  tuple of float
      ``(lo, hi)``, unbounded where it cannot be determined.
  """
  unbounded = (float("-inf"), float("inf"))
  support = getattr(obj, "support", None)
  if support is None:
    return unbounded
  if callable(support):
    try:
      lo, hi = support()
    except Exception:  # noqa: BLE001 - a foreign object may refuse; fall back
      return unbounded
    return (float(lo), float(hi))
  # A torch `Constraint`; only the interval-like ones carry usable bounds.
  lo = getattr(support, "lower_bound", None)
  hi = getattr(support, "upper_bound", None)
  return (
    float(lo) if lo is not None else float("-inf"),
    float(hi) if hi is not None else float("inf"),
  )


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
    var_type: str = "c",
    cdf_left: Optional[Callable[[Any], Any]] = None,
    support: Optional[tuple[float, float]] = None,
    family_name: Optional[str] = None,
  ) -> None:
    self.wrapped = obj
    self._pdf = pdf
    self._cdf = cdf
    self._icdf = icdf
    self._var_type = var_type
    self._cdf_left = cdf_left
    self._support = support if support is not None else _support_of(obj)
    self.family_name = family_name or type(obj).__name__

  @property
  def var_type(self) -> str:
    return self._var_type

  @property
  def support(self) -> tuple[float, float]:
    return self._support

  def pdf(self, x: Any) -> Any:
    return self._pdf(x)

  def cdf(self, x: Any) -> Any:
    return self._cdf(x)

  def icdf(self, p: Any) -> Any:
    return self._icdf(p)

  def cdf_left(self, x: Any) -> Any:
    if self._cdf_left is None:
      return super().cdf_left(x)
    return self._cdf_left(x)

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
      obj, pdf=obj.pdf, cdf=obj.cdf, icdf=obj.icdf, var_type="c"
    )
  return _WrappedMargin(
    obj,
    pdf=obj.pmf,
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
      cdf=obj.cdf,
      icdf=obj.ppf,
      var_type="c",
      family_name=name,
    )
  return _WrappedMargin(
    obj,
    pdf=obj.pmf,
    cdf=obj.cdf,
    icdf=obj.ppf,
    var_type="d",
    # Legacy `cdf` is a genuine step function, so `F(x - 1)` is exact on the
    # integer lattice and avoids canceling `F` against the mass in the tail.
    cdf_left=lambda x: obj.cdf(x - 1),
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
  them — so conformance cannot be inferred from member names here. A missing
  ``icdf`` falls back to inverting ``cdf`` numerically.
  """
  lo, hi = _support_of(obj)

  def _icdf(p: Any) -> Any:
    try:
      return obj.icdf(p)
    except NotImplementedError:
      from ..core._rootfind import solve_increasing

      return solve_increasing(obj.cdf, p, lo=lo, hi=hi)

  return _WrappedMargin(
    obj,
    pdf=lambda x: obj.log_prob(x).exp(),
    cdf=obj.cdf,
    icdf=_icdf,
    var_type="c",
    support=(lo, hi),
  )


def as_margin(obj: Any) -> MarginLike[Any]:
  """Present ``obj`` as a :class:`~pyvinecopulib.core.MarginLike`.

  Idempotent: anything this library produced is returned unchanged. Foreign
  objects are wrapped, because satisfying the contract's member *names* is not
  the same as satisfying its semantics — a modern SciPy discrete distribution
  has ``pdf`` / ``cdf`` / ``icdf`` and would pass an ``isinstance`` check while
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
