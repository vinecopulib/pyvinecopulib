"""Canonical partial implementation of the univariate-margin contract.

Written against the Array API (``array_api_compat.array_namespace``), so a
subclass works on NumPy or PyTorch without change. Arrays are typed ``Any``
per the ``pyvinecopulib.core`` typing policy (the Array API namespace is
untyped); the public signatures carry ``ArrayT`` so a concrete subclass
inherits precise return types.
"""

from __future__ import annotations

from abc import ABC
from typing import Any, Optional, cast

from array_api_compat import array_namespace

from ._rootfind import solve_increasing
from .protocols import _MARGIN_EXAMPLE, ArrayT, MarginLike

__all__ = ["MarginBase"]

#: Bisection steps behind the inherited :meth:`MarginBase.icdf`. Bisection halves
#: an *absolute* interval, so the accuracy it buys depends on where the answer
#: sits relative to the bracket: the h-inverses' 50 steps resolve the unit
#: interval to 1e-15, but on an unbounded support they leave a quantile of
#: magnitude 1e-16 quantized to multiples of 4.4e-16 -- a 32% error on a
#: `gamma(0.2)` left tail, and unbounded error further down. 110 steps carry the
#: same case to 5e-15 relative. The cost is one extra ``cdf`` call per step, and
#: it does not touch ``solve_increasing``'s own default, so the h-inverse path is
#: unchanged. A quantile whose magnitude is far below the bracket's scale is
#: still resolved only to that absolute floor; override ``icdf`` where the family
#: has a closed form.
_ICDF_ITER: int = 110

#: Variable types a margin may declare. ``"zi"`` (zero-inflated) is a margin-level
#: distinction; a vine sees it as ``"d"``, since what it needs is the left limit.


def support_of(obj: Any) -> tuple[float, float]:
  """Read ``(lo, hi)`` off an object that may spell its support three ways.

  SciPy exposes ``support()`` as a method, PyTorch a ``support`` property holding
  a ``Constraint``, and most objects nothing at all. A parameter-derived endpoint
  arrives as a tensor that may carry a gradient, so it is detached before being
  read as a scalar -- duck-typed, since ``core`` does not import torch.

  Parameters
  ----------
  obj : object
      The distribution to inspect.

  Returns
  -------
  tuple of float
      ``(lo, hi)``, unbounded on either side that cannot be determined.
  """
  unbounded = (float("-inf"), float("inf"))

  def scalar(bound: Any, fallback: float) -> float:
    if bound is None:
      return fallback
    return float(bound.detach() if hasattr(bound, "detach") else bound)

  support = getattr(obj, "support", None)
  if support is None:
    return unbounded
  if callable(support):
    try:
      lo, hi = support()
    except Exception:  # noqa: BLE001 - a foreign object may refuse; fall back
      return unbounded
    return (scalar(lo, unbounded[0]), scalar(hi, unbounded[1]))
  # A torch `Constraint`; only the interval-like ones carry usable bounds.
  return (
    scalar(getattr(support, "lower_bound", None), unbounded[0]),
    scalar(getattr(support, "upper_bound", None), unbounded[1]),
  )


def _reject_covariates(margin: Any, x: Optional[Any]) -> None:
  """Raise if ``x`` was supplied to a margin that fits unconditionally.

  Evaluation ignores covariates a margin does not read, since one distribution
  may mix conditional and unconditional margins and they all see the same ``x``.
  Fitting cannot: silently estimating ``f(y)`` when ``f(y | x)`` was asked for
  returns a different model than the caller believes they have.

  Raises
  ------
  ValueError
      If ``x`` is not ``None``.
  """
  if x is not None:
    raise ValueError(
      f"{type(margin).__name__} does not model covariates, so it cannot be "
      "fitted with x=; write a margin whose fit reads them and that declares "
      "supports_covariates, or drop x."
    )


def derive_cdf_left(
  margin: Any, y: Any, x: Optional[Any], var_type: str
) -> Any:
  """Left limit ``F(y^-)`` derived from a margin's ``cdf`` and ``var_type``.

  The one place the derivation lives, so a margin that declares atoms and
  supplies no ``cdf_left`` of its own gets a left limit that carries
  information rather than a copy of ``cdf``: a copy would make the pair-copula
  difference quotient vanish, and the marginal mass the copula divides out
  would then never cancel against the one the likelihood adds back.

  Parameters
  ----------
  margin : MarginLike
      The margin to read ``cdf`` (and, for ``"zi"``, ``pdf``) from.
  y : array, shape (n,), dtype float
      Observations on the original scale.
  x : array, shape (n, k), or None
      Exogenous covariates, forwarded only to a margin that reads them.
  var_type : {"c", "d", "zi"}
      The margin's variable type.

  Returns
  -------
  array, shape (n,), dtype float
      Left limits in ``[0, 1]``.

  Raises
  ------
  ValueError
      If ``var_type`` is ``"d"`` and ``y`` is not integer-valued, since the
      derivation steps back one lattice point.
  """
  if var_type == "c":
    return _margin_eval(margin, "cdf", y, x)

  ya: Any = y
  xp = array_namespace(ya)
  if var_type == "d":
    if not bool(xp.all(ya == xp.round(ya))):
      raise ValueError(
        f"{type(margin).__name__} declares var_type='d', so y must be "
        "integer-valued; give it a cdf_left for a support on another lattice."
      )
    return _margin_eval(margin, "cdf", ya - 1, x)

  # Zero-inflated: the only atom is at 0, and its mass is `pdf(0)`.
  upper: Any = _margin_eval(margin, "cdf", y, x)
  mass: Any = _margin_eval(margin, "pdf", y, x)
  return xp.where(ya == 0, xp.clip(upper - mass, 0.0, 1.0), upper)


def _margin_eval(margin: Any, name: str, values: Any, x: Optional[Any]) -> Any:
  """Evaluate a margin method, forwarding ``x`` only when it is read.

  Unlike the pair-copula gate, an unconditional margin *ignores* covariates
  rather than refusing them: one vine distribution may hold conditional and
  unconditional margins side by side, and the caller cannot be asked to know
  which is which. A margin opts in by declaring ``supports_covariates``, which
  is also what makes forgetting the declaration a visibly unconditional fit
  rather than a silently ignored keyword.
  """
  method = getattr(margin, name)
  if x is None or not getattr(margin, "supports_covariates", False):
    return method(values)
  return method(values, x=x)


class MarginBase(MarginLike[ArrayT], ABC):
  """Canonical partial implementation of :class:`MarginLike`.

  Subclasses implement ``pdf`` and ``cdf`` and inherit ``icdf`` (numerical
  inversion of ``cdf`` over :attr:`support`), ``logpdf``, ``cdf_left``,
  ``loglik``, ``sample`` and ``__repr__``. Override any of them with a
  closed form where one exists — in particular ``icdf`` and ``cdf_left``, whose
  defaults are correct but not the most accurate route for a family that knows
  its own quantile function or left limit.

  Every member takes the observations positionally and covariates as a
  keyword-only ``x``, so a subclass that models ``f(y)`` writes ``pdf(self, y)``
  and one that models ``f(y | x)`` writes ``pdf(self, y, *, x=None)`` and sets
  :attr:`supports_covariates`. The inherited members forward ``x`` only to a
  margin that declares it, so both shapes are complete.

  A margin is also its own estimator: hyperparameters go in ``__init__``, and
  :meth:`fit` estimates the remaining parameters in place and returns ``self``,
  as :class:`pyvinecopulib.core.Kde1d` and the scikit-learn estimators do. A
  margin whose parameters are fully specified at construction is already fitted
  and need not override :meth:`fit`.

  See Also
  --------
  MarginLike : The contract this fills in.
  pyvinecopulib.core.BicopBase : The pair-copula analog.
  """

  #: Whether :meth:`fit` accepts observation weights. Declared per family so a
  #: caller passing weights to a family that cannot use them gets an error
  #: rather than a silently unweighted fit.
  supports_weights: bool = False

  #: Which variable types this family can serve, as a superset of what any one
  #: instance currently declares through :attr:`var_type`. A selector uses it to
  #: judge admissibility before configuring a candidate, where ``var_type`` can
  #: only report what the candidate already is.
  supported_var_types: tuple[str, ...] = ("c",)

  #: Whether this margin reads the exogenous covariates ``x``, i.e. whether it
  #: models ``f(y | x)`` rather than ``f(y)``. Declared so consumers can omit
  #: ``x`` entirely for an unconditional margin, which is what lets a margin
  #: written against ``pdf(self, y)`` stay valid.
  supports_covariates: bool = False

  def __init_subclass__(cls, **kwargs: Any) -> None:
    """Reject an override of the pre-1.0 hook name.

    Parameters
    ----------
    **kwargs
        Forwarded to ``super().__init_subclass__``.
    """
    super().__init_subclass__(**kwargs)

  # --- optional capabilities, with continuous-correct defaults ------------- #

  @property
  def var_type(self) -> str:
    """Variable type: ``"c"``, ``"d"`` or ``"zi"``.

    Returns
    -------
    str
        ``"c"`` unless a subclass declares atoms.
    """
    return "c"

  @property
  def support(self) -> tuple[float, float]:
    """Closed bounds of the support, as ``(lo, hi)``.

    Infinite bounds are ``-inf`` / ``inf`` rather than ``None`` or ``NaN``, so
    comparisons against them are total. :meth:`icdf` brackets its search here.

    Returns
    -------
    tuple of float
        ``(-inf, inf)`` unless a subclass narrows it.
    """
    return (float("-inf"), float("inf"))

  @property
  def is_fitted(self) -> bool:
    """Whether the margin is ready to evaluate.

    Returns
    -------
    bool
        ``True`` unless a subclass defers parameters to :meth:`fit`.
    """
    return True

  def declare(
    self,
    *,
    var_type: Optional[str] = None,
    support: Optional[tuple[float, float]] = None,
  ) -> "MarginBase[ArrayT]":
    """Accept what the caller knows about the variable, before fitting it.

    A caller often knows the variable type and the declared support when the
    margin cannot infer them: an ordered categorical's levels, a column
    documented as a count, a rate bounded below by zero. Without this, a margin
    handed such a column re-infers both from the sample, which is strictly less
    information -- a count column whose smallest observation is 3 looks
    unbounded from below.

    The base implementation ignores both, so a margin whose type and support are
    fixed by construction needs nothing. An override must treat an explicit
    constructor argument as authoritative and only fill in what was left open,
    since the caller's schema is a default and not an instruction.

    Parameters
    ----------
    var_type : str or None, optional
        ``"c"``, ``"d"`` or ``"zi"``, or ``None`` when the caller does not know.
    support : tuple of float, or None, optional
        Declared bounds as ``(lo, hi)``, or ``None``. Infinite bounds are
        ``-inf`` / ``inf``.

    Returns
    -------
    MarginBase
        ``self``, so the call chains into :meth:`fit`.
    """
    return self

  def fit(
    self,
    y: ArrayT,
    /,
    *,
    x: Optional[ArrayT] = None,
    weights: Optional[ArrayT] = None,
  ) -> "MarginBase[ArrayT]":
    """Estimate the margin's parameters from data, in place.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation. A margin that reads them
        declares :attr:`supports_covariates`, and sets it when it does.
    weights : array, shape (n,), or None, optional
        Observation weights; only accepted when :attr:`supports_weights`.

    Returns
    -------
    MarginBase
        ``self``, so the call chains.

    Raises
    ------
    NotImplementedError
        If the subclass has no estimator; a margin specified entirely at
        construction is already fitted and does not need one.
    """
    raise NotImplementedError(
      f"{type(self).__name__}.fit is not defined; implement it to estimate "
      "this margin from data, or construct it with explicit parameters."
    )

  # --- derived evaluation surface ------------------------------------------ #

  def logpdf(self, y: ArrayT, /, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Log of :meth:`pdf`.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation.

    Returns
    -------
    array, shape (n,), dtype float
        Log-density, ``-inf`` where the density vanishes.
    """
    dens: Any = _margin_eval(self, "pdf", y, x)
    xp = array_namespace(dens)
    positive = dens > 0
    safe = xp.where(positive, dens, xp.ones_like(dens))
    return cast(
      ArrayT,
      xp.where(positive, xp.log(safe), xp.full_like(dens, float("-inf"))),
    )

  def cdf_left(self, y: ArrayT, /, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Left limit ``F(y^-)`` of the distribution function.

    A vine needs this wherever a margin has atoms: ``F(y)`` alone does not
    identify the copula there. The default is derived from :attr:`var_type` —
    it coincides with ``cdf`` for a continuous margin, steps back one lattice
    point for a discrete one, and removes the atom for a zero-inflated one.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation.

    Returns
    -------
    array, shape (n,), dtype float
        Left limits in ``[0, 1]``.

    Raises
    ------
    ValueError
        If :attr:`var_type` is ``"d"`` and ``y`` is not integer-valued, since
        the default steps back by one.
    """
    return cast(ArrayT, derive_cdf_left(self, y, x, self.var_type))

  def icdf(self, p: ArrayT, /, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse distribution function, by numerical inversion of ``cdf``.

    Bisects over :attr:`support`, widening an infinite bound outward first.
    Override with a closed form where the family has one — and necessarily so
    for a conditional margin whose support moves with ``x``, since
    :attr:`support` describes the margin as a whole.

    Parameters
    ----------
    p : array, shape (n,), dtype float
        Probabilities in ``[0, 1]``.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation.

    Returns
    -------
    array, shape (n,), dtype float
        Quantiles on the original scale.
    """
    lo, hi = self.support
    out: Any = solve_increasing(
      lambda v: _margin_eval(self, "cdf", v, x),
      p,
      lo=lo,
      hi=hi,
      n_iter=_ICDF_ITER,
    )
    if self.var_type == "d":
      # Snapped to the lattice the margin declares. Bisection converges to the
      # jump rather than landing on it, and an answer a few ulp below an integer
      # is one its own `cdf_left` would reject as non-integral.
      out = array_namespace(out).round(out)
    return cast(ArrayT, out)

  @property
  def _fitted_loglik(self) -> float:
    """Log-likelihood attained at :meth:`fit`.

    Returns
    -------
    float
        The stored value.

    Raises
    ------
    NotImplementedError
        Unless the subclass records one; a margin specified entirely at
        construction never attained a log-likelihood.
    """
    raise NotImplementedError(
      f"{type(self).__name__} records no fitted log-likelihood; pass data to "
      "loglik() to evaluate one."
    )

  def loglik(
    self,
    y: Optional[ArrayT] = None,
    /,
    *,
    x: Optional[ArrayT] = None,
    weights: Optional[ArrayT] = None,
  ) -> Any:
    """Log-likelihood of the observations under the margin.

    Called without data, returns the value attained when the margin was fitted,
    as ``Bicop.loglik`` and ``Vinecop.loglik`` do.

    Parameters
    ----------
    y : array, shape (n,), or None, optional
        Observations on the original scale; the fitted value when ``None``.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation.
    weights : array, shape (n,), or None, optional
        Observation weights; unweighted when ``None``.

    Returns
    -------
    array, shape (), dtype float, or float
        A 0-d array when evaluated, so it stays differentiable on autograd
        backends; a float when the fitted value is returned.

    Raises
    ------
    ValueError
        If ``weights`` is given without data.
    """
    if y is None:
      if weights is not None:
        raise ValueError("weights are only meaningful with data; pass y too")
      return self._fitted_loglik
    terms: Any = _margin_eval(self, "logpdf", y, x)
    xp = array_namespace(terms)
    if weights is not None:
      terms = terms * weights
    return cast(ArrayT, xp.sum(terms))

  def sample(
    self,
    n: int,
    *,
    x: Optional[ArrayT] = None,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT:
    """Draw ``n`` samples from the margin.

    Parameters
    ----------
    n : int
        Number of samples to draw. With covariates, ``x`` supplies one row per
        draw, so it must have ``n`` of them.
    x : array, shape (n, k), or None, optional
        Exogenous covariates to condition each draw on.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n,), dtype float
        Samples on the original scale.
    """
    base = self._sample_uniform(n, list(seeds) if seeds else [])
    return cast(ArrayT, _margin_eval(self, "icdf", base, x))

  def _sample_uniform(self, n: int, seeds: list[int]) -> ArrayT:
    """Draw ``n`` uniforms on the subclass's array namespace.

    Parameters
    ----------
    n : int
        Number of draws.
    seeds : list of int
        RNG seeds.

    Returns
    -------
    array, shape (n,), dtype float
        Uniforms in ``(0, 1)``.

    Raises
    ------
    NotImplementedError
        Always, unless overridden; the base class has no array namespace to
        draw from. Named after the exposed
        :func:`pyvinecopulib.utils.sample_uniform` free function.
    """
    raise NotImplementedError(
      f"{type(self).__name__}._sample_uniform is not defined; override it "
      "with the array namespace's RNG to enable sample()."
    )

  def __repr__(self) -> str:
    """Return a structural representation of the margin.

    Returns
    -------
    str
        The class name and its variable type.
    """
    return f"{type(self).__name__}(var_type={self.var_type!r})"


MarginBase.__doc__ = (MarginBase.__doc__ or "") + _MARGIN_EXAMPLE
