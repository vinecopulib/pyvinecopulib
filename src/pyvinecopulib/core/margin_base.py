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

from .._deprecations import _reject_renamed_hook
from ._rootfind import solve_increasing
from .protocols import _MARGIN_EXAMPLE, ArrayT, MarginLike

__all__ = ["MarginBase"]

#: Variable types a margin may declare. ``"zi"`` (zero-inflated) is a margin-level
#: distinction; a vine sees it as ``"d"``, since what it needs is the left limit.
VAR_TYPES: tuple[str, ...] = ("c", "d", "zi")


class MarginBase(MarginLike[ArrayT], ABC):
  """Canonical partial implementation of :class:`MarginLike`.

  Subclasses implement ``pdf`` and ``cdf`` and inherit ``icdf`` (numerical
  inversion of ``cdf`` over :attr:`support`), ``logpdf``, ``cdf_left``,
  ``loglik``, ``sample`` and ``__repr__``. Override any of them with a
  closed form where one exists — in particular ``icdf`` and ``cdf_left``, whose
  defaults are correct but not the most accurate route for a family that knows
  its own quantile function or left limit.

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

  def __init_subclass__(cls, **kwargs: Any) -> None:
    """Reject an override of the pre-1.0 hook name.

    Parameters
    ----------
    **kwargs
        Forwarded to ``super().__init_subclass__``.
    """
    super().__init_subclass__(**kwargs)
    _reject_renamed_hook(cls, "_simulate_uniform", "_sample_uniform")

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

  def fit(
    self, x: ArrayT, *, weights: Optional[ArrayT] = None
  ) -> "MarginBase[ArrayT]":
    """Estimate the margin's parameters from data, in place.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations on the original scale.
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

  def logpdf(self, x: ArrayT) -> ArrayT:
    """Log of :meth:`pdf`.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Log-density, ``-inf`` where the density vanishes.
    """
    dens: Any = self.pdf(x)
    xp = array_namespace(dens)
    positive = dens > 0
    safe = xp.where(positive, dens, xp.ones_like(dens))
    return cast(
      ArrayT,
      xp.where(positive, xp.log(safe), xp.full_like(dens, float("-inf"))),
    )

  def cdf_left(self, x: ArrayT) -> ArrayT:
    """Left limit ``F(x^-)`` of the distribution function.

    A vine needs this wherever a margin has atoms: ``F(x)`` alone does not
    identify the copula there. The default is derived from :attr:`var_type` —
    it coincides with ``cdf`` for a continuous margin, steps back one lattice
    point for a discrete one, and removes the atom for a zero-inflated one.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Left limits in ``[0, 1]``.

    Raises
    ------
    ValueError
        If :attr:`var_type` is ``"d"`` and ``x`` is not integer-valued, since
        the default steps back by one.
    """
    var_type = self.var_type
    if var_type == "c":
      return self.cdf(x)

    xa: Any = x
    xp = array_namespace(xa)
    if var_type == "d":
      if not bool(xp.all(xa == xp.round(xa))):
        raise ValueError(
          f"{type(self).__name__} declares var_type='d', so x must be "
          "integer-valued; override cdf_left for a support on another lattice."
        )
      return self.cdf(cast(ArrayT, xa - 1))

    # Zero-inflated: the only atom is at 0, and its mass is `pdf(0)`.
    upper: Any = self.cdf(x)
    at_atom = xa == 0
    return cast(
      ArrayT,
      xp.where(at_atom, xp.clip(upper - self.pdf(x), 0.0, 1.0), upper),
    )

  def icdf(self, p: ArrayT) -> ArrayT:
    """Inverse distribution function, by numerical inversion of ``cdf``.

    Bisects over :attr:`support`, widening an infinite bound outward first.
    Override with a closed form where the family has one.

    Parameters
    ----------
    p : array, shape (n,), dtype float
        Probabilities in ``[0, 1]``.

    Returns
    -------
    array, shape (n,), dtype float
        Quantiles on the original scale.
    """
    lo, hi = self.support
    return cast(
      ArrayT,
      solve_increasing(lambda v: self.cdf(cast(ArrayT, v)), p, lo=lo, hi=hi),
    )

  def loglik(self, x: ArrayT, *, weights: Optional[ArrayT] = None) -> ArrayT:
    """Log-likelihood of the observations under the margin.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations on the original scale.
    weights : array, shape (n,), or None, optional
        Observation weights; unweighted when ``None``.

    Returns
    -------
    array, shape (), dtype float
        A 0-d array, so it stays differentiable on autograd backends.
    """
    terms: Any = self.logpdf(x)
    xp = array_namespace(terms)
    if weights is not None:
      terms = terms * weights
    return cast(ArrayT, xp.sum(terms))

  def sample(self, n: int, *, seeds: Optional[list[int]] = None) -> ArrayT:
    """Draw ``n`` samples from the margin.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n,), dtype float
        Samples on the original scale.
    """
    base = self._sample_uniform(n, list(seeds) if seeds else [])
    return self.icdf(base)

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
