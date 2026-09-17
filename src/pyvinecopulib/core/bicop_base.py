"""Canonical partial implementation of the pair-copula contract.

``BicopBase`` is the array-agnostic (NumPy or PyTorch) base class for
``BicopLike``, and the short path to hosting a custom pair copula in a vine: a
subclass writes three ``_raw`` leaves -- ``_pdf_raw``, the density, and
``_hfunc1_raw`` / ``_hfunc2_raw``, the two h-functions -- and inherits the rest
of the evaluation surface together with the estimator surface ``fit`` /
``select`` / ``from_data``. The public ``pdf`` / ``hfunc1`` / ``hfunc2`` are
concrete *dispatchers* that prepare the argument and apply the pair's
``var_types`` before calling a leaf; overriding one orphans the leaf it hides.
The member inventory, and which members ship a raising default, is on
``BicopBase`` itself.

``x`` is keyword-only on every method, as it is on the contract: it carries the
covariates a conditional pair copula reads, row-aligned with ``u``. ``Bicop``
takes per-row ``parameters`` in that positional slot, so the two surfaces are
near-twins rather than one signature.

The rule for handing a pair copula its covariates is ``_covariates.pair_eval``
-- forward whenever there is one, so a pair accepting none, ``Bicop`` above
all, stays a valid host for a simplified vine -- kept beside the library's
other forwarding rule rather than here.

Array values are handled as ``Any`` inside the numeric bodies per the
``pyvinecopulib.core`` typing policy (``array_api_compat``, which resolves the
array namespace, is itself untyped); the generic ``ArrayT`` lives on the public
signatures.
"""

from __future__ import annotations

import copy
import math
from abc import ABC, abstractmethod
from collections.abc import Callable, Sequence
from typing import Any, ClassVar, Self, cast

from ._bicop_plot import (
  BICOP_PLOT_PARAMS,
  BICOP_PLOT_SUMMARY,
  bicop_plot,
)
from ._covariates import pair_eval, prepare_covariates
from ._loglik import safe_log, sum_loglik
from ._placement import PlacementMixin, QrngUniformMixin, to_numpy
from ._rootfind import solve_increasing
from ._trim import trim
from ._validation import check_var_types
from .margin_base import criteria
from .protocols import (
  _BICOP_EXAMPLE,
  ArrayT,
  BicopLike,
  BoolArray,
  ControlsLike,
  Namespace,
  array_namespace,
)

__all__ = ["BicopBase"]


#: Atom width below which a difference quotient is replaced by the derivative
#: at the atom's midpoint. ``AbstractBicop``'s own threshold, to the digit:
#: below it the numerator and denominator both vanish and the ratio is noise.
DELTA_MIN: float = 5e-5


def rect_prob_from_cdf(
  cdf: Callable[..., ArrayT],
  a1: ArrayT,
  b1: ArrayT,
  a2: ArrayT,
  b2: ArrayT,
  *,
  x: ArrayT | None = None,
) -> ArrayT:
  """``P((a1, b1] x (a2, b2])`` as the four-corner difference of ``_cdf_raw``.

  The generic route to a rectangle's probability, shared by
  :meth:`BicopBase.rect_prob` and by the fallback
  the mixed-discrete quotients take for a pair that overrides no
  ``rect_prob`` of its own, so the two cannot drift.

  Parameters
  ----------
  cdf : callable
      The pair's distribution function, taking an ``(n, 2)`` array.
  a1, b1 : array, shape (n,), dtype float
      Bounds in the first argument, in either order.
  a2, b2 : array, shape (n,), dtype float
      Bounds in the second argument, in either order.
  x : array, shape (n, p), or None, optional
      Exogenous covariates, forwarded to the pair.

  Returns
  -------
  array, shape (n,), dtype float
      Rectangle probabilities.
  """
  any_a1: Any = a1
  xp = array_namespace(any_a1)
  x0, x1 = xp.minimum(a1, b1), xp.maximum(a1, b1)
  y0, y1 = xp.minimum(a2, b2), xp.maximum(a2, b2)

  def at(p: ArrayT, q: ArrayT) -> ArrayT:
    val = pair_eval(cdf, xp.stack([p, q], axis=-1), x=x)
    # A bound of 0 is the distribution's own lower limit, so a corner on it
    # contributes nothing -- and `cdf` would have read a trimmed 1e-10 there.
    return xp.where((p <= 0.0) | (q <= 0.0), xp.zeros_like(val), val)

  # Summed in two pairs, the grouping the reference pair copula sums them in:
  # it is what makes the two agree to the last bit rather than to rounding.
  return (at(x1, y1) + at(x0, y0)) - (at(x0, y1) + at(x1, y0))


def cond_interval_prob_from_hfunc(
  hfunc: Callable[..., ArrayT],
  u_cond: ArrayT,
  lo: ArrayT,
  hi: ArrayT,
  cond_var: int,
  *,
  x: ArrayT | None = None,
) -> ArrayT:
  """``P(lo < U_free <= hi | U_cond = u_cond)`` as a difference of h-functions.

  The conditional counterpart of :func:`rect_prob_from_cdf`, and shared for the
  same reason. Each h-function value is clamped into the open unit interval, so
  on a narrow interval this difference is the less accurate of the two routes.

  Parameters
  ----------
  hfunc : callable
      ``hfunc1`` for ``cond_var=1`` and ``hfunc2`` for ``cond_var=2``.
  u_cond : array, shape (n,), dtype float
      The argument held fixed.
  lo, hi : array, shape (n,), dtype float
      Bounds in the free argument, in either order.
  cond_var : int
      ``1`` or ``2``, the argument held fixed.
  x : array, shape (n, p), or None, optional
      Exogenous covariates, forwarded to the pair.

  Returns
  -------
  array, shape (n,), dtype float
      Conditional probabilities.
  """
  any_lo: Any = lo
  xp = array_namespace(any_lo)
  a, b = xp.minimum(lo, hi), xp.maximum(lo, hi)

  def at(free: ArrayT) -> ArrayT:
    cols = [u_cond, free] if cond_var == 1 else [free, u_cond]
    return pair_eval(hfunc, xp.stack(cols, axis=-1), x=x)

  return at(b) - at(a)


class BicopBase(
  BicopLike[ArrayT], QrngUniformMixin[ArrayT], PlacementMixin, ABC
):
  """Canonical partial implementation of ``BicopLike``.

  A subclass writes three ``_raw`` leaves, and they are the only abstract
  members of the class:

  - ``_pdf_raw(u)``, the pair density.
  - ``_hfunc1_raw(u)`` / ``_hfunc2_raw(u)``, its two conditional distribution
    functions.

  ``_raw`` is ``AbstractBicop``'s own name for the primitive that ignores
  ``var_types``. The public ``pdf`` / ``cdf`` / ``hfunc1`` / ``hfunc2`` /
  ``hinv1`` / ``hinv2`` are concrete **dispatchers**: each places the argument
  on the pair's array namespace, checks its layout, clamps it into the open
  unit square, and applies the pair's ``var_types`` -- so a leaf always sees
  two already-prepared continuous columns and never calls ``_prep_args``
  itself.

  **Do not override a dispatcher.** Every internal caller goes through the
  public member -- :meth:`loglik` calls ``self.pdf``, the plots and the vine
  cascade go through ``pair_eval`` -- so an override leaves the mandatory leaf
  defined and never run. A subclass with per-call knobs puts them on its
  ``controls_class`` instead. Any pair that may sit in a *conditional* vine
  declares ``x`` on every leaf whether or not it reads one, because
  ``pair_eval`` forwards unconditionally.

  Four further leaves are optional, each backing a dispatcher that otherwise
  falls back: ``_cdf_raw`` (needed to declare the pair discrete),
  ``_hinv1_raw`` / ``_hinv2_raw`` (the inherited bisection costs
  ``O(iterations)`` in *model* calls, so a pair backed by a learned model or
  one that inverts in closed form should write them), and ``_flip_raw``.

  Everything else is inherited:

  - :meth:`hinv1` / :meth:`hinv2`, the h-function inverses. An h-function
    increases in the argument being inverted, so these need nothing beyond
    ``hfunc1`` / ``hfunc2``; override either where the family has a closed
    form.
  - :meth:`sample`, the pair's inverse Rosenblatt draw, which rests on
    :meth:`hinv1` and on ``_sample_uniform`` -- the array library's RNG, the
    one hook with no array-agnostic default.
  - :meth:`loglik`, :meth:`plot` and ``__repr__``.

  Three members make the class an estimator, so a vine can name it as its
  ``bicop_class`` and fit one pair per edge: :meth:`fit` estimates the current
  family's parameters in place, :meth:`select` chooses a family as well --
  defaulting to :meth:`fit` wherever there is nothing to choose -- and
  :meth:`from_data` constructs and selects in one call. :meth:`fit` raises
  here, since a pair copula handed its parameters at construction is already
  fitted.

  Two members raise by default, each needed to host the pair somewhere
  particular rather than to evaluate it:

  - :meth:`flip`, the pair with its arguments swapped, which structure
    selection needs to reorient a fitted pair onto its finalized slot;
    evaluation along a fixed structure never asks for it. Write
    ``_flip_raw``.
  - :meth:`cdf`, needed on a **discrete** edge, whose h-functions are
    difference quotients of the distribution function. Write ``_cdf_raw`` and
    declare the pair discrete with ``with_var_types``.

  Two more are supplied rather than raising, and exist to be overridden:
  :meth:`rect_prob` and :meth:`cond_interval_prob`, the two probabilities a
  discrete edge's difference quotients are built from. Both default to the
  difference of :meth:`cdf` or h-function values they replace, so a pair with
  only a ``cdf`` needs nothing; a pair that can measure an atom without that
  cancellation overrides them and is more accurate at a narrow atom.

  ``TorchTllBicop`` is the reference subclass: it supplies ``cdf``, both inverses,
  ``sample`` and ``flip`` natively, and fits its density grid in :meth:`fit`.

  See Also
  --------
  pyvinecopulib.core.BicopLike : The contract this implements.
  pyvinecopulib.torch.TorchTllBicop : A concrete (grid / TLL) subclass.
  """

  # The controls class this pair's fitter reads, or `None` where it reads none.
  # The one place that answers what `controls=None` means here, so a vine and
  # the lanes above read it instead of each naming a class of their own. A
  # plain comment, not a `#:` one: autosummary cannot page an attribute whose
  # value is a class.
  controls_class: ClassVar[type[ControlsLike] | None] = None

  # --- variable types --------------------------------------------------- #
  #: Mirrors ``AbstractBicop``'s in-class ``var_types_{"c", "c"}``: a class
  #: attribute rather than an ``__init__`` assignment, so a subclass that
  #: writes its own constructor without chaining is still continuous, and
  #: ``from_data``'s bare ``cls()`` keeps working.
  _var_types: tuple[str, ...] = ("c", "c")
  _d1: bool = False
  _d2: bool = False

  @property
  def var_types(self) -> list[str]:
    """The two variable types, ``"c"`` (continuous) or ``"d"`` (discrete).

    A pair copula declared discrete in either argument reads the four-column
    layout ``[u1, u2, u1^-, u2^-]`` and returns the mixed-discrete density and
    h-functions; a continuous one reads two columns. The types are state
    because they enter the *fit*: a grid fitted on a discrete edge is a
    different estimate, not the same copula viewed differently.

    Returns
    -------
    list of str
        The pair's ``[type1, type2]``.
    """
    return list(self._var_types)

  @var_types.setter
  def var_types(self, value: Sequence[str]) -> None:
    types = check_var_types(list(value), 2)
    self._var_types = types
    self._d1 = types[0] == "d"
    self._d2 = types[1] == "d"

  def with_var_types(self, var_types: Sequence[str] = ("c", "c")) -> Self:
    """The same pair copula under different variable types.

    Returns ``self`` when the types already match, so declaring a continuous
    pair continuous costs nothing; otherwise a shallow copy carrying the new
    types, which leaves the fitted parameters shared rather than duplicated.

    Parameters
    ----------
    var_types : sequence of str, default=("c", "c")
        The two types, ``"c"`` or ``"d"``.

    Returns
    -------
    BicopBase
        This pair copula, or a copy of it reading ``var_types``.
    """
    types = check_var_types(list(var_types), 2)
    if types == self._var_types:
      return self
    other = copy.copy(self)
    other.var_types = list(types)
    return other

  # --- the evaluation surface ------------------------------------------- #
  # Each public member dispatches on `var_types` and delegates to the `_*_raw`
  # leaf a subclass writes, which always sees two continuous columns --
  # `AbstractBicop`'s shape.
  def pdf(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Density with respect to each argument's own reference measure.

    A continuous argument contributes a derivative and a discrete one the
    probability of its atom, so a mixed pair gives a difference quotient and a
    pair with two discrete arguments the rectangle probability, each divided by
    the atom widths.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        ``[u1, u2]`` when both arguments are continuous, else
        ``[u1, u2, u1^-, u2^-]``.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        Density values.
    """
    xp, u1, u2, u1m, u2m, x = self._atoms(u, x)
    if self._d1 and self._d2:
      return self._pdf_d_d(xp, u1, u2, u1m, u2m, x)
    if self._d1 or self._d2:
      return self._pdf_mixed(
        xp, u1, u2, u1m, u2m, x, discrete=1 if self._d1 else 2
      )
    return self._pdf_c(xp, u1, u2, x)

  def hfunc1(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """``P(U2 <= u2 | U1)``, conditioning on the atom when ``U1`` is discrete.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values.
    """
    xp, u1, u2, u1m, _, x = self._atoms(u, x)
    if not self._d1:
      return self._h1_c(xp, u1, u2, x)
    # Conditioning on `u1^- < U1 <= u1` divides the rectangle probability by
    # the atom's width; the second argument enters at its value either way.
    return self._quotient(
      xp,
      self._strip(xp, u1m, u1, u2, x, axis=1),
      xp.abs(u1 - u1m),
      self._h1_c(xp, 0.5 * (u1 + u1m), u2, x),
    )

  def hfunc2(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """``P(U1 <= u1 | U2)``, conditioning on the atom when ``U2`` is discrete.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values.
    """
    xp, u1, u2, _, u2m, x = self._atoms(u, x)
    if not self._d2:
      return self._h2_c(xp, u1, u2, x)
    return self._quotient(
      xp,
      self._strip(xp, u2m, u2, u1, x, axis=2),
      xp.abs(u2 - u2m),
      self._h2_c(xp, u1, 0.5 * (u2 + u2m), x),
    )

  def logpdf(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Log-density of the pair copula at each observation.

    The other three bases answer this and `BicopBase` did not, so a custom
    pair was the one part a caller had to log themselves. Taken through
    `safe_log` rather than a bare logarithm: a copula density is legitimately
    zero off its support, where `log` warns or answers `nan`.

    Not on `BicopLike`: the compiled `Bicop` has no `logpdf`, and a contract
    member it lacks would put it outside its own contract -- `isinstance`
    compares member names.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation.

    Returns
    -------
    array, shape (n,), dtype float
        Log-density values, ``-inf`` where the density is zero.
    """
    return safe_log(self.pdf(u, x=x))

  def loglik(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Total log-likelihood ``sum(log c(u))`` of the pair at ``u``.

    An observation carrying a ``nan`` has no log-density and is left out of the
    total, as ``Bicop.loglik()`` leaves it out; a density of exactly zero is
    the model ruling an observation out and contributes ``-inf``.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (), dtype float
        The summed log-density, carrying gradients wherever the array library
        tracks them.
    """
    x = prepare_covariates(self, x, int(cast("Any", u).shape[0]))
    # `u` is left to the subclass's own `pdf`, which is where the two-column
    # layout is checked -- the base cannot know whether a pair is on a
    # discrete edge, whose argument is four columns wide.
    return sum_loglik(safe_log(self.pdf(u, x=x)))

  def aic(self, u: ArrayT, /) -> float:
    """Akaike information criterion at these observations, minimized.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Observations to evaluate the log-likelihood on. Required, unlike on
        ``Bicop`` / ``Vinecop``, which carry the value their own fit attained;
        a pair copula may host parts it did not fit, so there is no such value here.

    Returns
    -------
    float
        ``-2 loglik + 2 npars``.

    Raises
    ------
    NotImplementedError
        If ``npars`` reports none.

    See Also
    --------
    bic : The same, penalizing by ``log n`` per parameter.
    """
    npars = self.npars
    if not math.isfinite(npars):
      raise NotImplementedError(
        f"{type(self).__name__} reports no `npars`, so no criterion can "
        "penalize it. Implement `npars` to enable `aic` / `bic`."
      )
    return criteria(float(to_numpy(self.loglik(u))), npars, None)["aic"]

  def bic(self, u: ArrayT, /) -> float:
    """Bayesian information criterion at these observations, minimized.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Observations to evaluate the log-likelihood on; their row count is the
        ``n`` the penalty uses.

    Returns
    -------
    float
        ``-2 loglik + npars log n``.

    Raises
    ------
    NotImplementedError
        If ``npars`` reports none.
    """
    npars = self.npars
    if not math.isfinite(npars):
      raise NotImplementedError(
        f"{type(self).__name__} reports no `npars`, so no criterion can "
        "penalize it. Implement `npars` to enable `aic` / `bic`."
      )
    rows = float(self._prep(u).shape[0])
    return criteria(float(to_numpy(self.loglik(u))), npars, rows)["bic"]

  def hinv1(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Inverse of :meth:`hfunc1` in its second argument.

    Delegates to the continuous leaf when the conditioning argument is
    continuous, and bisects the mixed-discrete :meth:`hfunc1` otherwise.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        Column 0 is the conditioning value ``u1``; column 1 is the level to
        invert. See :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    xp, u1, p, u1m, _, x = self._atoms(u, x)
    if not self._d1:
      return pair_eval(self._hinv1_raw, xp.stack([u1, p], axis=-1), x=x)
    return solve_increasing(
      lambda v: self.hfunc1(xp.stack([u1, v, u1m, v], axis=-1), x=x),
      p,
    )

  def hinv2(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Inverse of :meth:`hfunc2` in its first argument.

    The counterpart of :meth:`hinv1`, on the other h-function.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        Column 0 is the level to invert; column 1 is the conditioning value
        ``u2``. See :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    xp, p, u2, _, u2m, x = self._atoms(u, x)
    if not self._d2:
      return pair_eval(self._hinv2_raw, xp.stack([p, u2], axis=-1), x=x)
    return solve_increasing(
      lambda v: self.hfunc2(xp.stack([v, u2, v, u2m], axis=-1), x=x),
      p,
    )

  def cdf(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Distribution function ``C(u)``, which the left limits do not enter.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values.

    Raises
    ------
    NotImplementedError
        If the subclass supplies no ``_cdf_raw``.
    """
    xp, u1, u2, _, _, x = self._atoms(u, x)
    return self._cdf_c(xp, u1, u2, x)

  def rect_prob(
    self,
    a1: ArrayT,
    b1: ArrayT,
    a2: ArrayT,
    b2: ArrayT,
    *,
    x: ArrayT | None = None,
  ) -> ArrayT:
    """Probability of the rectangle ``(a1, b1] x (a2, b2]``.

    The quantity a discrete argument's difference quotient is built from. This
    reads it as the four-corner difference of ``_cdf_raw``, the continuous
    leaf, which any pair copula with a distribution function can serve -- the
    leaf rather than :meth:`cdf` so the dispatcher cannot recur through it.
    The bounds are placed and not clamped: ``0`` and ``1`` are the
    distribution's own limits here, where a copula argument's ``1e-10`` would
    be. A pair that can evaluate the
    rectangle without that cancellation overrides this --
    :class:`~pyvinecopulib.torch.TorchTllBicop` does, reading the mass off its
    grid -- and gains accuracy at a narrow atom, where differencing amplifies
    an absolute error by ``4 / (w1 w2)`` in the atom widths.

    Parameters
    ----------
    a1, b1 : array, shape (n,), dtype float
        Bounds in the first argument, in either order.
    a2, b2 : array, shape (n,), dtype float
        Bounds in the second argument, in either order.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation.

    Returns
    -------
    array, shape (n,), dtype float
        Rectangle probabilities. A bound of ``0`` is the distribution's own
        lower limit, so a corner on it contributes nothing.

    See Also
    --------
    cond_interval_prob : The conditional counterpart, for a mixed edge.
    """
    return rect_prob_from_cdf(
      self._cdf_raw,
      self._prep(a1),
      self._prep(b1),
      self._prep(a2),
      self._prep(b2),
      x=x,
    )

  def cond_interval_prob(
    self,
    u_cond: ArrayT,
    lo: ArrayT,
    hi: ArrayT,
    cond_var: int,
    *,
    x: ArrayT | None = None,
  ) -> ArrayT:
    """Probability that the free argument falls in ``(lo, hi]``, given the other.

    What a mixed edge's density is built from, as :meth:`rect_prob` is what a
    doubly discrete one is built from. This reads it as the difference of two
    h-function values, each clamped into the open unit interval; a pair that
    can evaluate the mass itself overrides this, and is then not clamped. The
    arguments are placed, as :meth:`rect_prob`'s bounds are.

    Parameters
    ----------
    u_cond : array, shape (n,), dtype float
        The argument held fixed.
    lo, hi : array, shape (n,), dtype float
        Bounds in the free argument, in either order.
    cond_var : int
        ``1`` or ``2``, the argument held fixed.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional probabilities.
    """
    h = self._hfunc1_raw if cond_var == 1 else self._hfunc2_raw
    return cond_interval_prob_from_hfunc(
      h, self._prep(u_cond), self._prep(lo), self._prep(hi), cond_var, x=x
    )

  @classmethod
  def from_data(
    cls,
    u: ArrayT,
    /,
    controls: ControlsLike | None = None,
    *,
    var_types: list[str] | None = None,
    x: ArrayT | None = None,
  ) -> Self:
    """Construct a pair copula and select it from data.

    ``cls().select(u, ...)``, so a subclass that implements :meth:`fit` gets
    this for free -- which is what lets a vine name this class as its
    ``bicop_class`` and fit one pair per edge.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pseudo-observations in ``[0, 1]^2``.
    controls : ControlsLike, or None, optional
        Fit configuration, in whatever form the subclass accepts.
    var_types : list of str, or None, optional
        The two variable types of the edge this pair sits on, ``"c"``
        (continuous) or ``"d"`` (discrete) each. ``None`` means both are
        continuous.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, row-aligned with ``u``, for a conditional pair
        copula. Optional the way it is on the evaluation surface: a pair that
        models covariates reads them, and one that does not is never handed
        them, since a fit that quietly ignored them would return a different
        model than was asked for.

    Returns
    -------
    BicopBase
        The fitted pair copula.

    See Also
    --------
    select : Choose a family for an already-constructed pair copula, in place.
    fit : Estimate the current family's parameters, leaving the family alone.
    """
    return cls().select(u, controls, var_types=var_types, x=x)

  def fit(
    self,
    u: ArrayT,
    /,
    controls: ControlsLike | None = None,
    *,
    var_types: list[str] | None = None,
    x: ArrayT | None = None,
  ) -> Self:
    """Raise; override to estimate this pair copula from data, in place.

    The pair-copula analog of ``MarginBase.fit()``, and what makes a pair
    copula class usable as a fitter: a vine that names it as its
    ``bicop_class`` fits every edge through :meth:`from_data`. It stays
    optional, because a pair copula specified entirely at construction is
    already fitted, and a vine can always be given an explicit ``fit_edge``
    callback instead.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pseudo-observations in ``[0, 1]^2``.
    controls : ControlsLike, or None, optional
        Fit configuration, in whatever form the subclass accepts.
    var_types : list of str, or None, optional
        The two variable types of the edge this pair sits on. ``None`` means
        both are continuous.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, row-aligned with ``u``, for a conditional pair
        copula. Optional the way it is on the evaluation surface: a pair that
        models covariates reads them, and one that does not is never handed
        them, since a fit that quietly ignored them would return a different
        model than was asked for.

    Returns
    -------
    BicopBase
        ``self``, so the call chains -- only when a subclass overrides this
        method.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides an estimator.

    See Also
    --------
    select : Choose a family as well, where the pair copula has one to choose.
    from_data : Construct and fit in one call.
    """
    raise NotImplementedError(
      f"{type(self).__name__}.fit is not defined; implement it to estimate "
      "this pair copula from data, or construct it with explicit parameters."
    )

  def select(
    self,
    u: ArrayT,
    /,
    controls: ControlsLike | None = None,
    *,
    var_types: list[str] | None = None,
    x: ArrayT | None = None,
  ) -> Self:
    """Choose a family for this pair copula and estimate it, in place.

    Defaults to :meth:`fit`, which is the right answer whenever there is
    nothing to choose: a pair copula with a single fixed family, or a
    nonparametric one, is fully determined by its parameters. Override it in a
    subclass that searches a family set, and have it reset any state the search
    is meant to replace -- a rotation, a family tag -- before estimating.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pseudo-observations in ``[0, 1]^2``.
    controls : ControlsLike, or None, optional
        Fit configuration, in whatever form the subclass accepts.
    var_types : list of str, or None, optional
        The two variable types of the edge this pair sits on. ``None`` means
        both are continuous.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, row-aligned with ``u``, for a conditional pair
        copula. Forwarded to :meth:`fit` only when there is one.

    Returns
    -------
    BicopBase
        ``self``, so the call chains.

    See Also
    --------
    fit : Estimate the current family's parameters, leaving the family alone.
    from_data : Construct and select in one call.
    """
    # Forwarded one at a time, so a subclass whose `fit` declares only the
    # arguments it uses still works through `select`. All-or-nothing would
    # break such a subclass on `var_types` alone; this is the idiom
    # `MarginBase.select` uses, and it degrades the same way.
    passed: dict[str, Any] = {}
    if var_types is not None:
      passed["var_types"] = var_types
    if x is not None:
      passed["x"] = x
    if controls is None:
      return self.fit(u, **passed)
    return self.fit(u, controls, **passed)

  def flip(self) -> Self:
    """The pair copula with its two arguments swapped.

    The flipped copula satisfies ``c'(u1, u2) = c(u2, u1)`` with the two
    h-functions (and their inverses) exchanged, and with the variable types
    exchanged too -- a pair that carries them must swap them here, or a
    reoriented discrete slot would evaluate continuously. Required only to
    host the pair in structure *selection*, which reorients each selected pair
    onto its finalized slot (``VinecopBase.select()``); evaluation along a
    fixed structure never asks for it.

    Returns
    -------
    BicopBase
        The argument-swapped pair copula.

    Raises
    ------
    NotImplementedError
        If the subclass supplies no ``_flip_raw``.
    """
    # Through `with_var_types` rather than by assignment: a symmetric pair
    # returns `self` from `_flip_raw`, and writing the swapped types onto it
    # would mutate the copula being flipped.
    return self._flip_raw().with_var_types(
      (self._var_types[1], self._var_types[0])
    )

  def sample(
    self,
    n: int,
    *,
    x: ArrayT | None = None,
    qrng: bool = False,
    seeds: list[int] | None = None,
  ) -> ArrayT:
    """Draw ``n`` samples from the pair copula.

    Available on any subclass that supplies ``_sample_uniform``, the array
    library's RNG.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per sample, for a conditional draw.
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random ones.
    seeds : list of int, or None, optional
        RNG seeds for the draw.

    Returns
    -------
    array, shape (n, 2), dtype float
        Samples in the unit square.

    Raises
    ------
    NotImplementedError
        If the subclass supplies no ``_sample_uniform``.
    """
    x = prepare_covariates(self, x, n)
    base_u: Any = self._sample_uniform(n, qrng, list(seeds) if seeds else [])
    xp = array_namespace(base_u)
    u2: Any = self.hinv1(base_u, x=x)
    return cast("ArrayT", xp.stack([base_u[:, 0], u2], axis=-1))

  #: Whether this pair copula's leaves accept exogenous covariates. Declared
  #: because a bound class reports
  #: `(*args, **kwargs)`, so `pair_eval` cannot read the signature and asks
  #: the declaration before forwarding a matrix. A subclass whose leaves take
  #: `x` -- whether or not they read one -- sets this `True`.
  supports_covariates: bool = False

  #: Whether a fit of this pair copula honors `controls.weights`. ``False``
  #: here because `BicopBase.fit` raises; a subclass that fits declares it.
  supports_weights: bool = False

  # --- the leaves a subclass writes -------------------------------------- #
  # Each takes two continuous columns the dispatcher has already prepared, and
  # knows nothing about atoms. Every leaf declares `x` because `pair_eval`
  # forwards one unconditionally; whether the pair *reads* it is
  # `supports_covariates`.
  @abstractmethod
  def _pdf_raw(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Continuous pair-copula density at each observation.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation. ``pair_eval`` forwards
        one unconditionally whenever the caller supplies it, so every leaf
        declares it even where it reads none -- ``del x`` is what that looks
        like.

    Returns
    -------
    array, shape (n,), dtype float
        Density values.
    """

  @abstractmethod
  def _hfunc1_raw(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Continuous first h-function ``P(U2 <= u2 | U1 = u1)``.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, as on ``_pdf_raw``.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values.
    """

  @abstractmethod
  def _hfunc2_raw(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Continuous second h-function ``P(U1 <= u1 | U2 = u2)``.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, as on ``_pdf_raw``.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values.
    """

  def _cdf_raw(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Raise; override to give the pair copula a distribution ``C(u)``.

    Needed only to declare the pair discrete, whose h-functions are difference
    quotients of the distribution function. Nothing else asks for one -- a
    vine's own ``cdf`` is evaluated by Monte-Carlo simulation, which needs no
    per-pair distribution.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, as on ``_pdf_raw``.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values -- only when a subclass overrides this method.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides one.
    """
    del u, x
    raise NotImplementedError(
      f"{type(self).__name__} has no `cdf`; the vine cdf uses Monte-Carlo "
      "simulation and does not require a per-pair distribution. Implement "
      "`_cdf_raw` to declare this pair copula discrete, whose h-functions are "
      "difference quotients of the distribution function."
    )

  def _hinv1_raw(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Continuous inverse of ``_hfunc1_raw`` in its second argument.

    Solved by monotone bisection, so a subclass needs only ``_hfunc1_raw``;
    override where the family inverts in closed form.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the conditioning value; column 1 is the level to invert.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, carried into the h-function being inverted; a
        conditional pair inverts a different curve at every covariate value.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    ua: Any = u
    xp = array_namespace(ua)
    u1, p = ua[:, 0], ua[:, 1]
    return cast(
      "ArrayT",
      solve_increasing(
        lambda v: pair_eval(self._hfunc1_raw, xp.stack([u1, v], axis=-1), x=x),
        p,
      ),
    )

  def _hinv2_raw(self, u: ArrayT, *, x: ArrayT | None = None) -> ArrayT:
    """Continuous inverse of ``_hfunc2_raw`` in its first argument.

    The counterpart of :meth:`_hinv1_raw`, on the other h-function.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the level to invert; column 1 is the conditioning value.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, carried into the h-function being inverted; a
        conditional pair inverts a different curve at every covariate value.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    ua: Any = u
    xp = array_namespace(ua)
    p, u2 = ua[:, 0], ua[:, 1]
    return cast(
      "ArrayT",
      solve_increasing(
        lambda v: pair_eval(self._hfunc2_raw, xp.stack([v, u2], axis=-1), x=x),
        p,
      ),
    )

  def _flip_raw(self) -> Self:
    """Raise; override to return the pair with its arguments swapped.

    :meth:`flip` swaps the variable types on top of this, so an override need
    only exchange the two arguments of the copula itself.

    Returns
    -------
    BicopBase
        The argument-swapped pair copula -- only when a subclass overrides.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides one.
    """
    raise NotImplementedError(
      f"{type(self).__name__}.flip is not defined; implement `_flip_raw` "
      "(return the argument-swapped copula) to host this pair in structure "
      "selection."
    )

  # --- the mixed-discrete quotients -------------------------------------- #
  # On the pair copula, as `AbstractBicop` keeps them: they are what a discrete
  # declaration means, so no wrapper has to supply them.
  def _atoms(
    self, u: ArrayT, x: ArrayT | None
  ) -> tuple[Namespace[ArrayT], ArrayT, ArrayT, ArrayT, ArrayT, ArrayT | None]:
    """Namespace, the two values, their left limits, and the covariates."""
    ua: Any = self._prep_args(u)
    xp = array_namespace(ua)
    x = prepare_covariates(self, x, int(ua.shape[0]))
    u1, u2 = ua[:, 0], ua[:, 1]
    # A continuous argument's left limit is its own value
    # (``Bicop::format_data``), so the cascade's column for it is never read.
    return (
      xp,
      u1,
      u2,
      ua[:, 2] if self._d1 else u1,
      ua[:, 3] if self._d2 else u2,
      x,
    )

  def _pdf_c(
    self, xp: Namespace[ArrayT], a: ArrayT, b: ArrayT, x: ArrayT | None
  ) -> ArrayT:
    return pair_eval(self._pdf_raw, xp.stack([a, b], axis=-1), x=x)

  def _cdf_c(
    self, xp: Namespace[ArrayT], a: ArrayT, b: ArrayT, x: ArrayT | None
  ) -> ArrayT:
    return pair_eval(self._cdf_raw, xp.stack([a, b], axis=-1), x=x)

  def _h1_c(
    self, xp: Namespace[ArrayT], a: ArrayT, b: ArrayT, x: ArrayT | None
  ) -> ArrayT:
    return pair_eval(self._hfunc1_raw, xp.stack([a, b], axis=-1), x=x)

  def _h2_c(
    self, xp: Namespace[ArrayT], a: ArrayT, b: ArrayT, x: ArrayT | None
  ) -> ArrayT:
    return pair_eval(self._hfunc2_raw, xp.stack([a, b], axis=-1), x=x)

  @staticmethod
  def _take(value: ArrayT | None, mask: BoolArray) -> ArrayT | None:
    """Select rows from an optional conditioning matrix."""
    return None if value is None else value[mask]

  @staticmethod
  def _quotient(
    xp: Namespace[ArrayT], num: ArrayT, delta: ArrayT, fallback: ArrayT
  ) -> ArrayT:
    """``|num / delta|`` over a wide-enough atom, else ``|fallback|``."""
    wide = delta > DELTA_MIN
    safe = xp.where(wide, delta, xp.ones_like(delta))
    return xp.abs(xp.where(wide, num / safe, fallback))

  def _interval(
    self,
    u_cond: ArrayT,
    lo: ArrayT,
    hi: ArrayT,
    cond_var: int,
    x: ArrayT | None,
  ) -> ArrayT:
    """``P(lo < U_free <= hi | U_cond = u_cond)``, a mixed edge's numerator."""
    return pair_eval(self.cond_interval_prob, u_cond, lo, hi, cond_var, x=x)

  def _rect(
    self,
    a1: ArrayT,
    b1: ArrayT,
    a2: ArrayT,
    b2: ArrayT,
    x: ArrayT | None,
  ) -> ArrayT:
    """``P((a1, b1] x (a2, b2])``, through whichever route the pair declares."""
    return pair_eval(self.rect_prob, a1, b1, a2, b2, x=x)

  def _strip(
    self,
    xp: Namespace[ArrayT],
    a1: ArrayT,
    b1: ArrayT,
    b2: ArrayT,
    x: ArrayT | None,
    axis: int,
  ) -> ArrayT:
    """``P((a1, b1] x (0, b2])`` for ``axis=1``, transposed for ``axis=2``.

    The rectangle anchored at the origin, which an h-function's numerator is.
    A zero bound is the distribution's own lower limit, so the second pair of
    corners contributes nothing and the generic route collapses to the same
    two-term difference it always was.
    """
    zero = xp.zeros_like(b2)
    if axis == 1:
      return self._rect(a1, b1, zero, b2, x)
    return self._rect(zero, b2, a1, b1, x)

  def _pdf_mixed(
    self,
    xp: Namespace[ArrayT],
    u1: ArrayT,
    u2: ArrayT,
    u1m: ArrayT,
    u2m: ArrayT,
    x: ArrayT | None,
    *,
    discrete: int,
  ) -> ArrayT:
    """Evaluate only the quotient or derivative each row requires."""
    delta = xp.abs((u1 - u1m) if discrete == 1 else (u2 - u2m))
    wide = delta > DELTA_MIN
    out = xp.empty_like(delta)
    if bool(xp.any(wide)):
      x_wide = self._take(x, wide)
      if discrete == 1:
        # The discrete argument is integrated over its atom, the continuous
        # one is the coordinate conditioned on.
        num = self._interval(u2[wide], u1m[wide], u1[wide], 2, x_wide)
      else:
        num = self._interval(u1[wide], u2m[wide], u2[wide], 1, x_wide)
      out[wide] = num / delta[wide]
    narrow = ~wide
    if bool(xp.any(narrow)):
      out[narrow] = self._pdf_c(
        xp,
        0.5 * (u1[narrow] + u1m[narrow]),
        0.5 * (u2[narrow] + u2m[narrow]),
        self._take(x, narrow),
      )
    return xp.abs(out)

  def _pdf_d_d(
    self,
    xp: Namespace[ArrayT],
    u1: ArrayT,
    u2: ArrayT,
    u1m: ArrayT,
    u2m: ArrayT,
    x: ArrayT | None,
  ) -> ArrayT:
    """Rectangle probability per unit area, with the degenerate fallbacks."""
    d1, d2 = xp.abs(u1 - u1m), xp.abs(u2 - u2m)
    m1, m2 = 0.5 * (u1 + u1m), 0.5 * (u2 + u2m)
    narrow1, narrow2 = d1 < DELTA_MIN, d2 < DELTA_MIN
    both = xp.where(d1 > d2, d1, d2) < DELTA_MIN
    only1 = narrow1 & ~both
    only2 = narrow2 & ~both
    wide = ~(both | only1 | only2)
    out = xp.empty_like(d1)
    if bool(xp.any(wide)):
      out[wide] = self._rect(
        u1m[wide],
        u1[wide],
        u2m[wide],
        u2[wide],
        self._take(x, wide),
      ) / (d1[wide] * d2[wide])
    if bool(xp.any(only1)):
      x_only1 = self._take(x, only1)
      # A collapsed argument is held at the atom's midpoint in both terms.
      out[only1] = (
        self._h1_c(xp, m1[only1], u2[only1], x_only1)
        - self._h1_c(xp, m1[only1], u2m[only1], x_only1)
      ) / d2[only1]
    if bool(xp.any(only2)):
      x_only2 = self._take(x, only2)
      out[only2] = (
        self._h2_c(xp, u1[only2], m2[only2], x_only2)
        - self._h2_c(xp, u1m[only2], m2[only2], x_only2)
      ) / d1[only2]
    if bool(xp.any(both)):
      out[both] = self._pdf_c(xp, m1[both], m2[both], self._take(x, both))
    return xp.abs(out)

  def _prep_args(self, u: ArrayT) -> ArrayT:
    """Place ``u``, check its width, and clamp it into the unit square.

    The three steps a copula argument needs, in the one order that is correct:
    placement (``_prep``), then the two-column layout the contract
    specifies, then the domain clamp at the working precision. Covariates go
    through ``_prep`` alone, being reals rather than copula arguments.

    The admissible width follows the pair's own ``var_types``: two columns
    while both arguments are continuous, and the expanded four once either is
    declared discrete.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pseudo-observations to prepare.

    Returns
    -------
    array, shape (n, 2), dtype float
        ``u`` placed and clamped.

    Raises
    ------
    ValueError
        If ``u`` is not two-dimensional with exactly two columns.
    """
    ua: Any = self._layout(self._prep(u))
    return cast("ArrayT", trim(ua))

  def _layout(self, ua: ArrayT) -> ArrayT:
    """Require the layout this pair's variable types imply.

    Two columns for a continuous pair, and the expanded four
    ``[u1, u2, u1^-, u2^-]`` once either argument is declared discrete --
    which is what the vine's cascades hand a discrete slot. The compact
    ``(n, 2 + k)`` form ``Bicop`` also accepts is expanded by the caller
    rather than here, so a pair copula sees one width per declaration.
    """
    a: Any = ua
    expected = 4 if (self._d1 or self._d2) else 2
    if getattr(a, "ndim", None) != 2 or int(a.shape[1]) != expected:
      raise ValueError(
        f"u must have shape (n, {expected}) for var_types="
        f"{list(self._var_types)}; got {tuple(getattr(a, 'shape', ()))}"
      )
    return ua

  def plot(
    self,
    plot_type: str = "surface",
    margin_type: str = "unif",
    xylim: tuple[float, float] | None = None,
    grid_size: int | None = None,
    *,
    x: ArrayT | None = None,
  ) -> None:
    bicop_plot(
      self, plot_type, margin_type, xylim, grid_size, x=x, place=self._prep
    )

  def __repr__(self) -> str:
    return f"{type(self).__name__}()"


BicopBase.__doc__ = (BicopBase.__doc__ or "") + _BICOP_EXAMPLE

#: Composed, so the shared half is written once; `test_docs_examples.py`
#: runs the project's numpydoc checks over the result.
BicopBase.plot.__doc__ = (
  BICOP_PLOT_SUMMARY
  + """
    The evaluation grid is the one array this class manufactures from nothing,
    so it is placed through ``_prep`` before the pair copula sees it -- which
    means a pair copula on PyTorch plots without converting anything inside
    its own ``pdf``.

    Parameters
    ----------
"""
  + BICOP_PLOT_PARAMS
  + """    x : array, shape (p,) or (1, p), or None, optional
        One covariate row, for a conditional pair copula. Such a pair copula
        is a different surface at every covariate value, so a plot shows the
        one at this value; the row is repeated across the grid.

    Returns
    -------
    None
        The figure is drawn with matplotlib.

    Raises
    ------
    ValueError
        If ``x`` is not a single covariate row.
    TypeError
        If ``x`` is given and the pair copula's ``pdf`` takes none. A pair
        copula declares its covariates by its signature rather than by a
        flag, so the refusal is the argument binding's -- and it is a refusal,
        not an oversight: drawing the unconditional density under a
        conditional-looking call is the outcome it exists to prevent.

    See Also
    --------
    pyvinecopulib.core.Bicop.plot : The same plot on a fitted family.
"""
)
