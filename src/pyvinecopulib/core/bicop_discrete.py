"""A continuous pair copula evaluated on a discrete or mixed edge.

A discrete variable is described by two numbers, ``F(x)`` and its left limit
``F(x^-)``, and the copula quantities that are derivatives in a continuous
argument become difference quotients over the atom ``[F(x^-), F(x)]``.
:class:`DiscreteBicop` reads a four-column argument and supplies that surface
from the *continuous* ``pdf`` / ``cdf`` / ``hfunc1`` / ``hfunc2`` of the pair it
wraps, or forwards to the pair's own discrete surface where it has one.

Keeping the difference quotients here rather than on the pair copulas is what
lets :class:`~pyvinecopulib.core.BicopLike` stay a purely continuous, two-column
contract: a custom pair copula gains discrete support by implementing ``cdf``,
and needs to know nothing else.
"""

from __future__ import annotations

from types import ModuleType
from typing import Any, Optional, Protocol, cast

from array_api_compat import array_namespace

from ._covariates import pair_eval, prepare_covariates
from ._rootfind import solve_increasing
from ._trim import trim
from ._vinecop_discrete import check_var_types
from .bicop_base import (
  BicopBase,
  cond_interval_prob_from_hfunc,
  flip_of,
  rect_prob_from_cdf,
)
from .protocols import ArrayT, BicopLike

__all__ = ["DiscreteBicop", "continuous_view"]


class _ContinuousBicop(Protocol[ArrayT]):
  """The unconditional evaluations :class:`DiscreteBicop` calls on its pair.

  **Why this is not** :class:`~pyvinecopulib.core.BicopLike`. That contract
  declares a keyword-only ``x`` on every method, which ``Bicop`` does not have
  -- its own optional argument is a per-row ``parameters`` matrix -- so
  ``Bicop`` satisfies it only *nominally*, by method name, and a type checker
  rejects it. Since ``Bicop`` is the pair this class most often wraps, the
  annotation has to be a surface ``Bicop`` satisfies structurally. ``cdf`` is
  required here and optional there, for the same reason it is needed at all: an
  atom's probability is built from it wherever the pair offers nothing better.

  ``sample`` is absent although it *is* called, and the omission is forced: the
  call has two shapes, with and without ``x``, and only a pair that declares one
  can be typed against the second. Declaring it here without ``x`` makes the
  conditional call a type error, and declaring it with one puts ``Bicop`` back
  outside the protocol -- so the one call site casts instead. A conditioning
  matrix is forwarded dynamically either way (see ``pair_eval``), which is what
  makes a pair that cannot accept one fail loudly rather than silently.
  """

  def pdf(self, u: ArrayT) -> ArrayT: ...

  def cdf(self, u: ArrayT) -> ArrayT: ...

  def hfunc1(self, u: ArrayT) -> ArrayT: ...

  def hfunc2(self, u: ArrayT) -> ArrayT: ...

  def hinv1(self, u: ArrayT) -> ArrayT: ...

  def hinv2(self, u: ArrayT) -> ArrayT: ...


#: Atom width below which a difference quotient is numerically unstable and the
#: derivative is used instead (``AbstractBicop``'s threshold).
DELTA_MIN: float = 5e-5


def continuous_view(pair: object) -> _ContinuousBicop[ArrayT]:
  """Return ``pair`` evaluated as a continuous copula, when it can be.

  A pair copula may carry its own variable types -- ``Bicop`` does -- in which
  case its two-column evaluation surface is unavailable. Such a pair is used
  through its optional ``with_var_types`` capability, whose default argument is
  the continuous declaration; one that does not advertise it is already
  continuous-only and is returned as is.

  Parameters
  ----------
  pair : BicopLike
      The pair copula to view.

  Returns
  -------
  BicopLike
      The continuous view, or ``pair`` itself.

  Raises
  ------
  ValueError
      If ``pair`` declares discrete variables but offers no continuous view.
  """
  view = getattr(pair, "with_var_types", None)
  if view is not None:
    return cast("_ContinuousBicop[ArrayT]", view())
  types = getattr(pair, "var_types", None)
  if types is not None and any(t != "c" for t in types):
    raise ValueError(
      f"{type(pair).__name__} declares var_types={list(types)} but has no "
      "with_var_types(); the inverse Rosenblatt cascade evaluates every pair "
      "copula as continuous, since it produces the values a left limit would "
      "be taken of. Host the continuous copula in a DiscreteBicop instead, or "
      "add with_var_types()."
    )
  return cast("_ContinuousBicop[ArrayT]", pair)


class DiscreteBicop(BicopBase[ArrayT]):
  """A continuous pair copula evaluated as a mixed-discrete one.

  Reads the four-column layout ``[F(u1), F(u2), F(u1^-), F(u2^-)]`` and returns
  the mixed-discrete density and h-functions, built from the wrapped pair's
  continuous ``pdf`` / ``cdf`` / ``hfunc1`` / ``hfunc2``: in a discrete argument
  a derivative becomes a difference quotient over the atom, and the derivative
  itself is used where the atom is narrower than ``5e-5``. A continuous
  argument's left-limit column is ignored -- it equals its value.

  This is how a pair copula that knows nothing about discreteness is hosted on a
  discrete edge of a :class:`~pyvinecopulib.core.VinecopBase`: wrap it in
  ``get_pair_copula`` with the types the vine derives for that slot,
  :meth:`~pyvinecopulib.core.VinecopBase.pair_var_types`. The h-functions of a
  discrete argument, and the density of a pair with two, are built from the
  copula's distribution function, so a pair hosted on a discrete edge must
  implement ``cdf`` -- the cascades evaluate h-functions at every edge.

  Parameters
  ----------
  pair : BicopLike
      The pair copula to wrap; only its two-column ``pdf`` / ``cdf`` / ``hfunc1``
      / ``hfunc2`` are used, through ``with_var_types()`` when it advertises one.
  var_types : tuple of str
      The edge's two variable types, ``"c"`` or ``"d"``.

  See Also
  --------
  pyvinecopulib.core.VinecopBase.pair_var_types : The types to wrap a slot with.
  pyvinecopulib.core.Bicop : The reference pair copula, discrete-aware itself.

  Notes
  -----
  ``examples/10_extending_pyvinecopulib.ipynb`` hosts a custom pair copula on a
  discrete edge and checks it against ``Vinecop``.

  The output matches ``Bicop``'s own discrete surface bit-for-bit, for every
  family including the nonparametric ``tll``, and satisfies
  ``sum_atoms c * (u1 - u1^-) == 1`` exactly. That identity is what pins the
  quotients: they telescope over the atoms, so no tolerance argument is needed
  to say which of two candidate surfaces is right.
  """

  #: The batched grid surface has no distribution-function lookup, which the
  #: discrete h-functions need, so a wrapped pair never takes that fast path.
  supports_batched: bool = False

  def __init__(
    self, pair: _ContinuousBicop[ArrayT], var_types: tuple[str, str]
  ) -> None:
    self._pair = continuous_view(pair)
    self.var_types = list(var_types)

  @property
  def var_types(self) -> list[str]:
    """The two variable types, ``"c"`` (continuous) or ``"d"`` (discrete).

    Returns
    -------
    list of str
        The wrapped edge's ``[type1, type2]``.
    """
    return list(self._var_types)

  @var_types.setter
  def var_types(self, value: list[str]) -> None:
    types = check_var_types(list(value), 2)
    self._var_types = types
    self._d1 = types[0] == "d"
    self._d2 = types[1] == "d"
    # A pair that models atoms itself needs none of the quotients below, so it
    # is asked for its own surface and this class forwards to it. Resolved here
    # rather than per call because the view depends on the types.
    view = getattr(self._pair, "with_var_types", None)
    self._native = view(list(types)) if view is not None else None
    # The two rectangle routines, read once, for a pair that has no whole
    # discrete surface but can still measure a rectangle without differencing
    # -- `TorchTllBicop` above all. `BicopBase` supplies both, so the `getattr`
    # is for a foreign object implementing `BicopLike` directly.
    self._rect_prob = getattr(self._pair, "rect_prob", None)
    self._cond_prob = getattr(self._pair, "cond_interval_prob", None)

  def with_var_types(
    self, var_types: tuple[str, str] = ("c", "c")
  ) -> BicopLike[ArrayT]:
    """The same wrapped pair under different variable types.

    The continuous declaration is the default, and gives back the pair copula
    this wraps rather than a wrapper around it: with no atoms to integrate over
    there is nothing for the quotients to do.

    Parameters
    ----------
    var_types : tuple of str, default=("c", "c")
        The two types, ``"c"`` or ``"d"``.

    Returns
    -------
    BicopLike
        The wrapped pair copula itself when both types are continuous, else a
        wrapper around it reading ``var_types``.
    """
    if all(t == "c" for t in var_types):
      return cast("BicopLike[ArrayT]", self._pair)
    return DiscreteBicop(self._pair, (var_types[0], var_types[1]))

  def flip(self) -> "DiscreteBicop[ArrayT]":
    """The argument-swapped copula, with its variable types swapped too.

    Returns
    -------
    DiscreteBicop
        A wrapper around the flipped pair, reading ``[u2, u1, u2^-, u1^-]``.
    """
    flipped = flip_of(self._pair)
    return DiscreteBicop(flipped, (self._var_types[1], self._var_types[0]))

  def sample(
    self,
    n: int,
    *,
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT:
    """Draw from the wrapped continuous copula.

    Discreteness changes how a copula density is evaluated at atoms, not the
    continuous latent uniforms the copula samples. The wrapped pair therefore
    owns the RNG and inverse Rosenblatt transform.

    Parameters
    ----------
    n : int
        Number of samples.
    x : array, shape (n, p), or None, optional
        External covariates for a conditional draw.
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random ones.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n, 2), dtype float
        Samples in the unit square.
    """
    x = prepare_covariates(self, x, n)
    method = cast("BicopLike[ArrayT]", self._pair).sample
    draw_seeds = list(seeds) if seeds else []
    if x is None:
      return method(n, qrng=qrng, seeds=draw_seeds)
    return method(n, x=x, qrng=qrng, seeds=draw_seeds)

  def __repr__(self) -> str:
    return f"DiscreteBicop({self._pair!r}, var_types={list(self._var_types)})"

  # --- argument handling ------------------------------------------------ #
  # The helpers below compute on the columns rather than forwarding them --
  # differences, quotients, comparisons -- so they hold their arguments as
  # `Any`, which is what `protocols.py` says an unbounded `ArrayT` requires of
  # a body that does arithmetic. Every one is private to this class; the
  # module-level functions above, which only forward, are typed `ArrayT`.
  def _layout(self, ua: ArrayT) -> ArrayT:
    """Require the layout this edge's variable types imply.

    Applied on both routes, so the message names this class whichever one the
    wrapped pair puts the evaluation on.
    """
    any_u: Any = ua
    expected = 4 if (self._d1 or self._d2) else 2
    if any_u.ndim != 2 or int(any_u.shape[1]) != expected:
      raise ValueError(
        f"u must have shape (n, {expected}) for var_types="
        f"{list(self._var_types)}; got {tuple(any_u.shape)}"
      )
    return ua

  def _split(self, u: ArrayT) -> tuple[ModuleType, Any, Any, Any, Any]:
    """Namespace plus the two values and their left limits."""
    ua: Any = self._layout(u)
    xp = array_namespace(ua)
    # Trimmed before anything is subtracted, as ``Bicop::prep_for_abstract``
    # does: an h-function feeding the next tree may land exactly on 0 or 1, and
    # the atom width in the denominator has to be the width of the trimmed atom.
    ut: Any = trim(u, xp)
    u1, u2 = ut[:, 0], ut[:, 1]
    # A continuous argument's left limit is its own value
    # (``Bicop::format_data``), so the cascade's column for it is never read.
    return (
      xp,
      u1,
      u2,
      ut[:, 2] if self._d1 else u1,
      ut[:, 3] if self._d2 else u2,
    )

  def _pdf(self, xp: ModuleType, a: Any, b: Any, x: Optional[ArrayT]) -> Any:
    return pair_eval(self._pair.pdf, xp.stack([a, b], axis=-1), x=x)

  def _cdf(self, xp: ModuleType, a: Any, b: Any, x: Optional[ArrayT]) -> Any:
    return pair_eval(self._pair.cdf, xp.stack([a, b], axis=-1), x=x)

  def _h1(self, xp: ModuleType, a: Any, b: Any, x: Optional[ArrayT]) -> Any:
    return pair_eval(self._pair.hfunc1, xp.stack([a, b], axis=-1), x=x)

  def _h2(self, xp: ModuleType, a: Any, b: Any, x: Optional[ArrayT]) -> Any:
    return pair_eval(self._pair.hfunc2, xp.stack([a, b], axis=-1), x=x)

  @staticmethod
  def _take(value: Optional[Any], mask: Any) -> Optional[Any]:
    """Select rows from an optional conditioning matrix."""
    return None if value is None else value[mask]

  @staticmethod
  def _quotient(xp: ModuleType, num: Any, delta: Any, fallback: Any) -> Any:
    """``|num / delta|`` over a wide-enough atom, else ``|fallback|``."""
    wide = delta > DELTA_MIN
    safe = xp.where(wide, delta, xp.ones_like(delta))
    return xp.abs(xp.where(wide, num / safe, fallback))

  def _pdf_mixed(
    self,
    xp: ModuleType,
    u1: Any,
    u2: Any,
    u1m: Any,
    u2m: Any,
    x: Optional[ArrayT],
    *,
    discrete: int,
  ) -> Any:
    """Evaluate only the quotient or derivative each row requires."""
    delta = xp.abs((u1 - u1m) if discrete == 1 else (u2 - u2m))
    wide = delta > DELTA_MIN
    out = xp.empty_like(delta)
    if bool(xp.any(wide)):
      x_wide = self._take(x, wide)
      if discrete == 1:
        # The discrete argument is integrated over its atom, the continuous one
        # is the coordinate conditioned on.
        num = self._interval(u2[wide], u1m[wide], u1[wide], 2, x_wide)
      else:
        num = self._interval(u1[wide], u2m[wide], u2[wide], 1, x_wide)
      out[wide] = num / delta[wide]
    narrow = ~wide
    if bool(xp.any(narrow)):
      out[narrow] = self._pdf(
        xp,
        0.5 * (u1[narrow] + u1m[narrow]),
        0.5 * (u2[narrow] + u2m[narrow]),
        self._take(x, narrow),
      )
    return xp.abs(out)

  def _interval(
    self,
    u_cond: Any,
    lo: Any,
    hi: Any,
    cond_var: int,
    x: Optional[ArrayT],
  ) -> Any:
    """``P(lo < U_free <= hi | U_cond = u_cond)``, the mixed edge's numerator.

    The conditional counterpart of :meth:`_rect`, and the same choice of route:
    the pair's own ``cond_interval_prob`` where it declares one, and otherwise
    the difference of two h-function values. The second clamps each value into
    the open unit interval and the mass does not, which is most of the
    difference between them at a narrow atom.
    """
    if self._cond_prob is not None:
      return pair_eval(self._cond_prob, u_cond, lo, hi, cond_var, x=x)
    h = self._pair.hfunc1 if cond_var == 1 else self._pair.hfunc2
    return cond_interval_prob_from_hfunc(h, u_cond, lo, hi, cond_var, x=x)

  def _rect(
    self,
    xp: ModuleType,
    a1: Any,
    b1: Any,
    a2: Any,
    b2: Any,
    x: Optional[ArrayT],
  ) -> Any:
    """``P((a1, b1] x (a2, b2])``, through the pair's own route where it has one.

    A pair copula that can measure a rectangle without the cancellation a
    four-corner difference carries declares ``rect_prob``, and the density is
    then more accurate at exactly the atom widths a vine's inner trees reach:
    differencing amplifies an absolute error by ``4 / (w1 w2)``, reading the
    mass by ``1 / w2`` alone. The reference pair copula makes the same choice
    per family, so taking it here is what keeps the two in step.
    """
    del xp
    if self._rect_prob is not None:
      return pair_eval(self._rect_prob, a1, b1, a2, b2, x=x)
    return rect_prob_from_cdf(self._pair.cdf, a1, b1, a2, b2, x=x)

  def _strip(
    self,
    xp: ModuleType,
    a1: Any,
    b1: Any,
    b2: Any,
    x: Optional[ArrayT],
    axis: int,
  ) -> Any:
    """``P((a1, b1] x (0, b2])`` for ``axis=1``, transposed for ``axis=2``.

    The rectangle anchored at the origin, which an h-function's numerator is.
    A zero bound is the distribution's own lower limit, so the second pair of
    corners contributes nothing and the generic route collapses to the same
    two-term difference it always was. See :meth:`_rect` for the choice of
    route.
    """
    zero = xp.zeros_like(b2)
    if axis == 1:
      return self._rect(xp, a1, b1, zero, b2, x)
    return self._rect(xp, zero, b2, a1, b1, x)

  # --- the mixed-discrete surface --------------------------------------- #
  def pdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Density with respect to each argument's own reference measure.

    A continuous argument contributes a derivative, a discrete one the
    probability of its atom, so a mixed pair gives a difference quotient and two
    discrete arguments the rectangle probability, each divided by the atom
    widths.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        ``[u1, u2]`` when both arguments are continuous, else
        ``[u1, u2, u1^-, u2^-]``.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        Density values.
    """
    # A pair that models atoms itself answers for the whole surface: the
    # quotients below exist to build what it already has, so building them
    # again from its own `cdf` would be both slower and less accurate.
    if self._native is not None:
      return cast("ArrayT", pair_eval(self._native.pdf, self._layout(u), x=x))
    xp, u1, u2, u1m, u2m = self._split(u)
    if self._d1 and self._d2:
      return cast("ArrayT", self._pdf_d_d(xp, u1, u2, u1m, u2m, x))
    if self._d1:
      return cast(
        "ArrayT", self._pdf_mixed(xp, u1, u2, u1m, u2m, x, discrete=1)
      )
    if self._d2:
      return cast(
        "ArrayT", self._pdf_mixed(xp, u1, u2, u1m, u2m, x, discrete=2)
      )
    return cast("ArrayT", self._pdf(xp, u1, u2, x))

  def _pdf_d_d(
    self,
    xp: ModuleType,
    u1: Any,
    u2: Any,
    u1m: Any,
    u2m: Any,
    x: Optional[ArrayT],
  ) -> Any:
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
        xp,
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
        self._h1(xp, m1[only1], u2[only1], x_only1)
        - self._h1(xp, m1[only1], u2m[only1], x_only1)
      ) / d2[only1]
    if bool(xp.any(only2)):
      x_only2 = self._take(x, only2)
      out[only2] = (
        self._h2(xp, u1[only2], m2[only2], x_only2)
        - self._h2(xp, u1m[only2], m2[only2], x_only2)
      ) / d1[only2]
    if bool(xp.any(both)):
      out[both] = self._pdf(xp, m1[both], m2[both], self._take(x, both))
    return xp.abs(out)

  def hfunc1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """``P(U2 <= u2 | U1)``, conditioning on the atom when ``U1`` is discrete.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf`.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values.
    """
    if self._native is not None:
      return cast(
        "ArrayT", pair_eval(self._native.hfunc1, self._layout(u), x=x)
      )
    xp, u1, u2, u1m, _ = self._split(u)
    if not self._d1:
      return cast("ArrayT", self._h1(xp, u1, u2, x))
    # Conditioning on `u1^- < U1 <= u1` divides the rectangle probability by the
    # atom's width; the second argument enters at its value either way.
    return cast(
      "ArrayT",
      self._quotient(
        xp,
        self._strip(xp, u1m, u1, u2, x, axis=1),
        xp.abs(u1 - u1m),
        self._h1(xp, 0.5 * (u1 + u1m), u2, x),
      ),
    )

  def hfunc2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """``P(U1 <= u1 | U2)``, conditioning on the atom when ``U2`` is discrete.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf`.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values.
    """
    if self._native is not None:
      return cast(
        "ArrayT", pair_eval(self._native.hfunc2, self._layout(u), x=x)
      )
    xp, u1, u2, _, u2m = self._split(u)
    if not self._d2:
      return cast("ArrayT", self._h2(xp, u1, u2, x))
    return cast(
      "ArrayT",
      self._quotient(
        xp,
        self._strip(xp, u2m, u2, u1, x, axis=2),
        xp.abs(u2 - u2m),
        self._h2(xp, u1, 0.5 * (u2 + u2m), x),
      ),
    )

  def cdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Distribution function, which the left limits do not enter.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf`.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values.
    """
    xp, u1, u2, _, _ = self._split(u)
    return cast("ArrayT", self._cdf(xp, u1, u2, x))

  def hinv1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of :meth:`hfunc1` in its second argument.

    Analytic through the wrapped pair when the conditioning argument is
    continuous, and a monotone bisection of the mixed-discrete ``hfunc1``
    otherwise.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        Column ``0`` is the conditioning value and column ``1`` the level; see
        :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        The inverse values.
    """
    xp, u1, p, u1m, _ = self._split(u)
    if not self._d1:
      return cast(
        "ArrayT", pair_eval(self._pair.hinv1, xp.stack([u1, p], axis=-1), x=x)
      )
    return cast(
      "ArrayT",
      solve_increasing(
        lambda v: self.hfunc1(
          cast("ArrayT", xp.stack([u1, v, u1m, v], axis=-1)), x=x
        ),
        p,
      ),
    )

  def hinv2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of :meth:`hfunc2` in its first argument.

    Analytic through the wrapped pair when the conditioning argument is
    continuous, and a monotone bisection of the mixed-discrete ``hfunc2``
    otherwise.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        Column ``0`` is the level and column ``1`` the conditioning value; see
        :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        The inverse values.
    """
    xp, p, u2, _, u2m = self._split(u)
    if not self._d2:
      return cast(
        "ArrayT", pair_eval(self._pair.hinv2, xp.stack([p, u2], axis=-1), x=x)
      )
    return cast(
      "ArrayT",
      solve_increasing(
        lambda v: self.hfunc2(
          cast("ArrayT", xp.stack([v, u2, v, u2m], axis=-1)), x=x
        ),
        p,
      ),
    )
