"""Canonical partial implementation of the pair-copula contract.

``BicopBase`` is the array-agnostic (NumPy or PyTorch) base class for
``BicopLike``, and the short path to hosting a custom pair copula in a vine: a
subclass writes the density ``pdf`` and the two h-functions ``hfunc1`` /
``hfunc2``, and inherits the rest of the evaluation surface together with the
estimator surface ``fit`` / ``select`` / ``from_data``. The member inventory,
and which members ship a raising default, is on ``BicopBase`` itself.

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

from abc import ABC
from typing import Any, Optional, Self, TypeVar, cast

from array_api_compat import array_namespace

from ._bicop_plot import bicop_plot
from ._covariates import prepare
from ._placement import PlacementMixin, QrngUniformMixin
from ._trim import trim
from ._rootfind import solve_increasing

from .protocols import ArrayT, BicopLike, ControlsLike, _BICOP_EXAMPLE

__all__ = ["BicopBase"]

# A pair copula's own type, which `flip_of` hands back: every `flip` in the
# library returns a pair of the same kind as the one it was asked of, and a
# caller that established more about its pair than the contract requires --
# a `cdf`, for a pair on a discrete edge -- still has it afterwards. Unbounded
# on purpose: a foreign object that satisfies nothing at all is exactly what
# the raise below is for, so this cannot demand `BicopLike`.
_PairT = TypeVar("_PairT")


def flip_of(pair: _PairT) -> _PairT:
  """The argument-swapped pair, for a caller that has established it has one.

  ``flip`` is an optional capability on :class:`BicopLike` -- needed only to
  host a pair in structure *selection* or in a relabeling, never to evaluate
  one -- so the contract does not require it and the type checker will not let
  it be called unguarded. This is the single place that reads it, so the guard
  each caller relies on is named once rather than cast away four times:

  - ``VinecopBase.select`` refuses a pair class without one up front
    (``_check_selectable``), and probes the first fitted pair behind an opaque
    ``fit_edge``;
  - ``reorient`` and the reoriented view only reach slots of a vine that was
    selected or built with flippable pairs;
  - ``DiscreteBicop.flip`` delegates to the continuous pair it wraps.

  Parameters
  ----------
  pair : BicopLike
      The pair copula to flip.

  Returns
  -------
  BicopLike
      The pair with its two arguments swapped.

  Raises
  ------
  NotImplementedError
      If the pair has no ``flip``, which is what
      :class:`~pyvinecopulib.core.BicopBase` raises and what the guards above
      are there to turn into an earlier, clearer failure.
  """
  method = getattr(pair, "flip", None)
  if method is None:
    raise NotImplementedError(
      f"{type(pair).__name__} has no `flip` (the argument-swapped copula), "
      "which structure selection and relabeling need to reorient a pair onto "
      "its slot. Implement it, or supply a structure and fit along it."
    )
  return cast("_PairT", method())


class BicopBase(
  BicopLike[ArrayT], QrngUniformMixin[ArrayT], PlacementMixin, ABC
):
  """Canonical partial implementation of ``BicopLike``.

  A subclass writes three methods -- ``pdf``, the pair density, and ``hfunc1``
  / ``hfunc2``, its two conditional distributions -- and inherits the rest of
  the evaluation surface:

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
    evaluation along a fixed structure never asks for it.
  - :meth:`cdf`, needed on a **discrete** edge, whose h-functions are
    difference quotients of the distribution function. Add one and wrap the
    pair in :class:`~pyvinecopulib.core.DiscreteBicop` to sit on such an edge.

  ``TorchTllBicop`` is the reference subclass: it supplies ``cdf``, both inverses,
  ``sample`` and ``flip`` natively, and fits its density grid in :meth:`fit`.

  See Also
  --------
  pyvinecopulib.core.BicopLike : The contract this implements.
  pyvinecopulib.torch.TorchTllBicop : A concrete (grid / TLL) subclass.
  """

  def loglik(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Total log-likelihood ``sum(log c(u))`` of the pair at ``u``.

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
    x = prepare(self, x, int(cast("Any", u).shape[0]))
    # `u` is left to the subclass's own `pdf`, which is where the two-column
    # layout is checked -- the base cannot know whether a pair is on a
    # discrete edge, whose argument is four columns wide.
    dens: Any = self.pdf(u, x=x)
    xp = array_namespace(dens)
    return cast("ArrayT", xp.sum(xp.log(dens)))

  def hinv1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of ``hfunc1`` in its second argument.

    Solved numerically, so a subclass needs only ``hfunc1``; override it where
    the family inverts in closed form.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the conditioning value ``u1``; column 1 is the level to
        invert.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    ua: Any = self._prep_args(u)
    x = prepare(self, x, int(ua.shape[0]))
    xp = array_namespace(ua)
    u1, p = ua[:, 0], ua[:, 1]
    return cast(
      "ArrayT",
      solve_increasing(
        lambda v: self.hfunc1(xp.stack([u1, v], axis=-1), x=x), p
      ),
    )

  def hinv2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of ``hfunc2`` in its first argument.

    The counterpart of :meth:`hinv1`, on the other h-function.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the level to invert; column 1 is the conditioning value
        ``u2``.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    ua: Any = self._prep_args(u)
    x = prepare(self, x, int(ua.shape[0]))
    xp = array_namespace(ua)
    p, u2 = ua[:, 0], ua[:, 1]
    return cast(
      "ArrayT",
      solve_increasing(
        lambda v: self.hfunc2(xp.stack([v, u2], axis=-1), x=x), p
      ),
    )

  def cdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Raise; override to give the pair copula a distribution ``C(u)``.

    Needed only to host the pair on a discrete edge, whose h-functions are
    difference quotients of the distribution function: add a ``cdf``, then
    wrap the pair in :class:`~pyvinecopulib.core.DiscreteBicop`. Nothing else
    asks for one -- a vine's own ``cdf`` is evaluated by Monte-Carlo
    simulation, which needs no per-pair distribution.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values in ``[0, 1]`` -- only when a subclass overrides
        this method.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides a ``cdf``.
    """
    del u, x
    raise NotImplementedError(
      f"{type(self).__name__}.cdf is not defined; the vine cdf uses "
      "Monte-Carlo simulation and does not require a per-pair cdf. Implement it "
      "to host this pair copula on a discrete edge, whose h-functions are "
      "difference quotients of the distribution function."
    )

  @classmethod
  def from_data(
    cls,
    u: ArrayT,
    /,
    controls: Optional[ControlsLike] = None,
    *,
    var_types: Optional[list[str]] = None,
    x: Optional[ArrayT] = None,
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
    controls: Optional[ControlsLike] = None,
    *,
    var_types: Optional[list[str]] = None,
    x: Optional[ArrayT] = None,
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
    controls: Optional[ControlsLike] = None,
    *,
    var_types: Optional[list[str]] = None,
    x: Optional[ArrayT] = None,
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

  def flip(self) -> "BicopBase[ArrayT]":
    """Raise; override to return the pair with its arguments swapped.

    The flipped copula satisfies ``c'(u1, u2) = c(u2, u1)`` with the two
    h-functions (and their inverses) exchanged. It is required only to host the
    pair in structure *selection*, which reorients each selected pair onto its
    finalized slot (``VinecopBase.select()``); evaluation along a fixed
    structure never asks for it.

    Returns
    -------
    BicopBase
        The argument-swapped pair copula -- only when a subclass overrides this
        method.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides a ``flip``.
    """
    raise NotImplementedError(
      f"{type(self).__name__}.flip is not defined; implement it (return the "
      "argument-swapped copula) to host this pair in structure selection."
    )

  def sample(
    self,
    n: int,
    *,
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
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
    x = prepare(self, x, n)
    base_u: Any = self._sample_uniform(n, qrng, list(seeds) if seeds else [])
    xp = array_namespace(base_u)
    u2: Any = self.hinv1(base_u, x=x)
    return cast("ArrayT", xp.stack([base_u[:, 0], u2], axis=-1))

  #: Whether a vine may stack this pair into its grid-batched
  #: cascade. ``False`` here because the fast path reads an interpolation grid
  #: off each pair and a pair copula in general has none -- so a subclass opts
  #: in only if it exposes one. Declared rather than discovered, and declared
  #: on the base rather than left to a ``getattr`` default, so that the third
  #: state (declared ``False``) is distinguishable from "never heard of it"
  #: and a subclass author can find the flag without tripping its error.
  supports_batched: bool = False

  def _prep_args(self, u: ArrayT) -> ArrayT:
    """Place ``u``, check its width, and clamp it into the unit square.

    The three steps a copula argument needs, in the one order that is correct:
    placement (``_prep``), then the two-column layout the contract
    specifies, then the domain clamp at the working precision. Covariates go
    through ``_prep`` alone, being reals rather than copula arguments.

    A discrete edge is reached through
    :class:`~pyvinecopulib.core.DiscreteBicop`, which owns the four-column
    layout and hands each wrapped pair two columns at a time -- so this stays
    the continuous two-column contract.

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
    return cast("ArrayT", trim(array_namespace(ua), ua))

  def _layout(self, ua: ArrayT) -> ArrayT:
    """Check the two-column layout the pair-copula contract specifies."""
    a: Any = ua
    if getattr(a, "ndim", None) != 2 or int(a.shape[1]) != 2:
      raise ValueError(
        f"u must have shape (n, 2); got {tuple(getattr(a, 'shape', ()))}"
      )
    return ua

  def plot(
    self,
    plot_type: str = "surface",
    margin_type: str = "unif",
    xylim: Optional[tuple[float, float]] = None,
    grid_size: Optional[int] = None,
    *,
    x: Optional[ArrayT] = None,
  ) -> None:
    """Plot the pair-copula density, as a contour or a 3-D surface.

    Mirrors ``Bicop.plot()``, and adds ``x`` for a conditional pair copula.

    The evaluation grid is the one place this class manufactures an array from
    nothing, so it is the one place a subclass could be handed the wrong array
    type. It is placed through ``_prep`` first, which means a pair copula
    on PyTorch plots without converting anything inside its own ``pdf``.

    Parameters
    ----------
    plot_type : str, default="surface"
        ``"surface"`` for a 3-D surface, ``"contour"`` for a contour plot.
    margin_type : str, default="unif"
        Margins the density is shown on: ``"unif"``, ``"norm"`` or ``"exp"``.
    xylim : tuple of float, or None, optional
        Axis limits; ``None`` uses a default per ``margin_type``.
    grid_size : int, or None, optional
        Number of grid points per axis; ``None`` uses a default per
        ``plot_type``.
    x : array, shape (p,) or (1, p), or None, optional
        One covariate row, for a conditional pair copula. The density is a
        different surface at every covariate value, so a plot shows the slice
        at this one; the row is repeated across the grid. A pair copula that
        reads no covariates refuses it, rather than drawing an unconditional
        surface under a conditional-looking call.

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
    """
    bicop_plot(
      self, plot_type, margin_type, xylim, grid_size, x=x, place=self._prep
    )

  def __repr__(self) -> str:
    return f"{type(self).__name__}()"


BicopBase.__doc__ = (BicopBase.__doc__ or "") + _BICOP_EXAMPLE
