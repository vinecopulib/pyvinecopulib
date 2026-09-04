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

``_pair_eval`` lives here too, shared with the vine cascades: the rule that a
pair-copula method is handed ``x`` only when there is one, so a pair copula
accepting no covariates -- ``Bicop`` above all -- stays a valid host for a
simplified vine.

Array values are handled as ``Any`` inside the numeric bodies per the
``pyvinecopulib.core`` typing policy (``array_api_compat``, which resolves the
array namespace, is itself untyped); the generic ``ArrayT`` lives on the public
signatures.
"""

from __future__ import annotations

from abc import ABC
from typing import Any, Callable, Optional, Self, cast

from array_api_compat import array_namespace

from ._rootfind import solve_increasing
from ._validation import validate_covariates
from .protocols import ArrayT, BicopLike, _BICOP_EXAMPLE

__all__ = ["BicopBase"]


def _pair_eval(method: Callable[..., Any], u: Any, x: Optional[Any]) -> Any:
  """Evaluate a pair-copula method, forwarding ``x`` only when there is one.

  A pair copula that takes no conditioning argument -- ``Bicop`` above all --
  is a valid host for a simplified vine, so an absent context must not reach
  it as a keyword. When there *is* a conditioning matrix, it is passed by
  keyword, which is what makes a pair copula that cannot accept one fail
  loudly instead of binding it to some unrelated positional parameter.
  """
  return method(u) if x is None else method(u, x=x)


class BicopBase(BicopLike[ArrayT], ABC):
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
    pair in :class:`~pyvinecopulib.core.DiscretePair` to sit on such an edge.

  ``TorchBicop`` is the reference subclass: it supplies ``cdf``, both inverses,
  ``sample`` and ``flip`` natively, and fits its density grid in :meth:`fit`.

  See Also
  --------
  pyvinecopulib.core.BicopLike : The contract this implements.
  pyvinecopulib.torch.TorchBicop : A concrete (grid / TLL) subclass.
  """

  def loglik(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Total log-likelihood ``sum(log c(u))`` of the pair at ``u``.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (), dtype float
        The summed log-density, carrying gradients wherever the array library
        tracks them.
    """
    validate_covariates(x, int(cast(Any, u).shape[0]))
    dens: Any = self.pdf(u, x=x)
    xp = array_namespace(dens)
    return cast(ArrayT, xp.sum(xp.log(dens)))

  def hinv1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of ``hfunc1`` in its second argument.

    Solved numerically, so a subclass needs only ``hfunc1``; override it where
    the family inverts in closed form.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the conditioning value ``u1``; column 1 is the level to
        invert.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    ua: Any = u
    validate_covariates(x, int(ua.shape[0]))
    xp = array_namespace(ua)
    u1, p = ua[:, 0], ua[:, 1]
    return cast(
      ArrayT,
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
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation; ignored by an
        unconditional pair copula.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    ua: Any = u
    validate_covariates(x, int(ua.shape[0]))
    xp = array_namespace(ua)
    p, u2 = ua[:, 0], ua[:, 1]
    return cast(
      ArrayT,
      solve_increasing(
        lambda v: self.hfunc2(xp.stack([v, u2], axis=-1), x=x), p
      ),
    )

  def cdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Raise; override to give the pair copula a distribution ``C(u)``.

    Needed only to host the pair on a discrete edge, whose h-functions are
    difference quotients of the distribution function: add a ``cdf``, then
    wrap the pair in :class:`~pyvinecopulib.core.DiscretePair`. Nothing else
    asks for one -- a vine's own ``cdf`` is evaluated by Monte-Carlo
    simulation, which needs no per-pair distribution.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, k), or None, optional
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
    validate_covariates(x, int(cast(Any, u).shape[0]))
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
    controls: Optional[Any] = None,
    var_types: Optional[list[str]] = None,
  ) -> Self:
    """Construct a pair copula and select it from data.

    ``cls().select(u, ...)``, so a subclass that implements :meth:`fit` gets
    this for free -- which is what lets a vine name this class as its
    ``bicop_class`` and fit one pair per edge.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pseudo-observations in ``[0, 1]^2``.
    controls : object, or None, optional
        Fit configuration, in whatever form the subclass accepts.
    var_types : list of str, or None, optional
        The two variable types of the edge this pair sits on, ``"c"``
        (continuous) or ``"d"`` (discrete) each. ``None`` means both are
        continuous.

    Returns
    -------
    BicopBase
        The fitted pair copula.

    See Also
    --------
    select : Choose a family for an already-constructed pair copula, in place.
    fit : Estimate the current family's parameters, leaving the family alone.
    """
    return cls().select(u, controls, var_types)

  def fit(
    self,
    u: ArrayT,
    /,
    controls: Optional[Any] = None,
    var_types: Optional[list[str]] = None,
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
    controls : object, or None, optional
        Fit configuration, in whatever form the subclass accepts.
    var_types : list of str, or None, optional
        The two variable types of the edge this pair sits on. ``None`` means
        both are continuous.

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
    controls: Optional[Any] = None,
    var_types: Optional[list[str]] = None,
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
    controls : object, or None, optional
        Fit configuration, in whatever form the subclass accepts.
    var_types : list of str, or None, optional
        The two variable types of the edge this pair sits on. ``None`` means
        both are continuous.

    Returns
    -------
    BicopBase
        ``self``, so the call chains.

    See Also
    --------
    fit : Estimate the current family's parameters, leaving the family alone.
    from_data : Construct and select in one call.
    """
    if controls is None and var_types is None:
      # Forwarded only when there is something to forward, so a subclass whose
      # `fit` takes neither still works through `select`.
      return self.fit(u)
    return self.fit(u, controls, var_types)

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
    x : array, shape (n, k), or None, optional
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
    validate_covariates(x, n)
    base_u: Any = self._sample_uniform(n, qrng, list(seeds) if seeds else [])
    xp = array_namespace(base_u)
    u2: Any = self.hinv1(base_u, x=x)
    return cast(ArrayT, xp.stack([base_u[:, 0], u2], axis=-1))

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> ArrayT:
    """Draw the ``(n, 2)`` base uniforms ``sample`` transforms.

    Raises unless a subclass overrides it: NumPy and PyTorch differ on RNG, so
    this is the one hook with no array-agnostic default, and overriding it is
    all ``sample`` needs. Named after the ``sample_uniform`` free function in
    ``pyvinecopulib.utils``.
    """
    raise NotImplementedError(
      f"{type(self).__name__} does not implement _sample_uniform; override it "
      "to enable sample()."
    )

  def plot(
    self,
    plot_type: str = "surface",
    margin_type: str = "unif",
    xylim: Optional[tuple[float, float]] = None,
    grid_size: Optional[int] = None,
  ) -> None:
    """Plot the pair-copula density, as a contour or a 3-D surface.

    Mirrors ``Bicop.plot()``.

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

    Returns
    -------
    None
        The figure is drawn with matplotlib.
    """
    from .._python_helpers.bicop import bicop_plot

    bicop_plot(self, plot_type, margin_type, xylim, grid_size)

  def __repr__(self) -> str:
    return f"{type(self).__name__}()"


BicopBase.__doc__ = (BicopBase.__doc__ or "") + _BICOP_EXAMPLE
