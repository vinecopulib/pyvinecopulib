"""Canonical partial implementation of the pair-copula contract.

:class:`BicopBase` is the array-agnostic (NumPy / PyTorch) base class for
:class:`~pyvinecopulib.core.BicopLike`. A subclass supplies only the three
primitives ``pdf`` / ``hfunc1`` / ``hfunc2`` and inherits ``hinv1`` / ``hinv2``
(numerical inversion of the h-functions), ``simulate`` (inverse Rosenblatt of
the pair), ``loglik``, ``plot`` and ``__repr__`` for free; each is overridable
when a native exact form exists (as :class:`~pyvinecopulib.torch.TorchBicop`
does for ``cdf`` / ``hinv`` / ``simulate``). ``cdf`` raises unless overridden —
the vine CDF is Monte-Carlo, so a per-pair ``cdf`` is rarely needed.

Written against the Array API (:func:`array_api_compat.array_namespace`) so the
same code runs on numpy and torch (via ``pdf`` / ``hfunc`` outputs); the numeric
bodies handle arrays as ``Any`` (the generic ``ArrayT`` lives on the public
signatures). RNG for ``simulate`` depends on the array namespace, so it is
delegated to a ``_simulate_uniform`` hook.
"""

from __future__ import annotations

from abc import ABC
from typing import Any, Optional, cast

from array_api_compat import array_namespace

from ._rootfind import solve_increasing
from .protocols import ArrayT, BicopLike, _BICOP_EXAMPLE

__all__ = ["BicopBase"]


class BicopBase(BicopLike[ArrayT], ABC):
  """Canonical partial implementation of :class:`~pyvinecopulib.core.BicopLike`.

  Subclasses implement ``pdf`` / ``hfunc1`` / ``hfunc2`` and inherit ``hinv1`` /
  ``hinv2`` (bisection of the h-functions), ``simulate``, ``loglik``, ``plot``
  and ``__repr__``. ``cdf`` raises unless overridden. To enable ``simulate``,
  override ``_simulate_uniform`` with the array namespace's RNG.

  See Also
  --------
  pyvinecopulib.core.BicopLike : The contract this implements.
  pyvinecopulib.torch.TorchBicop : A concrete (grid / TLL) subclass.
  """

  def loglik(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT:
    """Total log-likelihood ``sum(log c(u))`` of the pair at ``u``.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pseudo-observations in the unit square.
    x : array, shape (n, k), or None, optional
        Conditioning variables forwarded to :meth:`pdf`.

    Returns
    -------
    array, shape (), dtype float
        The summed log-density (a differentiable scalar under autograd, e.g.
        PyTorch). The per-observation density is floored at ``1e-20`` before the
        log.
    """
    dens: Any = self.pdf(u, x)
    xp = array_namespace(dens)
    return cast(ArrayT, xp.sum(xp.log(xp.clip(dens, 1e-20))))

  def hinv1(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT:
    """Numerically invert :meth:`hfunc1` in its second argument.

    Solves ``hfunc1([u1, .], x) = u2`` for the second argument by bisection
    (``solve_increasing``), since ``hfunc1``
    is increasing in it.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the conditioning value ``u1``; column 1 is the level to
        invert.
    x : array, shape (n, k), or None, optional
        Conditioning variables forwarded to :meth:`hfunc1`.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    ua: Any = u
    xp = array_namespace(ua)
    u1, p = ua[:, 0], ua[:, 1]
    return cast(
      ArrayT,
      solve_increasing(lambda v: self.hfunc1(xp.stack([u1, v], axis=-1), x), p),
    )

  def hinv2(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT:
    """Numerically invert :meth:`hfunc2` in its first argument.

    Solves ``hfunc2([., u2], x) = u1`` for the first argument by bisection
    (``solve_increasing``), since ``hfunc2``
    is increasing in it.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the level to invert; column 1 is the conditioning value
        ``u2``.
    x : array, shape (n, k), or None, optional
        Conditioning variables forwarded to :meth:`hfunc2`.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """
    ua: Any = u
    xp = array_namespace(ua)
    p, u2 = ua[:, 0], ua[:, 1]
    return cast(
      ArrayT,
      solve_increasing(lambda v: self.hfunc2(xp.stack([v, u2], axis=-1), x), p),
    )

  def cdf(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT:
    """Raise; the vine CDF is Monte-Carlo, so a per-pair ``cdf`` is optional.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations (unused in the raising default).
    x : array, shape (n, k), or None, optional
        Conditioning variables (unused in the raising default).

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values — only when a subclass overrides this method.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides a native ``cdf``.
    """
    raise NotImplementedError(
      f"{type(self).__name__}.cdf is not defined; the vine cdf uses "
      "Monte-Carlo simulation and does not require a per-pair cdf."
    )

  def flip(self) -> "BicopBase[ArrayT]":
    """Raise; override to return the pair with its arguments swapped.

    The flipped copula satisfies ``c'(u1, u2) = c(u2, u1)`` with the two
    h-functions (and inverses) exchanged. It is only required for hosting the
    pair in structure *selection*
    (:meth:`~pyvinecopulib.core.VinecopBase.select` reorients selection-time
    pairs onto their finalized slots with it); evaluation along a fixed
    structure never calls it.

    Returns
    -------
    BicopBase
        The argument-swapped pair copula — only when a subclass overrides this
        method.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides a native ``flip``.
    """
    raise NotImplementedError(
      f"{type(self).__name__}.flip is not defined; implement it (return the "
      "argument-swapped copula) to host this pair in structure selection."
    )

  def simulate(
    self,
    n: int,
    *,
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT:
    """Draw ``n`` samples via the pair's inverse Rosenblatt transform.

    Draws two independent uniforms ``(w1, w2)`` from ``_simulate_uniform``
    and returns ``(w1, hinv1([w1, w2], x))``, so the pair carries its fitted
    dependence.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    x : array, shape (n, k), or None, optional
        Conditioning variables (one row per sample) for a conditional draw.
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random ones.
    seeds : list of int or None, optional
        RNG seeds forwarded to ``_simulate_uniform``.

    Returns
    -------
    array, shape (n, 2), dtype float
        Samples in the unit square.
    """
    base_u: Any = self._simulate_uniform(n, qrng, list(seeds) if seeds else [])
    xp = array_namespace(base_u)
    u2: Any = self.hinv1(base_u, x)
    return cast(ArrayT, xp.stack([base_u[:, 0], u2], axis=-1))

  def _simulate_uniform(self, n: int, qrng: bool, seeds: list[int]) -> ArrayT:
    """Draw ``(n, 2)`` base uniforms for :meth:`simulate` (namespace-dependent RNG).

    Raises unless a subclass overrides it (numpy / torch differ on RNG);
    overriding it is all that is needed to enable :meth:`simulate`. Named after
    the exposed :func:`pyvinecopulib.utils.simulate_uniform` free function.
    """
    raise NotImplementedError(
      f"{type(self).__name__} does not implement _simulate_uniform; override it "
      "to enable simulate()."
    )

  def plot(
    self,
    plot_type: str = "surface",
    margin_type: str = "unif",
    xylim: Optional[tuple[float, float]] = None,
    grid_size: Optional[int] = None,
  ) -> None:
    """Plot the pair-copula density (contour or 3-D surface).

    Evaluates :meth:`pdf` on a grid and renders it with matplotlib, mirroring
    :meth:`pyvinecopulib.core.Bicop.plot`.

    Parameters
    ----------
    plot_type : str, default="surface"
        ``"surface"`` for a 3-D surface, ``"contour"`` for a contour plot.
    margin_type : str, default="unif"
        Margins the density is shown on: ``"unif"``, ``"norm"`` or ``"exp"``.
    xylim : tuple of float, or None, optional
        Axis limits; a sensible default per ``margin_type`` is used when ``None``.
    grid_size : int, or None, optional
        Number of grid points per axis; a default per ``plot_type`` is used when
        ``None``.

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
