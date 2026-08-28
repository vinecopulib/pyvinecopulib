"""The independence copula, on any array namespace.

What :meth:`VinecopBase.select` puts on an edge whose dependence falls below
``threshold``. Upstream leaves such an edge holding a default-constructed
``vinecopulib::Bicop``, which is the independence copula, and never fits it
(``tools_select.ipp`` ``fit_or_reuse_pair_copula``); this is the port's
equivalent.

Written against the Array API, so one class serves a NumPy and a PyTorch vine
alike. A subclass that stores its pairs in a container of its own -- as
:class:`~pyvinecopulib.torch.TorchVinecop` stores ``nn.Module`` grids -- is free
to substitute its own representation of independence when it takes delivery of
``select``'s pairs.
"""

from __future__ import annotations

from typing import Any, Optional, cast

from array_api_compat import array_namespace

from .bicop_base import BicopBase
from .protocols import ArrayT


class IndependencePair(BicopBase[ArrayT]):
  """The independence copula ``C(u1, u2) = u1 * u2``.

  The evaluation surface is the closed form rather than an inherited
  numerical one, so nothing here bisects and nothing depends on a grid.
  ``sample`` and ``loglik`` come from :class:`BicopBase` unchanged, the
  closed form buying them nothing.

  See Also
  --------
  pyvinecopulib.core.VinecopBase.select : Places this on a thresholded edge.
  """

  def pdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Density, which is one everywhere.

    Parameters
    ----------
    u : array, shape (n, 2) or wider, dtype float
        Copula-scale observations. Only the shape is read.
    x : array, optional
        Ignored; the pair is unconditional.

    Returns
    -------
    array, shape (n,), dtype float
        Ones.
    """
    del x
    cols = cast(Any, u)
    xp = array_namespace(cols)
    return cast(ArrayT, xp.ones_like(cols[:, 0]))

  def cdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Distribution function ``u1 * u2``.

    Parameters
    ----------
    u : array, shape (n, 2) or wider, dtype float
        Copula-scale observations.
    x : array, optional
        Ignored; the pair is unconditional.

    Returns
    -------
    array, shape (n,), dtype float
        The product of the two arguments.
    """
    del x
    cols = cast(Any, u)
    return cast(ArrayT, cols[:, 0] * cols[:, 1])

  def hfunc1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """``P(U2 <= u2 | U1) = u2``.

    Parameters
    ----------
    u : array, shape (n, 2) or wider, dtype float
        Copula-scale observations.
    x : array, optional
        Ignored; the pair is unconditional.

    Returns
    -------
    array, shape (n,), dtype float
        The second argument.
    """
    del x
    return cast(ArrayT, cast(Any, u)[:, 1])

  def hfunc2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """``P(U1 <= u1 | U2) = u1``.

    Parameters
    ----------
    u : array, shape (n, 2) or wider, dtype float
        Copula-scale observations.
    x : array, optional
        Ignored; the pair is unconditional.

    Returns
    -------
    array, shape (n,), dtype float
        The first argument.
    """
    del x
    return cast(ArrayT, cast(Any, u)[:, 0])

  def hinv1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of :meth:`hfunc1` in its second argument, which is identity.

    Parameters
    ----------
    u : array, shape (n, 2) or wider, dtype float
        Conditioning value in column 0, level in column 1.
    x : array, optional
        Ignored; the pair is unconditional.

    Returns
    -------
    array, shape (n,), dtype float
        The second argument.
    """
    del x
    return cast(ArrayT, cast(Any, u)[:, 1])

  def hinv2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of :meth:`hfunc2` in its first argument, which is identity.

    Parameters
    ----------
    u : array, shape (n, 2) or wider, dtype float
        Level in column 0, conditioning value in column 1.
    x : array, optional
        Ignored; the pair is unconditional.

    Returns
    -------
    array, shape (n,), dtype float
        The first argument.
    """
    del x
    return cast(ArrayT, cast(Any, u)[:, 0])

  def flip(self) -> "IndependencePair[ArrayT]":
    """The argument-swapped copula, which is this one.

    Returns
    -------
    IndependencePair
        ``self``; the independence copula is symmetric.
    """
    return self

  def __repr__(self) -> str:
    """Short representation.

    Returns
    -------
    str
        The class name.
    """
    return "IndependencePair()"
