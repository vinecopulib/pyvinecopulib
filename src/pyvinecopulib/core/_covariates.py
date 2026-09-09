"""Forwarding exogenous covariates ``x`` to a part being evaluated.

Every evaluation method in the library takes an optional ``x``, and every layer
that composes parts has to decide, per call, whether to hand it on. There are
exactly **two** rules, and the split is the point of this module: they look
alike and are not interchangeable.

- :func:`pair_eval` -- **forward whenever there is one.** A conforming pair
  copula declares ``x`` in its signature, which ``ty`` enforces on every
  ``BicopBase`` subclass, so the signature *is* the declaration and nothing
  else has to be maintained alongside it. Forwarding unconditionally is also
  what makes a pair that cannot take one fail loudly -- ``Bicop`` raises on the
  keyword rather than quietly modeling something else.
- :func:`declared_eval` -- **forward only to a part that declares
  ``supports_covariates``.** Margins and whole copulas are reached through
  structural protocols that foreign objects satisfy (a SciPy distribution, an
  adapter, ``Vinecop`` itself), so their signatures answer nothing: ``Vinecop``
  reports every method as ``(*args, **kwargs)`` under introspection. And one
  distribution may legitimately hold conditional and unconditional parts side
  by side -- a covariate-driven margin next to a plain one is a model the
  caller chose per column, not an accident -- so an undeclared part must be
  *skipped*, not refused.

What keeps the second rule from hiding a silent downgrade is that the refusal
happens one level up, at the object the covariates were handed to:
``VinedistBase._check_covariates`` raises when *nothing* reads them, which is
the case where silence would be indistinguishable from a conditional answer.
The fit-time member of the same family is ``_validation.reject_covariates``,
which refuses outright -- fitting cannot skip an ``x`` and stay correct, because
estimating ``f(y)`` when ``f(y | x)`` was asked for returns a different model.
"""

from __future__ import annotations

from typing import Any, Callable, Optional, Union, cast

from ._placement import place
from ._validation import validate_covariates
from .protocols import ArrayT

__all__ = ["declared_eval", "pair_eval", "prepare_covariates"]


def pair_eval(
  method: Callable[..., ArrayT], u: ArrayT, x: Optional[ArrayT]
) -> ArrayT:
  """Evaluate a pair-copula method, forwarding ``x`` whenever there is one.

  Parameters
  ----------
  method : callable
      The bound method to call, e.g. ``pair.hfunc1``.
  u : array
      The pair-copula argument, passed positionally.
  x : array, shape (n, p), or None, optional
      The edge's conditioning matrix, passed by keyword when not ``None``.

  Returns
  -------
  array
      Whatever ``method`` returns.
  """
  return method(u) if x is None else method(u, x=x)


def declared_eval(
  part: object,
  name: str,
  # The forwarded method's first positional argument: an array for `pdf` /
  # `cdf` / `icdf`, a draw count for `sample`.
  values: Union[ArrayT, int],
  x: Optional[ArrayT],
  # Forwarded to a method on a foreign object, so the keyword set is open --
  # and the answer may legitimately be another array type than `values`, which
  # is why the callers read `xp` off what came back rather than off the input.
  **kwargs: Any,  # noqa: ANN401
) -> Any:  # noqa: ANN401
  """Call ``part.name(values)``, forwarding ``x`` only if ``part`` reads one.

  Parameters
  ----------
  part : object
      The margin or copula to evaluate; its ``supports_covariates`` attribute
      decides, and an absent one means ``False``.
  name : str
      The method to call.
  values : object
      The first argument, passed positionally -- observations for a margin,
      copula-scale data or a sample size for a copula. Positional because a
      margin may name it whatever its own ecosystem does.
  x : array, shape (n, p), or None, optional
      The covariates, passed by keyword to a part that declares them.
  **kwargs : Any
      Further keyword arguments, forwarded either way.

  Returns
  -------
  object
      Whatever the method returns.
  """
  method = getattr(part, name)
  if x is None or not getattr(part, "supports_covariates", False):
    return method(values, **kwargs)
  return method(values, x=x, **kwargs)


def prepare_covariates(
  onto: object, x: Optional[ArrayT], n: int
) -> Optional[ArrayT]:
  """Validate covariates and place them where the numerics run.

  The two steps ``x`` needs and the third it must not get: it is checked for
  the row-aligned two-dimensional layout, and placed onto the namespace, dtype
  and device the values it will be combined with live on -- but never trimmed,
  being arbitrary reals rather than copula arguments.

  Placing is not cosmetic. A non-simplified vine concatenates ``x`` with the
  conditioning columns it gathered from the observations, so a NumPy ``x``
  handed to a PyTorch vine has to be brought across before they can meet.

  Parameters
  ----------
  onto : object
      The part whose placement to match, or the array to match directly --
      the observations, for a static fit engine that holds nothing itself.
  x : array, shape (n, p), or None, optional
      The covariates the caller supplied.
  n : int
      Number of observations they must align with.

  Returns
  -------
  array, shape (n, p), or None
      The placed covariates, or ``None`` when there were none.

  Raises
  ------
  ValueError
      If ``x`` is not two-dimensional or not row-aligned with ``n``.
  """
  if x is None:
    return None
  validate_covariates(x, n)
  return cast("ArrayT", place(onto, x))
