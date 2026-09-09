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

from array_api_compat import array_namespace

from ._placement import place
from ._validation import validate_covariates
from .protocols import ArrayT

__all__ = [
  "covariate_column",
  "covariate_row",
  "declared_eval",
  "pair_eval",
  "prepare_covariates",
]


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

  The layout it requires is not a subclass's to widen at its own entry
  points: the bases call this directly from ``logpdf`` / ``cdf_left`` /
  ``loglik`` / ``sample`` and from the vine and pair cascades, so a subclass
  that accepts a different ``x`` shape on the methods it wrote still meets
  this check on every method it inherited. Changing the layout means changing
  it where the composite runs, not at the entry points.

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
  # Through the `_prep` hook where there is one, so that an object whose
  # placement is *declared* rather than inferable is honored here as it is on
  # the argument path -- `_prep`'s own docstring promises it is "equally
  # correct for exogenous covariates", which routing around it made false.
  # Guarded because `onto` is not always an object: the static fit engines
  # pass an *array* as its own placement reference, and an array has no hook.
  hook = getattr(onto, "_prep", None)
  placed = hook(x) if callable(hook) else place(onto, x)
  return cast("ArrayT", placed)


def covariate_row(x: ArrayT, name: str = "x") -> ArrayT:
  """Check that ``x`` is one covariate row, and shape it ``(1, p)``.

  One of the two ways to resolve what :func:`prepare_covariates` refuses.
  That one rejects a one-dimensional ``x`` because ``(n,)`` is ambiguous --
  ``n`` observations of one covariate, or one observation of ``n`` of them --
  and it is row-aligned with the data, so guessing would silently align the
  wrong values. Naming a function is how the caller says which reading it
  means: this one reads ``(p,)`` as a single row, and
  :func:`covariate_column` reads ``(n,)`` as a single column.

  This direction needs no observation count to be safe. A single row is one
  row whatever ``n`` is -- a conditional object is a different surface at
  every covariate value, so drawing one shows the slice at one value, and the
  caller tiles it. ``covariate_column`` does need one, because there the
  length *is* the observation count and an unchecked length is exactly the
  silent misalignment this pair exists to prevent.

  Placement is not applied here. This shapes the row; the caller places it,
  and tiles it to ``(n, p)`` where it needs one row per observation.

  Parameters
  ----------
  x : array, shape (p,) or (1, p), dtype float
      One covariate row.
  name : str, default="x"
      Name to use in the error message.

  Returns
  -------
  array, shape (1, p), dtype float
      The same values, with a leading axis of length one.

  Raises
  ------
  ValueError
      If ``x`` has more than one row, or more than two axes.
  """
  a: Any = x
  xp = array_namespace(a)
  if getattr(a, "ndim", None) == 1:
    a = xp.reshape(a, (1, -1))
  if getattr(a, "ndim", None) != 2 or int(a.shape[0]) != 1:
    raise ValueError(
      f"{name} must be a single covariate row, shape (p,) or (1, p); got "
      f"{tuple(getattr(x, 'shape', ()))}. For one covariate per observation, "
      "which is the other reading of a one-dimensional x, use "
      "`covariate_column(x, n)`."
    )
  return cast("ArrayT", a)


def covariate_column(x: ArrayT, n: int, name: str = "x") -> ArrayT:
  """Read a one-dimensional ``x`` as one covariate per observation: ``(n, 1)``.

  The other resolution of the ambiguity :func:`prepare_covariates` refuses,
  and the sibling of :func:`covariate_row`. A conditional model over a single
  covariate produces exactly this shape -- a ``linspace`` over one variable,
  a single column read from a frame -- and every consumer of one otherwise
  writes the same reshape in front of ``prepare_covariates``.

  ``n`` is required, and is the whole of what makes this safe. The shape
  ``(n,)`` carries no statement about which axis is which, so the caller
  supplies that statement by choosing this function, and ``n`` is what lets
  the values be *checked* against it rather than trusted: an ``x`` of the
  wrong length is a misalignment, and one silently reshaped to ``(len(x), 1)``
  would pair each observation with some other observation's covariate. That
  is the failure the two-name split exists to prevent, so this does not
  default ``n`` and does not infer it.

  Already-shaped input passes through: ``(n, 1)`` states the same thing
  explicitly and is returned unchanged, after the same row check. A second
  column is refused rather than accepted, since the caller asked for one
  covariate.

  Placement is not applied here, as in :func:`covariate_row`. Compose it:
  ``prepare_covariates(onto, covariate_column(x, n), n)`` validates the
  layout again and places the result.

  Parameters
  ----------
  x : array, shape (n,) or (n, 1), dtype float
      One covariate, one value per observation.
  n : int
      Number of observations the values must align with.
  name : str, default="x"
      Name to use in the error message.

  Returns
  -------
  array, shape (n, 1), dtype float
      The same values, as a single column.

  Raises
  ------
  ValueError
      If ``x`` does not have exactly ``n`` values, if it has more than one
      column, or if it has more than two axes. The message names the other
      reading, since choosing the wrong sibling is the likely cause.
  """
  a: Any = x
  xp = array_namespace(a)
  shape = tuple(getattr(a, "shape", ()))
  if getattr(a, "ndim", None) == 1:
    a = xp.reshape(a, (-1, 1))
  if getattr(a, "ndim", None) != 2 or int(a.shape[1]) != 1:
    raise ValueError(
      f"{name} must be one covariate per observation, shape (n,) or (n, 1); "
      f"got {shape}. For a single covariate row shared across observations, "
      "which is the other reading of a one-dimensional x, use "
      "`covariate_row(x)`."
    )
  rows = int(a.shape[0])
  if rows != n:
    raise ValueError(
      f"{name} has {rows} values but {n} observations were declared; a "
      "one-dimensional x is read here as one value per observation, so the "
      "two have to agree. For a single covariate row of that length, use "
      "`covariate_row(x)`."
    )
  return cast("ArrayT", a)
