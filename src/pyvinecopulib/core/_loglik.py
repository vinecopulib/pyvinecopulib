"""Turning densities into log-densities, and log-densities into a total.

Two operations the four levels share, kept together because they are the two
places a likelihood can lose information it need not lose. :func:`safe_log`
maps a density of exactly zero to ``-inf`` rather than to a warning or a
``nan``; :func:`sum_loglik` leaves an observation with no log-density out of
the total rather than propagating its ``nan`` over the whole sample.

Internal: reached through the bases' own ``logpdf`` / ``loglik``.
"""

from __future__ import annotations

from typing import Any, cast

from .protocols import array_namespace

from .protocols import ArrayT

__all__ = ["safe_log", "sum_loglik"]


def safe_log(dens: ArrayT) -> ArrayT:
  """Log of a density, with a zero mapped to ``-inf`` rather than a warning.

  A density is legitimately zero off its support, and ``log(0)`` there is the
  right answer -- but computing it directly warns, and on some namespaces
  returns a ``nan``. The mask is applied before the log rather than after, so
  nothing invalid is evaluated.

  Parameters
  ----------
  dens : array
      Density or mass values, nonnegative.

  Returns
  -------
  array
      ``log(dens)``, and ``-inf`` wherever ``dens`` is not positive.
  """
  d: Any = dens
  xp = array_namespace(d)
  positive = d > 0
  safe = xp.where(positive, d, xp.ones_like(d))
  return cast(
    "ArrayT", xp.where(positive, xp.log(safe), xp.full_like(d, float("-inf")))
  )


def sum_loglik(logdens: ArrayT) -> ArrayT:
  """Total of a vector of log-densities, over the observations that have one.

  An observation carrying a ``nan`` has a ``nan`` log-density and is left out,
  so one missing value costs its own row rather than the whole sample. A
  log-density of ``-inf`` is not missing -- it is the model ruling an
  observation out -- and is kept.

  Written as a masked sum rather than a boolean gather, so the shape is static
  and the total stays differentiable.

  Parameters
  ----------
  logdens : array, shape (n,), dtype float
      Per-observation log-densities.

  Returns
  -------
  array
      A 0-d array holding the total.
  """
  d: Any = logdens
  xp = array_namespace(d)
  return cast("ArrayT", xp.sum(xp.where(xp.isnan(d), xp.zeros_like(d), d)))
