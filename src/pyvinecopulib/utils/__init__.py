"""Utility primitives used alongside the core copula classes.

This subpackage groups the supporting functionality that copula
modeling typically needs but that doesn't fit inside the
:class:`pyvinecopulib.core.Bicop` / :class:`pyvinecopulib.core.Vinecop`
classes themselves:

- **Univariate kernel density estimation** — :class:`Kde1d` is a
  boundary-corrected 1d KDE with built-in support for continuous and
  discrete margins. It is the marginal estimator used internally by
  :class:`pyvinecopulib.sklearn.VineDensity` /
  :class:`pyvinecopulib.sklearn.VineRegressor` and is reusable
  standalone for any 1d density problem.
- **Pseudo-observations** — :func:`to_pseudo_obs` rank-transforms a
  data matrix into the unit hypercube, the canonical input shape for
  fitting copulas.
- **Dependence measures** — :func:`wdm` computes weighted versions of
  Kendall's :math:`\\tau`, Spearman's :math:`\\rho`, Pearson,
  Chatterjee's :math:`\\xi`, etc.
- **Latent samples for discrete data** — :func:`find_latent_sample`
  recovers a continuous sample from interval-censored copula data,
  the transform a nonparametric fit on discrete margins runs on.
- **Low-discrepancy sequences** — :func:`sobol`, :func:`ghalton`, and
  :func:`sample_uniform` (the high-level driver) produce
  quasi-random uniform points used by Monte-Carlo evaluation of
  copula CDFs and by random vine-structure generation.
- **Plotting helper** — :func:`pairs_copula_data` produces a
  pair-plot of a copula sample or a fitted ``Vinecop``.
- **Benchmarking** — :func:`benchmark` runs a quick comparison of
  available pair-copula families on synthetic data; convenient for
  smoke-testing a fresh install.

Most users reach for these via :mod:`pyvinecopulib.sklearn` rather
than calling them directly. The ``examples/07_kde1d.ipynb`` and
``examples/06_weighted_dependence_measures.ipynb`` notebooks demo
:class:`Kde1d` and :func:`wdm` in isolation.

Notes
-----
The :doc:`concepts page </concepts>` introduces
:ref:`pseudo-observations <concepts-sklar>` (the canonical input
shape for copula fits) and dependence measures alongside the
vine-copula factorization.
"""

import warnings
from typing import Any

from ..pyvinecopulib_ext import (
  Kde1d,
  benchmark,
  find_latent_sample,
  ghalton,
  sample_uniform,
  sobol,
  to_pseudo_obs,
  wdm,
)
from .._deprecations import _method_alias
from ._pair_plots import pairs_copula_data

__all__ = [
  "Kde1d",
  "benchmark",
  "find_latent_sample",
  "ghalton",
  "pairs_copula_data",
  "sample_uniform",
  "sobol",
  "to_pseudo_obs",
  "wdm",
]


# `Kde1d.simulate` shipped in 0.7.6; the canonical name is now `sample`.
Kde1d.simulate = _method_alias(Kde1d.sample, "simulate", "utils.Kde1d")

# `simulate_uniform` shipped in 0.7.6 under that name. Served from `__getattr__`
# rather than assigned, so it stays out of `__all__`, the generated stubs and the
# docs while still resolving for existing callers.
_DEPRECATED_FUNCTIONS = {"simulate_uniform": "sample_uniform"}


def __getattr__(name: str) -> Any:
  """Resolve a deprecated function name.

  Parameters
  ----------
  name : str
      Attribute being looked up.

  Returns
  -------
  object
      The current binding.

  Raises
  ------
  AttributeError
      If the name is neither current nor deprecated.
  """
  if name in _DEPRECATED_FUNCTIONS:
    new = _DEPRECATED_FUNCTIONS[name]
    warnings.warn(
      f"`pyvinecopulib.utils.{name}` is deprecated; use "
      f"`pyvinecopulib.utils.{new}` instead.",
      DeprecationWarning,
      stacklevel=2,
    )
    return globals()[new]
  raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
  """List the module's public and deprecated names.

  Returns
  -------
  list of str
      Sorted names.
  """
  return sorted(set(__all__) | set(_DEPRECATED_FUNCTIONS))
