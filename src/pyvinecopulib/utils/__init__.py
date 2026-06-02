"""Utility primitives used alongside the core copula classes.

This subpackage groups the supporting functionality that copula
modelling typically needs but that doesn't fit inside the
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
  Kendall's :math:`\\tau`, Spearman's :math:`\\rho`, Pearson, etc.
- **Low-discrepancy sequences** — :func:`sobol`, :func:`ghalton`, and
  :func:`simulate_uniform` (the high-level driver) produce
  quasi-random uniform points used by Monte-Carlo evaluation of
  copula CDFs and by random vine-structure generation.
- **Plotting helper** — :func:`pairs_copula_data` produces a
  pair-plot of a copula sample or a fitted :class:`Vinecop`.
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
vine-copula factorisation.
"""

from ..pyvinecopulib_ext import (
  Kde1d,
  benchmark,
  ghalton,
  simulate_uniform,
  sobol,
  to_pseudo_obs,
  wdm,
)
from ._pair_plots import pairs_copula_data

__all__ = [
  "Kde1d",
  "benchmark",
  "ghalton",
  "pairs_copula_data",
  "simulate_uniform",
  "sobol",
  "to_pseudo_obs",
  "wdm",
]
