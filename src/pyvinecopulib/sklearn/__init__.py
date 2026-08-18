"""Scikit-learn-compatible vine-copula estimators.

This subpackage wraps the core pyvinecopulib machinery behind the
standard sklearn ``BaseEstimator`` / ``fit`` / ``predict`` interface.
Two estimators ship:

- :class:`VineDensity` — joint-density estimator. Fits univariate
  margins, then a vine copula on the resulting pseudo-observations.
  Exposes ``score_samples`` / ``pdf`` / ``cdf`` / ``sample``.
- :class:`VineRegressor` — non-parametric conditional mean / quantile
  regressor built from a vine copula over ``(Y, X)``. Predictions are
  weighted statistics of the training responses.

If you have not used vine copulas before, the
:doc:`concepts page </concepts>` introduces pair copulas, R-vines,
and the default *Transformed Local Likelihood* (TLL) pair-copula
family in ~5 minutes.

Requires scikit-learn and pandas. Install with
``pip install pyvinecopulib[sklearn]``.

Notes
-----
**Margins.** The marginal half of the model is configured with
``margins=``, in any form :func:`pyvinecopulib.margins.resolve_margins`
accepts: an alias (``"kde"``, the default, ``"empirical"``,
``"parametric"``), one margin broadcast to every column, a per-column
sequence, a mapping keyed by feature name, or a callable. Fitting
assembles both halves into a :class:`pyvinecopulib.core.Vinedist`,
published as ``distribution_``; ``selection_report_`` carries the
per-candidate table of any margin that chose its own family.

**Backends.** By default the estimators run on ``Vinecop`` (the
default backend), so the sklearn module
**does not require PyTorch**. Pass a configured
:class:`~pyvinecopulib.sklearn.backends.TorchVinecopBackend` via the
``backend=`` kwarg to route through the torch backend instead — see
:mod:`pyvinecopulib.sklearn.backends` for a comparison and examples.

**DataFrame input.** Every estimator accepts both NumPy arrays and
pandas DataFrames. DataFrames may mix numeric, ordered-categorical,
and unordered-categorical columns; the latter are expanded to ordered
``{0, 1}`` dummies before fitting, and the same expansion is
re-applied at predict time.

**Low-level knobs.** Pair family, threading, structure-selection
algorithm, etc. are passed through to
:class:`pyvinecopulib.core.FitControlsVinecop` and
:class:`pyvinecopulib.core.RVineStructure` (carried inside the
backend object). Reach for those directly whenever you need control
beyond the sklearn convenience layer. See each class docstring for
the full methodology and references.
"""

try:
  import sklearn  # noqa: F401
except ImportError as e:
  raise ImportError(
    "pyvinecopulib.sklearn requires scikit-learn. "
    "Install it with `pip install pyvinecopulib[sklearn]`."
  ) from e

from . import backends
from .backends import (
  TorchVinecopBackend,
  VinecopBackend,
  resolve_backend,
)
from .density import VineDensity
from .regressor import VineRegressor

__all__ = [
  "TorchVinecopBackend",
  "VineDensity",
  "VineRegressor",
  "VinecopBackend",
  "backends",
  "resolve_backend",
]
