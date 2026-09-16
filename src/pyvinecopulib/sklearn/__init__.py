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
**Margins.** Which margin class each variable gets is the distribution's:
``distribution=`` names a :class:`pyvinecopulib.core.VinedistBase` subclass,
and its ``margin_class`` is what an estimator fits per column. *How* each
one is fitted is ``margin_controls=``, in any form
:func:`pyvinecopulib.margins.resolve_margin_controls` accepts: one controls
object broadcast to every column, a per-column sequence, or a mapping keyed
by feature name. Fitting runs the whole two-step estimator on the
distribution, published as ``distribution_``, and ``margin_summary_``
describes the margin each variable ended up with.

**Which lane.** By default the estimators fit a
:class:`pyvinecopulib.core.Vinedist` --- ``Vinecop`` paired with ``Kde1d``
--- so the sklearn module **does not require PyTorch**. Pass
``distribution=TorchVinedist`` to route the same pipeline through the
PyTorch evaluator (GPU placement, autograd); importing that class to name
it is the explicit opt-in.

**DataFrame input.** Every estimator accepts both NumPy arrays and
pandas DataFrames. DataFrames may mix numeric, ordered-categorical,
and unordered-categorical columns; the latter are expanded to ordered
``{0, 1}`` dummies before fitting, and the same expansion is
re-applied at predict time.

**Low-level knobs.** Pair family, threading and the structure-selection
algorithm are :class:`pyvinecopulib.core.FitControlsVinecop`'s, passed as
``controls=``; a pre-specified :class:`pyvinecopulib.core.RVineStructure`
is passed as ``structure=``. Reach for those directly whenever you need
control beyond the sklearn convenience layer. See each class docstring for
the full methodology and references.
"""

from ..core._validation import extra_required

with extra_required(  # noqa: RUF067
  extra="sklearn",
  requirement="pyvinecopulib.sklearn requires scikit-learn.",
):
  import sklearn  # noqa: F401

from .density import VineDensity
from .regressor import VineRegressor

__all__ = [
  "VineDensity",
  "VineRegressor",
]
