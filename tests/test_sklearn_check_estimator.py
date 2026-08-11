"""sklearn ``check_estimator`` compliance tests.

Runs :func:`sklearn.utils.estimator_checks.parametrize_with_checks` on
both sklearn estimators and documents which standard checks are
opted out of. The skip list captures two categories:

- **Genuine opt-outs**: checks that don't apply to a joint-density /
  vine-regressor archetype (sparse inputs, 1-D / 1-feature degenerate
  cases, complex-valued data, etc.).
- **Known WIP**: checks that surface real sklearn-compliance gaps the
  initial backend-refactor PR did not fix. Listed with ``TODO``
  comments so future PRs can address them and shrink this list.

The goal is to (a) prove that the basic sklearn-developer-guide
contract holds on the easy checks (constructor / clone / fit returns
self / etc.), and (b) make every unintended deviation explicit and
trackable.
"""

from __future__ import annotations

import inspect

import pytest

pytest.importorskip("sklearn")

from sklearn.utils.estimator_checks import parametrize_with_checks  # noqa: E402

# ``xfail_strict`` was added in sklearn 1.7; on older versions we fall
# back to the default (strict) and rely on the ``expected_failed_checks``
# list staying accurate for the available checks.
_PWC_KWARGS: dict = {}
if "xfail_strict" in inspect.signature(parametrize_with_checks).parameters:
  _PWC_KWARGS["xfail_strict"] = False

from pyvinecopulib.sklearn import VineDensity, VineRegressor  # noqa: E402


# ---------------------------------------------------------------------------
# Per-estimator opt-outs
# ---------------------------------------------------------------------------


_GENUINE_OPT_OUTS = {
  "check_estimator_sparse_array": "vines require dense U-scale data",
  "check_estimator_sparse_matrix": "vines require dense U-scale data",
  "check_estimator_sparse_tag": "vines require dense U-scale data",
  "check_fit2d_1feature": "1-D vine is degenerate (no pair copulas)",
  "check_fit2d_1sample": "Kde1d requires >= 2 samples to estimate a bandwidth",
  "check_fit1d": "vines require 2-D U-scale data",
  "check_fit2d_predict1d": "predict requires 2-D inputs",
  "check_estimators_nan_inf": "NaN inputs not supported (documented)",
  "check_complex_data": "complex inputs not supported",
  "check_estimators_dtypes": "internal float64 promotion is intentional",
  "check_dtype_object": "internal float64 promotion is intentional",
  "check_estimators_empty_data_messages": (
    "empty data error messages — TODO: surface a sklearn-shaped error early"
  ),
  "check_n_features_in_after_fitting": (
    "TODO: n_features_in_ is set, but the sklearn check inspects "
    "the post-fit error path on wrong-shape input"
  ),
  "check_mixin_order": (
    "TODO: swap to (DensityMixin, VineBase) / (RegressorMixin, VineBase) — "
    "sklearn 1.8 enforces Mixin-before-BaseEstimator order"
  ),
  "check_pipeline_consistency": (
    "density score is mean log-likelihood (negative); standard "
    "Pipeline.score consistency check assumes positive scores"
  ),
}


# Checks that actually pass on each estimator despite living in the
# baseline skip list. Listed here so they get popped before
# ``parametrize_with_checks`` sees the dict — otherwise pytest emits
# XPASS noise. Shrink the per-estimator list as the underlying check
# semantics shift across sklearn versions; the matching baseline
# entries can move into ``_GENUINE_OPT_OUTS`` once the pass is
# universal.
_KNOWN_PASSING_PER_ESTIMATOR: dict[str, tuple[str, ...]] = {
  "VineDensity": (
    "check_estimators_dtypes",
    "check_fit2d_1feature",
    "check_fit2d_predict1d",
    "check_pipeline_consistency",
  ),
  "VineRegressor": (
    "check_estimators_dtypes",
    "check_fit2d_1feature",
  ),
}


def _prune(skips: dict[str, str], cls_name: str) -> dict[str, str]:
  for k in _KNOWN_PASSING_PER_ESTIMATOR.get(cls_name, ()):
    skips.pop(k, None)
  return skips


def _density_opt_outs() -> dict[str, str]:
  return _prune(dict(_GENUINE_OPT_OUTS), "VineDensity")


def _regressor_opt_outs() -> dict[str, str]:
  d = dict(_GENUINE_OPT_OUTS)
  d.update(
    {
      # Quantile-stacked output is a per-call shape, not a fitted
      # property; sklearn's multi-output checks require
      # __sklearn_tags__.target_tags.multi_output = True which is
      # context-sensitive here.
      "check_regressors_train": "quantile-stacked output is not strict multi-output",
      "check_regressor_data_not_an_array": "quantile-stacked output shape",
      "check_fit_score_takes_y": (
        "TODO: regressor with quantiles returns multi-output predictions; "
        "sklearn's R^2 score check expects single output"
      ),
    }
  )
  return _prune(d, "VineRegressor")


# ---------------------------------------------------------------------------
# Parametrized run
# ---------------------------------------------------------------------------


_ESTIMATORS = [
  (VineDensity(), _density_opt_outs()),
  (VineRegressor(quantiles=[0.5]), _regressor_opt_outs()),
]


def _expected_failed_checks(estimator):
  for est, skips in _ESTIMATORS:
    if type(est) is type(estimator):
      return skips
  return {}


@parametrize_with_checks(
  [est for est, _ in _ESTIMATORS],
  expected_failed_checks=_expected_failed_checks,
  **_PWC_KWARGS,
)
def test_sklearn_compliance(estimator, check):
  check(estimator)
