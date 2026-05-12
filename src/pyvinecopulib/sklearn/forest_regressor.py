from __future__ import annotations

import numpy as np
from joblib import Parallel, delayed
from sklearn.base import RegressorMixin

from ._base import (
  _DOC_DISCRETE,
  _DOC_FACTORIZATION,
  _DOC_FOREST,
  _DOC_PIPELINE,
  _DOC_REFERENCES,
)
from ._forest_base import VineForestBase
from .regressor import VineRegressor


class VineForestRegressor(VineForestBase, RegressorMixin):
  def __init__(
    self,
    base_params: dict | None = None,
    n_vines: int = 100,
    vines_sampling: str = "uniform",
    bootstrap: bool = True,
    val_fraction: float = 0.25,
    best_only: bool = False,
    method: str | None = "da_mcs_marg",
    alpha: float = 0.05,
    add_dissmann: bool = True,
    seed: int = 42,
    n_jobs: int = 1,
    verbose: bool = False,
  ) -> None:
    """Ensemble of vine-copula regressors with random structures.

    Builds ``n_vines`` :class:`VineRegressor` base learners on
    randomly sampled vine structures, prunes them via the model
    confidence set (MCS), and averages the conditional weights across
    survivors before computing predictions. See the class docstring
    for the full methodology.

    Parameters
    ----------
    base_params : dict, optional
        Keyword arguments forwarded to each :class:`VineRegressor`
        ``__init__``. Example:
        ``{"quantiles": [0.1, 0.5, 0.9], "batch_size": 200}``.
    n_vines : int, default=100
        Number of random base estimators (before MCS pruning).
    vines_sampling : {"uniform", "local"}, default="uniform"
        Random-structure generator. See :class:`VineForestBase`.
    bootstrap : bool, default=True
        Bootstrap-resample the training set for each base estimator.
    val_fraction : float, default=0.25
        Held-out fraction used for MCS survivor selection. ``0``
        disables validation.
    best_only : bool, default=False
        Keep only the single best survivor rather than the full MCS.
    method : {"da_mcs_marg", "da_mcs_unif", None}, default="da_mcs_marg"
        Survivor-selection method. ``None`` keeps anything strictly
        better than the Dissmann baseline.
    alpha : float, default=0.05
        Significance level for the MCS selector.
    add_dissmann : bool, default=True
        Include the Dissmann-structure baseline among the candidates.
    seed : int, default=42
        Random seed for reproducibility.
    n_jobs : int, default=1
        Number of joblib workers used during fit and predict.
    verbose : bool, default=False
        Warn if no random estimator beats the default.
    """
    super().__init__(
      base_class=VineRegressor,
      base_params=base_params,
      n_vines=n_vines,
      vines_sampling=vines_sampling,
      bootstrap=bootstrap,
      val_fraction=val_fraction,
      best_only=best_only,
      method=method,
      alpha=alpha,
      add_dissmann=add_dissmann,
      seed=seed,
      n_jobs=n_jobs,
      verbose=verbose,
    )

  def _create_base_estimator(self) -> VineRegressor:
    # The base learner emits raw (un-normalised) weights so the
    # forest can average across trees and normalise once at the
    # ensemble level. `_normalize_weights` is no longer an __init__
    # parameter on VineRegressor; set it via attribute access.
    est = VineRegressor(**self.base_params)
    est._normalize_weights = False
    return est

  def _loglik_estimator(self, estimator, X, y=None):
    """Conditional log-density log f(y | x) per row.

    Uses the same three-term identity as the single-vine regressor's
    private hooks:

    .. math::

       \\log f_{Y \\mid X}(y \\mid x)
         = \\log c_{Y, X}(u_y, u_x) - \\log c_X(u_x) + \\log f_Y(y).

    Y-only inputs (``y is None``) aren't supported on the regressor
    forest — by construction ``fit(X, y)`` always supplies a target.
    """
    if y is None:
      raise ValueError(
        "VineForestRegressor._loglik_estimator requires y; the regressor "
        "forest is fit on (X, y) pairs."
      )

    log_c_joint = estimator._pdf_samples(X, y, log=True, copula_only=True)
    log_c_x = estimator._copula_marginal_density(X, log=True)
    eps = 1e-10
    log_f_y = np.log(
      np.clip(estimator._y_kde1d.pdf(np.asarray(y).ravel()), eps, None)
    )
    return log_c_joint - log_c_x + log_f_y

  def fit(self, X, y):
    """Fit the ensemble of vine regressors.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Training covariates.
    y : ndarray of shape (n_samples,)
        Training responses.

    Returns
    -------
    self : VineForestRegressor
    """
    return self._fit_ensemble(X, y)

  def _iter_weights(self, X):
    """Yield batched ensemble-averaged conditional weights.

    Iterates over batches of ``X`` and, for each batch, averages the
    raw weights produced by every surviving :class:`VineRegressor`
    (recall that the forest sets ``_normalize_weights = False`` on
    each base learner so the unnormalised copula densities can be
    averaged sensibly). The averaged weights are then row-normalised
    once at the ensemble level.

    Parameters
    ----------
    X : ndarray of shape (n_test, n_features)
        Test covariates on the original (un-transformed) scale.

    Yields
    ------
    weights : ndarray of shape (batch, n_train_or_grid)
        Row-normalised ensemble weights for the current batch.
    start : int
        Index of the first row in this batch.
    end : int
        Index one past the last row in this batch.
    """
    n_test = X.shape[0]
    batch_size = self._estimators[0].batch_size
    outer_n_jobs = self._adjust_estimators_num_threads(self._estimators)

    with Parallel(n_jobs=outer_n_jobs) as parallel:
      for start in range(0, n_test, batch_size):
        end = min(start + batch_size, n_test)

        def get_one_weight(estimator, start=start, end=end):
          weights, _, _ = next(estimator._iter_weights(X[start:end]))
          return weights

        weights = parallel(
          delayed(get_one_weight)(est) for est in self._estimators
        )

        avg = np.mean(weights, axis=0)
        avg /= avg.sum(axis=1, keepdims=True)
        yield avg, start, end

  def predict(self, X):
    """Predict conditional mean and/or quantiles from the ensemble.

    For each test row :math:`x`, computes the ensemble-averaged
    weights via :meth:`_iter_weights` and applies the same weighted
    statistics as :meth:`VineRegressor.predict`: weighted average of
    training responses for the mean, weighted quantile (inverted CDF)
    for each requested level.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test covariates.

    Returns
    -------
    ndarray
        Predictions of shape ``(n_samples,)`` if a single output is
        requested (mean *or* a single quantile), or
        ``(n_samples, n_outputs)`` otherwise. Output columns are
        ordered: mean (if enabled), then quantiles in the order
        configured on the base estimator.
    """
    X = self._prepare_prediction_data(X)
    return self._estimators[0]._predict_from_iter(X, self._iter_weights)


VineForestRegressor.__doc__ = f"""Forest of vine-copula regressors.

An ensemble of :class:`VineRegressor` base learners fitted on
randomly sampled vine structures. Survivors are selected via a model
confidence set (MCS) on a held-out validation split, and the
conditional weights :math:`w_i(x)` from the single-vine estimating
equation framework are averaged across survivors before being used
to compute conditional means or quantiles.

{_DOC_FOREST}
{_DOC_PIPELINE}
{_DOC_FACTORIZATION}
{_DOC_DISCRETE}

Examples
--------
>>> import numpy as np
>>> from pyvinecopulib.sklearn import VineForestRegressor
>>> rng = np.random.default_rng(0)
>>> X = rng.standard_normal((300, 3))
>>> y = X @ [1.5, -0.8, 0.4] + 0.2 * rng.standard_normal(300)
>>> forest = VineForestRegressor(
...     base_params={{"quantiles": [0.1, 0.5, 0.9]}},
...     n_vines=10, n_jobs=1,
... ).fit(X[:200], y[:200])
>>> forest.predict(X[200:205])           # columns: mean, q10, q50, q90

{_DOC_REFERENCES}
"""
