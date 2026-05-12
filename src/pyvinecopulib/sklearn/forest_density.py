from __future__ import annotations

import numpy as np
from joblib import Parallel, delayed
from sklearn.base import DensityMixin

from ._base import (
  _DOC_DISCRETE,
  _DOC_FACTORIZATION,
  _DOC_FOREST,
  _DOC_PIPELINE,
  _DOC_REFERENCES,
)
from ._forest_base import VineForestBase
from .density import VineDensity


class VineForestDensity(VineForestBase, DensityMixin):
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
    """Ensemble of vine-copula density estimators with random structures.

    Builds ``n_vines`` :class:`VineDensity` base learners on randomly
    sampled vine structures, prunes them via the model confidence set
    (MCS), and averages the per-sample density values at predict
    time. See the class docstring for the full methodology.

    Parameters
    ----------
    base_params : dict, optional
        Keyword arguments forwarded to each :class:`VineDensity`
        ``__init__``. Example: ``{"batch_size": 200}``.
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
        Number of joblib workers used during fit and pdf evaluation.
    verbose : bool, default=False
        Warn if no random estimator beats the default.
    """
    super().__init__(
      base_class=VineDensity,
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

  def _create_base_estimator(self) -> VineDensity:
    return VineDensity(**self.base_params)

  def _loglik_estimator(self, estimator, X, y=None):
    return estimator.score_samples(X)

  def fit(self, X):
    """Fit the ensemble of vine density estimators.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Training data.

    Returns
    -------
    self : VineForestDensity
    """
    return self._fit_ensemble(X)

  def score_samples(self, X):
    """Log of the ensemble-averaged density at each sample.

    Returns :math:`\\log \\hat f_{\\text{ens}}(x)` where the ensemble
    density is the simple average over surviving base estimators,
    :math:`\\hat f_{\\text{ens}}(x) = M^{-1} \\sum_{m=1}^M
    \\hat f_m(x)`.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test samples.

    Returns
    -------
    ndarray of shape (n_samples,)
        Log-density values.
    """
    return np.log(self.pdf(X))

  def score(self, X, y=None):
    """Mean log-density across samples (sklearn ``score`` convention).

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test samples.
    y : ignored
        Present for API compatibility with ``DensityMixin.score``.

    Returns
    -------
    float
        Mean of :meth:`score_samples` over the rows of ``X``.
    """
    return float(self.score_samples(X).mean())

  def pdf(self, X):
    """Ensemble-averaged joint density.

    Evaluates :class:`VineDensity.pdf` for each surviving base
    estimator in parallel and averages the resulting per-sample
    density values.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test samples.

    Returns
    -------
    ndarray of shape (n_samples,)
        Density values.
    """
    X = self._prepare_prediction_data(X)

    outer_n_jobs = self._adjust_estimators_num_threads(self._estimators)
    pdf_list = Parallel(n_jobs=outer_n_jobs)(
      delayed(estimator.pdf)(X) for estimator in self._estimators
    )

    return np.mean(pdf_list, axis=0)

  def sample(self, n_samples=1):
    """Draw samples from the ensemble density.

    Samples are generated by randomly selecting one surviving
    estimator per draw (uniformly over the ensemble) and drawing one
    sample from it. This is equivalent to sampling from the
    mixture-of-densities distribution :math:`M^{-1} \\sum_m
    \\hat f_m(x)`.

    Parameters
    ----------
    n_samples : int, default=1
        Number of samples to generate.

    Returns
    -------
    ndarray of shape (n_samples, n_features)
        Generated samples in the original feature scale.
    """
    self._check_fitted()
    rng = np.random.default_rng(self.seed)
    selected_estimators = rng.choice(
      self._estimators, size=n_samples, replace=True
    )
    samples = [estimator.sample(1)[0] for estimator in selected_estimators]
    return np.array(samples)

  def sample_from_estimator(self, estimator_idx, n_samples=1):
    """Draw samples from a specific surviving estimator.

    Parameters
    ----------
    estimator_idx : int
        Index into ``self._estimators``.
    n_samples : int, default=1
        Number of samples to generate.

    Returns
    -------
    ndarray of shape (n_samples, n_features)
        Generated samples in the original feature scale.
    """
    self._check_fitted()
    if estimator_idx >= len(self._estimators):
      raise ValueError(
        f"estimator_idx {estimator_idx} is out of range. "
        f"Ensemble has {len(self._estimators)} estimators."
      )
    return self._estimators[estimator_idx].sample(n_samples)


VineForestDensity.__doc__ = f"""Forest of vine-copula density estimators.

An ensemble of :class:`VineDensity` base learners fitted on randomly
sampled vine structures. Survivors are selected via a model
confidence set (MCS) on a held-out validation split, and predictions
are averaged across survivors. Targets the random-search-plus-MCS
methodology of the underlying paper.

{_DOC_FOREST}
{_DOC_PIPELINE}
{_DOC_FACTORIZATION}
{_DOC_DISCRETE}

Examples
--------
>>> import numpy as np
>>> from pyvinecopulib.sklearn import VineForestDensity
>>> rng = np.random.default_rng(0)
>>> X = rng.standard_normal((300, 3))
>>> forest = VineForestDensity(n_vines=10, n_jobs=1).fit(X[:200])
>>> forest.score_samples(X[200:205])      # log-density per row
>>> forest.pdf(X[200:205])                # ensemble-averaged density
>>> forest.sample(n_samples=5)            # mixture sampling

{_DOC_REFERENCES}
"""
