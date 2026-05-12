from typing import Any

import numpy as np
import pandas as pd
import pyvinecopulib as pv
from sklearn.base import DensityMixin

from ._base import VineBase


class VineDensity(VineBase, DensityMixin):
  """
  Vine-copula based density estimator.

  This estimator fits a vine copula to the joint distribution of input features
  and provides methods for density evaluation and sampling.
  """

  def __init__(
    self,
    controls: object | None = None,
    structure: object | None = None,
    batch_size: int = 100,
    schema: dict[str, list[str]] | None = None,
  ) -> None:
    """
    Vine-copula based density estimator.

    Parameters
    ----------
    controls : pv.FitControlsVinecop, optional
        Controls for vinecop fitting. If None, defaults to tll family with 1 thread.
    structure : pv.RVineStructure, optional
        Vine structure. If None, structure will be selected automatically.
    batch_size : int, default=100
        Number of test points to process per batch when evaluating density.
    schema : dict | None, default=None
        Pre-specified metadata about the input. If None, it will be inferred from the training data.
        Currently, only _kde1d_types is used.
        Example: {"kde1d_types": ["continuous", "discrete", "continuous", ...]}
        Supported types are "continuous" and "discrete" (for ordered variables).
    """
    if controls is not None and not isinstance(controls, pv.FitControlsVinecop):
      raise TypeError(
        f"controls must be pv.FitControlsVinecop or None, got {type(controls).__name__}"
      )
    if structure is not None and not isinstance(structure, pv.RVineStructure):
      raise TypeError(
        f"structure must be pv.RVineStructure or None, got {type(structure).__name__}"
      )
    super().__init__(
      controls=controls,
      structure=structure,
      batch_size=batch_size,
      schema=schema,
    )

  def fit(self, X: np.ndarray | pd.DataFrame) -> "VineDensity":
    """
    Fit the density estimator to the input data.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Training data.

    Returns
    -------
    self : VineDensity
        Fitted estimator.
    """
    result = self._check_and_expand_fit(X)
    assert isinstance(result, np.ndarray)  # y is None ⇒ scalar X return
    X = result

    self._fit_marginals(X)
    U = self._to_u_scale(X)

    assert self.schema is not None  # Guaranteed after _check_and_expand_fit
    var_types = [x[0] for x in self.schema["kde1d_types"]]
    self._fit_vine(U, var_types=var_types)
    return self

  def score_samples(self, X: np.ndarray | pd.DataFrame) -> np.ndarray:
    """
    Evaluate the log-likelihood of each sample under the fitted density.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Test samples.

    Returns
    -------
    log_likelihood : ndarray of shape (n_samples,)
        Log-likelihood of each sample.
    """
    return self._pdf_samples(X, y=None, log=True, copula_only=False)

  def score(
    self,
    X: Any,  # MatrixLike per DensityMixin; runtime accepts ndarray/DataFrame.
    y: Any = None,
  ) -> float:
    """
    Compute the average log-likelihood of the samples.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Test samples. We support np.ndarray and pd.DataFrame.
        Sparse matrices are not supported.
    y : array-like, optional
        Not used, present for API consistency with DensityMixin.

    Returns
    -------
    score : float
        Average log-likelihood.

    Raises
    ------
    TypeError
        If X is a sparse matrix or other unsupported array type.
    """

    if isinstance(X, pd.DataFrame):
      X = X.to_numpy()
    elif not isinstance(X, np.ndarray):
      try:
        X = np.asarray(X)
      except Exception as e:
        raise TypeError(f"Unsupported array type {type(X)} for scoring") from e

    return float(self.score_samples(X).mean())

  def sample(self, n_samples: int = 1, seeds: list[int] = [42]) -> np.ndarray:
    """
    Generate random samples from the fitted density.

    Parameters
    ----------
    n_samples : int, default=1
        Number of samples to generate.
    seeds : list of int, default=[42]
        List of seeds for reproducibility.

    Returns
    -------
    X_sampled : ndarray of shape (n_samples, n_features)
        Generated samples.
    """
    if not hasattr(self, "_vine"):
      raise RuntimeError("Model not fitted yet.")

    U_sampled = np.asarray(self._vine.simulate(n_samples, seeds=seeds))
    X_sampled = np.empty((n_samples, self.n_features_in_))
    for j in range(self.n_features_in_):
      X_sampled[:, j] = self._x_kde1d[j].quantile(U_sampled[:, j])
    return X_sampled

  def pdf(self, X: np.ndarray, copula_only: bool = False) -> np.ndarray:
    """
    Evaluate the probability density function.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Test samples.
    copula_only : bool, default=True
        If True, return only the copula density. If False, return the full joint density

    Returns
    -------
    density : ndarray of shape (n_samples,)
        Density values.
    """
    return self._pdf_samples(X, y=None, log=False, copula_only=False)
