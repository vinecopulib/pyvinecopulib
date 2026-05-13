from typing import Any

import numpy as np
import pandas as pd
import pyvinecopulib as pv
from sklearn.base import DensityMixin

from ._base import (
  _DOC_DISCRETE,
  _DOC_FACTORIZATION,
  _DOC_PIPELINE,
  _DOC_REFERENCES,
  VineBase,
)


class VineDensity(VineBase, DensityMixin):
  def __init__(
    self,
    controls: pv.FitControlsVinecop | None = None,
    structure: pv.RVineStructure | None = None,
    batch_size: int = 100,
  ) -> None:
    """Vine-copula based density estimator.

    Parameters
    ----------
    controls : pyvinecopulib.core.FitControlsVinecop, optional
        Controls for the vine copula fit. If ``None``, defaults to the
        nonparametric ``tll`` pair family with a single thread.
    structure : pyvinecopulib.core.RVineStructure, optional
        Pre-specified vine structure. If ``None``, the structure is
        selected automatically via :meth:`Vinecop.from_data`.
    batch_size : int, default=100
        Number of test points processed per batch when evaluating the
        density. Higher values trade memory for throughput.
    """
    super().__init__(
      controls=controls,
      structure=structure,
      batch_size=batch_size,
    )

  def fit(self, X: np.ndarray | pd.DataFrame) -> "VineDensity":
    """Fit the density estimator to the training data.

    Runs the shared three-step pipeline (marginal KDEs →
    pseudo-observations → vine copula) on ``X``; see the class
    docstring for the full description.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Training data. ``DataFrame`` input may mix numeric, ordered
        categorical, and unordered categorical columns; the latter are
        expanded to ordered ``{0, 1}`` dummies before fitting.

    Returns
    -------
    self : VineDensity
        Fitted estimator. ``self._vine``, ``self._x_kde1d``,
        ``self._schema`` and ``self.n_features_in_`` are set.
    """
    result = self._check_and_expand_fit(X)
    assert isinstance(result, np.ndarray)  # y is None ⇒ scalar X return
    X = result

    self._fit_marginals(X)
    U = self._to_u_scale(X)

    assert self._schema is not None  # Guaranteed after _check_and_expand_fit
    var_types = [x[0] for x in self._schema["kde1d_types"]]
    self._fit_vine(U, var_types=var_types)
    return self

  def score_samples(self, X: np.ndarray | pd.DataFrame) -> np.ndarray:
    """Evaluate the per-sample log-likelihood under the fitted density.

    Equivalent to ``np.log(self.pdf(X, copula_only=False))`` but
    computed in log-space throughout for numerical stability.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test samples. Must have the same columns / dtypes / order as
        the training data.

    Returns
    -------
    ndarray of shape (n_samples,)
        Log-density values :math:`\\log \\hat f(\\mathbf{x}_i)`.
    """
    return self._pdf_samples(X, y=None, log=True, copula_only=False)

  def score(
    self,
    X: Any,  # MatrixLike per DensityMixin; runtime accepts ndarray/DataFrame.
    y: Any = None,
  ) -> float:
    """Mean log-likelihood over a sample (sklearn ``score`` convention).

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test samples. We support ``np.ndarray`` and ``pd.DataFrame``.
        Sparse matrices are not supported.
    y : ignored
        Present for API compatibility with ``DensityMixin.score``.

    Returns
    -------
    float
        Mean of ``score_samples(X)``.

    Raises
    ------
    TypeError
        If ``X`` is a sparse matrix or another unsupported array type.
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
    """Draw samples from the fitted joint density.

    Draws :math:`U \\sim C` via :meth:`Vinecop.simulate` and pushes
    each component back through the inverse marginal CDF
    :math:`F_j^{-1}` to obtain a sample in the original feature space.

    Parameters
    ----------
    n_samples : int, default=1
        Number of samples to generate.
    seeds : list of int, default=[42]
        Seeds forwarded to :meth:`Vinecop.simulate` for reproducibility.

    Returns
    -------
    ndarray of shape (n_samples, n_features)
        Generated samples in the original feature scale.
    """
    if not hasattr(self, "_vine"):
      raise RuntimeError("Model not fitted yet.")

    U_sampled = np.asarray(self._vine.simulate(n_samples, seeds=seeds))
    X_sampled = np.empty((n_samples, self.n_features_in_))
    for j in range(self.n_features_in_):
      X_sampled[:, j] = self._x_kde1d[j].quantile(U_sampled[:, j])
    return X_sampled

  def pdf(self, X: np.ndarray, copula_only: bool = False) -> np.ndarray:
    """Evaluate the joint density at the given samples.

    By default returns the full joint density
    :math:`\\hat f(\\mathbf{x}) = \\hat c(\\hat F_1(x_1), \\ldots,
    \\hat F_d(x_d)) \\prod_j \\hat f_j(x_j)`. Set ``copula_only=True``
    to skip the marginal-density product and return only the copula
    factor :math:`\\hat c(\\hat F_1(x_1), \\ldots, \\hat F_d(x_d))`.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test samples.
    copula_only : bool, default=False
        If ``True``, return the copula density evaluated at the
        pseudo-observations only (skip the marginal-density product).

    Returns
    -------
    ndarray of shape (n_samples,)
        Density values.
    """
    return self._pdf_samples(X, y=None, log=False, copula_only=copula_only)

  def cdf(
    self,
    X: np.ndarray | pd.DataFrame,
    N: int = 10000,
    seeds: list[int] | None = None,
  ) -> np.ndarray:
    """Evaluate the joint CDF at the given samples.

    Returns :math:`\\hat F(\\mathbf{x}) = \\hat C\\bigl(\\hat F_1(x_1),
    \\ldots, \\hat F_d(x_d)\\bigr)` by applying the marginal CDFs to
    each column to obtain pseudo-observations, then evaluating the
    fitted copula CDF via :meth:`pyvinecopulib.core.Vinecop.cdf`.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test samples.
    N : int, default=10000
        Number of quasi-random points used by :meth:`Vinecop.cdf` for
        the Monte-Carlo integration; larger ``N`` gives more accurate
        CDF values at the cost of more compute.
    seeds : list of int, optional
        Seeds forwarded to the underlying quasi-random generator. If
        ``None``, the generator is seeded randomly (results then
        differ from one call to the next).

    Returns
    -------
    ndarray of shape (n_samples,)
        Joint CDF values in :math:`[0, 1]`.

    Notes
    -----
    Because the underlying copula CDF is approximated by Monte-Carlo
    quasi-random integration, values are stochastic for fixed inputs
    unless ``seeds`` is provided. Use a larger ``N`` if the noise
    floor is significant for your application.
    """
    if not hasattr(self, "_vine"):
      raise RuntimeError("Model not fitted yet.")
    X = self._check_and_expand_predict(X)
    U = self._to_u_scale(X)
    return np.asarray(
      self._vine.cdf(U, N=N, seeds=seeds if seeds is not None else [])
    )


VineDensity.__doc__ = f"""Vine-copula based density estimator.

A scikit-learn-compatible non-parametric joint density estimator. The
joint density is factorized into one-dimensional marginals (estimated
with kernel density) and a vine copula capturing the dependence
structure. After fitting, :meth:`pdf` / :meth:`score_samples` evaluate
the density at new points and :meth:`sample` draws from the fitted
joint distribution.

{_DOC_PIPELINE}
{_DOC_FACTORIZATION}
{_DOC_DISCRETE}

Examples
--------
>>> import numpy as np
>>> from pyvinecopulib.sklearn import VineDensity
>>> rng = np.random.default_rng(0)
>>> X = rng.standard_normal((500, 3))
>>> density = VineDensity().fit(X)
>>> density.score_samples(X[:3])          # log-density at first three rows
>>> density.pdf(X[:3])                    # density on the natural scale
>>> density.cdf(X[:3])                    # joint CDF (MC-approximated)
>>> samples = density.sample(n_samples=100)

{_DOC_REFERENCES}
"""
