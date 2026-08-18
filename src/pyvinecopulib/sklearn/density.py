from typing import Any

import numpy as np
import pandas as pd
from sklearn.base import DensityMixin
from sklearn.utils.validation import check_is_fitted

from ._base import (
  _DOC_DISCRETE,
  _DOC_FACTORIZATION,
  _DOC_PIPELINE,
  _DOC_REFERENCES,
  VineBase,
)


class VineDensity(DensityMixin, VineBase):
  def __init__(
    self,
    backend=None,
    margins=None,
    batch_size: int = 100,
    random_state=None,
  ) -> None:
    """Vine-copula based density estimator.

    Parameters
    ----------
    backend : VinecopBackend or compatible, default=None
        Backend instance bundling fit-time controls and an optional
        pre-specified structure. `None` resolves to a default
        ``VinecopBackend`` at fit time, which calls ``Vinecop.from_data()``
        with the nonparametric ``tll`` pair family. Pass
        ``TorchVinecopBackend`` for the PyTorch backend.
    margins : object, default=None
        The marginal half of the model, in any form
        :func:`pyvinecopulib.margins.resolve_margins` accepts. `None`
        fits a ``Kde1dMargin`` per column with the variable type
        inferred from the input.
    batch_size : int, default=100
        Number of test points processed per batch when evaluating
        the density. Higher values trade memory for throughput.
    random_state : int, RandomState instance or None, default=None
        Seeds the RNG used by stochastic operations (`cdf` quasi-MC,
        `sample`). Resolved via `sklearn.utils.check_random_state`
        inside `fit`.
    """
    super().__init__(
      backend=backend,
      margins=margins,
      batch_size=batch_size,
      random_state=random_state,
    )

  def fit(self, X, y=None) -> "VineDensity":
    """Fits the joint density to the training data.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype float, or DataFrame
        Training data. DataFrame columns may mix numeric,
        ordered-categorical, and unordered-categorical dtypes;
        unordered categoricals are expanded to ordered ``{0, 1}``
        dummies before fitting.
    y : None
        Ignored; present for sklearn API compatibility.

    Returns
    -------
    self : VineDensity
        The fitted estimator.
    """
    self._validate_params()
    X = self._validate_input(X, reset=True)
    self._resolve_runtime_state()
    self._fit_marginals(X)
    self._fit_vine(self._to_u_scale(X))
    self._bind_distribution(self._x_margins)
    return self

  def score_samples(self, X: np.ndarray | pd.DataFrame) -> np.ndarray:
    """Evaluates the per-sample log-likelihood under the fitted density.

    Equivalent to ``np.log(self.pdf(X, copula_only=False))`` but
    computed in log-space throughout for numerical stability.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype float, or DataFrame
        Test samples. Columns must match the training schema.

    Returns
    -------
    ndarray, shape (n_samples,), dtype float
        Log-density values.
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
    X : ndarray, shape (n_samples, n_features), dtype float, or DataFrame
        Test samples. Sparse matrices are not supported.
    y : None
        Ignored; present for `DensityMixin.score` compatibility.

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

  def sample(self, n_samples: int = 1, random_state=None) -> np.ndarray:
    """Draws samples from the fitted joint density.

    Samples :math:`U \\sim C` from the fitted copula and pushes each
    component back through the inverse marginal CDF :math:`F_j^{-1}`
    to obtain a sample in the original feature space.

    Parameters
    ----------
    n_samples : int, default=1
        Number of samples to generate.
    random_state : int, RandomState instance or None, default=None
        Resolved fresh via `sklearn.utils.check_random_state` for
        this call. `None` reuses the RNG resolved at fit time.

    Returns
    -------
    ndarray, shape (n_samples, n_features), dtype float
        Generated samples, in the *expanded* feature space -- the one
        `fit` modeled, with ``n_features_in_`` columns. Dummies of an
        expanded unordered categorical come back as separate columns
        rather than as the original level, since a vine draw can put a
        ``1`` in two of them at once and choosing between those is a
        modeling decision, not an inverse.
    """
    check_is_fitted(self, attributes=["distribution_"])
    if random_state is None:
      rng = self.random_state_
    else:
      from sklearn.utils.validation import check_random_state

      rng = check_random_state(random_state)
    seeds = [int(x) for x in rng.randint(0, 2**31 - 1, size=5)]
    return np.asarray(self.distribution_.simulate(n_samples, seeds=seeds))

  def pdf(self, X: np.ndarray, copula_only: bool = False) -> np.ndarray:
    """Evaluates the joint density at the given samples.

    Returns the full joint density
    :math:`\\hat f(\\mathbf{x}) = \\hat c(\\hat F_1(x_1), \\ldots,
    \\hat F_d(x_d)) \\prod_j \\hat f_j(x_j)` by default. Set
    ``copula_only=True`` to skip the marginal-density product and
    return only the copula factor.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype float, or DataFrame
        Test samples.
    copula_only : bool, default=False
        If ``True``, return only the copula density evaluated at
        the pseudo-observations.

    Returns
    -------
    ndarray, shape (n_samples,), dtype float
        Density values.
    """
    return self._pdf_samples(X, y=None, log=False, copula_only=copula_only)

  def cdf(
    self,
    X,
    N: int = 10000,
    random_state=None,
  ) -> np.ndarray:
    """Evaluates the joint CDF at the given samples.

    Returns
    :math:`\\hat F(\\mathbf{x}) = \\hat C(\\hat F_1(x_1), \\ldots,
    \\hat F_d(x_d))` by applying the marginal CDFs to obtain
    pseudo-observations and evaluating the fitted copula CDF via
    the backend's quasi-Monte-Carlo routine.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype float, or DataFrame
        Test samples.
    N : int, default=10000
        Number of quasi-random points used for the Monte-Carlo
        integration. Larger `N` gives more accurate CDF values at
        the cost of more compute.
    random_state : int, RandomState instance or None, default=None
        Resolved fresh via `sklearn.utils.check_random_state` for
        this call. `None` reuses the RNG resolved at fit time.

    Returns
    -------
    ndarray, shape (n_samples,), dtype float
        Joint CDF values in :math:`[0, 1]`. Values are stochastic
        across calls with different `random_state`; use a larger
        `N` if the MC noise floor is significant for your
        application.
    """
    check_is_fitted(self, attributes=["distribution_"])
    X = self._validate_input(X, reset=False)
    if random_state is None:
      rng = self.random_state_
    else:
      from sklearn.utils.validation import check_random_state

      rng = check_random_state(random_state)
    seeds = [int(x) for x in rng.randint(0, 2**31 - 1, size=5)]
    return np.asarray(
      self.distribution_.cdf(np.asarray(X, dtype=float), N=N, seeds=seeds)
    )


VineDensity.__doc__ = f"""Vine-copula based density estimator.

A scikit-learn-compatible non-parametric joint density estimator.
The joint density is factorized into one-dimensional marginals
(estimated with kernel density) and a vine copula capturing the
dependence structure. After fitting, `pdf` / `score_samples`
evaluate the density at new points and `sample` draws from the
fitted joint distribution.

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
>>> density.score_samples(X[:3])

{_DOC_REFERENCES}
"""
