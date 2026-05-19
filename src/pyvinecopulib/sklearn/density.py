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
  _DOC_WRAPPER,
  VineBase,
)


class VineDensity(VineBase, DensityMixin):
  def __init__(
    self,
    backend=None,
    batch_size: int = 100,
    random_state=None,
  ) -> None:
    """Vine-copula based density estimator.

    Parameters
    ----------
    backend : :class:`~pyvinecopulib.sklearn.backends.VinecopBackend` or compatible, optional
        Backend strategy bundling fit-time controls and an optional
        pre-specified structure. ``None`` (default) resolves to
        :class:`~pyvinecopulib.sklearn.backends.VinecopBackend` at fit
        time, which uses :meth:`pyvinecopulib.core.Vinecop.from_data`
        with the nonparametric ``tll`` pair family. Pass
        :class:`~pyvinecopulib.sklearn.backends.TorchVinecopBackend`
        for the PyTorch backend.
    batch_size : int, default=100
        Number of test points processed per batch when evaluating the
        density. Higher values trade memory for throughput.
    random_state : int, RandomState, Generator or None
        Seeds the RNG used by stochastic operations (``cdf`` quasi-MC,
        ``sample``). Stored as-is and resolved via
        :func:`sklearn.utils.check_random_state` inside ``fit``.
    """
    super().__init__(
      backend=backend, batch_size=batch_size, random_state=random_state
    )

  def fit(self, X, y=None) -> "VineDensity":
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
    y : ignored
        Present for sklearn API compatibility.

    Returns
    -------
    self : VineDensity
        Fitted estimator. Sets ``self._vine``, ``self._x_kde1d``,
        ``self.schema_``, ``self.structure_``, ``self.backend_``,
        ``self.random_state_``, ``self.n_features_in_``, and
        ``self.feature_names_in_`` (for DataFrame input).
    """
    self._validate_params()
    X = self._validate_input(X, reset=True)
    self._resolve_runtime_state()
    self._fit_marginals(X)
    U = self._to_u_scale(X)
    var_types = [x[0] for x in self.schema_["kde1d_types"]]
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

  def sample(self, n_samples: int = 1, random_state=None) -> np.ndarray:
    """Draw samples from the fitted joint density.

    Draws :math:`U \\sim C` via :meth:`Vinecop.simulate` and pushes
    each component back through the inverse marginal CDF
    :math:`F_j^{-1}` to obtain a sample in the original feature space.

    Parameters
    ----------
    n_samples : int, default=1
        Number of samples to generate.
    random_state : int, RandomState, Generator or None
        If ``None`` (default), reuses the RNG resolved at ``fit`` time
        (``self.random_state_``). Otherwise resolved fresh via
        :func:`sklearn.utils.check_random_state` for this call only —
        useful for getting independent draws without re-fitting.

    Returns
    -------
    ndarray of shape (n_samples, n_features)
        Generated samples in the original feature scale.
    """
    check_is_fitted(self, attributes=["_vine"])
    if random_state is None:
      rng = self.random_state_
    else:
      from sklearn.utils.validation import check_random_state

      rng = check_random_state(random_state)
    seeds = [int(x) for x in rng.randint(0, 2**31 - 1, size=5)]
    U_sampled = np.asarray(
      self.backend_.simulate(self._vine, n_samples, seeds=seeds)
    )
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
    X,
    N: int = 10000,
    random_state=None,
  ) -> np.ndarray:
    """Evaluate the joint CDF at the given samples.

    Returns :math:`\\hat F(\\mathbf{x}) = \\hat C\\bigl(\\hat F_1(x_1),
    \\ldots, \\hat F_d(x_d)\\bigr)` by applying the marginal CDFs to
    each column to obtain pseudo-observations, then evaluating the
    fitted copula CDF via the configured backend's quasi-Monte-Carlo
    routine.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test samples.
    N : int, default=10000
        Number of quasi-random points used by the Monte-Carlo
        integration; larger ``N`` gives more accurate CDF values at the
        cost of more compute.
    random_state : int, RandomState, Generator or None
        If ``None`` (default), reuses the RNG resolved at ``fit`` time
        (``self.random_state_``). Otherwise resolved fresh via
        :func:`sklearn.utils.check_random_state` for this call.

    Returns
    -------
    ndarray of shape (n_samples,)
        Joint CDF values in :math:`[0, 1]`.

    Notes
    -----
    Because the underlying copula CDF is approximated by Monte-Carlo
    quasi-random integration, values are stochastic across calls with
    different ``random_state`` values. Use a larger ``N`` if the noise
    floor is significant for your application.
    """
    check_is_fitted(self, attributes=["_vine"])
    if not self.backend_.supports_cdf:
      raise NotImplementedError(
        f"backend '{self.backend_.name}' does not support cdf()."
      )
    X = self._validate_input(X, reset=False)
    U = self._to_u_scale(X)
    if random_state is None:
      rng = self.random_state_
    else:
      from sklearn.utils.validation import check_random_state

      rng = check_random_state(random_state)
    seeds = [int(x) for x in rng.randint(0, 2**31 - 1, size=5)]
    return np.asarray(self.backend_.cdf(self._vine, U, N=N, seeds=seeds))


VineDensity.__doc__ = f"""Vine-copula based density estimator.

A scikit-learn-compatible non-parametric joint density estimator. The
joint density is factorized into one-dimensional marginals (estimated
with kernel density) and a vine copula capturing the dependence
structure. After fitting, :meth:`pdf` / :meth:`score_samples` evaluate
the density at new points and :meth:`sample` draws from the fitted
joint distribution.

{_DOC_WRAPPER}
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

Use the PyTorch backend for GPU placement / autograd:

>>> from pyvinecopulib.sklearn.backends import TorchVinecopBackend
>>> density_gpu = VineDensity(backend=TorchVinecopBackend()).fit(X)

See also
--------

* :class:`VineRegressor` — sister regressor on the same primitives.
* :class:`VineForestDensity` — ensemble variant.
* :class:`pyvinecopulib.sklearn.backends.VinecopBackend`,
  :class:`pyvinecopulib.sklearn.backends.TorchVinecopBackend` — the
  two backend choices.
* :class:`pyvinecopulib.core.Vinecop` — the underlying vine-copula
  class if you need control beyond the sklearn convenience layer.

{_DOC_REFERENCES}
"""
