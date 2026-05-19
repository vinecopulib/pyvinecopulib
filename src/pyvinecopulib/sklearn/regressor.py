import numpy as np
from sklearn.base import RegressorMixin

import pyvinecopulib as pv

from ._base import (
  _DOC_DISCRETE,
  _DOC_FACTORIZATION,
  _DOC_PIPELINE,
  _DOC_REFERENCES,
  _DOC_WRAPPER,
  VineBase,
)


class VineRegressor(VineBase, RegressorMixin):
  def __init__(
    self,
    mean: bool = True,
    quantiles: list[float] | np.ndarray | None = None,
    controls: pv.FitControlsVinecop | None = None,
    structure: pv.RVineStructure | None = None,
    batch_size: int = 100,
    use_grid: bool = True,
  ) -> None:
    """Sklearn-compatible vine-copula regressor.

    Predicts the conditional mean :math:`\\hat{\\mathbb{E}}[Y \\mid X = x]`
    and/or conditional quantiles using the weighted-sample estimator
    derived in the class docstring.

    Parameters
    ----------
    mean : bool, default=True
        If ``True``, predict the conditional mean. Set to ``False`` to
        get quantile-only predictions (``quantiles`` must then be set).
    quantiles : list of float or ndarray, optional
        Quantile levels in ``(0, 1)`` to predict. If ``None``, quantile
        prediction is disabled.
    controls : :class:`pyvinecopulib.core.FitControlsVinecop`, optional
        Controls forwarded to the underlying vine copula fit
        (:meth:`pyvinecopulib.core.Vinecop.from_data`). If ``None``,
        defaults to the ``tll`` nonparametric pair family with a
        single thread.
    structure : :class:`pyvinecopulib.core.RVineStructure`, optional
        Pre-specified vine structure on ``(Y, X_1, ..., X_d)``. ``Y``
        always occupies the first dimension. If ``None``, structure
        selection is delegated to
        :meth:`pyvinecopulib.core.Vinecop.from_data` (Dissmann
        algorithm by default).
    batch_size : int, default=100
        Number of test points processed per batch in :meth:`predict`.
        ``1`` minimises memory at the cost of speed; ``n_test`` is the
        opposite extreme; intermediate values trade off memory and
        throughput.
    use_grid : bool, default=True
        Controls how training responses are represented for the
        weighted-sample predictor:

        - ``False`` (importance weighting over training rows): weights
          are :math:`w_i(x) \\propto c_{Y,X}(\\hat F_Y(y_i),
          \\hat F_X(x))` for :math:`i = 1, \\dots, n_{\\text{train}}`.
        - ``True`` (importance weighting over a marginal grid): the
          training response set is replaced by the Kde1d grid points
          and the weights pick up an extra :math:`\\hat f_Y(y_g)`
          factor. Cheaper when ``n_train`` is large because the grid
          size stays fixed (independent of ``n_train``); see
          Nagler & Vatter (2024), and the supplement to Vatter & Nagler
          (2026) for the grid-variant derivation.

    Notes
    -----
    Advanced / ensemble use: ``self._normalize_weights`` (default
    ``True``) can be assigned between construction and :meth:`fit` to
    disable the row-wise sum-to-one normalisation of weights inside
    :meth:`_iter_weights`. Forest wrappers set it to ``False`` so they
    can average raw weights across trees and normalise once at the
    ensemble level. Non-default values are not preserved by
    :func:`sklearn.base.clone`.
    """
    if not isinstance(mean, bool):
      raise TypeError(f"mean must be bool, got {type(mean).__name__}")
    if not isinstance(use_grid, bool):
      raise TypeError(f"use_grid must be bool, got {type(use_grid).__name__}")

    super().__init__(
      controls=controls,
      structure=structure,
      batch_size=batch_size,
    )

    self.mean = mean
    if quantiles is None:
      self.quantiles = None
    else:
      q_arr = np.atleast_1d(np.asarray(quantiles, dtype=float))
      if q_arr.ndim != 1 or q_arr.size == 0:
        raise ValueError(
          f"quantiles must be a non-empty 1d sequence, got {quantiles!r}"
        )
      if np.any(q_arr <= 0) or np.any(q_arr >= 1):
        raise ValueError(f"quantiles must lie in (0, 1), got {quantiles!r}")
      self.quantiles = q_arr
    if (not self.mean) and (self.quantiles is None):
      raise ValueError("At least one of mean or quantiles must be enabled.")
    self.use_grid = use_grid
    self._normalize_weights: bool = True

  def fit(self, X: np.ndarray, y: np.ndarray) -> "VineRegressor":
    """Fit a vine copula to the joint distribution of ``(Y, X)``.

    Runs the shared three-step pipeline (marginal KDEs →
    pseudo-observations → vine copula) on the joint vector
    :math:`(Y, X_1, \\dots, X_d)`, with :math:`Y` in the first
    position and treated as continuous. See the class docstring for
    the full mathematical setup.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Training covariates.
    y : ndarray of shape (n_samples,)
        Training responses (continuous).

    Returns
    -------
    self : VineRegressor
        Fitted estimator. ``self._vine``, ``self._x_kde1d``,
        ``self._y_kde1d``, ``self._y_train``, ``self._uy_train``,
        ``self._schema`` and ``self.n_features_in_`` are set.
    """
    X, y = self._check_and_expand_fit(X, y)
    self._fit_marginals(X, y)

    uy_train = self._to_u_scale(y, is_y=True)
    ux = self._to_u_scale(X)

    assert self._schema is not None  # Guaranteed after _check_and_expand_fit
    var_types = ["c"] + [x[0] for x in self._schema["kde1d_types"]]
    self._fit_vine(np.column_stack([uy_train, ux]), var_types=var_types)

    if not self.use_grid:
      self._y_train = y
      self._uy_train = uy_train
    else:
      self._y_train = self._y_kde1d.grid_points
      uy_cdf = self._y_kde1d.cdf(self._y_train)
      self._uy_train = np.asarray(uy_cdf).reshape(-1, 1)
      self._y_density = self._y_kde1d.values[np.newaxis, :]
    return self

  def _copula_marginal_density(
    self, X: np.ndarray, log: bool = False, n_grid: int = 101
  ) -> np.ndarray:
    """Numerical approximation of :math:`c_X(u_X)`.

    Computes

    .. math::

       c_X(u_X) = \\int_0^1 c_{Y, X}(u_Y, u_X)\\, du_Y

    by Simpson's rule on a uniform grid in :math:`[0, 1]`. Used
    internally by forest-style ensembles to evaluate the conditional
    log-likelihood
    :math:`\\log f_{Y \\mid X}(y \\mid x) = \\log c_{Y,X}(u_Y, u_X) -
    \\log c_X(u_X) + \\log f_Y(y)`; not normally called by users.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Conditioning covariates.
    log : bool, default=False
        Whether to return the log-density.
    n_grid : int, default=101
        Number of integration nodes. Simpson's rule requires an odd
        count; an even value is silently incremented.

    Returns
    -------
    ndarray of shape (n_samples,)
        Marginal copula density :math:`c_X(u_X)` (or its log).
    """
    if not hasattr(self, "_vine"):
      raise RuntimeError("Model not fitted yet.")

    X = np.asarray(X)
    ux = self._to_u_scale(X)
    n_test = ux.shape[0]

    if n_grid % 2 == 0:
      n_grid += 1

    eps = 1e-10
    uy_nodes = np.linspace(eps, 1 - eps, n_grid).reshape(-1, 1)

    # Simpson weights: 1,4,2,4,...,4,1
    w = np.ones(n_grid)
    w[1:-1:2] = 4
    w[2:-1:2] = 2
    simpson_factor = 1.0 / (n_grid - 1) / 3.0

    out = np.empty(n_test)

    for start in range(0, n_test, self.batch_size):
      end = min(start + self.batch_size, n_test)
      ux_batch = np.repeat(ux[start:end], n_grid, axis=0)
      uy_rep = np.tile(uy_nodes, (end - start, 1))
      u = np.column_stack([uy_rep, ux_batch])

      vals = self._vine.pdf(u, num_threads=self.controls.num_threads)
      vals = np.asarray(vals).reshape(end - start, n_grid)
      out[start:end] = simpson_factor * (vals * w[None, :]).sum(axis=1)

    return np.log(np.clip(out, eps, None)) if log else out

  # Generator form (`yield (w, start, end)` per batch) is preserved so a
  # future forest PR can interleave batches across trees in `_predict_from_iter`.
  def _iter_weights(self, X):
    """Yield batched conditional weights for the weighted-sample predictor.

    For each batch of test rows, computes weights

    .. math::

       w_i(x) \\propto
       \\begin{cases}
         c_{Y, X}\\bigl(\\hat F_Y(y_i), \\hat F_X(x)\\bigr) &
            \\text{if } \\texttt{use\\_grid = False}, \\\\
         c_{Y, X}\\bigl(\\hat F_Y(y_g), \\hat F_X(x)\\bigr)\\,
         \\hat f_Y(y_g) &
            \\text{if } \\texttt{use\\_grid = True},
       \\end{cases}

    over the training responses (``use_grid=False``) or the Kde1d
    grid points (``use_grid=True``). Weights are normalised row-wise
    when ``normalize_weights=True``. The generator shape lets a
    future forest implementation interleave batches across trees.

    Parameters
    ----------
    X : ndarray of shape (n_test, n_features)
        Test covariates on the original (un-transformed) scale.

    Yields
    ------
    w : ndarray of shape (batch, n_train_or_grid)
        Weights for the current batch.
    start : int
        Index of the first row in this batch.
    end : int
        Index one past the last row in this batch.
    """
    X = np.asarray(X)
    n_test = X.shape[0]
    n_train = self._y_train.shape[0]
    ux_test = self._to_u_scale(X)

    for start in range(0, n_test, self.batch_size):
      end = min(start + self.batch_size, n_test)
      ux_batch = np.repeat(ux_test[start:end], n_train, axis=0)
      uy_rep = np.tile(self._uy_train, (end - start, 1)).reshape(-1, 1)
      u_test = np.column_stack([uy_rep, ux_batch])

      w = self._vine.pdf(u_test, num_threads=self.controls.num_threads)
      w = np.asarray(w).reshape(end - start, n_train)
      if self.use_grid:
        w *= self._y_density
      if self._normalize_weights:
        w /= np.sum(w, axis=1, keepdims=True)

      yield w, start, end

  def _predict_from_iter(self, X, iter_weights):
    """Combine batched weights with training responses to form predictions.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Test covariates (already validated / expanded).
    iter_weights : callable
        Generator returning ``(weights, start, end)`` triples. Usually
        :meth:`_iter_weights`, but a forest ensemble can pass an
        averaged-across-trees variant.

    Returns
    -------
    ndarray
        Predictions of shape ``(n_samples,)`` when a single output is
        produced (mean or single quantile) or
        ``(n_samples, n_outputs)`` otherwise. Column order is mean
        first (if enabled), then quantiles in the order requested.
    """
    n_test = X.shape[0]
    n_outputs = (1 if self.mean else 0) + (
      len(self.quantiles) if self.quantiles is not None else 0
    )
    y_pred = np.empty((n_test, n_outputs))

    for w, start, end in iter_weights(X):
      col = 0
      if self.mean:
        y_pred[start:end, col] = w @ self._y_train
        col += 1
      if self.quantiles is not None:
        batch_preds = [
          np.quantile(
            a=self._y_train,
            q=self.quantiles,
            weights=row_w,
            method="inverted_cdf",
          )
          for row_w in w
        ]
        y_pred[start:end, col : col + len(self.quantiles)] = np.vstack(
          batch_preds
        )

    return y_pred.squeeze()

  def predict(self, X):
    """Predict the conditional mean and/or quantiles of ``Y`` given ``X``.

    For each test row :math:`x`, computes the weights
    :math:`w_i(x)` described in :meth:`_iter_weights` and returns the
    weighted statistics:

    - **Mean**: :math:`\\hat{\\mathbb{E}}[Y \\mid X = x] = \\sum_i
      w_i(x)\\, y_i` (closed-form solution of the estimating equation
      :math:`\\int (y - \\beta) \\hat f(y \\mid x)\\, dy = 0`).
    - **Quantile** :math:`\\tau`: weighted quantile of ``self._y_train``
      with the same weights, computed via
      :func:`numpy.quantile` with ``method="inverted_cdf"``.

    Parameters
    ----------
    X : ndarray or DataFrame of shape (n_samples, n_features)
        Test covariates. Must match the training schema.

    Returns
    -------
    ndarray
        Predictions of shape ``(n_samples,)`` if only one output is
        requested (mean *or* a single quantile), or
        ``(n_samples, n_outputs)`` otherwise. Output columns are
        ordered: mean (if ``self.mean``), then quantiles in
        ``self.quantiles`` order.
    """
    X = self._check_and_expand_predict(X)
    return self._predict_from_iter(X, self._iter_weights)


VineRegressor.__doc__ = f"""Vine-copula based regressor (mean and quantile).

A scikit-learn-compatible non-parametric regressor that predicts the
conditional mean :math:`\\mathbb{{E}}[Y \\mid X = x]` and/or
conditional :math:`\\tau`-quantiles of ``Y`` given covariates ``X``,
by fitting a vine copula to the joint distribution of :math:`(Y, X)`
and reducing prediction to a weighted statistic of the training
responses.

**Estimating-equation framework.** Following Nagler & Vatter (2024),
for a target functional :math:`\\beta(x)` characterised by
:math:`\\mathbb{{E}}[\\psi_\\beta(Y) \\mid X = x] = 0` we plug in the
fitted conditional density :math:`\\hat f_{{Y \\mid X}}(\\cdot \\mid x)`
and solve

.. math::

   \\int \\psi_\\beta(y)\\, \\hat f_{{Y \\mid X}}(y \\mid x)\\, dy = 0.

Two standard choices of :math:`\\psi_\\beta` give

- :math:`\\psi_\\beta(y) = y - \\beta` → conditional mean,
- :math:`\\psi_\\beta(y) = \\mathbf{{1}}\\{{y < \\beta\\}} - \\tau` →
  conditional :math:`\\tau`-quantile.

Both reduce in practice to weighted statistics of the training
responses. With weights :math:`w_i(x) \\propto c_{{Y, X}}(\\hat F_Y(y_i),
\\hat F_X(x))` derived from the copula density (and optionally
reweighted by :math:`\\hat f_Y` on a fixed grid; see ``use_grid``),
the conditional mean is :math:`\\sum_i w_i(x)\\, y_i` (closed form)
and the conditional quantile is the weighted quantile of
:math:`\\{{y_i\\}}_i` via :func:`numpy.quantile` with
``method="inverted_cdf"``.

{_DOC_WRAPPER}
{_DOC_PIPELINE}
{_DOC_FACTORIZATION}
{_DOC_DISCRETE}

Examples
--------
>>> import numpy as np
>>> from pyvinecopulib.sklearn import VineRegressor
>>> rng = np.random.default_rng(0)
>>> X = rng.standard_normal((200, 3))
>>> y = X @ [1.5, -0.8, 0.4] + 0.2 * rng.standard_normal(200)
>>> est = VineRegressor(quantiles=[0.1, 0.5, 0.9]).fit(X, y)
>>> est.predict(X[:5])          # columns: mean, q10, q50, q90

{_DOC_REFERENCES}
"""
