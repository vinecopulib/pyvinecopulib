import math
from numbers import Integral

import numpy as np
from sklearn.base import RegressorMixin
from sklearn.metrics import r2_score
from sklearn.utils._param_validation import Interval
from sklearn.utils.validation import check_is_fitted

from ..core import Vinedist
from ._base import (
  _DOC_DISCRETE,
  _DOC_FACTORIZATION,
  _DOC_PIPELINE,
  _DOC_REFERENCES,
  VineBase,
  _as_ndarray,
)

# Half-width, in standard deviations, of the probit substitution behind the
# quadrature nodes: the outermost node sits at Phi(-a), so this is how far into
# the response's tails the rule reaches. Five leaves 3e-7 of the probability
# mass outside the grid, worth about 1e-5 response standard deviations on the
# conditional mean; three leaves enough out to cost 3e-2 of one.
_PROBIT_HALF_WIDTH = 5.0


class VineRegressor(RegressorMixin, VineBase):
  # Inherits VineBase._parameter_constraints; extend with regressor knobs.
  _parameter_constraints: dict = {
    **VineBase._parameter_constraints,
    "mean": ["boolean"],
    "quantiles": ["array-like", None],
    "use_grid": ["boolean"],
    "n_nodes": [Interval(Integral, 2, None, closed="left")],
    "normalize_weights": ["boolean"],
  }

  def __init__(
    self,
    mean: bool = True,
    quantiles=None,
    backend=None,
    margins=None,
    batch_size: int = 100,
    use_grid: bool = True,
    n_nodes: int = 401,
    normalize_weights: bool = True,
    random_state=None,
    n_jobs=None,
  ) -> None:
    """Sklearn-compatible vine-copula regressor.

    Predicts the conditional mean
    :math:`\\hat{\\mathbb{E}}[Y \\mid X = x]` and/or conditional
    quantiles using the weighted-sample estimator derived in the
    class docstring.

    Parameters
    ----------
    mean : bool, default=True
        If ``True``, predict the conditional mean. Set to ``False``
        to get quantile-only predictions (``quantiles`` must then be
        set).
    quantiles : array-like of float, shape (n_quantiles,), default=None
        Quantile levels in ``(0, 1)`` to predict. ``None`` disables
        quantile prediction.
    backend : VinecopBackend or compatible, default=None
        Backend instance bundling fit-time controls and an optional
        pre-specified structure on ``(Y, X_1, ..., X_d)`` (`Y`
        always in the first dimension). `None` resolves to a default
        ``VinecopBackend`` with the ``tll`` pair family at fit time.
    margins : object, default=None
        The marginal half of the model, in any form
        :func:`pyvinecopulib.margins.resolve_margins` accepts. `None`
        fits a ``Kde1d`` per column. The specification addresses
        the covariates; one that broadcasts (an alias, a single margin,
        a callable) applies to the response as well. The response
        margin must be continuous, and any such margin works with
        either `use_grid` setting.
    batch_size : int, default=100
        Number of test points processed per batch in `predict`.
    use_grid : bool, default=True
        Where the response values entering the weighted statistic
        come from. ``True`` (default) evaluates
        :math:`\\mathbb{E}[Y \\mid X = x] = \\int_0^1 \\hat
        F_Y^{-1}(p)\\, c_{Y,X}(p, \\hat F_X(x))\\, dp` on a fixed
        grid of `n_nodes` probability levels, which costs the same
        whatever the sample size. ``False`` sums over the training
        rows instead, with
        :math:`w_i(x) \\propto c_{Y,X}(\\hat F_Y(y_i), \\hat
        F_X(x))` -- the same integral approximated by Monte Carlo,
        which in exchange cannot reach past the largest observed
        response.
    n_nodes : int, default=401
        Number of quadrature nodes when ``use_grid=True``, ignored
        otherwise. The default keeps the quadrature error on the
        conditional mean near 1e-5 response standard deviations, well
        under the model's own error; predicted quantiles are resolved
        to about one node spacing, so raise it when they are the
        output of interest.
    normalize_weights : bool, default=True
        If ``True`` (default), per-row weights produced by
        `_iter_weights` are rescaled to sum to one. Set to
        ``False`` to get the raw copula weights instead -- useful
        when a caller combines the weights of several fitted
        estimators and wants to rescale once, after combining.
    random_state : int, RandomState instance or None, default=None
        Seeds the RNG used by stochastic operations. Resolved via
        `sklearn.utils.check_random_state` inside `fit`.
    n_jobs : int or None, default=None
        Threads the vine may use, for fitting and for every evaluation
        (`pdf`, `cdf`, `sample`, and the prediction paths built on them).
        `None` means one thread and `-1` every processor, following the
        scikit-learn convention. Results never depend on it: the fitted
        structure, the fitted pair copulas and every evaluated value are
        bit-identical at any thread count.

        `None` is deliberate: a caller that parallelizes *over* vines owns
        the parallelism, and nesting it would oversubscribe the machine. Set
        it when a single vine is the whole job.
    """
    super().__init__(
      backend=backend,
      margins=margins,
      batch_size=batch_size,
      random_state=random_state,
      n_jobs=n_jobs,
    )
    self.mean = mean
    self.quantiles = quantiles
    self.use_grid = use_grid
    self.n_nodes = n_nodes
    self.normalize_weights = normalize_weights

  def fit(self, X: np.ndarray, y: np.ndarray) -> "VineRegressor":
    """Fits a vine copula to the joint distribution of ``(Y, X)``.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype float
        Training covariates.
    y : ndarray, shape (n_samples,), dtype float
        Training responses (continuous).

    Returns
    -------
    self : VineRegressor
        The fitted estimator.
    """
    self._validate_params()
    if y is None:
      raise ValueError(
        "VineRegressor requires y to be passed, but the target y is None"
      )
    X, y = self._validate_input(X, y, reset=True)
    self._resolve_runtime_state()
    if self.quantiles is None and not self.mean:
      raise ValueError("At least one of mean or quantiles must be enabled.")
    if self.quantiles is not None:
      q_arr = np.atleast_1d(np.asarray(self.quantiles, dtype=float))
      if q_arr.ndim != 1 or q_arr.size == 0:
        raise ValueError(
          f"quantiles must be a non-empty 1d sequence, got {self.quantiles!r}"
        )
      if np.any(q_arr <= 0) or np.any(q_arr >= 1):
        raise ValueError(
          f"quantiles must lie in (0, 1), got {self.quantiles!r}"
        )
      self.quantiles_ = q_arr
    else:
      self.quantiles_ = None
    self._fit_marginals(X, y)
    if getattr(self._y_margin, "supports_covariates", False):
      raise ValueError(
        "VineRegressor cannot use a conditional response margin "
        f"({type(self._y_margin).__name__} declares supports_covariates): the "
        "quadrature inverts one shared probability grid, which a margin whose "
        "quantiles move with the covariates turns into one grid per test row."
      )

    uy_train = self._to_u_scale(y, is_y=True)
    ux = self._to_u_scale(X)

    # The response leads the joint model, and being continuous it adds no
    # left-limit column, so its `u` column simply prepends to the covariates'
    # own layout.
    margins = (self._y_margin, *self._x_margins)
    var_types = Vinedist.copula_var_types(margins)
    self._fit_vine(np.column_stack([uy_train, ux]), var_types=var_types)
    self._bind_distribution(margins)

    if not self.use_grid:
      self._u_nodes = uy_train
      self._y_nodes = y
      self._node_weights = None
    else:
      p_nodes, node_weights = self._probability_grid()
      self._u_nodes = p_nodes.reshape(-1, 1)
      self._y_nodes = _as_ndarray(self._y_margin.icdf(p_nodes)).ravel()
      self._node_weights = node_weights
    return self

  def _probability_grid(self) -> tuple[np.ndarray, np.ndarray]:
    """Quadrature nodes on ``(0, 1)``, and their weights.

    The rule integrates
    :math:`\\int_0^1 F_Y^{-1}(p)\\, c(p, u_x)\\, dp` under the
    substitution :math:`p = \\Phi(z)` on a uniform ``z`` grid, which
    spends nodes on the tails of the response rather than on the bulk
    of its probability scale.

    Returns
    -------
    p : ndarray, shape (n_nodes,), dtype float
        Increasing probability levels, all strictly inside ``(0, 1)``.
    dp : ndarray, shape (1, n_nodes), dtype float
        Node weights :math:`dp/dz`, up to the constant factor the
        row normalization absorbs.
    """
    z = np.linspace(-_PROBIT_HALF_WIDTH, _PROBIT_HALF_WIDTH, int(self.n_nodes))
    root_two = math.sqrt(2.0)
    p = np.array([0.5 * (1.0 + math.erf(zi / root_two)) for zi in z])
    return p, np.exp(-0.5 * z**2)[np.newaxis, :]

  def _copula_marginal_density(
    self, X: np.ndarray, log: bool = False, n_grid: int = 101
  ) -> np.ndarray:
    """Computes :math:`c_X(u_X) = \\int_0^1 c_{Y, X}(u_Y, u_X)\\, du_Y`.

    Numerical approximation via Simpson's rule. This is the term
    that turns the joint copula density into the conditional
    log-likelihood
    :math:`\\log f_{Y \\mid X}(y \\mid x) = \\log c_{Y,X}(u_Y, u_X)
    - \\log c_X(u_X) + \\log f_Y(y)`.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype float
        Conditioning covariates.
    log : bool, default=False
        If ``True``, return the log-density.
    n_grid : int, default=101
        Number of Simpson nodes. An even value is silently
        incremented.

    Returns
    -------
    ndarray, shape (n_samples,), dtype float
        Marginal copula densities :math:`c_X(u_X)` (or their log).
    """
    check_is_fitted(self, attributes=["_vine"])

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

      vals = self.backend_.pdf(self._vine, u)
      vals = np.asarray(vals).reshape(end - start, n_grid)
      out[start:end] = simpson_factor * (vals * w[None, :]).sum(axis=1)

    return np.log(np.clip(out, eps, None)) if log else out

  def _weights_for_batch(self, X_batch):
    """Conditional copula weights for one batch of test rows.

    Single source of truth for the weight math: `_iter_weights` is
    the batched generator over it and `_predict_from_iter` the
    consumer. Kept as a separate, directly callable seam so external
    code can reuse the exact weight definition. Each weight pairs one
    node :math:`y_k` -- a training response when ``use_grid=False``,
    else :math:`\\hat F_Y^{-1}(p_k)` -- with

    .. math::

       w_k(x) \\propto
       c_{Y, X}\\bigl(u_k, \\hat F_X(x)\\bigr) \\cdot \\Delta_k,

    where :math:`u_k` is :math:`\\hat F_Y(y_k)` over training rows
    and :math:`p_k` on the probability grid, and :math:`\\Delta_k` is
    the node spacing there (constant, hence absent, over training
    rows). Row-normalized when ``normalize_weights=True``.

    Parameters
    ----------
    X_batch : ndarray, shape (batch, n_features), dtype float
        Test covariates on the original (un-transformed) scale.

    Returns
    -------
    w : ndarray, shape (batch, n_nodes), dtype float
        Per-row weights for the batch.
    """
    ux_batch_rows = self._to_u_scale(np.asarray(X_batch))
    m = ux_batch_rows.shape[0]
    n_nodes = self._y_nodes.shape[0]
    ux_batch = np.repeat(ux_batch_rows, n_nodes, axis=0)
    uy_rep = np.tile(self._u_nodes, (m, 1)).reshape(-1, 1)
    u_test = np.column_stack([uy_rep, ux_batch])

    w = np.asarray(self.backend_.pdf(self._vine, u_test)).reshape(m, n_nodes)
    if self._node_weights is not None:
      w = w * self._node_weights
    if self.normalize_weights:
      w /= np.sum(w, axis=1, keepdims=True)
    return w

  def _iter_weights(self, X):
    """Yields ``(weights, start, end)`` per batch for `_predict_from_iter`.

    Thin batching wrapper over `_weights_for_batch` (which defines the
    weights). The generator form is what `_predict_from_iter` consumes.

    Parameters
    ----------
    X : ndarray, shape (n_test, n_features), dtype float
        Test covariates on the original (un-transformed) scale.

    Yields
    ------
    w : ndarray, shape (batch, n_nodes), dtype float
        Per-row weights for the current batch.
    start : int
        Index of the first row in the batch.
    end : int
        One past the index of the last row in the batch.
    """
    X = np.asarray(X)
    n_test = X.shape[0]
    for start in range(0, n_test, self.batch_size):
      end = min(start + self.batch_size, n_test)
      yield self._weights_for_batch(X[start:end]), start, end

  def _predict_from_iter(self, X, iter_weights):
    """Combines batched weights with the response nodes into predictions.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype float
        Test covariates (already validated / expanded).
    iter_weights : Callable
        Generator factory yielding ``(weights, start, end)``
        triples. Usually `_iter_weights`; taking it as an argument
        keeps the prediction step reusable with any other weight
        source that follows the same batching contract.

    Returns
    -------
    ndarray, shape (n_samples,) or (n_samples, n_outputs), dtype float
        Predictions. Shape ``(n_samples,)`` when a single output is
        produced (mean or single quantile), ``(n_samples, n_outputs)``
        otherwise. Column order is mean first (if enabled), then
        quantiles in the order requested.
    """
    n_test = X.shape[0]
    quantiles = self.quantiles_
    n_outputs = (1 if self.mean else 0) + (
      len(quantiles) if quantiles is not None else 0
    )
    y_pred = np.empty((n_test, n_outputs))

    for w, start, end in iter_weights(X):
      col = 0
      if self.mean:
        y_pred[start:end, col] = w @ self._y_nodes
        col += 1
      if quantiles is not None:
        batch_preds = [
          np.quantile(
            a=self._y_nodes,
            q=quantiles,
            weights=row_w,
            method="inverted_cdf",
          )
          for row_w in w
        ]
        y_pred[start:end, col : col + len(quantiles)] = np.vstack(batch_preds)

    # Drop the output axis for a single output, never the sample axis: a bare
    # `squeeze()` turns a one-row prediction into a scalar.
    return y_pred[:, 0] if y_pred.shape[1] == 1 else y_pred

  def predict(self, X):
    """Predicts the conditional mean and/or quantiles of ``Y`` given ``X``.

    Computes weights :math:`w_k(x)` from the fitted copula
    (`_iter_weights`) and returns the weighted statistics over the
    response nodes :math:`y_k`:
    :math:`\\hat{\\mathbb{E}}[Y \\mid X = x] = \\sum_k w_k(x)\\, y_k`
    for the mean (closed-form solution of the estimating equation
    :math:`\\int (y - \\beta) \\hat f(y \\mid x)\\, dy = 0`) and the
    weighted quantile via :func:`numpy.quantile` with
    ``method="inverted_cdf"`` for each requested level.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype float, or DataFrame
        Test covariates. Must match the training schema.

    Returns
    -------
    ndarray, shape (n_samples,) or (n_samples, n_outputs), dtype float
        Predictions. Shape ``(n_samples,)`` if only one output is
        requested (mean or a single quantile), otherwise
        ``(n_samples, n_outputs)``. Output columns are ordered:
        mean (if `self.mean`), then quantiles in `self.quantiles`
        order.
    """
    check_is_fitted(self, attributes=["_vine"])
    X = self._validate_input(X, reset=False)
    return self._predict_from_iter(X, self._iter_weights)

  def score(self, X, y, sample_weight=None) -> float:
    """Return :math:`R^2` for the fitted conditional mean.

    Quantile columns are supplementary prediction outputs, not independent
    response targets. When the conditional mean is disabled there is no
    scalar regressor score, so callers must choose a quantile-specific metric.

    Parameters
    ----------
    X : array-like or pandas.DataFrame
        Covariates to predict.
    y : array-like
        Observed response values.
    sample_weight : array-like, optional
        Per-observation weights for the coefficient of determination.

    Returns
    -------
    float
        Coefficient of determination for the conditional mean.
    """
    if not self.mean:
      raise ValueError(
        "VineRegressor.score is defined for mean predictions only; "
        "use a quantile-specific metric when mean=False."
      )
    prediction = np.asarray(self.predict(X))
    if prediction.ndim == 2:
      prediction = prediction[:, 0]
    return float(r2_score(y, prediction, sample_weight=sample_weight))


VineRegressor.__doc__ = f"""Vine-copula based regressor (mean and quantile).

A scikit-learn-compatible non-parametric regressor that predicts the
conditional mean :math:`\\mathbb{{E}}[Y \\mid X = x]` and/or
conditional :math:`\\tau`-quantiles of ``Y`` given covariates ``X``
by fitting a vine copula to the joint distribution of :math:`(Y, X)`
and reducing prediction to a weighted statistic of the training
responses.

For a target functional :math:`\\beta(x)` characterized by
:math:`\\mathbb{{E}}[\\psi_\\beta(Y) \\mid X = x] = 0` (Nagler &
Vatter, 2024), the fitted conditional density
:math:`\\hat f_{{Y \\mid X}}(\\cdot \\mid x)` solves

.. math::

   \\int \\psi_\\beta(y)\\, \\hat f_{{Y \\mid X}}(y \\mid x)\\, dy
   = 0.

Setting :math:`\\psi_\\beta(y) = y - \\beta` recovers the
conditional mean (a closed-form weighted average of the training
responses); setting :math:`\\psi_\\beta(y) =
\\mathbf{{1}}\\{{y < \\beta\\}} - \\tau` recovers the conditional
:math:`\\tau`-quantile (a weighted quantile via
:func:`numpy.quantile` with ``method="inverted_cdf"``).

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
>>> est.predict(X[:5])

{_DOC_REFERENCES}
"""
