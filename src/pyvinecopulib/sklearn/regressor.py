import numpy as np
import pyvinecopulib as pv
from sklearn.base import RegressorMixin

from ._base import VineBase


class VineRegressor(VineBase, RegressorMixin):
  def __init__(
    self,
    mean: bool = True,
    quantiles: list[float] | np.ndarray | None = None,
    controls: object | None = None,
    structure: object | None = None,
    batch_size: int = 100,
    use_grid: bool = True,
    normalize_weights: bool = True,
    schema: dict[str, list[str]] | None = None,
  ) -> None:
    """
    Sklearn-compatible vine-copula based regressor for conditional expectations/quantiles.

    Parameters
    ----------
    mean : bool, default=True
        Whether to enable mean prediction.
    quantiles : list[float] | None, default=None
        List of quantiles to enable quantile prediction. If None, quantile prediction is disabled.
    controls : pv.FitControlsVinecop, optional Controls for vinecop fitting.
        If None, defaults to tll family with 1 thread.
    structure : pv.RVineStructure, optional Vine structure.
    batch_size : int, default=100
        Number of test points to process per batch when predicting.
        - 1 = "loop" mode (minimal memory, slowest)
        - n_test = "stack" mode (maximal memory, fastest)
        - in between = trade-off
    use_grid : bool, default=False
        Whether to use the grid points of the target for predictions
        instead of the raw values.
    normalize_weights : bool, default=True
        Whether to normalize the conditional weights to sum to 1.
    schema : dict | None, default=None
        Pre-specified metadata about the input. If None, it will be inferred from the training data. Currently, only _kde1d_types is used.
        Example: {"kde1d_types": ["continuous", "discrete", "continuous", ...]}
        Supported types are "continuous" and "discrete" (for ordered variables).
    """
    # Validate types for pyvinecopulib compatibility
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

    # Regressor-specific parameters
    self.mean = mean
    if quantiles is None:
      self.quantiles = None
    else:
      self.quantiles = np.atleast_1d(quantiles)
    if (not self.mean) and (self.quantiles is None):
      raise ValueError("At least one of mean or quantiles must be enabled.")
    self.use_grid = use_grid
    self.normalize_weights = normalize_weights

  def fit(self, X: np.ndarray, y: np.ndarray) -> "VineRegressor":
    """Fit the regressor, a vine copula to joint distribution of (X, y)."""

    # Check and possibly expand inputs
    X, y = self._check_and_expand_fit(X, y)

    # Fit marginal distributions
    self._fit_marginals(X, y)

    # Convert to pseudo-observations
    uy_train = self._to_u_scale(y, is_y=True)
    ux = self._to_u_scale(X)

    # Fit vine copula to (U_Y, U_X)
    assert self.schema is not None  # Guaranteed after _check_and_expand_fit
    var_types = ["c"] + [x[0] for x in self.schema["kde1d_types"]]
    self._fit_vine(np.column_stack([uy_train, ux]), var_types=var_types)

    # Store training data for predictions
    if not self.use_grid:
      self._y_train = y
      self._uy_train = uy_train
    else:
      self._y_train = self._y_kde1d.grid_points
      uy_cdf = self._y_kde1d.cdf(self._y_train)
      self._uy_train = np.asarray(uy_cdf).reshape(-1, 1)
      self._y_density = self._y_kde1d.values[np.newaxis, :]
    return self

  def pdf(self, X: np.ndarray, y: np.ndarray, log: bool = False) -> np.ndarray:
    return self._pdf_samples(X, y=y, log=log, copula_only=True)

  def _copula_marginal_density(
    self, X: np.ndarray, log: bool = False, n_grid: int = 101
  ) -> np.ndarray:
    """
    Compute c_X(u_X) = ∫ c_{Y,X}(u_Y, u_X) du_Y
    using Simpson's rule on a uniform grid in [0, 1].

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Conditioning covariates.
    log : bool, default=False
        Whether to return log-density.
    n_grid : int, default=31
        Number of grid points for integration over u_Y.
        Simpson's rule requires an odd number; if even is given,
        it is increased by 1.

    Returns
    -------
    density : ndarray of shape (n_samples,)
        Copula marginal density values c_X(u_X).
    """
    if not hasattr(self, "_vine"):
      raise RuntimeError("Model not fitted yet.")

    X = np.asarray(X)
    ux = self._to_u_scale(X)
    n_test = ux.shape[0]

    # Ensure odd number of grid points
    if n_grid % 2 == 0:
      n_grid += 1

    eps = 1e-10
    uy_nodes = np.linspace(eps, 1 - eps, n_grid).reshape(-1, 1)

    # Simpson weights: 1,4,2,4,...,4,1
    w = np.ones(n_grid)
    w[1:-1:2] = 4
    w[2:-1:2] = 2
    h = 1.0 / (n_grid - 1)  # grid spacing
    simpson_factor = h / 3.0

    out = np.empty(n_test)

    for start in range(0, n_test, self.batch_size):
      end = min(start + self.batch_size, n_test)
      ux_batch = np.repeat(ux[start:end], n_grid, axis=0)
      uy_rep = np.tile(uy_nodes, (end - start, 1))
      u = np.column_stack([uy_rep, ux_batch])

      vals = self._vine.pdf(u, num_threads=self.controls.num_threads)
      vals = np.asarray(vals).reshape(end - start, n_grid)

      # Simpson's integration
      c_x = simpson_factor * (vals * w[None, :]).sum(axis=1)
      out[start:end] = c_x

    return np.log(np.clip(out, eps, None)) if log else out

  # Generator form (`yield (w, start, end)` per batch) is preserved so a
  # future forest PR can interleave batches across trees in `_predict_from_iter`.
  def _iter_weights(self, X):
    """
    Compute conditional weights for each test point:
        w_i(x) ∝ c_{X,Y}(U_x, U_{y_i}), for i in {1, ..., n_train}, normalized to sum to 1.
    Shapes:
        X: (n_test, d)
        returns: (n_test, n_train)
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
      if self.normalize_weights:
        w /= np.sum(w, axis=1, keepdims=True)

      yield w, start, end

  def _predict_from_iter(self, X, iter_weights):
    """
    Core prediction logic.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Test covariates.
    iter_weights : callable
        Function taking X and yielding (weights, start, end).

    Returns
    -------
    y_pred : ndarray of shape (n_samples, n_outputs)
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
    """Predict conditional mean E[Y|X] and/or quantiles."""
    X = self._check_and_expand_predict(X)
    return self._predict_from_iter(X, self._iter_weights)
