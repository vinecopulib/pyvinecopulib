import numpy as np
import pandas as pd
import pyvinecopulib as pv
from sklearn.base import BaseEstimator


def expand_factors(df: pd.DataFrame) -> pd.DataFrame:
  """
  - Keep numeric and ordered categoricals unchanged.
  - Expand unordered categoricals into 0/1 dummies, dropping the first level.
  - Dummies are cast to ordered categorical with levels {0,1}.
  """
  out_parts: list[pd.Series | pd.DataFrame] = []

  for colname, x in df.items():
    # Numeric: leave as-is
    if pd.api.types.is_numeric_dtype(x):
      out_parts.append(x)

    # Categorical ordered: leave as-is
    elif isinstance(x.dtype, pd.CategoricalDtype) and x.dtype.ordered:
      out_parts.append(x)

    # Categorical unordered: expand
    elif isinstance(x.dtype, pd.CategoricalDtype) and not x.dtype.ordered:
      dummies = pd.get_dummies(x, drop_first=True)
      # cast each dummy to ordered categorical {0,1}
      dummies_int = dummies.astype("int")
      # Convert each column individually to ordered categorical
      for col in dummies_int.columns:
        dummies_int[col] = pd.Categorical(
          dummies_int[col], categories=[0, 1], ordered=True
        )
      # rename columns to match R's levels[-1]
      new_names = [f"{colname}_{cat}" for cat in x.cat.categories[1:]]
      dummies_int.columns = new_names
      out_parts.append(dummies_int)

    else:
      raise ValueError(f"Unsupported dtype for column {colname}: {x.dtype}")

  return pd.concat(out_parts, axis=1)


class VineBase(BaseEstimator):
  """
  Base class for vine-copula based estimators.

  Provides common functionality for:
  - Marginal distribution fitting (Kde1d)
  - Data preprocessing and validation
  - Pseudo-observation transformation
  - Vine copula fitting
  - Batched operations
  """

  def __init__(
    self,
    controls: pv.FitControlsVinecop | None = None,
    structure: pv.RVineStructure | None = None,
    batch_size: int = 100,
    schema: dict[str, list[str]] | None = None,
  ) -> None:
    """
    Base vine copula estimator.

    Parameters
    ----------
    controls : pv.FitControlsVinecop, optional
        Controls for vinecop fitting. If None, defaults to tll family with 1 thread.
    structure : pv.RVineStructure, optional
        Vine structure. If None, structure will be selected automatically.
    batch_size : int, default=100
        Number of test points to process per batch when making predictions.
        - 1 = "loop" mode (minimal memory, slowest)
        - n_test = "stack" mode (maximal memory, fastest)
        - in between = trade-off
    schema : dict | None, default=None
        Pre-specified metadata about the input. If None, it will be inferred from the training data.
        Currently, only _kde1d_types is used.
        Example: {"kde1d_types": ["continuous", "discrete", "continuous", ...]}
        Supported types are "continuous" and "discrete" (for ordered variables).
    """
    if controls is None:
      controls = pv.FitControlsVinecop(
        family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
      )
    self.controls = controls
    self.structure = structure
    self.batch_size = batch_size
    self.schema = schema

  def _check_and_expand_fit(
    self, X: np.ndarray | pd.DataFrame, y: np.ndarray | None = None
  ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """
    Check and expand input data for fitting.

    Parameters
    ----------
    X : array-like
        Input features.
    y : array-like, optional
        Target values. If None, only X is processed (for density estimation).

    Returns
    -------
    X : ndarray
        Processed and expanded X.
    y : ndarray or None
        Processed y if provided, None otherwise.
    """
    # Check inputs
    if not isinstance(X, np.ndarray) and not isinstance(X, pd.DataFrame):
      raise ValueError("X must be a numpy array or pandas DataFrame")

    if y is not None:
      if not isinstance(y, np.ndarray):
        raise ValueError("y must be a numpy array")
      if not X.shape[0] == y.shape[0]:
        raise ValueError("X and y must have the same number of samples.")

    # Discrete variables handling when the input is a DataFrame
    if isinstance(X, pd.DataFrame):
      if self.schema is not None:
        raise ValueError("When schema is already set, X must be a numpy array.")
      self._used_columns = list(X.columns)
      self._dtypes = X.dtypes.to_dict()
      self._categories = {
        c: X[c].cat.categories
        for c in X.select_dtypes(include="category").columns
      }
      X = expand_factors(X)
      self._expanded_columns = list(X.columns)
      kde1d_types = [
        "discrete" if isinstance(dtype, pd.CategoricalDtype) else "continuous"
        for dtype in X.dtypes
      ]
      X = X.to_numpy()
    elif isinstance(X, np.ndarray):
      kde1d_types = (
        ["continuous"] * X.shape[1]
        if (self.schema is None) or ("kde1d_types" not in self.schema)
        else self.schema["kde1d_types"]
      )
      if len(kde1d_types) != X.shape[1]:
        raise ValueError(
          "schema['kde1d_types'] length does not match number of features in X."
        )

    # Either way, we now know the number of features after expansion
    self.n_features_in_ = X.shape[1]

    if self.schema is None:
      self.schema = {
        "kde1d_types": kde1d_types,
      }

    return X if y is None else (X, y)

  def _check_and_expand_predict(
    self, X: np.ndarray | pd.DataFrame
  ) -> np.ndarray:
    """
    Check and expand input data for prediction.

    Parameters
    ----------
    X : array-like
        Input features.

    Returns
    -------
    X : ndarray
        Processed and expanded X.
    """
    if isinstance(X, np.ndarray):
      # Numeric array, just check shape
      if X.shape[1] != self.n_features_in_:
        raise ValueError("X has wrong number of features.")
      return X

    elif isinstance(X, pd.DataFrame):
      # Check original columns
      if list(X.columns) != self._used_columns:
        raise ValueError("Column names/order do not match training data.")

      # Enforce dtypes/categories
      for col in self._used_columns:
        dtype_expected = self._dtypes[col]
        if isinstance(dtype_expected, pd.CategoricalDtype):
          if not isinstance(X[col].dtype, pd.CategoricalDtype):
            raise ValueError(f"Column {col} must be categorical.")
          # align categories
          X[col] = X[col].cat.set_categories(
            dtype_expected.categories, ordered=dtype_expected.ordered
          )
        # numeric check is lighter (float/int are ok)

      # Expand factors
      X_exp = expand_factors(X)

      # Reorder to match training
      X_exp = X_exp[self._expanded_columns]

      return X_exp.to_numpy()

    else:
      raise ValueError("X must be a numpy array or pandas DataFrame")

  def _fit_marginals(
    self, X: np.ndarray, y: np.ndarray | None = None
  ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """
    Fit marginal distributions (Kde1d) for features and optionally target.

    Parameters
    ----------
    X : ndarray
        Input features.
    y : ndarray, optional
        Target values. If None, only X marginals are fitted.

    Returns
    -------
    X : ndarray
        Input features (unchanged).
    y : ndarray or None
        Target values (unchanged) if provided, None otherwise.
    """
    # Fit 1d KDEs for X marginals
    self._x_kde1d = []
    assert self.schema is not None  # Guaranteed after _check_and_expand_fit
    for j in range(self.n_features_in_):
      kde = pv.utils.Kde1d(type=self.schema["kde1d_types"][j])
      kde.fit(X[:, j])
      self._x_kde1d.append(kde)

    # Fit 1d KDE for y marginal if provided
    if y is not None:
      self._y_kde1d = pv.utils.Kde1d()
      self._y_kde1d.fit(y)
      return X, y
    else:
      return X

  def _to_u_scale(self, Z: np.ndarray, is_y: bool = False) -> np.ndarray:
    """
    Compute marginal pseudo-observations for new data, handling discrete vars.

    Parameters
    ----------
    Z : array-like
        If is_y=False: shape (n_samples, d).
        If is_y=True : shape (n_samples,).
    is_y : bool, default=False
        If True, transform response with y_kde1d.
        If False, transform covariates with x_kde1d list.

    Returns
    -------
    U : ndarray
        If is_y=False: shape (n_samples, d [+ optional discrete sub-columns]).
        If is_y=True : shape (n_samples, 1).
    """
    if not hasattr(self, "_x_kde1d"):
      raise RuntimeError("Model not fitted yet.")

    if is_y and not hasattr(self, "_y_kde1d"):
      raise RuntimeError(
        "Target marginal not fitted (y was not provided during fit)."
      )

    Z = np.asarray(Z)

    if is_y:
      # For the response we always treat it continuous -> plain PIT
      result = self._y_kde1d.cdf(Z.squeeze(), check_fitted=False)
      return np.asarray(result).reshape(-1, 1)

    n_samples, n_features = Z.shape
    u_cols = []
    u_sub_cols = []

    for j in range(n_features):
      kde = self._x_kde1d[j]
      zj = Z[:, j]

      # always get the main PIT
      u = kde.cdf(zj, check_fitted=False)
      u_cols.append(u)

      if kde.type == "discrete":
        # compute sub-CDF: F(x^-)
        u_sub = kde.cdf(zj - 1, check_fitted=False)
        u_sub_cols.append(u_sub)

    # Clamp to (0,1) to avoid numerical issues
    eps = 1e-10
    u_cols = [np.clip(u, eps, 1 - eps) for u in u_cols]
    if u_sub_cols:
      u_sub_cols = [np.clip(u, eps, 1 - eps) for u in u_sub_cols]
      return np.column_stack([*u_cols, *u_sub_cols])
    else:
      return np.column_stack(u_cols)

  def _fit_vine(
    self, U: np.ndarray, var_types: list[str] | None = None
  ) -> "VineBase":
    """
    Fit vine copula to pseudo-observations.

    Parameters
    ----------
    U : ndarray
        Pseudo-observations in [0,1]^d.
    var_types : list, optional
        Variable types for vine fitting. If None, inferred from schema.

    Returns
    -------
    self
    """
    if var_types is None:
      assert self.schema is not None  # Guaranteed after _check_and_expand_fit
      var_types = [x[0] for x in self.schema["kde1d_types"]]

    self._vine = pv.Vinecop.from_data(
      data=U,
      structure=self.structure,
      var_types=var_types,
      controls=self.controls,
    )
    self.structure = self._vine.structure
    return self

  # `copula_only` is unused by the standalone estimators; it's kept so a
  # future forest PR can combine pure copula densities across trees.
  def _pdf_samples(
    self,
    X: np.ndarray | pd.DataFrame,
    y: np.ndarray | None = None,
    log: bool = False,
    copula_only: bool = False,
  ) -> np.ndarray:
    """
    Unified method to compute density/likelihood of samples.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Input features.
    y : array-like of shape (n_samples,), optional
        Target values. If provided, computes joint density of (X, y).
        If None, computes density of X only (for density estimation).
    log : bool, default=False
        Whether to return log-density.
    copula_only : bool, default=False
        If True, returns only the copula density component.
        If False, returns full density (copula + marginals).

    Returns
    -------
    density : ndarray of shape (n_samples,)
        Density or log-density values.
    """
    if not hasattr(self, "_vine"):
      raise RuntimeError("Model not fitted yet.")

    # Check inputs
    X = self._check_and_expand_predict(X)
    if y is not None:
      y = np.asarray(y).reshape(-1, 1)
      if X.shape[0] != y.shape[0]:
        raise ValueError("X and y must have the same number of samples.")

    # Convert to pseudo-observations
    if y is not None:
      uy = self._to_u_scale(y, is_y=True)
      ux = self._to_u_scale(X)
      U = np.column_stack([uy, ux])
    else:
      U = self._to_u_scale(X)

    # Compute copula density
    pdf_vals = np.asarray(
      self._vine.pdf(U, num_threads=self.controls.num_threads)
    )
    log_c = np.asarray(np.log(pdf_vals))

    # Return only copula part if requested
    if copula_only:
      return np.asarray(log_c if log else np.exp(log_c))

    # Add marginal densities for full joint density
    # Clamp to >1e-10 to avoid log(0)
    eps = 1e-10
    f_x = np.array(
      [
        np.clip(self._x_kde1d[j].pdf(X[:, j]), eps, None)
        for j in range(self.n_features_in_)
      ]
    )
    log_f_x = np.sum(
      np.log(f_x),
      axis=0,
    )
    log_density = log_c + log_f_x
    if y is not None:
      log_f_y = np.log(self._y_kde1d.pdf(y.squeeze()))
      log_density += log_f_y

    return np.asarray(log_density if log else np.exp(log_density))
