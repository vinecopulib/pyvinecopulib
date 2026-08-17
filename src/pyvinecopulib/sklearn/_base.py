from numbers import Integral

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.utils._param_validation import Interval
from sklearn.utils.validation import check_is_fitted, check_random_state

import pyvinecopulib as pv

from .backends import resolve_backend

# Shared docstring fragments interpolated into VineDensity / VineRegressor
# class docstrings via f-strings. Defined once here, used by both subclasses
# so changes propagate without copy-paste.

_DOC_PIPELINE = r"""The estimator follows the standard pyvinecopulib
two-step pipeline: a univariate kernel density estimator
(``Kde1d``) is fit to each column, the marginal CDFs transform the data
to pseudo-observations
:math:`U_j = \hat F_j(X_j) \in [0, 1]`, and a vine copula is fit on
the pseudo-observations. For discrete columns the left limit
:math:`\hat F_j(X_j^-)` is also stacked so the vine sees a continuous
proxy. Unordered categoricals are first expanded to ordered
``{0, 1}`` dummies via `expand_factors`.

Fit-time configuration is bundled in a backend object passed via
``backend=``. The default ``VinecopBackend`` wraps ``Vinecop`` and has
no extra dependencies; ``TorchVinecopBackend`` routes the same
pipeline through the PyTorch evaluator (GPU / autograd). See the
:doc:`concepts page </concepts>` for the underlying vine-copula
construction.
"""

_DOC_FACTORIZATION = r"""By Sklar's theorem the joint density
factorizes as

.. math::

   f(\mathbf{x})
   = c\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr)\,
     \prod_{j=1}^{d} f_j(x_j),

with :math:`c` further decomposed into pair copulas indexed by a
vine structure (Bedford & Cooke, 2002; Aas et al., 2009). Passing
``copula_only=True`` to `pdf` returns the copula factor
:math:`c(\mathbf{u})` alone.
"""

_DOC_DISCRETE = r"""Discrete (or expanded unordered-categorical)
columns are handled via ``Kde1d``'s ``type="discrete"`` mode:
pseudo-observations stack :math:`\hat F_j(X_j)` and
:math:`\hat F_j(X_j^-)` so the vine evaluation sees the appropriate
continuous proxy. Handled transparently by `fit` and `pdf`.
"""

_DOC_REFERENCES = r"""References
----------
.. [1] Bedford, T. and Cooke, R. M. (2002).
       *Vines--a new graphical model for dependent random variables.*
       The Annals of Statistics, 30(4), 1031--1068.
.. [2] Aas, K., Czado, C., Frigessi, A. and Bakken, H. (2009).
       *Pair-copula constructions of multiple dependence.*
       Insurance: Mathematics and Economics, 44(2), 182--198.
.. [3] Nagler, T. and Vatter, T. (2024).
       *Solving Estimating Equations With Copulas.*
       Journal of the American Statistical Association, 119(546),
       1168--1180.
"""


def expand_factors(df: pd.DataFrame) -> pd.DataFrame:
  """Expand unordered categoricals to ordered ``{0, 1}`` dummies.

  Unordered categorical columns are not directly usable by the vine
  copula (which assumes orderable marginals). This helper one-hot
  encodes them, drops the first level (avoiding collinearity), and
  re-casts each dummy to an ordered categorical with levels
  ``[0, 1]`` so they pass through Kde1d's ``"discrete"`` mode
  unchanged. Numeric and already-ordered categorical columns are
  passed through.

  Parameters
  ----------
  df : pandas.DataFrame
      Input DataFrame. Columns may be numeric, ordered categorical, or
      unordered categorical. Any other dtype raises ``ValueError``.

  Returns
  -------
  pandas.DataFrame
      DataFrame with unordered categoricals expanded. Dummy columns are
      named ``"<original>_<category>"`` (matching R's
      ``model.matrix`` convention, dropping the first level), and
      typed as ordered categorical over ``[0, 1]``.
  """
  out_parts: list[pd.Series | pd.DataFrame] = []

  for colname, x in df.items():
    if pd.api.types.is_numeric_dtype(x):
      out_parts.append(x)
    elif isinstance(x.dtype, pd.CategoricalDtype) and x.dtype.ordered:
      out_parts.append(x)
    elif isinstance(x.dtype, pd.CategoricalDtype) and not x.dtype.ordered:
      dummies_int = pd.get_dummies(x, drop_first=True).astype("int")
      for col in dummies_int.columns:
        dummies_int[col] = pd.Categorical(
          dummies_int[col], categories=[0, 1], ordered=True
        )
      # Names match R's levels[-1].
      dummies_int.columns = [f"{colname}_{cat}" for cat in x.cat.categories[1:]]
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
  - Vine copula fitting (via a backend strategy object)
  - Batched operations

  Per the scikit-learn developer guide, ``__init__`` only stores
  parameters as-is; validation runs lazily in ``fit()`` via
  ``_validate_params()``. Fitted attributes use the trailing-underscore
  convention.
  """

  _parameter_constraints: dict = {
    "backend": [object, None],
    "batch_size": [Interval(Integral, 1, None, closed="left")],
    "random_state": ["random_state"],
  }

  def __init__(
    self,
    backend=None,
    batch_size: int = 100,
    random_state=None,
  ) -> None:
    """Base vine copula estimator.

    Parameters
    ----------
    backend : VinecopBackend or compatible, default=None
        Backend strategy that holds fit-time controls (a
        ``FitControlsVinecop`` for the default backend or a
        ``FitControlsTorchVinecop`` for the torch backend) and an
        optional structure. `None` resolves to a default
        ``VinecopBackend`` at fit time.
    batch_size : int, default=100
        Number of test points to process per batch when making
        predictions. ``1`` minimizes memory at the cost of speed;
        ``n_test`` is the opposite extreme; intermediate values
        trade off memory and throughput.
    random_state : int, RandomState instance or None, default=None
        Seeds the RNG used by stochastic operations (e.g.
        `simulate`, `cdf` quasi-MC, structure simulation). Stored
        as-is; resolved via `sklearn.utils.check_random_state`
        inside `fit`.
    """
    self.backend = backend
    self.batch_size = batch_size
    self.random_state = random_state

  def _validate_input(
    self,
    X,
    y=None,
    *,
    reset: bool,
  ):
    """Validate the ``X`` (and optional ``y``) input.

    For DataFrames, captures the canonical
    :attr:`feature_names_in_`, expands unordered categoricals via
    :func:`expand_factors`, and infers ``kde1d_types``. For ndarrays,
    captures :attr:`n_features_in_` and assumes all-continuous unless
    a previously-set ``schema_["kde1d_types"]`` says otherwise.

    Parameters
    ----------
    X : ndarray or DataFrame
    y : array-like, optional
    reset : bool
        If ``True`` (called from ``fit``), set fitted attributes; if
        ``False`` (called from ``predict``-style methods), validate
        against the previously-set ones.
    """
    if not isinstance(X, (np.ndarray, pd.DataFrame)):
      raise ValueError("X must be a numpy array or pandas DataFrame")
    if y is not None:
      y = np.asarray(y)
      if X.shape[0] != y.shape[0]:
        raise ValueError("X and y must have the same number of samples.")

    if isinstance(X, pd.DataFrame):
      if reset:
        self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        self._dtypes = X.dtypes.to_dict()
        self._categories = {
          c: X[c].cat.categories
          for c in X.select_dtypes(include="category").columns
        }
        X_exp = expand_factors(X)
        self._expanded_columns = list(X_exp.columns)
        kde1d_types = [
          "discrete" if isinstance(dtype, pd.CategoricalDtype) else "continuous"
          for dtype in X_exp.dtypes
        ]
        self.schema_ = {"kde1d_types": kde1d_types}
        self.n_features_in_ = X_exp.shape[1]
        X_arr = X_exp.to_numpy()
      else:
        check_is_fitted(self, attributes=["feature_names_in_"])
        if list(X.columns) != list(self.feature_names_in_):
          raise ValueError("Column names/order do not match training data.")
        for col in self.feature_names_in_:
          dtype_expected = self._dtypes[col]
          if isinstance(dtype_expected, pd.CategoricalDtype):
            if not isinstance(X[col].dtype, pd.CategoricalDtype):
              raise ValueError(f"Column {col} must be categorical.")
            X[col] = X[col].cat.set_categories(
              dtype_expected.categories, ordered=dtype_expected.ordered
            )
        X_arr = expand_factors(X)[self._expanded_columns].to_numpy()
    else:
      if reset:
        # ndarray input carries no dtype information: respect a
        # ``schema_`` set by the caller (or by a subclass) before
        # ``fit`` to declare per-column kde1d types, and otherwise
        # treat every column as continuous.
        existing = getattr(self, "schema_", None)
        if existing is not None and "kde1d_types" in existing:
          kde1d_types = existing["kde1d_types"]
          if len(kde1d_types) != X.shape[1]:
            raise ValueError(
              "schema_['kde1d_types'] length does not match number of "
              "features in X."
            )
        else:
          kde1d_types = ["continuous"] * X.shape[1]
        self.schema_ = {"kde1d_types": kde1d_types}
        self.n_features_in_ = X.shape[1]
      else:
        check_is_fitted(self, attributes=["n_features_in_"])
        if X.shape[1] != self.n_features_in_:
          raise ValueError("X has wrong number of features.")
      X_arr = X

    if y is not None:
      return X_arr, y
    return X_arr

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
    self._x_kde1d = []
    assert (
      self.schema_ is not None
    )  # Guaranteed after _validate_input(reset=True)
    for j in range(self.n_features_in_):
      kde = pv.utils.Kde1d(type=self.schema_["kde1d_types"][j])
      kde.fit(X[:, j])
      self._x_kde1d.append(kde)

    if y is not None:
      self._y_kde1d = pv.utils.Kde1d()
      self._y_kde1d.fit(y)
      return X, y
    else:
      return X

  def _resolve_runtime_state(self) -> None:
    """Resolve the random-state and backend at fit time. Sets
    ``self.random_state_`` and ``self.backend_`` so subclasses can reuse
    them throughout ``fit`` and post-fit methods.
    """
    self.random_state_ = check_random_state(self.random_state)
    self.backend_ = resolve_backend(self.backend)

  def _draw_seeds(self, size: int = 5) -> list[int]:
    """Derive a list of ints suitable for C++ ``seeds=[...]`` kwargs
    from the resolved RNG. Reproducible iff ``random_state_`` is."""
    return [int(x) for x in self.random_state_.randint(0, 2**31 - 1, size=size)]

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
    check_is_fitted(self, attributes=["_x_kde1d"])
    if is_y and not hasattr(self, "_y_kde1d"):
      raise RuntimeError(
        "Target marginal not fitted (y was not provided during fit)."
      )

    Z = np.asarray(Z)

    if is_y:
      # Response is always treated continuous -> plain PIT.
      result = self._y_kde1d.cdf(Z.squeeze(), check_fitted=False)
      return np.asarray(result).reshape(-1, 1)

    n_samples, n_features = Z.shape
    u_cols = []
    u_sub_cols = []

    for j in range(n_features):
      kde = self._x_kde1d[j]
      zj = Z[:, j]
      u_cols.append(kde.cdf(zj, check_fitted=False))
      if kde.type == "discrete":
        # Sub-CDF F(x^-) for the discrete component.
        u_sub_cols.append(kde.cdf(zj - 1, check_fitted=False))

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
      assert (
        self.schema_ is not None
      )  # Guaranteed after _validate_input(reset=True)
      var_types = [x[0] for x in self.schema_["kde1d_types"]]

    backend = self.backend_
    self._vine = backend.fit_vine(U, var_types=var_types)
    self.structure_ = backend.structure_of(self._vine)
    return self

  # `copula_only=True` skips the marginal-density product and returns
  # the copula factor c(u) alone; `VineDensity.pdf` exposes it directly.
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
    check_is_fitted(self, attributes=["_vine"])

    X = self._validate_input(X, reset=False)
    if y is not None:
      y = np.asarray(y).reshape(-1, 1)
      if X.shape[0] != y.shape[0]:
        raise ValueError("X and y must have the same number of samples.")
      U = np.column_stack([self._to_u_scale(y, is_y=True), self._to_u_scale(X)])
    else:
      U = self._to_u_scale(X)

    pdf_vals = np.asarray(self.backend_.pdf(self._vine, U))
    log_c = np.asarray(np.log(pdf_vals))

    if copula_only:
      return np.asarray(log_c if log else np.exp(log_c))

    # Marginal densities; clamp >0 to avoid log(0).
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
