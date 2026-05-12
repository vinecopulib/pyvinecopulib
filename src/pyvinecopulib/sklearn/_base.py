import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator

import pyvinecopulib as pv

# Shared docstring fragments interpolated into VineDensity / VineRegressor
# class docstrings via f-strings. Defined once here, used by both subclasses
# so changes propagate without copy-paste.

_DOC_WRAPPER = r"""**Relation to the core API.** This estimator is a
thin scikit-learn wrapper around the lower-level pyvinecopulib
machinery. Marginals are estimated with
:class:`pyvinecopulib.utils.Kde1d`; the joint copula is fit by
:meth:`pyvinecopulib.core.Vinecop.from_data`, accepting an optional
:class:`pyvinecopulib.core.RVineStructure` (``structure`` argument) and
a :class:`pyvinecopulib.core.FitControlsVinecop` (``controls`` argument)
for the family set, threading, structure-selection algorithm, etc.
The sklearn class on top adds the standard ``fit`` / ``predict``-style
interface, ``DataFrame`` input handling with auto-inferred
``continuous`` / ``discrete`` schema, and batched evaluation; reach
for the underlying classes whenever you need control beyond the
sklearn convenience layer.
"""

_DOC_PIPELINE = r"""**Estimation pipeline.** The estimator follows a
three-step pipeline shared with :class:`VineDensity` /
:class:`VineRegressor`:

1. **Marginals.** A univariate kernel density estimator
   (:class:`pyvinecopulib.utils.Kde1d`) is fitted to each feature
   column. Continuous and ordered-discrete dtypes are inferred from the
   input (or read from a user-supplied ``schema``); unordered
   categoricals are first expanded into ordered ``{0, 1}`` dummies by
   :func:`pyvinecopulib.sklearn._base.expand_factors`.

2. **Pseudo-observations.** Each marginal CDF is applied to its column
   to produce pseudo-observations :math:`U_j = \hat F_j(X_j) \in [0, 1]`.
   For discrete columns the left limit :math:`\hat F_j(X_j^-)` is also
   computed and stacked, so the vine sees a continuous proxy for the
   discrete margin.

3. **Vine copula.** A vine copula is fitted to the pseudo-observations
   with :meth:`pyvinecopulib.core.Vinecop.from_data`, using the
   ``structure``, ``controls``, and (where relevant) variable-type
   tags. The default ``controls`` use the nonparametric ``tll`` pair
   family for every edge.
"""

_DOC_FACTORIZATION = r"""**Joint-density factorization.** Both
estimators rely on Sklar's theorem, which expresses the joint density
as a product of marginals and a copula density,

.. math::

   f(\mathbf{x})
   = c\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr)\,
     \prod_{j=1}^{d} f_j(x_j),

and on the pair-copula construction of Bedford & Cooke (2002) and
Aas et al. (2009), which writes :math:`c` as a product of bivariate
(pair-)copulas indexed by a vine structure. The ``copula_only=True`` flag on
:meth:`pdf` returns the copula factor :math:`c(\mathbf{u})` alone,
without the marginal product.
"""

_DOC_DISCRETE = r"""**Discrete variables.** Discrete (or expanded
unordered-categorical) columns are handled via the Kde1d
``type="discrete"`` mode: pseudo-observations stack
:math:`\hat F_j(X_j)` and :math:`\hat F_j(X_j^-)` so that the vine
copula evaluation sees the appropriate (continuous) proxy. This is
done transparently by :meth:`fit` and :meth:`pdf`.
"""

_DOC_FOREST = r"""**Forest construction.** The ensemble is built via
hold-out random search with model-confidence-set (MCS) survivor
selection:

1. **Validation split.** A fraction ``val_fraction`` of the training
   data is held out for survivor selection.
2. **Random search.** ``n_vines`` candidate vine structures are
   sampled (uniformly via Joe's algorithm when
   ``vines_sampling="uniform"``, or by a weighted random walk on the
   structure space when ``"local"``). Each candidate is fit on a
   bootstrap resample of the training portion.
3. **Survivor selection.** Each candidate's per-sample log-likelihood
   is evaluated on the validation set; the dual-split DA test of
   Hansen, Lunde & Nason (2011) — refined by Kim, Olsen, Nagler &
   Vatter (2025) — is applied to the resulting loss matrix to compute
   a model confidence set (marginal or uniform error control,
   controlled by ``method``). ``add_dissmann=True`` adds the
   Dissmann-structure baseline as an extra candidate.
4. **Prediction averaging.** Survivors are refit on the full training
   data; predictions are averaged across survivors at evaluation time.
   For the regressor, the conditional weights derived from
   :math:`c_{Y, X}` are averaged then row-normalised once at the
   ensemble level.
"""

_DOC_REFERENCES = r"""References
----------
- Bedford, T. and Cooke, R. M. (2002).
  *Vines--a new graphical model for dependent random variables.*
  The Annals of Statistics, 30(4), 1031--1068.
- Aas, K., Czado, C., Frigessi, A. and Bakken, H. (2009).
  *Pair-copula constructions of multiple dependence.*
  Insurance: Mathematics and Economics, 44(2), 182--198.
- Nagler, T. and Vatter, T. (2024).
  *Solving Estimating Equations With Copulas.*
  Journal of the American Statistical Association, 119(546), 1168--1180.
- Vatter, T. and Nagler, T. (2026).
  *Throwing Vines at the Wall: Structure Learning via Random Search.*
  arXiv preprint arXiv:2510.20035.
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
  - Vine copula fitting
  - Batched operations
  """

  def __init__(
    self,
    controls: pv.FitControlsVinecop | None = None,
    structure: pv.RVineStructure | None = None,
    batch_size: int = 100,
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

    Notes
    -----
    Advanced / ensemble use: ``self._schema`` (a dict with key
    ``"kde1d_types"``) can be assigned between construction and
    :meth:`fit` to override the auto-inferred marginal types
    (``"continuous"`` / ``"discrete"``). It is not exposed as an
    ``__init__`` parameter because casual users should never need it,
    and exposing it would clutter the sklearn-facing signature.
    Non-default values are not preserved by :func:`sklearn.base.clone`.
    """
    if controls is not None and not isinstance(controls, pv.FitControlsVinecop):
      raise TypeError(
        "controls must be pv.FitControlsVinecop or None, "
        f"got {type(controls).__name__}"
      )
    if structure is not None and not isinstance(structure, pv.RVineStructure):
      raise TypeError(
        "structure must be pv.RVineStructure or None, "
        f"got {type(structure).__name__}"
      )
    if not isinstance(batch_size, int) or isinstance(batch_size, bool):
      raise TypeError(
        f"batch_size must be int, got {type(batch_size).__name__}"
      )
    if batch_size < 1:
      raise ValueError(f"batch_size must be >= 1, got {batch_size}")

    if controls is None:
      controls = pv.FitControlsVinecop(
        family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
      )
    self.controls = controls
    self.structure = structure
    self.batch_size = batch_size
    self._schema: dict[str, list[str]] | None = None

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
    if not isinstance(X, np.ndarray) and not isinstance(X, pd.DataFrame):
      raise ValueError("X must be a numpy array or pandas DataFrame")

    if y is not None:
      if not isinstance(y, np.ndarray):
        raise ValueError("y must be a numpy array")
      if not X.shape[0] == y.shape[0]:
        raise ValueError("X and y must have the same number of samples.")

    if isinstance(X, pd.DataFrame):
      if self._schema is not None:
        raise ValueError(
          "When self._schema is already set, X must be a numpy array."
        )
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
        if (self._schema is None) or ("kde1d_types" not in self._schema)
        else self._schema["kde1d_types"]
      )
      if len(kde1d_types) != X.shape[1]:
        raise ValueError(
          "schema['kde1d_types'] length does not match number of features in X."
        )

    self.n_features_in_ = X.shape[1]

    if self._schema is None:
      self._schema = {
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
      if X.shape[1] != self.n_features_in_:
        raise ValueError("X has wrong number of features.")
      return X

    elif isinstance(X, pd.DataFrame):
      if list(X.columns) != self._used_columns:
        raise ValueError("Column names/order do not match training data.")

      for col in self._used_columns:
        dtype_expected = self._dtypes[col]
        if isinstance(dtype_expected, pd.CategoricalDtype):
          if not isinstance(X[col].dtype, pd.CategoricalDtype):
            raise ValueError(f"Column {col} must be categorical.")
          X[col] = X[col].cat.set_categories(
            dtype_expected.categories, ordered=dtype_expected.ordered
          )

      X_exp = expand_factors(X)[self._expanded_columns]
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
    self._x_kde1d = []
    assert self._schema is not None  # Guaranteed after _check_and_expand_fit
    for j in range(self.n_features_in_):
      kde = pv.utils.Kde1d(type=self._schema["kde1d_types"][j])
      kde.fit(X[:, j])
      self._x_kde1d.append(kde)

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
      assert self._schema is not None  # Guaranteed after _check_and_expand_fit
      var_types = [x[0] for x in self._schema["kde1d_types"]]

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

    X = self._check_and_expand_predict(X)
    if y is not None:
      y = np.asarray(y).reshape(-1, 1)
      if X.shape[0] != y.shape[0]:
        raise ValueError("X and y must have the same number of samples.")
      U = np.column_stack([self._to_u_scale(y, is_y=True), self._to_u_scale(X)])
    else:
      U = self._to_u_scale(X)

    pdf_vals = np.asarray(
      self._vine.pdf(U, num_threads=self.controls.num_threads)
    )
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
