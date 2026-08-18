import copy
import warnings
from numbers import Integral
from typing import Any, Optional

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.exceptions import DataConversionWarning
from sklearn.utils._param_validation import Interval
from sklearn.utils.validation import (
  assert_all_finite,
  check_is_fitted,
  check_random_state,
)

from ..core import Vinedist
from ..margins import Kde1dMargin, resolve_margins
from ..margins._resolve import fit_margin
from .backends import _BackendVinecop, resolve_backend

# Shared docstring fragments interpolated into VineDensity / VineRegressor
# class docstrings via f-strings. Defined once here, used by both subclasses
# so changes propagate without copy-paste.

_DOC_PIPELINE = r"""The estimator follows the standard pyvinecopulib
two-step pipeline: a univariate margin is fit to each column, the
marginal CDFs transform the data to pseudo-observations
:math:`U_j = \hat F_j(X_j) \in [0, 1]`, and a vine copula is fit on
the pseudo-observations. For columns with atoms the left limit
:math:`\hat F_j(X_j^-)` is also stacked so the vine sees a continuous
proxy. Unordered categoricals are first expanded to ordered
``{0, 1}`` dummies via `expand_factors`.

The two halves are configured separately and symmetrically. Margins
come from ``margins=``, which defaults to a boundary-corrected kernel
density (``Kde1dMargin``) per column and accepts anything
:func:`pyvinecopulib.margins.resolve_margins` understands --- an alias
such as ``"empirical"`` or ``"parametric"``, one margin broadcast to
every column, a per-column sequence, or a mapping keyed by feature
name. The copula comes from ``backend=``: the default
``VinecopBackend`` wraps ``Vinecop`` and has no extra dependencies,
while ``TorchVinecopBackend`` routes the same pipeline through the
PyTorch evaluator (GPU / autograd).

Fitting assembles a ``Vinedist`` --- the copula and its margins as one
object --- and every post-fit method evaluates through it. It is
published as ``distribution_``, so the fitted joint distribution is
usable outside the estimator; ``selection_report_`` carries the
per-candidate table of any margin that selected its own family. See
the :doc:`concepts page </concepts>` for the underlying vine-copula
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
columns are handled by the margin's own variable type --- for the
default ``Kde1dMargin`` that is ``type="discrete"``:
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
  - Marginal distribution fitting (via ``margins=``)
  - Data preprocessing and validation
  - Pseudo-observation transformation
  - Vine copula fitting (via a backend strategy object)
  - Assembling both halves into a fitted ``Vinedist``
  - Batched operations

  Per the scikit-learn developer guide, ``__init__`` only stores
  parameters as-is; validation runs lazily in ``fit()`` via
  ``_validate_params()``. Fitted attributes use the trailing-underscore
  convention.
  """

  # `schema_` is heterogeneous by design -- per-column variable types, and the
  # bounds of the columns whose support the input pinned down. Declared here so
  # the type checker knows that; `_validate_input(reset=True)` sets it, and
  # `getattr(self, "schema_", None)` is still how a caller-set one is read.
  schema_: dict[str, Any]

  _parameter_constraints: dict = {
    "backend": [object, None],
    "margins": [object, None],
    "batch_size": [Interval(Integral, 1, None, closed="left")],
    "random_state": ["random_state"],
  }

  def __init__(
    self,
    backend=None,
    margins=None,
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
    margins : object, default=None
        What to fit to each column, in any form
        :func:`pyvinecopulib.margins.resolve_margins` accepts: an alias
        (``"kde"``, ``"empirical"``, ``"parametric"``), one margin
        broadcast to every column, a sequence of length
        ``n_features_in_``, a mapping keyed by feature name or
        position, or a callable taking a column and returning a
        margin. `None` fits a ``Kde1dMargin`` per column, carrying the
        variable type inferred from the input and, for the ``{0, 1}``
        dummies of an expanded unordered categorical, the bounds of its
        support. Stored as-is and never mutated: every specification is
        fitted on a copy.
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
    self.margins = margins
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
      # A column vector is the shape callers reach for after slicing a
      # DataFrame; sklearn's convention is to ravel it and say so.
      if y.ndim == 2 and y.shape[1] == 1:
        warnings.warn(
          "A column-vector y was passed when a 1d array was expected. "
          "Please change the shape of y to (n_samples,), for example using "
          "ravel().",
          DataConversionWarning,
          stacklevel=2,
        )
        y = y.ravel()
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
        # `expand_factors` only ever adds columns, so an expanded name absent
        # from the input is a {0, 1} dummy: its support is known exactly, and
        # saying so keeps the marginal fit off values that cannot occur. The
        # bounds of any other discrete column are a modeling assumption the
        # caller has to make, so they stay unset.
        original = set(X.columns)
        bounds: list[Optional[tuple[float, float]]] = [
          None if name in original else (0.0, 1.0)
          for name in self._expanded_columns
        ]
        self.schema_ = {"kde1d_types": kde1d_types, "bounds": bounds}
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
        existing = getattr(self, "schema_", None) or {}
        if "kde1d_types" in existing:
          kde1d_types = existing["kde1d_types"]
          if len(kde1d_types) != X.shape[1]:
            raise ValueError(
              "schema_['kde1d_types'] length does not match number of "
              "features in X."
            )
        else:
          kde1d_types = ["continuous"] * X.shape[1]
        bounds = existing.get("bounds") or [None] * X.shape[1]
        if len(bounds) != X.shape[1]:
          raise ValueError(
            "schema_['bounds'] length does not match number of features in X."
          )
        self.schema_ = {"kde1d_types": kde1d_types, "bounds": bounds}
        self.n_features_in_ = X.shape[1]
      else:
        check_is_fitted(self, attributes=["n_features_in_"])
        if X.shape[1] != self.n_features_in_:
          raise ValueError("X has wrong number of features.")
      X_arr = X

    # Margins are fitted on float columns, so a complex input would otherwise
    # be cast and lose its imaginary part without saying so.
    if np.iscomplexobj(X_arr) or (y is not None and np.iscomplexobj(y)):
      raise ValueError("Complex data not supported")
    # The marginal estimator is C++ and reads NaN / inf as a segmentation
    # fault rather than an error, so nothing downstream can report this.
    assert_all_finite(X_arr, input_name="X")
    if y is not None:
      assert_all_finite(y, input_name="y")
      return X_arr, y
    return X_arr

  def _column_name(self, j: int) -> str:
    """Name of expanded column ``j``, for reports and messages.

    Parameters
    ----------
    j : int
        Column index.

    Returns
    -------
    str
        The feature name for DataFrame input, else a positional stand-in.
    """
    names = getattr(self, "_expanded_columns", None)
    if names is not None and j < len(names):
      return str(names[j])
    return f"x{j}"

  def _default_margin_specs(self) -> list[Kde1dMargin]:
    """One unfitted kernel-density margin per column, from the schema.

    This is what ``margins=None`` means, and what a column a ``margins=``
    mapping does not address falls back to.

    Returns
    -------
    list of Kde1dMargin
        One margin per feature, carrying the variable type and whatever bounds
        the input told us about.
    """
    types = self.schema_["kde1d_types"]
    bounds = self.schema_.get("bounds") or [None] * len(types)
    specs = []
    for type_, bound in zip(types, bounds):
      lo, hi = (
        (None, None) if bound is None else (float(bound[0]), float(bound[1]))
      )
      specs.append(Kde1dMargin(type=type_, xmin=lo, xmax=hi))
    return specs

  def _fit_one_margin(self, spec: Any, column: np.ndarray, name: str) -> Any:
    """Fit one column's margin.

    Parameters
    ----------
    spec : object
        A specification from :func:`pyvinecopulib.margins.resolve_margins`.
    column : ndarray, shape (n_samples,), dtype float
        The column to fit.
    name : str
        The variable's name, for a selector's report.

    Returns
    -------
    MarginLike
        The fitted margin.
    """
    # Fitting a margin mutates it, and `self.margins` is a constructor argument
    # that has to survive `fit` untouched so `clone` reproduces the estimator;
    # hence every specification is fitted on a copy of itself.
    spec = copy.deepcopy(spec)
    # A selector records which variable each candidate was fitted to; supply
    # the name when the caller has not.
    if hasattr(spec, "report_") and getattr(spec, "name", "") is None:
      spec.name = name
    return fit_margin(spec, np.asarray(column, dtype=float))

  def _response_margin_spec(self) -> Any:
    """The specification for the response margin.

    ``margins`` addresses the features, so a per-variable form -- a sequence, or
    a mapping over the feature names -- says nothing about the response, which
    keeps the default. Anything that broadcasts (an alias, one margin, a
    callable) applies to the response too, so ``margins="parametric"`` means
    what it looks like.

    Returns
    -------
    object
        One specification, not yet fitted.
    """
    if isinstance(self.margins, (list, tuple, dict)):
      return Kde1dMargin()
    return resolve_margins(self.margins, 1)[0]

  def _fit_marginals(
    self, X: np.ndarray, y: np.ndarray | None = None
  ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """
    Fit one margin per feature column, and the response's when given.

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
    specs = resolve_margins(
      self.margins,
      self.n_features_in_,
      names=getattr(self, "_expanded_columns", None),
      default=self._default_margin_specs(),
    )
    self._x_margins = tuple(
      self._fit_one_margin(specs[j], X[:, j], self._column_name(j))
      for j in range(self.n_features_in_)
    )
    fitted = list(self._x_margins)

    if y is not None:
      self._y_margin = self._fit_one_margin(
        self._response_margin_spec(), np.asarray(y, dtype=float), "y"
      )
      # The joint model orders the response first and gives it no left-limit
      # column, which the prediction paths rely on.
      if Vinedist.copula_var_types([self._y_margin])[0] != "c":
        raise ValueError(
          "The response margin must be continuous, but "
          f"{type(self._y_margin).__name__} declares "
          f"var_type={getattr(self._y_margin, 'var_type', 'c')!r}."
        )
      fitted.append(self._y_margin)

    self.selection_report_ = [
      dict(row) for margin in fitted for row in getattr(margin, "report_", ())
    ]
    if y is not None:
      return X, y
    return X

  def _bind_distribution(self, margins: Any) -> None:
    """Publish the fitted vine and its margins as one ``Vinedist``.

    The copula is wrapped so that it evaluates through the backend, which is
    what keeps ``distribution_`` numerically identical to the estimator's own
    methods -- including on the torch backend, where the raw vine returns
    tensors on its own device.

    Parameters
    ----------
    margins : sequence of MarginLike
        The fitted margins, in the vine's variable order.

    Returns
    -------
    None
    """
    self.distribution_ = Vinedist(
      _BackendVinecop(self.backend_, self._vine), list(margins)
    )

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

    A thin delegation to ``Vinedist.copula_data``, which owns the layout: one
    column per variable, then one left-limit column per variable with atoms.

    Parameters
    ----------
    Z : array-like
        If is_y=False: shape (n_samples, d).
        If is_y=True : shape (n_samples,).
    is_y : bool, default=False
        If True, transform the response with its own margin.
        If False, transform the covariates with theirs.

    Returns
    -------
    U : ndarray
        If is_y=False: shape (n_samples, d [+ optional discrete sub-columns]).
        If is_y=True : shape (n_samples, 1).
    """
    check_is_fitted(self, attributes=["_x_margins"])
    if is_y and not hasattr(self, "_y_margin"):
      raise RuntimeError(
        "Target marginal not fitted (y was not provided during fit)."
      )

    Z = np.asarray(Z, dtype=float)
    if is_y:
      return np.asarray(
        Vinedist.copula_data([self._y_margin], Z.reshape(-1, 1))
      )
    return np.asarray(Vinedist.copula_data(self._x_margins, Z))

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
      var_types = Vinedist.copula_var_types(self._x_margins)

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
    check_is_fitted(self, attributes=["distribution_"])

    X = self._validate_input(X, reset=False)
    if y is not None:
      y = np.asarray(y, dtype=float).reshape(-1, 1)
      if X.shape[0] != y.shape[0]:
        raise ValueError("X and y must have the same number of samples.")
      # The joint model orders the response first.
      Z = np.column_stack([y, X])
    else:
      Z = np.asarray(X, dtype=float)

    dist = self.distribution_
    if copula_only:
      u = Vinedist.copula_data(dist.margins, Z)
      out = np.log(np.asarray(dist.copula.pdf(u)))
    else:
      out = np.asarray(dist.logpdf(Z))

    return np.asarray(out if log else np.exp(out))
