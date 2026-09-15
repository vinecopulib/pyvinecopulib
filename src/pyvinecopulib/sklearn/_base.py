import copy
import os
import warnings
from numbers import Integral
from typing import Any, Optional, Sequence, Union, overload

import numpy as np
import pandas as pd
from pandas.api.extensions import ExtensionDtype
from sklearn.base import BaseEstimator
from sklearn.exceptions import DataConversionWarning
from sklearn.utils._param_validation import Interval, Options
from sklearn.utils.validation import (
  assert_all_finite,
  check_is_fitted,
  check_random_state,
)

import pyvinecopulib as pv

from ..core import (
  ControlsLike,
  FitControlsMargin,
  MarginLike,
  Vinedist,
  VinedistBase,
)
from ..margins import resolve_margins
from ..core._loglik import safe_log
from ..core._margins import MarginSpec, fit_margin
from ..core.extend import to_numpy

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
density (``Kde1d``) per column and accepts anything
:func:`pyvinecopulib.margins.resolve_margins` understands --- an alias
such as ``"parametric"``, one margin broadcast to
every column, a per-column sequence, or a mapping keyed by feature
name. Both halves come from ``distribution=``: the default
``Vinedist`` pairs ``Vinecop`` with ``Kde1d`` and has no extra
dependencies, while :class:`pyvinecopulib.torch.TorchVinedist` routes the
same pipeline through the PyTorch evaluator (GPU / autograd).

Fitting assembles a ``Vinedist`` --- the copula and its margins as one
object --- and every post-fit method evaluates through it. It is
published as ``distribution_``, so the fitted joint distribution is
usable outside the estimator; under ``TorchVinedist`` that object is a
``TorchVinedist``, which stays on its device and differentiable even
though the estimator's own methods return arrays; ``margin_summary_``
describes the margin each variable ended up with. See the
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
columns are handled by the margin's own variable type --- for the
default ``Kde1d`` that is ``type="discrete"``:
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


#: `schema_["kde1d_types"]` carries `Kde1d`'s spellings; a margin declares the
#: contract's. One map, so the two never drift apart.
_VAR_TYPE_OF = {"continuous": "c", "discrete": "d", "zero-inflated": "zi"}

#: What an estimator accepts for ``X``: the two shapes the pipeline models
#: directly, or any nested sequence of rows -- sklearn's convention is that
#: an array-like is valid input, and `VineBase._validate_input` coerces one.
_XLike = Union[np.ndarray, pd.DataFrame, Sequence[Sequence[object]]]

#: What it accepts for ``y``, and for per-observation weights: one column of
#: numbers, as an array or a sequence.
_YLike = Union[np.ndarray, Sequence[float]]

#: What ``random_state=`` accepts, per the scikit-learn convention that
#: `sklearn.utils.check_random_state` implements.
_RandomStateLike = Union[int, np.random.RandomState, None]


def _named_for(name: str, exc: BaseException) -> BaseException:
  """``exc`` with the column it came from named, keeping its type where it can.

  A margin sees one array and cannot say which column it was, so the estimator
  says it. The type is preserved by rebuilding the exception from the new
  message, which most exceptions accept; the ones that do not -- anything whose
  constructor takes more than a message -- become a ``ValueError`` rather than
  a confusing ``TypeError`` raised from inside the handler.

  Parameters
  ----------
  name : str
      The column's name.
  exc : BaseException
      What the margin raised.

  Returns
  -------
  BaseException
      The exception to raise, to be chained from ``exc``.
  """
  message = f"margin for {name!r}: {exc}"
  try:
    return type(exc)(message)
  except Exception:
    return ValueError(message)


def _categorical_bounds(
  dtype: Union[np.dtype, ExtensionDtype],
) -> Optional[tuple[float, float]]:
  """Exact support of an ordered categorical column, when it states one.

  Parameters
  ----------
  dtype : numpy.dtype or pandas.api.extensions.ExtensionDtype
      A column's dtype.

  Returns
  -------
  tuple of float, or None
      ``(min, max)`` over the declared categories, or ``None`` when the dtype is
      not categorical or its categories are not numeric.
  """
  if not isinstance(dtype, pd.CategoricalDtype):
    return None
  try:
    levels = np.asarray(dtype.categories, dtype=float)
  except (TypeError, ValueError):
    return None
  if levels.size == 0:
    return None
  return (float(levels.min()), float(levels.max()))


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
  - Vine copula fitting (via ``distribution``'s own vine class)
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

  #: Whether this estimator reports a density on the original scale, and so
  #: cannot accept a margin that has none. `VineRegressor` does not: its
  #: quadrature reads the copula density and the response margin's `icdf` only.
  _needs_marginal_density: bool = False

  _parameter_constraints: dict[str, list[object]] = {
    "distribution": [type, None],
    "controls": [object, None],
    "structure": [object, None],
    "margins": [object, None],
    "batch_size": [Interval(Integral, 1, None, closed="left")],
    "random_state": ["random_state"],
    "n_jobs": [
      Interval(Integral, 1, None, closed="left"),
      Options(Integral, {-1}),
      None,
    ],
  }

  def __init__(
    self,
    distribution: Optional[type[VinedistBase[Any]]] = None,
    controls: Optional[ControlsLike] = None,
    structure: Optional[pv.RVineStructure] = None,
    margins: object = None,
    batch_size: int = 100,
    random_state: _RandomStateLike = None,
    n_jobs: Optional[int] = None,
  ) -> None:
    """Base vine copula estimator.

    Parameters
    ----------
    distribution : type, or None, optional
        The ``VinedistBase`` subclass to fit --- ``Vinedist`` (the default)
        or :class:`pyvinecopulib.torch.TorchVinedist`. It names both halves
        of the model, so naming it is the whole lane choice: which vine class
        is fitted, which margin class an unaddressed column gets, and which
        array namespace the published ``distribution_`` lives on. Importing
        ``TorchVinedist`` to name it here is the explicit opt-in to PyTorch.
    controls : ControlsLike, or None, optional
        Fit-time controls for the copula half: a ``FitControlsVinecop``, or a
        :class:`pyvinecopulib.torch.FitControlsTorchVinecop` alongside
        ``TorchVinedist``. `None` fits nonparametric (TLL) pair copulas
        truncated at depth 20. Stored as-is and never mutated --- the
        estimator copies before writing its own seeds in.
    structure : RVineStructure, or None, optional
        A pre-specified vine structure; when given, structure selection is
        skipped and the vine is fitted on it. `None` selects one.
    margins : object, or None, optional
        What to fit to each column, in any form
        :func:`pyvinecopulib.margins.resolve_margins` accepts: an alias
        (``"kde"``, ``"parametric"``), one margin
        broadcast to every column, a sequence of length
        ``n_features_in_``, a mapping keyed by feature name or
        position, or a callable taking a column and returning a
        margin. `None` fits the distribution's own kernel-density margin
        per column --- ``Kde1d``, or ``TorchKde1d`` under
        ``TorchVinedist`` --- carrying the variable type inferred from the input
        and, for the ``{0, 1}`` dummies of an expanded unordered
        categorical, the bounds of its support. Stored as-is and never
        mutated: every specification is fitted on a copy.
    batch_size : int, default=100
        Number of test points to process per batch when making
        predictions. ``1`` minimizes memory at the cost of speed;
        ``n_test`` is the opposite extreme; intermediate values
        trade off memory and throughput.
    random_state : int, RandomState instance, or None, optional
        Seeds the RNG used by stochastic operations (e.g.
        `sample`, `cdf` quasi-MC, structure simulation). Stored
        as-is; resolved via `sklearn.utils.check_random_state`
        inside `fit`.
    n_jobs : int, or None, optional
        Threads the vine may use, for fitting and for every evaluation
        (`pdf`, `cdf`, `sample`, and the prediction paths built on them).
        `None` means one thread and `-1` means every processor, following
        the scikit-learn convention. Results never depend on it: the fitted
        structure, the fitted pair copulas and every evaluated value are
        bit-identical at any thread count.

        `None` has a reason: a caller that parallelizes *over* vines owns
        the parallelism, and nesting would oversubscribe the machine. Set it
        when one vine is the whole job.
    """
    self.distribution = distribution
    self.controls = controls
    self.structure = structure
    self.margins = margins
    self.batch_size = batch_size
    self.random_state = random_state
    self.n_jobs = n_jobs

  def _reset_fitted_schema(self) -> None:
    """Drop the state a previous ``fit`` derived, before deriving it again.

    A caller may pre-set ``schema_`` to declare per-column types for array
    input, so only a schema this estimator derived itself is discarded --
    otherwise a second ``fit`` would validate the new data against the
    previous fit's schema, which sklearn's contract forbids.
    """
    if getattr(self, "_schema_from_fit", False):
      self.__dict__.pop("schema_", None)
    for name in (
      "_schema_from_fit",
      "feature_names_in_",
      "_dtypes",
      "_categories",
      "_expanded_columns",
      "n_model_features_",
    ):
      self.__dict__.pop(name, None)

  @overload
  def _validate_input(
    self,
    X: _XLike,
    y: None = None,
    *,
    reset: bool,
  ) -> np.ndarray: ...

  @overload
  def _validate_input(
    self,
    X: _XLike,
    y: _YLike,
    *,
    reset: bool,
  ) -> tuple[np.ndarray, np.ndarray]: ...

  def _validate_input(
    self,
    X: _XLike,
    y: Optional[_YLike] = None,
    *,
    reset: bool,
  ) -> Union[np.ndarray, tuple[np.ndarray, np.ndarray]]:
    """Validate the ``X`` (and optional ``y``) input.

    For DataFrames, captures the canonical
    :attr:`feature_names_in_`, expands unordered categoricals via
    :func:`expand_factors`, and infers ``kde1d_types``. For ndarrays,
    captures :attr:`n_features_in_` and assumes all-continuous unless
    a previously-set ``schema_["kde1d_types"]`` says otherwise.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features), or DataFrame
    y : array-like, or None, optional
    reset : bool
        If ``True`` (called from ``fit``), set fitted attributes; if
        ``False`` (called from ``predict``-style methods), validate
        against the previously-set ones.
    """
    if not isinstance(X, (np.ndarray, pd.DataFrame)):
      # `.format` plus `.nnz` identifies a scipy sparse matrix; `.format`
      # alone also matches `str` and anything else carrying that attribute.
      if hasattr(X, "format") and hasattr(X, "nnz"):
        raise TypeError("Sparse input is not supported; pass a dense array.")
      # sklearn's convention is that any array-like is acceptable input, so a
      # list of rows must not be refused.
      try:
        X = np.asarray(X, dtype=float)
      except (TypeError, ValueError) as error:
        raise ValueError(
          "X must be an array-like of floats, a numpy array or a pandas "
          f"DataFrame; got {type(X).__name__}"
        ) from error
    if isinstance(X, np.ndarray):
      if X.ndim != 2:
        raise ValueError(f"Expected 2D array, got {X.ndim}D array instead.")
      if X.shape[0] == 0:
        raise ValueError(
          "Found array with 0 sample(s) while a minimum of 1 is required."
        )
      if X.shape[1] == 0:
        raise ValueError(
          "0 feature(s) (shape=(%d, 0)) while a minimum of 1 is required."
          % X.shape[0]
        )
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

    if reset:
      self._reset_fitted_schema()
    elif (
      isinstance(X, pd.DataFrame)
      and hasattr(self, "n_features_in_")
      and not hasattr(self, "feature_names_in_")
    ):
      # Fitted on an array, so there are no names to match against. sklearn's
      # convention is to accept a frame positionally rather than refuse it --
      # but a categorical column was never modeled as one, and silently
      # reading its codes as a numeric column would be a wrong answer.
      categorical = [
        c for c in X.columns if isinstance(X[c].dtype, pd.CategoricalDtype)
      ]
      if categorical:
        raise ValueError(
          f"{type(self).__name__} was fitted without feature names, so the "
          f"categorical column(s) {categorical} cannot be expanded the way "
          "they would have been at fit time; refit on a DataFrame."
        )
      X = X.to_numpy()

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
        # Bounds are passed wherever the input states them and left unset
        # wherever they would be a guess: a categorical declares its levels, a
        # continuous column states nothing. Without them the kernel-density
        # grid is padded past the data and puts mass on values that cannot
        # occur -- a count column picks up density below zero.
        original = set(X.columns)
        bounds: list[Optional[tuple[float, float]]] = []
        for name, dtype in zip(self._expanded_columns, X_exp.dtypes):
          if name not in original:
            bounds.append((0.0, 1.0))
          else:
            bounds.append(_categorical_bounds(dtype))
        self.schema_ = {"kde1d_types": kde1d_types, "bounds": bounds}
        self._schema_from_fit = True
        self.n_features_in_ = X.shape[1]
        self.n_model_features_ = X_exp.shape[1]
        X_arr = X_exp.to_numpy()
      else:
        check_is_fitted(self, attributes=["n_features_in_"])
        if list(X.columns) != list(self.feature_names_in_):
          raise ValueError("Column names/order do not match training data.")
        X_for_expansion = X.copy(deep=True)
        for col in self.feature_names_in_:
          dtype_expected = self._dtypes[col]
          if isinstance(dtype_expected, pd.CategoricalDtype):
            if not isinstance(X_for_expansion[col].dtype, pd.CategoricalDtype):
              raise ValueError(f"Column {col} must be categorical.")
            recoded = X_for_expansion[col].cat.set_categories(
              dtype_expected.categories, ordered=dtype_expected.ordered
            )
            # `set_categories` maps a level the fit never saw to NaN, whose
            # dummy expansion is an all-zero row -- indistinguishable from the
            # reference level, so it has to be caught rather than evaluated.
            unseen = recoded.isna() & X_for_expansion[col].notna()
            if bool(unseen.any()):
              levels = sorted(
                {str(v) for v in X_for_expansion.loc[unseen, col].unique()}
              )
              raise ValueError(
                f"Column {col} contains categories not seen during fit: "
                f"{levels}"
              )
            X_for_expansion[col] = recoded
        X_arr = expand_factors(X_for_expansion)[
          self._expanded_columns
        ].to_numpy()
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
        # Only a schema this estimator derived on its own may be discarded on
        # the next `fit`. A caller who pre-set one is declaring something the
        # array cannot show -- which columns are discrete, where a variable is
        # bounded -- and re-inferring it would silently change the model on the
        # second fit. That is what an ensembling wrapper does to every survivor.
        self._schema_from_fit = not existing
        self.n_features_in_ = X.shape[1]
        self.n_model_features_ = X.shape[1]
      else:
        check_is_fitted(self, attributes=["n_features_in_"])
        # `sample` emits the modeled layout, which is wider than the public
        # one whenever a categorical was expanded, so the estimator's own
        # output has to be a legal input to its own density.
        accepted = {self.n_features_in_, self.n_model_features_}
        if X.shape[1] not in accepted:
          expected = " or ".join(str(n) for n in sorted(accepted))
          raise ValueError(
            f"X has {X.shape[1]} features, but {type(self).__name__} is "
            f"expecting {expected} features as input"
          )
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

  def _default_margin_specs(self) -> list[MarginLike[Any]]:
    """One unfitted kernel-density margin per column, from the schema.

    This is what ``margins=None`` means, and what a column a ``margins=``
    mapping does not address falls back to. The *class* comes from
    ``distribution``, so a torch copula gets torch margins: fitting NumPy ones
    onto it would put the whole distribution on two array namespaces and stop
    every gradient at the marginal transform.

    Returns
    -------
    list of MarginLike
        One margin per feature, carrying the variable type and whatever bounds
        the input told us about.
    """
    types = self.schema_["kde1d_types"]
    bounds = self.schema_.get("bounds") or [None] * len(types)
    specs = []
    for j, (type_, bound) in enumerate(zip(types, bounds)):
      pair = None if bound is None else (float(bound[0]), float(bound[1]))
      try:
        specs.append(self._default_margin(type_, pair))
      except (ValueError, RuntimeError) as exc:
        # The bounds come from the column's own dtype, so a margin that refuses
        # them is a statement about that column -- most often an ordered
        # categorical whose levels are not integers, which cannot be a discrete
        # `Kde1d`. The margin cannot name the column; we can.
        raise _named_for(self._column_name(j), exc) from exc
    return specs

  def _default_margin(
    self, var_type: str, bounds: Optional[tuple[float, float]]
  ) -> MarginLike[Any]:
    """One unfitted margin of the distribution's own margin class.

    Routed through ``_default_margins``, which is where each lane already
    states how its margin class takes a variable type and bounds at
    construction --- a kernel density fitted unbounded has padded past the
    data by the time anything could tell it otherwise.

    Resolved rather than read off a fitted attribute: this is an internal a
    caller may reach before ``fit`` has pinned one, and the resolution is
    cheap and idempotent.

    Parameters
    ----------
    var_type : str
        ``Kde1d``'s spelling of the variable type, as ``schema_`` carries it.
    bounds : tuple of float, or None, optional
        Declared support, or ``None`` where the input states none.

    Returns
    -------
    MarginLike
        An unfitted margin, on the distribution's own array namespace.
    """
    cls: type[VinedistBase[Any]] = getattr(
      self, "distribution_class_", None
    ) or (self.distribution or Vinedist)
    controls = getattr(self, "controls_", None) or self.controls
    declared = FitControlsMargin(
      var_type=_VAR_TYPE_OF.get(var_type, var_type), support=bounds
    )
    margins = cls._default_margins(1, controls, [declared])
    if margins is None:
      raise TypeError(
        f"{cls.__name__} names no `margin_class`, so it cannot supply the "
        "default margin an unaddressed column needs. Pass `margins=` "
        "explicitly, or name a `margin_class` on the distribution."
      )
    return margins[0]

  def _declared_for(
    self, index: Optional[int]
  ) -> tuple[Optional[str], Optional[tuple[float, float]]]:
    """What :attr:`schema_` knows about one feature.

    ``_default_margin_specs`` only reaches the columns a ``margins=`` argument
    leaves unaddressed, so without this a `margins="parametric"` on a frame with
    an ordered categorical would have the margin re-infer both the type and the
    bounds from the sample -- less than the input already stated.

    Parameters
    ----------
    index : int, or None, optional
        The feature's position, or ``None`` for the response.

    Returns
    -------
    tuple
        ``(var_type, support)``, either entry ``None`` when unknown.
    """
    if index is None:
      return None, None
    types = self.schema_["kde1d_types"]
    var_type = _VAR_TYPE_OF.get(types[index])
    bounds = (self.schema_.get("bounds") or [None] * len(types))[index]
    support = None if bounds is None else (float(bounds[0]), float(bounds[1]))
    return var_type, support

  def _fit_one_margin(
    self,
    spec: MarginSpec,
    column: np.ndarray,
    name: str,
    *,
    index: Optional[int] = None,
  ) -> MarginLike[Any]:
    """Fit one column's margin.

    Parameters
    ----------
    spec : MarginSpec
        A specification from :func:`pyvinecopulib.margins.resolve_margins`.
    column : ndarray, shape (n_samples,), dtype float
        The column to fit.
    name : str
        The variable's name, for a selector's report.
    index : int, or None, optional
        Which feature this is, used to hand the margin the variable type and
        bounds from :attr:`schema_`. ``None`` for the response, whose type the
        estimator already constrains to continuous.

    Returns
    -------
    MarginLike
        The fitted margin.
    """
    # Fitting a margin mutates it, and `self.margins` is a constructor argument
    # that has to survive `fit` untouched so `clone` reproduces the estimator;
    # hence every specification is fitted on a copy of itself.
    spec = copy.deepcopy(spec)
    var_type, support = self._declared_for(index)
    try:
      margin = fit_margin(
        spec,
        np.asarray(column, dtype=float),
        var_type=var_type,
        support=support,
      )
    except (ValueError, RuntimeError) as exc:
      # Whatever went wrong, the caller needs to know which column it was: the
      # margin only ever sees one array and cannot name it. `Kde1d` refusing a
      # non-integer discrete bound or observation is the common case, and it is
      # reachable from an ordered categorical whose levels are not integers.
      raise _named_for(name, exc) from exc
    if self._needs_marginal_density:
      self._require_density(margin, name, column)
    return margin

  def _require_density(
    self, margin: MarginLike[Any], name: str, column: np.ndarray
  ) -> None:
    """Refuse a margin with no density, at fit time rather than at first score.

    A margin whose ``pdf`` is undefined -- an atomic distribution, whose mass
    is not a density -- cannot serve an estimator that reports one. Probed
    rather than declared, so a margin from anywhere is caught.

    Parameters
    ----------
    margin : MarginLike
        The fitted margin.
    name : str
        The variable's name, for the message.
    column : ndarray, shape (n_samples,), dtype float
        The data it was fitted to, to probe at.

    Raises
    ------
    ValueError
        If the margin has no density.
    """
    probe = np.asarray(column, dtype=float)[:1]
    density: Any = getattr(margin, "logpdf", None) or margin.pdf
    try:
      density(probe)
    except NotImplementedError as exc:
      raise ValueError(
        f"{type(self).__name__} needs a density on the original scale, but the "
        f"margin for {name!r} ({type(margin).__name__}) has none: {exc}. Pass "
        'margins="kde" (the default), or any margin that reports a density, '
        "for that column."
      ) from exc

  def _response_margin_spec(self) -> MarginSpec:
    """The specification for the response margin.

    ``margins`` addresses the features, so a per-variable form -- a sequence, or
    a mapping over the feature names -- says nothing about the response, which
    keeps the default. Anything that broadcasts (an alias, one margin, a
    callable) applies to the response too, so ``margins="parametric"`` means
    what it looks like.

    Returns
    -------
    MarginSpec
        One specification, not yet fitted.
    """
    if isinstance(self.margins, (list, tuple, dict)):
      # Same hook as `_default_margin_specs`: without it the torch lane would
      # give the covariates torch margins and the response a NumPy one. Under
      # `Vinedist` this *is* `Kde1d()`.
      return self._default_margin("continuous", None)
    return resolve_margins(self.margins, 1)[0]

  @staticmethod
  def _check_response_is_continuous(margin: object) -> None:
    """Refuse a response margin with atoms.

    The joint model orders the response first and gives it no left-limit
    column, which the prediction paths rely on. Called twice: once on the
    specification, so a margin that declares its type up front is refused
    before anything is fitted, and once on the fitted margin, for the
    specifications that can only declare one afterwards.

    Parameters
    ----------
    margin : object
        A margin specification or a fitted margin.

    Raises
    ------
    ValueError
        If the margin declares a discrete or zero-inflated variable type.
    """
    try:
      kind = Vinedist.copula_var_types([margin])[0]
    except (RuntimeError, TypeError):
      # Not every specification can answer yet: a selector has no variable type
      # until it has chosen one, and a callable is not a margin until it has
      # been called. Nothing to check here; the call on the fitted margin binds.
      return
    if kind == "c":
      return
    raise ValueError(
      "The response margin must be continuous, but "
      f"{type(margin).__name__} declares "
      f"var_type={getattr(margin, 'var_type', 'c')!r}."
    )

  def _fit_marginals(
    self, X: np.ndarray, y: np.ndarray | None = None
  ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """
    Fit one margin per feature column, and the response's when given.

    Parameters
    ----------
    X : ndarray
        Input features.
    y : ndarray, or None, optional
        Target values. If None, only X marginals are fitted.

    Returns
    -------
    X : ndarray
        Input features (unchanged).
    y : ndarray, or None
        Target values (unchanged) if provided, None otherwise.
    """
    specs = resolve_margins(
      self.margins,
      self.n_model_features_,
      names=getattr(self, "_expanded_columns", None),
      # Lazily: building it constructs one margin per column, which can refuse
      # a column outright, and a specification naming every column never uses
      # it.
      default=self._default_margin_specs,
    )
    y_spec = None if y is None else self._response_margin_spec()
    if y_spec is not None:
      self._check_response_is_continuous(y_spec)

    self._x_margins = tuple(
      self._fit_one_margin(specs[j], X[:, j], self._column_name(j), index=j)
      for j in range(self.n_model_features_)
    )
    fitted = list(self._x_margins)

    if y_spec is not None:
      self._y_margin = self._fit_one_margin(
        y_spec, np.asarray(y, dtype=float), "y"
      )
      self._check_response_is_continuous(self._y_margin)
      fitted.append(self._y_margin)

    if y is not None:
      return X, y
    return X

  def _bind_distribution(self, margins: Sequence[MarginLike[Any]]) -> None:
    """Publish the fitted vine and its margins as one distribution.

    The class is ``distribution``'s, so the object is on the same array
    namespace as the vine that was fitted --- a ``TorchVinedist`` stays
    differentiable and movable. A class that lifts a margin hands back the
    lifted one, so the margins are re-read from the result rather than kept in
    two places.

    Parameters
    ----------
    margins : sequence of MarginLike
        The fitted margins, in the vine's variable order.

    Returns
    -------
    None
    """
    self.distribution_ = self.distribution_class_(self._vine, list(margins))
    bound = self.distribution_.margins
    if hasattr(self, "_y_margin"):
      self._y_margin = bound[0]
      self._x_margins = tuple(bound[1:])
    else:
      self._x_margins = tuple(bound)
    # An optional capability on `VinedistLike`, read the way the contract
    # says -- but this estimator publishes `margin_summary_`, so a
    # distribution that has none is named rather than silently summaryless.
    summary = getattr(self.distribution_, "margin_summary", None)
    if summary is None:
      raise TypeError(
        f"{type(self.distribution_).__name__} has no `margin_summary`, which "
        f"{type(self).__name__} publishes as `margin_summary_`. Return a "
        "`VinedistBase` subclass as `distribution`, or add a "
        "`margin_summary()` to the one you pass."
      )
    self.margin_summary_ = summary()

  def _resolve_runtime_state(self) -> None:
    """Resolve the random state, the distribution class and the controls.

    Sets ``random_state_``, ``distribution_class_`` and ``controls_`` so
    subclasses can reuse them throughout ``fit`` and the post-fit methods.
    ``controls_`` is a copy whenever this writes to it, so the caller's own
    ``controls`` object is never mutated.
    """
    self.random_state_ = check_random_state(self.random_state)
    self.distribution_class_ = self.distribution or Vinedist
    controls = self.controls
    if controls is None:
      controls = self._default_copula_controls()
    if self.n_jobs is not None and hasattr(controls, "num_threads"):
      threads = os.cpu_count() or 1 if self.n_jobs == -1 else int(self.n_jobs)
      # A local `Any`: `ControlsLike` names only `to_dict`, and `num_threads`
      # is one lane's field -- the `hasattr` above is the check.
      threaded: Any = copy.copy(controls)
      threaded.num_threads = threads
      controls = threaded
    self.controls_ = controls

  def _default_copula_controls(self) -> Optional[ControlsLike]:
    """The copula controls when the caller named none.

    Nonparametric (TLL) pair copulas truncated at depth 20 --- a density
    estimator's default, rather than the copula library's parametric search.
    A lane whose controls this cannot describe gets ``None``, which lets the
    vine class resolve its own default.

    Returns
    -------
    ControlsLike, or None
        The controls to fit with.
    """
    if self.distribution_class_.vinecop_class is pv.Vinecop:
      return pv.FitControlsVinecop(
        family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
      )
    return None

  @property
  def _num_threads(self) -> int:
    """Threads the copula evaluates on; a lane naming none evaluates on one.

    Returns
    -------
    int
        The thread count.
    """
    return int(getattr(self.controls_, "num_threads", 1))

  def _draw_seeds(self, size: int = 5) -> list[int]:
    """Derive a list of ints suitable for ``FitControlsVinecop.seeds``
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
      return to_numpy(
        Vinedist.copula_data([self._y_margin], Z.reshape(-1, 1)), dtype=float
      )
    return to_numpy(Vinedist.copula_data(self._x_margins, Z), dtype=float)

  def _fit_vine(
    self, U: np.ndarray, var_types: list[str] | None = None
  ) -> "VineBase":
    """
    Fit vine copula to pseudo-observations.

    Parameters
    ----------
    U : ndarray
        Pseudo-observations in [0,1]^d.
    var_types : list, or None, optional
        Variable types for vine fitting. If None, inferred from schema.

    Returns
    -------
    self
    """
    if var_types is None:
      var_types = Vinedist.copula_var_types(self._x_margins)

    controls = self.controls_
    # Seeds the caller named win; `random_state` only fills in for a random
    # tree search that would otherwise draw from entropy.
    if (
      self.structure is None
      and not getattr(controls, "seeds", None)
      and getattr(controls, "tree_algorithm", "").startswith("random")
    ):
      seeded: Any = copy.copy(controls)
      seeded.seeds = self._draw_seeds()
      controls = seeded
    vinecop_class: Any = self.distribution_class_.vinecop_class
    self._vine = vinecop_class.from_data(
      U, controls, structure=self.structure, var_types=var_types
    )
    self.structure_ = self._vine.structure
    return self

  # `copula_only=True` skips the marginal-density product and returns
  # the copula factor c(u) alone; `VineDensity.pdf` exposes it directly.
  def _pdf_samples(
    self,
    X: _XLike,
    y: np.ndarray | None = None,
    log: bool = False,
    copula_only: bool = False,
  ) -> np.ndarray:
    """
    Unified method to compute density/likelihood of samples.

    Parameters
    ----------
    X : array-like of float, shape (n_samples, n_features), or DataFrame
        Input features.
    y : array-like of shape (n_samples,), or None, optional
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
      u = dist.copula_layout(Z)
      # The copula's own log-density where it has one, as `logpdf` reads it in
      # the branch below: the density is a product over edges and underflows
      # on a deep or strongly dependent vine.
      # A local `Any`: `logpdf` is an optional capability the `getattr` below
      # checks for, and `num_threads` is a hint every `VinecopBase` subclass
      # and the core `Vinecop` accept -- neither belongs on the protocol,
      # which stays the narrow evaluation contract a foreign vine must meet.
      vine: Any = dist.vinecop
      nt = self._num_threads
      if getattr(vine, "logpdf", None) is not None:
        out = to_numpy(vine.logpdf(u, num_threads=nt), dtype=float)
      else:
        out = safe_log(to_numpy(vine.pdf(u, num_threads=nt), dtype=float))
    else:
      out = to_numpy(dist.logpdf(Z), dtype=float)

    return np.asarray(out if log else np.exp(out))
