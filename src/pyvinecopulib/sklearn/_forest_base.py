"""Private base class for the vine-copula forest estimators.

Implements the shared pipeline (hold-out validation split, random
structure generation, optional bootstrap, MCS-based survivor
selection, parallel fitting, ensemble re-fit on full data) for
:class:`VineForestDensity` and :class:`VineForestRegressor`. The two
public subclasses only override the abstract hooks
``_create_base_estimator`` and ``_loglik_estimator``.
"""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from numbers import Integral, Real
from typing import Any

import numpy as np
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator
from sklearn.model_selection import train_test_split
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.validation import check_is_fitted, check_random_state

from ._mcs import da_mcs_marg, da_mcs_unif

_VALID_METHODS = {"da_mcs_marg", "da_mcs_unif", None}

# Forest-specific docstring fragments, interpolated into VineForestDensity
# / VineForestRegressor class docstrings via f-string. Kept here (not in
# `_base.py`) so the single-vine class docs stay free of forest concepts;
# `_DOC_FOREST_REFERENCES` is appended to the shared single-vine
# `_DOC_REFERENCES` for the forest classes.

_DOC_FOREST = r"""**Forest construction.** The ensemble is built via
hold-out random search with model-confidence-set (MCS) survivor
selection:

1. **Validation split.** A fraction ``val_fraction`` of the training
   data is held out for survivor selection.
2. **Random search.** ``n_vines`` candidate vine structures are
   sampled. When ``vines_sampling="uniform"``, the structure is drawn
   uniformly over R-vines via Joe's generation algorithm (Joe,
   Cooke & Kurowicka 2011; Algorithm 13 of Joe 2014; implemented in
   :meth:`pyvinecopulib.core.RVineStructure.simulate`). When
   ``"local"``, each tree :math:`T` is drawn from the
   Kendall's-:math:`\tau`-weighted distribution
   :math:`P(T) \propto \prod_{(j, k) \in T} |\hat\tau_{j, k}|`
   (with the Dissmann maximum spanning tree as its mode) via Wilson's
   loop-erased random walk (Wilson 1996); pyvinecopulib exposes this
   as the ``"random_weighted"`` tree algorithm in
   :class:`pyvinecopulib.core.FitControlsVinecop`. Each candidate is
   fit on a bootstrap resample of the training portion.
3. **Survivor selection.** Each candidate's per-sample log-likelihood
   is evaluated on the validation set. The resulting :math:`n_{\text{val}}
   \times M` loss matrix is fed to a dual-split DA test — adapted by
   Vatter & Nagler (2026) from the discrete-argmin-inference framework
   of Kim & Ramdas (2025) — to compute the model confidence set (Hansen,
   Lunde & Nason 2011 for the foundational definition). ``method``
   selects between marginal per-model coverage (``"da_mcs_marg"``) and
   uniform / familywise coverage (``"da_mcs_unif"``).
   ``add_dissmann=True`` adds the Dissmann-structure baseline as an
   extra candidate.
4. **Prediction averaging.** Survivors are refit on the full training
   data; predictions are averaged across survivors at evaluation time.
   For the regressor, the conditional weights derived from
   :math:`c_{Y, X}` are averaged then row-normalised once at the
   ensemble level.
"""

_DOC_FOREST_REFERENCES = r"""- Joe, H., Cooke, R. M. and Kurowicka, D. (2011).
  *Regular vines: generation algorithm and number of equivalence
  classes.* In Dependence Modeling: Vine Copula Handbook, 219--231.
  World Scientific.
- Joe, H. (2014). *Dependence Modeling with Copulas.* CRC Press.
- Wilson, D. B. (1996). *Generating random spanning trees more
  quickly than the cover time.* Proceedings of STOC '96, 296--303.
- Hansen, P. R., Lunde, A. and Nason, J. M. (2011).
  *The Model Confidence Set.* Econometrica, 79(2), 453--497.
- Kim, I. and Ramdas, A. (2025). *Locally minimax optimal and
  dimension-agnostic discrete argmin inference.* arXiv:2503.21639.
"""


class VineForestBase(BaseEstimator, ABC):
  """Base class for vine-copula ensemble estimators.

  Centralises random-search structure generation, bootstrap
  resampling, validation split, MCS-based survivor selection, and
  parallel fitting / prediction across the ensemble. Subclasses bind
  the base-learner type (``VineDensity`` or ``VineRegressor``) and
  the per-estimator log-likelihood used to score candidates against
  the validation set.

  Per the scikit-learn developer guide, ``__init__`` performs no
  validation; ``_parameter_constraints`` + ``_validate_params()`` run
  at fit time. Fitted attributes use the trailing-underscore
  convention.
  """

  _parameter_constraints: dict = {
    "base_class": [type],
    "base_params": [dict, None],
    "n_vines": [Interval(Integral, 1, None, closed="left")],
    "vines_sampling": [StrOptions({"uniform", "local"})],
    "bootstrap": ["boolean"],
    "val_fraction": [Interval(Real, 0.0, 1.0, closed="left")],
    "best_only": ["boolean"],
    "method": [StrOptions({"da_mcs_marg", "da_mcs_unif"}), None],
    "alpha": [Interval(Real, 0.0, 1.0, closed="neither")],
    "add_dissmann": ["boolean"],
    "random_state": ["random_state"],
    "n_jobs": [Integral, None],
    "verbose": ["boolean"],
  }

  def __init__(
    self,
    base_class: type[Any],
    base_params: dict[str, Any] | None = None,
    n_vines: int = 100,
    vines_sampling: str = "uniform",
    bootstrap: bool = True,
    val_fraction: float = 0.25,
    best_only: bool = False,
    method: str | None = "da_mcs_marg",
    alpha: float = 0.05,
    add_dissmann: bool = True,
    random_state=None,
    n_jobs: int = 1,
    verbose: bool = False,
  ) -> None:
    """
    Base ensemble estimator for vine copulas.

    Parameters
    ----------
    base_class : class
        The base estimator class (e.g. :class:`VineRegressor`,
        :class:`VineDensity`).
    base_params : dict, optional
        Keyword arguments forwarded to each base estimator's
        ``__init__``. Example: ``{"batch_size": 200}``.
    n_vines : int, default=100
        Number of base estimators in the ensemble (before MCS-based
        pruning).
    vines_sampling : {"uniform", "local"}, default="uniform"
        How random structures are generated.

        - ``"uniform"`` draws uniformly over all R-vine structures
          via Joe's generation algorithm (Joe, Cooke & Kurowicka 2011;
          Algorithm 13 of Joe 2014); see
          :meth:`pyvinecopulib.core.RVineStructure.simulate`.
        - ``"local"`` draws each tree :math:`T` from the
          Kendall's-:math:`\\tau`-weighted distribution
          :math:`P(T) \\propto \\prod_{(j, k) \\in T}
          |\\hat\\tau_{j, k}|` (with the Dissmann MST as its mode)
          via Wilson's loop-erased random walk (Wilson 1996), exposed
          as the ``"random_weighted"`` tree algorithm in
          :class:`pyvinecopulib.core.FitControlsVinecop`.
    bootstrap : bool, default=True
        If ``True``, each random estimator is fit on a bootstrap
        resample of the training set.
    val_fraction : float, default=0.25
        Fraction of training data held out for the MCS-based
        selection step. ``0`` skips validation (every random
        estimator survives).
    best_only : bool, default=False
        Keep only the single best survivor instead of the entire MCS.
    method : {"da_mcs_marg", "da_mcs_unif", None}, default="da_mcs_marg"
        Survivor-selection method. ``"da_mcs_marg"`` / ``"da_mcs_unif"``
        run the dual-split DA test with marginal / uniform error
        control. ``None`` keeps every random estimator strictly
        better than the Dissmann baseline.
    alpha : float, default=0.05
        Significance level for the MCS selectors.
    add_dissmann : bool, default=True
        Add the Dissmann-structure baseline as an extra candidate.
    random_state : int, RandomState, Generator or None
        Seeds reproducibility across random structure generation,
        bootstrap resampling, and MCS randomisation. Stored as-is;
        resolved via :func:`sklearn.utils.check_random_state` in
        ``_fit_ensemble``.
    n_jobs : int, default=1
        Number of joblib workers for parallel fit / predict.
    verbose : bool, default=False
        Emit a warning if no random estimator beats the default.
    """
    self.base_class = base_class
    self.base_params = base_params
    self.n_vines = n_vines
    self.vines_sampling = vines_sampling
    self.bootstrap = bootstrap
    self.val_fraction = val_fraction
    self.best_only = best_only
    self.random_state = random_state
    self.n_jobs = n_jobs
    self.verbose = verbose
    self.method = method
    self.alpha = alpha
    self.add_dissmann = add_dissmann

  def _resolved_base_params(self) -> dict[str, Any]:
    """Return a defensive copy of ``base_params`` for forwarding to the
    base estimator. Defers validation to fit time."""
    return {} if self.base_params is None else dict(self.base_params)

  @abstractmethod
  def _create_base_estimator(self):
    """Build a fresh base estimator with any forest-specific overrides
    applied. After the sklearn-guide refactor, ``normalize_weights`` is
    a real ``__init__`` parameter on :class:`VineRegressor`, so this
    method is now equivalent to ``self.base_class(**base_params)`` and
    a follow-up :func:`sklearn.base.clone` round-trip preserves all
    state."""

  @abstractmethod
  def _loglik_estimator(self, estimator, X, y=None):
    """Per-sample log-likelihood used to score candidate estimators."""

  def _fit_estimator(
    self, estimator, X: np.ndarray, y: np.ndarray | None = None
  ):
    return estimator.fit(X) if y is None else estimator.fit(X, y)

  def _split_data(
    self, X: np.ndarray, y: np.ndarray | None = None
  ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None]:
    if self.val_fraction > 0:
      if y is not None:
        X_train, X_val, y_train, y_val = train_test_split(
          X,
          y,
          test_size=self.val_fraction,
          random_state=self._split_rng_seed,
        )
      else:
        X_train, X_val = train_test_split(
          X,
          test_size=self.val_fraction,
          random_state=self._split_rng_seed,
        )
        y_train, y_val = None, None
    else:
      X_train, y_train = X, y
      X_val, y_val = X, y

    return X_train, X_val, y_train, y_val

  def _fit_random_estimator(
    self,
    seed: int,
    X_train: np.ndarray,
    y_train: np.ndarray | None,
    X_val: np.ndarray,
    y_val: np.ndarray | None,
  ):
    """Fit a single random estimator and score it on the validation set."""
    # ``_create_base_estimator()`` returns a fresh, fully-cloneable
    # instance built from ``base_params``; since
    # ``normalize_weights`` is now a real init parameter, no
    # post-init attribute mutation is needed (and ``clone()`` would
    # preserve everything anyway).
    estimator = self._create_base_estimator()

    # Propagate the schema inferred on the full input (in particular,
    # discrete-vs-continuous flags from DataFrame columns) — by this
    # point X_train is a pre-expanded numpy array, so the base
    # estimator's own inference would default everything to continuous.
    if getattr(self._base_estimator, "schema_", None) is not None:
      estimator.schema_ = dict(self._base_estimator.schema_)

    local_rng = np.random.default_rng(seed)
    local_seeds = [int(x) for x in local_rng.integers(0, 2**31 - 1, size=5)]

    from .backends import resolve_backend

    parent_backend = resolve_backend(estimator.backend)
    if self.vines_sampling == "uniform":
      structure_size = X_train.shape[1] + int(y_train is not None)
      estimator.backend = parent_backend.with_random_structure(
        structure_size, seeds=local_seeds
      )
    elif self.vines_sampling == "local":
      estimator.backend = parent_backend.with_local_random(local_seeds)
    else:
      raise ValueError(f"Unknown vines_sampling method: {self.vines_sampling}")

    if self.bootstrap:
      idx = local_rng.integers(0, X_train.shape[0], size=X_train.shape[0])
      self._fit_estimator(
        estimator, X_train[idx], y_train[idx] if y_train is not None else None
      )
    else:
      self._fit_estimator(estimator, X_train, y_train)

    loglik = self._loglik_estimator(estimator, X_val, y_val)
    return estimator, loglik

  def _fit_ensemble(
    self, X: np.ndarray, y: np.ndarray | None = None
  ) -> "VineForestBase":
    """Core ensemble fitting pipeline."""
    self._validate_params()
    self.random_state_ = check_random_state(self.random_state)
    # `train_test_split` wants a plain seed/RandomState, not a Generator.
    # Derive a stable int from the resolved RNG so `_split_data` reuses
    # it across the default-estimator + survivor splits.
    self._split_rng_seed = int(self.random_state_.randint(0, 2**31 - 1))
    self._base_estimator = self._create_base_estimator()

    if y is not None:
      X, _ = self._base_estimator._validate_input(X, y, reset=True)
    else:
      X = self._base_estimator._validate_input(X, reset=True)
    self.n_features_in_ = X.shape[1]

    X_train, X_val, y_train, y_val = self._split_data(X, y)

    # Default estimator (Dissmann-structure baseline). Propagate the
    # inferred schema so fitting on the already-expanded numpy doesn't
    # default everything back to continuous.
    default = self._create_base_estimator()
    if getattr(self._base_estimator, "schema_", None) is not None:
      default.schema_ = dict(self._base_estimator.schema_)
    if y_train is not None:
      default.fit(X_train, y_train)
    else:
      default.fit(X_train)
    default_loglik = self._loglik_estimator(default, X_val, y_val)
    default_loglik_val = default_loglik.mean()

    seeds = self.random_state_.randint(0, 2**31 - 1, size=self.n_vines)

    def _fit_one(seed):
      return self._fit_random_estimator(seed, X_train, y_train, X_val, y_val)

    results = Parallel(n_jobs=self.n_jobs)(delayed(_fit_one)(s) for s in seeds)

    if self.add_dissmann:
      results.append((default, default_loglik))

    # Survivor selection.
    random_logliks_val = np.array([loglik.mean() for _, loglik in results])
    if self.method is None:
      survivors = random_logliks_val > default_loglik_val
    else:
      # Negative logliks → losses; selector minimises losses.
      nll = -np.array([loglik for _, loglik in results]).T
      mcs_rng = np.random.default_rng(self._split_rng_seed)
      mcs_fn = da_mcs_unif if self.method == "da_mcs_unif" else da_mcs_marg
      mcs = mcs_fn(nll, alpha=self.alpha, randomize=True, rng=mcs_rng)
      survivors = mcs["decision"]

    if sum(survivors) == 0:
      self._estimators = [default]
      self._estimators_logliks = np.array(default_loglik_val).reshape(1, -1)
      self._default_selected = True
      self._selection_rate = 0.0
      if self.verbose:
        warnings.warn(
          "No random estimator improved over the default; "
          f"default loglik {default_loglik_val:.4f}, "
          f"max random loglik {max(random_logliks_val):.4f}.",
          UserWarning,
          stacklevel=2,
        )
    else:
      self._estimators = [
        estimator
        for (estimator, _), selected in zip(results, survivors)
        if selected
      ]
      self._estimators_logliks = random_logliks_val[survivors].reshape(-1, 1)
      self._selection_rate = sum(survivors) / len(survivors)
      self._default_selected = (
        bool(survivors[-1]) if self.add_dissmann else False
      )
      self._default_loglik = default_loglik_val
      self._default_estimator = default

    if self.best_only and (len(self._estimators) > 1):
      best_idx = int(np.argmax(self._estimators_logliks))
      self._estimators = [self._estimators[best_idx]]
      self._estimators_logliks = self._estimators_logliks[best_idx].reshape(
        1, -1
      )

    # Refit survivors on the full training data.
    outer_n_jobs = self._adjust_estimators_num_threads(self._estimators)
    self._estimators = Parallel(n_jobs=outer_n_jobs)(
      delayed(lambda est: self._fit_estimator(est, X, y))(est)
      for est in self._estimators
    )
    self._seeds = seeds
    return self

  def _adjust_estimators_num_threads(self, estimators: list) -> int:
    """Allocate ``n_jobs`` across the outer (estimator) and inner (thread)
    loops. Routes through each estimator's backend via
    ``with_num_threads``: the C++ backend rebinds
    ``controls.num_threads``; the torch backend treats it as a no-op."""
    n_estimators = len(estimators)
    outer_n_jobs = min(self.n_jobs, n_estimators)
    if n_estimators < self.n_jobs:
      num_threads = max(1, self.n_jobs // n_estimators)
      from .backends import resolve_backend

      for estimator in estimators:
        target = (
          estimator.backend_
          if hasattr(estimator, "backend_")
          else resolve_backend(estimator.backend)
        )
        estimator.backend = target.with_num_threads(num_threads)
        if hasattr(estimator, "backend_"):
          # Keep the pinned post-fit backend in sync so subsequent
          # pdf/cdf calls observe the new thread count.
          estimator.backend_ = resolve_backend(estimator.backend)
    return outer_n_jobs

  def _prepare_prediction_data(self, X: np.ndarray) -> np.ndarray:
    check_is_fitted(self, attributes=["_estimators"])
    result = self._base_estimator._validate_input(X, reset=False)
    return np.asarray(result)
