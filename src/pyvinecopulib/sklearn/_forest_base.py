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
from typing import Any

import numpy as np
import pyvinecopulib as pv
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator
from sklearn.model_selection import train_test_split

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
   the forest paper from the discrete-argmin-inference framework of
   Kim & Ramdas (2025) — to compute the model confidence set (Hansen,
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
  """

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
    seed: int = 42,
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
    seed : int, default=42
        Random seed for reproducibility.
    n_jobs : int, default=1
        Number of joblib workers for parallel fit / predict.
    verbose : bool, default=False
        Emit a warning if no random estimator beats the default.
    """
    if not isinstance(base_class, type):
      raise TypeError(
        f"base_class must be a class, got {type(base_class).__name__}"
      )
    if base_params is not None and not isinstance(base_params, dict):
      raise TypeError(
        f"base_params must be dict or None, got {type(base_params).__name__}"
      )
    if not isinstance(n_vines, int) or isinstance(n_vines, bool):
      raise TypeError(f"n_vines must be int, got {type(n_vines).__name__}")
    if n_vines < 1:
      raise ValueError(f"n_vines must be >= 1, got {n_vines}")
    if vines_sampling not in {"uniform", "local"}:
      raise ValueError(
        f"vines_sampling must be 'uniform' or 'local', got {vines_sampling!r}"
      )
    if not isinstance(bootstrap, bool):
      raise TypeError(f"bootstrap must be bool, got {type(bootstrap).__name__}")
    if not isinstance(val_fraction, (int, float)) or isinstance(
      val_fraction, bool
    ):
      raise TypeError(
        f"val_fraction must be float, got {type(val_fraction).__name__}"
      )
    if not (0.0 <= val_fraction < 1.0):
      raise ValueError(f"val_fraction must be in [0, 1), got {val_fraction}")
    if not isinstance(best_only, bool):
      raise TypeError(f"best_only must be bool, got {type(best_only).__name__}")
    if method not in _VALID_METHODS:
      raise ValueError(
        f"method must be one of {sorted(str(m) for m in _VALID_METHODS)}, "
        f"got {method!r}"
      )
    if not isinstance(alpha, (int, float)) or isinstance(alpha, bool):
      raise TypeError(f"alpha must be float, got {type(alpha).__name__}")
    if not (0.0 < alpha < 1.0):
      raise ValueError(f"alpha must be in (0, 1), got {alpha}")
    if not isinstance(add_dissmann, bool):
      raise TypeError(
        f"add_dissmann must be bool, got {type(add_dissmann).__name__}"
      )
    if not isinstance(seed, int) or isinstance(seed, bool):
      raise TypeError(f"seed must be int, got {type(seed).__name__}")
    if not isinstance(n_jobs, int) or isinstance(n_jobs, bool):
      raise TypeError(f"n_jobs must be int, got {type(n_jobs).__name__}")
    if n_jobs == 0 or n_jobs < -1:
      raise ValueError(f"n_jobs must be -1 or >= 1, got {n_jobs}")
    if not isinstance(verbose, bool):
      raise TypeError(f"verbose must be bool, got {type(verbose).__name__}")

    self.base_class = base_class
    self.base_params = {} if base_params is None else base_params
    self.n_vines = n_vines
    self.vines_sampling = vines_sampling
    self.bootstrap = bootstrap
    self.val_fraction = val_fraction
    self.best_only = best_only
    self.seed = seed
    self.n_jobs = n_jobs
    self.verbose = verbose
    self.method = method
    self.alpha = alpha
    self.add_dissmann = add_dissmann

  @abstractmethod
  def _create_base_estimator(self):
    """Build a fresh base estimator with any forest-specific overrides applied.

    Subclasses use this in place of :func:`sklearn.base.clone` so
    that attribute-set overrides (e.g. ``_normalize_weights = False``
    on the regressor) survive the round-trip — ``clone`` would
    re-instantiate from ``__init__`` parameters only and reset such
    attributes to their defaults.
    """

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
          X, y, test_size=self.val_fraction, random_state=self.seed
        )
      else:
        X_train, X_val = train_test_split(
          X, test_size=self.val_fraction, random_state=self.seed
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
    # `_create_base_estimator()` (not clone) so any attribute-set
    # overrides like `_normalize_weights = False` are re-applied.
    estimator = self._create_base_estimator()

    # Propagate the schema inferred on the full input (in particular,
    # discrete-vs-continuous flags from DataFrame columns) — by this
    # point X_train is a pre-expanded numpy array, so the base
    # estimator's own inference would default everything to continuous.
    if self._base_estimator._schema is not None:
      estimator._schema = dict(self._base_estimator._schema)

    local_rng = np.random.default_rng(seed)
    local_seeds = [int(x) for x in local_rng.integers(0, 2**31 - 1, size=5)]

    if self.vines_sampling == "uniform":
      structure_size = X_train.shape[1] + int(y_train is not None)
      estimator.structure = pv.RVineStructure.simulate(
        structure_size, seeds=local_seeds
      )
    elif self.vines_sampling == "local":
      estimator.controls.tree_algorithm = "random_weighted"
      estimator.controls.seeds = local_seeds
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
    self._base_estimator = self._create_base_estimator()

    if y is not None:
      X, _ = self._base_estimator._check_and_expand_fit(X, y)
    else:
      X = self._base_estimator._check_and_expand_fit(X)
    self.n_features_in_ = X.shape[1]

    X_train, X_val, y_train, y_val = self._split_data(X, y)

    # Default estimator (Dissmann-structure baseline). Propagate the
    # inferred schema so fitting on the already-expanded numpy doesn't
    # default everything back to continuous.
    default = self._create_base_estimator()
    if self._base_estimator._schema is not None:
      default._schema = dict(self._base_estimator._schema)
    if y_train is not None:
      default.fit(X_train, y_train)
    else:
      default.fit(X_train)
    default_loglik = self._loglik_estimator(default, X_val, y_val)
    default_loglik_val = default_loglik.mean()

    rng = np.random.default_rng(self.seed)
    seeds = rng.integers(0, 2**31 - 1, size=self.n_vines)

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
      mcs_rng = np.random.default_rng(self.seed)
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
    """Allocate `n_jobs` across the outer (estimator) and inner (thread) loops."""
    n_estimators = len(estimators)
    outer_n_jobs = min(self.n_jobs, n_estimators)
    if n_estimators < self.n_jobs:
      num_threads = max(1, self.n_jobs // n_estimators)
      for estimator in estimators:
        estimator.controls.num_threads = num_threads
    return outer_n_jobs

  def _check_fitted(self) -> None:
    if not hasattr(self, "_estimators"):
      raise ValueError(f"{self.__class__.__name__} not fitted yet.")

  def _prepare_prediction_data(self, X: np.ndarray) -> np.ndarray:
    self._check_fitted()
    result = self._base_estimator._check_and_expand_predict(X)
    return np.asarray(result)
