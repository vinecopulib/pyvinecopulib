"""OpenTURNS distributions as univariate margins.

OpenTURNS spells the univariate surface its own way -- ``computePDF`` /
``computeCDF`` / ``computeQuantile``, over ``Point`` and ``Sample`` rather than
NumPy arrays -- so it needs an adapter of its own rather than one of the SciPy
ones. Two things ship here: that adapter, registered on import so
:func:`as_margin` accepts an OpenTURNS distribution; and
:class:`OpenTURNSMargin`, one family fitted by a ``DistributionFactory`` --
or, left unnamed, whichever family of the registry an information criterion
picks.

OpenTURNS is an optional dependency: it is imported when one of these is
constructed, not when this module is imported, and the adapter's predicate
never imports it at all, so :mod:`pyvinecopulib.margins` keeps working without
it.
"""

from __future__ import annotations

import warnings
from typing import Any, Optional

import numpy as np

from ..core import MarginBase, MarginLike
from ..core.margin_base import _reject_covariates
from ..core._validation import validate_univariate
from ._adapters import register_margin_adapter
from .controls import CRITERIA, FitControlsMargin

__all__ = ["OpenTURNSMargin"]

#: The ``openturns.FittingTest`` entry point behind each criterion name. All
#: three are minimized, and all three are reported per observation.
_FITTING_TEST: dict[str, str] = {"aic": "AIC", "bic": "BIC", "aicc": "AICC"}


def _openturns() -> Any:
  """Import ``openturns``, or explain which extra provides it.

  Returns
  -------
  module
      The ``openturns`` module.

  Raises
  ------
  ImportError
      If OpenTURNS is not installed.
  """
  try:
    import openturns
  except ImportError as e:  # pragma: no cover - exercised in a subprocess
    raise ImportError(
      "pyvinecopulib's OpenTURNS margins require OpenTURNS. Install it with "
      "`pip install pyvinecopulib[openturns]`."
    ) from e
  return openturns


# --- marshaling ----------------------------------------------------------- #


def _at_points(method: Any, x: Any) -> np.ndarray:
  """Evaluate a per-point OpenTURNS method over a NumPy array.

  Parameters
  ----------
  method : callable
      A bound ``computePDF``, ``computeCDF`` or ``computeLogPDF``.
  x : array, shape (n,), dtype float
      Points on the original scale.

  Returns
  -------
  array, dtype float
      One value per entry of ``x``, in ``x``'s shape.
  """
  openturns = _openturns()
  values = np.asarray(x, dtype=float)
  if values.size == 0:
    # An empty `Sample` loses its dimension, and OpenTURNS then reads the
    # argument as a 0-dimensional point rather than an empty batch.
    return np.empty(values.shape, dtype=float)
  # A univariate argument crosses as an (n, 1) `Sample`: a flat length-n array
  # is read as one n-dimensional point instead, and raises.
  sample = openturns.Sample(np.atleast_1d(values).reshape(-1, 1))
  out = np.asarray(method(sample), dtype=float).reshape(-1)
  return out.reshape(values.shape)


def _at_probabilities(distribution: Any, p: Any) -> np.ndarray:
  """Evaluate ``computeQuantile`` over a NumPy array of probabilities.

  Parameters
  ----------
  distribution : openturns.Distribution
      The distribution to invert.
  p : array, shape (n,), dtype float
      Probabilities in ``[0, 1]``.

  Returns
  -------
  array, dtype float
      One quantile per entry of ``p``, in ``p``'s shape.
  """
  openturns = _openturns()
  values = np.asarray(p, dtype=float)
  if values.size == 0:
    return np.empty(values.shape, dtype=float)
  # The shape convention inverts here. A probability is a scalar, so the vector
  # form takes a flat `Point` of length n and rejects the (n, 1) array that
  # `computePDF` requires; the result comes back as an (n, 1) `Sample` either
  # way.
  point = openturns.Point(np.atleast_1d(values).reshape(-1))
  out = np.asarray(distribution.computeQuantile(point), dtype=float).reshape(-1)
  return out.reshape(values.shape)


def _left_limit(distribution: Any, x: Any) -> np.ndarray:
  """Left limit ``F(x^-)`` of a discrete OpenTURNS distribution.

  Parameters
  ----------
  distribution : openturns.Distribution
      A discrete distribution.
  x : array, shape (n,), dtype float
      Points on the original scale.

  Returns
  -------
  array, dtype float
      Left limits in ``[0, 1]``.
  """
  # `computeCDF` is a genuine step function and `computePDF` is the mass, so
  # removing one from the other is exact at an atom and a no-op between atoms.
  # Unlike stepping back a lattice point, it does not assume the support is the
  # integers, which a `FiniteDiscreteDistribution`'s is not.
  upper = _at_points(distribution.computeCDF, x)
  return np.clip(upper - _at_points(distribution.computePDF, x), 0.0, 1.0)


def _support(distribution: Any) -> tuple[float, float]:
  """Read a distribution's support off its range.

  Parameters
  ----------
  distribution : openturns.Distribution
      The distribution.

  Returns
  -------
  tuple of float
      ``(lo, hi)``, infinite wherever the range says the bound is.
  """
  interval = distribution.getRange()
  # A bound OpenTURNS flags as infinite still carries a number -- the
  # quantile-based cutoff it truncates plots at -- which is not a support bound.
  lo = (
    float(interval.getLowerBound()[0])
    if interval.getFiniteLowerBound()[0]
    else float("-inf")
  )
  hi = (
    float(interval.getUpperBound()[0])
    if interval.getFiniteUpperBound()[0]
    else float("inf")
  )
  return (lo, hi)


def _univariate(distribution: Any) -> Any:
  """Check that an OpenTURNS distribution can stand in for a margin.

  Parameters
  ----------
  distribution : openturns.Distribution
      The distribution to check.

  Returns
  -------
  openturns.Distribution
      ``distribution`` itself.

  Raises
  ------
  ValueError
      If it is not univariate.
  """
  dimension = int(distribution.getDimension())
  if dimension != 1:
    raise ValueError(
      "a margin is univariate, but this OpenTURNS distribution has "
      f"dimension {dimension}; take one of its marginals with getMarginal(j)"
    )
  return distribution


def _resolve_factory(spec: Any) -> Any:
  """Resolve a factory specification into an OpenTURNS factory.

  Parameters
  ----------
  spec : str or openturns.DistributionFactory
      A factory, or the name of one. ``"Gamma"`` and ``"GammaFactory"`` both
      resolve, since OpenTURNS registers factories under the suffixed name.

  Returns
  -------
  openturns.DistributionFactory
      The factory.

  Raises
  ------
  ValueError
      If ``spec`` names no registered factory.
  """
  openturns = _openturns()
  if not isinstance(spec, str):
    return spec
  name = spec if spec.endswith("Factory") else spec + "Factory"
  try:
    return openturns.DistributionFactory.GetByName(name)
  except TypeError as e:
    raise ValueError(
      f"unknown OpenTURNS factory {spec!r}; expected a distribution name such "
      "as 'Gamma', or an openturns.DistributionFactory"
    ) from e


def _factory_name(factory: Any) -> str:
  """Name the family a factory estimates.

  Parameters
  ----------
  factory : openturns.DistributionFactory
      The factory.

  Returns
  -------
  str
      The family name, with OpenTURNS' ``Factory`` suffix removed.
  """
  # `DistributionFactory` is an interface object wrapping the concrete factory,
  # so its own class name is the wrapper's; a concrete factory has no
  # implementation to unwrap.
  implementation = getattr(factory, "getImplementation", None)
  name = str(
    implementation().getClassName()
    if implementation is not None
    else factory.getClassName()
  )
  return name.removesuffix("Factory")


# --- one family ------------------------------------------------------------ #


class OpenTURNSMargin(MarginBase[np.ndarray]):
  """A univariate margin from one OpenTURNS family.

  Construct it with the family to estimate -- ``OpenTURNSMargin("Gamma")`` --
  and call :meth:`fit`; or wrap a distribution that already carries its
  parameters with :meth:`from_distribution`, which is what :func:`as_margin`
  does for an OpenTURNS object.

  Parameters
  ----------
  factory : str or openturns.DistributionFactory, or None, optional
      The family to estimate, as a factory or the name of one (``"Gamma"`` and
      ``"GammaFactory"`` both resolve). Omit it only alongside
      ``distribution``.
  distribution : openturns.Distribution, or None, optional
      A distribution that already carries its parameters. Given here, the
      margin is already fitted and :attr:`n_parameters` is 0, since nothing was
      estimated from data. :meth:`from_distribution` is the readable spelling.

  Raises
  ------
  ImportError
      If OpenTURNS is not installed.
  ValueError
      If neither or both of ``factory`` and ``distribution`` are given, if
      ``factory`` names no registered factory, or if ``distribution`` is not
      univariate.

  See Also
  --------
  pyvinecopulib.margins.SciPyMargin : The SciPy counterpart.
  pyvinecopulib.margins.SciPyMargin : The same thing over SciPy.

  Notes
  -----
  **A discrete family's left limit does not assume a lattice.** OpenTURNS'
  ``computeCDF`` is a genuine step function and ``computePDF`` is the mass, so
  :meth:`cdf_left` removes one from the other. That is exact at an atom, a
  no-op between atoms, and right for a ``FiniteDiscreteDistribution`` whose
  support is not the integers.

  **Observation weights are not supported.** An OpenTURNS ``Sample`` carries
  none, and a ``DistributionFactory`` estimates from nothing else, so
  :meth:`fit` raises rather than silently returning an unweighted fit.

  Examples
  --------
  ::

      import numpy as np
      import openturns
      from pyvinecopulib.margins import OpenTURNSMargin, as_margin

      x = np.random.default_rng(0).gamma(2.5, 1.5, size=500)
      m = OpenTURNSMargin("Gamma").fit(x)
      m.parameter_names            # -> ('k', 'lambda', 'gamma')
      as_margin(openturns.Normal(0.0, 1.0)).cdf(np.zeros(1))   # -> array([0.5])
  """

  supports_weights: bool = False

  def __init__(
    self,
    factory: Optional[Any] = None,
    *,
    distribution: Optional[Any] = None,
  ) -> None:
    _openturns()
    self._declared_var_type: Optional[str] = None
    if factory is None and distribution is None:
      # Neither: the family is chosen by `select`, from the registry for
      # whichever variable type the data read as.
      self._factory = None
      self._distribution = None
      self._var_type = "c"
      self._n_free = 0
      self._loglik: Optional[float] = None
      self._nobs: Optional[int] = None
      self._unnamed = True
      return
    if factory is not None and distribution is not None:
      raise ValueError(
        "OpenTURNSMargin takes either a factory to estimate a family with or "
        "a distribution that already carries its parameters, not both"
      )
    self._unnamed = False
    # Branch on the arguments rather than on the attributes: the check above
    # guarantees exactly one is set, but that invariant is not expressible in
    # the attribute types.
    known: Any
    if distribution is not None:
      resolved_distribution: Any = _univariate(distribution)
      resolved_factory: Any = None
      known = resolved_distribution
    else:
      resolved_distribution = None
      resolved_factory = _resolve_factory(factory)
      # Whether a family is discrete is a property of the family, not of the
      # fit, so a factory answers it from the default distribution it builds
      # without data.
      known = resolved_factory.build()
    self._factory = resolved_factory
    self._distribution = resolved_distribution
    self._var_type = "d" if known.isDiscrete() else "c"
    self._n_free = 0
    self._loglik = None
    self._nobs = None

  @staticmethod
  def from_distribution(distribution: Any) -> "OpenTURNSMargin":
    """Wrap a distribution that already carries its parameters.

    Parameters
    ----------
    distribution : openturns.Distribution
        A univariate OpenTURNS distribution.

    Returns
    -------
    OpenTURNSMargin
        A fitted margin, with :attr:`n_parameters` 0.
    """
    return OpenTURNSMargin(distribution=distribution)

  # --- identity ------------------------------------------------------------ #

  @property
  def distribution(self) -> Any:
    """The underlying OpenTURNS distribution.

    Returns
    -------
    openturns.Distribution
        The one given at construction, or the one :meth:`fit` built.

    Raises
    ------
    RuntimeError
        If the margin has neither been fitted nor been given a distribution.
    """
    if self._distribution is None:
      raise RuntimeError(
        f"OpenTURNSMargin({self.family_name!r}) is not fitted; call fit(y) or "
        "pass distribution="
      )
    return self._distribution

  @property
  def family_name(self) -> str:
    """Name of the OpenTURNS family.

    Returns
    -------
    str
        The distribution's own name once fitted, and the factory's family name
        before that.

    Raises
    ------
    RuntimeError
        If the family is neither named nor selected yet.
    """
    if self._distribution is not None:
      return str(self._distribution.getName())
    if self._factory is None:
      raise RuntimeError(
        "OpenTURNSMargin() has no family yet; name a factory at construction "
        "or call select(y) to choose from the registry"
      )
    return _factory_name(self._factory)

  @property
  def nobs(self) -> Optional[int]:
    """Number of observations the fit used.

    Returns
    -------
    int or None
        The sample size, or ``None`` on a margin built from a distribution
        that already carried its parameters. It is what penalizes ``bic`` and
        ``aicc``.
    """
    return self._nobs

  @property
  def parameter_names(self) -> tuple[str, ...]:
    """Parameter names in OpenTURNS' order.

    Returns
    -------
    tuple of str
        One name per entry of :attr:`parameters`.
    """
    return tuple(
      str(name) for name in self.distribution.getParameterDescription()
    )

  @property
  def parameters(self) -> tuple[float, ...]:
    """Parameter values in OpenTURNS' order.

    Returns
    -------
    tuple of float
        One value per entry of :attr:`parameter_names`.
    """
    return tuple(float(value) for value in self.distribution.getParameter())

  @property
  def n_parameters(self) -> float:
    """Number of freely estimated parameters.

    Returns
    -------
    float
        What :meth:`fit` estimated; 0 for a margin built from a distribution
        that already carried its parameters.
    """
    return float(self._n_free)

  @property
  def _fitted_loglik(self) -> float:
    """Log-likelihood attained at the fit.

    Returns
    -------
    float
        The maximized log-likelihood.

    Raises
    ------
    RuntimeError
        If the margin was never fitted.
    """
    if self._loglik is None:
      raise RuntimeError(
        f"OpenTURNSMargin({self.family_name!r}).loglik() is only defined after "
        "fit(y); pass a sample to evaluate one instead"
      )
    return self._loglik

  @property
  def is_fitted(self) -> bool:
    return self._distribution is not None

  @property
  def var_type(self) -> str:
    return self._var_type

  @property
  def support(self) -> tuple[float, float]:
    """Closed bounds of the fitted support.

    Returns
    -------
    tuple of float
        ``(lo, hi)`` at the current parameters, and ``(-inf, inf)`` while the
        margin is unfitted, since the support of most families depends on them.
    """
    if self._distribution is None:
      return (float("-inf"), float("inf"))
    return _support(self._distribution)

  # --- estimation ---------------------------------------------------------- #

  def fit(
    self,
    y: Any,
    /,
    controls: Optional[Any] = None,
    *,
    x: Optional[Any] = None,
    weights: Optional[Any] = None,
  ) -> "OpenTURNSMargin":
    """Estimate the family's parameters with its OpenTURNS factory.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations; NaNs are dropped.
    controls : object, or None, optional
        Unused; this margin's family is named at construction, so there is
        nothing left to configure.
    x : array, shape (n, k), or None, optional
        Not supported; passing covariates raises rather than silently
        fitting an unconditional margin.
    weights : array, shape (n,), or None, optional
        Not supported; see the Notes on the class.

    Returns
    -------
    OpenTURNSMargin
        ``self``.

    Raises
    ------
    TypeError
        If ``weights`` is given, or if the margin was built from a
        distribution that already carries its parameters.
    ValueError
        If no observation survives.
    """
    _reject_covariates(self, x)
    if weights is not None:
      raise TypeError(
        f"OpenTURNSMargin({self.family_name!r}) cannot use observation "
        "weights: an OpenTURNS Sample carries none. Pass margins='kde' for a "
        "weighted fit, or drop weights="
      )
    if self._factory is None:
      raise TypeError(
        f"OpenTURNSMargin({self.family_name!r}) was built from a distribution "
        "that already carries its parameters; construct it with a factory to "
        "estimate one from data"
      )
    openturns = _openturns()
    data = validate_univariate(np.asarray(y, dtype=float))
    data = data[~np.isnan(data)]
    if data.size == 0:
      raise ValueError(
        f"OpenTURNSMargin({self.family_name!r}).fit got no usable observation"
      )
    self._distribution = _univariate(
      self._factory.build(openturns.Sample(data.reshape(-1, 1)))
    )
    self._n_free = len(self._distribution.getParameter())
    self._loglik = float(np.sum(self.logpdf(data)))
    self._nobs = int(data.size)
    return self

  # --- selection ----------------------------------------------------------- #

  def declare(
    self,
    *,
    var_type: Optional[str] = None,
    support: Optional[tuple[float, float]] = None,
  ) -> "OpenTURNSMargin":
    """Accept what the caller knows, to be honored by :meth:`select`.

    A named factory already fixes the variable type, so this only steers a
    search: it decides which of OpenTURNS' two registries is searched.

    Parameters
    ----------
    var_type : str or None, optional
        ``"c"``, ``"d"`` or ``"zi"``. ``"zi"`` reduces to ``"d"``, the
        partition OpenTURNS exposes for families with atoms.
    support : tuple of float, or None, optional
        Ignored. The registries are partitioned by variable type rather than
        by support, so bounds do not narrow them.

    Returns
    -------
    OpenTURNSMargin
        ``self``, so the call chains into :meth:`select`.
    """
    del support
    if var_type is not None:
      self._declared_var_type = "d" if var_type == "zi" else var_type
    return self

  def select(
    self,
    y: Any,
    /,
    controls: Optional[Any] = None,
    *,
    x: Optional[Any] = None,
    weights: Optional[Any] = None,
  ) -> "OpenTURNSMargin":
    """Choose a family from OpenTURNS' registry, fit it, and become it.

    Every candidate the registry offers for the variable's type is built and
    scored, and the best on the criterion wins; this margin then *is* that
    family. Which registry is searched follows the data, or
    ``controls.var_type`` when the caller states it: OpenTURNS partitions its
    univariate families into continuous and discrete, and a probability mass
    is not comparable with a Lebesgue density on one criterion.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations; NaNs are dropped.
    controls : FitControlsMargin, or None, optional
        Bounds the search: ``family_set`` names the candidate factories,
        ``selection_criterion`` scores them, and ``var_type`` says which
        registry to search. ``support`` is ignored, as in :meth:`declare`.
    x : array, shape (n, k), or None, optional
        Not supported; passing covariates raises rather than silently
        selecting an unconditional margin.
    weights : array, shape (n,), or None, optional
        Not supported; an OpenTURNS ``Sample`` carries none.

    Returns
    -------
    OpenTURNSMargin
        ``self``, now carrying the chosen family.

    Raises
    ------
    TypeError
        If ``weights`` is given.
    ValueError
        If no observation survives, or if every candidate was refused.

    See Also
    --------
    fit : Estimate a named family, leaving the family alone.

    Examples
    --------
    ::

        import numpy as np
        from pyvinecopulib.margins import FitControlsMargin, OpenTURNSMargin

        y = np.random.default_rng(0).gamma(2.5, 1.5, size=500)
        chosen = OpenTURNSMargin().select(
          y, FitControlsMargin(selection_criterion="bic")
        )
        chosen.family_name
    """
    _reject_covariates(self, x)
    if weights is not None:
      raise TypeError(
        "OpenTURNSMargin cannot use observation weights: an OpenTURNS Sample "
        "carries none. Pass margins='kde' for a weighted fit, or drop weights="
      )
    settings = controls if controls is not None else FitControlsMargin()
    criterion = getattr(settings, "selection_criterion", "aic")
    self.declare(var_type=getattr(settings, "var_type", None))

    openturns = _openturns()
    data = validate_univariate(np.asarray(y, dtype=float))
    data = data[~np.isnan(data)]
    if data.size == 0:
      raise ValueError("OpenTURNSMargin.select got no usable observation")
    sample = openturns.Sample(data.reshape(-1, 1))

    discrete = self._reads_as_discrete(data)
    scored: list[tuple[float, OpenTURNSMargin]] = []
    refused: list[str] = []
    for factory in self._candidate_factories(
      discrete=discrete, family_set=getattr(settings, "family_set", None)
    ):
      candidate, reason = _try_factory(factory, data, sample, discrete=discrete)
      if candidate is None:
        refused.append(f"{_factory_name(factory)}: {reason}")
        continue
      score = _openturns_criteria(
        sample,
        candidate.distribution,
        int(candidate.n_parameters),
        candidate.loglik(),
      )[criterion]
      scored.append((score, candidate))

    if not scored:
      detail = "\n  ".join(refused) or "(the registry offered no candidate)"
      raise ValueError(
        "no OpenTURNS family fits this variable; every candidate was "
        f"refused:\n  {detail}\nName a factory directly, or narrow "
        "family_set."
      )
    scored.sort(key=lambda pair: pair[0])
    return self._adopt(scored[0][1])

  def _reads_as_discrete(self, data: np.ndarray) -> bool:
    """Whether the discrete registry is the candidate set.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations.

    Returns
    -------
    bool
        ``True`` for the discrete group. Integer-valued data qualify whatever
        their sign, since OpenTURNS' discrete families include ``Skellam`` on
        all of the integers.
    """
    if self._declared_var_type is not None:
      return self._declared_var_type == "d"
    return bool(np.all(data == np.round(data)))

  @staticmethod
  def _candidate_factories(
    *, discrete: bool, family_set: Optional[Any]
  ) -> list[Any]:
    """Resolve the factories to try.

    Parameters
    ----------
    discrete : bool
        Whether the discrete group applies.
    family_set : sequence, or None
        Family names or OpenTURNS factories, or ``None`` for the whole
        registry of the applicable type.

    Returns
    -------
    list
        OpenTURNS factories.
    """
    if family_set is None:
      registry = _openturns().DistributionFactory
      return list(
        registry.GetDiscreteUniVariateFactories()
        if discrete
        else registry.GetContinuousUniVariateFactories()
      )
    return [_resolve_factory(entry) for entry in family_set]

  def _adopt(self, winner: "OpenTURNSMargin") -> "OpenTURNSMargin":
    """Become the selected candidate.

    Parameters
    ----------
    winner : OpenTURNSMargin
        The candidate that won.

    Returns
    -------
    OpenTURNSMargin
        ``self``.
    """
    self._factory = winner._factory
    self._distribution = winner._distribution
    self._var_type = winner._var_type
    self._n_free = winner._n_free
    self._loglik = winner._loglik
    self._nobs = winner._nobs
    self._unnamed = False
    return self

  # --- evaluation ---------------------------------------------------------- #

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return _at_points(self.distribution.computePDF, y)

  def logpdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    """Log of :meth:`pdf`, from OpenTURNS' own log-density.

    Overrides the inherited ``log(pdf(y))``, which loses the tails to underflow
    exactly where a log-likelihood is most sensitive to them.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Ignored; this margin is unconditional.

    Returns
    -------
    array, shape (n,), dtype float
        Log-density, ``-inf`` outside the support.
    """
    return _at_points(self.distribution.computeLogPDF, y)

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return _at_points(self.distribution.computeCDF, y)

  def icdf(self, p: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return _at_probabilities(self.distribution, p)

  def cdf_left(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    if self._var_type == "c":
      return self.cdf(y)
    return _left_limit(self.distribution, y)

  def _sample_uniform(self, n: int, seeds: list[int]) -> np.ndarray:
    return np.random.default_rng(seeds or None).random(n)

  def __repr__(self) -> str:
    if self._factory is None and self._distribution is None:
      return "OpenTURNSMargin(no family yet)"
    if self._distribution is None:
      return f"OpenTURNSMargin({self.family_name!r}, unfitted)"
    shown = ", ".join(
      f"{name}={value:.4g}"
      for name, value in zip(self.parameter_names, self.parameters)
    )
    return f"OpenTURNSMargin({self.family_name!r}, {shown})"


# --- family selection ------------------------------------------------------ #


def _try_factory(
  factory: Any,
  data: np.ndarray,
  sample: Any,
  *,
  discrete: bool,
) -> tuple[Optional["OpenTURNSMargin"], Optional[str]]:
  """Build one candidate and report why it is inadmissible, or ``None``.

  Parameters
  ----------
  factory : openturns.DistributionFactory
      The candidate family.
  data : array, shape (n,), dtype float
      The observations.
  sample : openturns.Sample
      The same observations, shape ``(n, 1)``.
  discrete : bool
      Whether the data read as discrete.

  Returns
  -------
  tuple
      The fitted candidate and ``None``, or ``None`` and a one-line reason.
  """
  del sample
  # A factory's optimizer chatter says nothing a caller can act on -- the
  # outcome is admissible or not either way -- so it is collected rather than
  # emitted once per family tried.
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
      candidate = OpenTURNSMargin(factory).fit(data)
    except Exception as e:  # noqa: BLE001 - any build failure is a rejection
      return None, f"{type(e).__name__}: {e}"
  if (candidate.var_type == "d") != discrete:
    return None, (
      f"a {candidate.var_type!r} family on data read as "
      f"{'d' if discrete else 'c'}; masses and densities are not comparable"
    )
  loglik = candidate.loglik()
  if not np.isfinite(loglik):
    return None, f"non-finite log-likelihood ({loglik})"
  return candidate, None


def _openturns_criteria(
  sample: Any, distribution: Any, k: int, loglik: float
) -> dict[str, float]:
  """Evaluate every criterion at one fit, through ``openturns.FittingTest``.

  Parameters
  ----------
  sample : openturns.Sample
      The observations, shape ``(n, 1)``.
  distribution : openturns.Distribution
      The fitted candidate.
  k : int
      Number of freely estimated parameters.
  loglik : float
      Log-likelihood attained, used only to reject a degenerate fit.

  Returns
  -------
  dict
      One entry per criterion name; ``inf`` where the log-likelihood is not
      finite, so a candidate can never win by being undefined.
  """
  if not np.isfinite(loglik):
    return dict.fromkeys(CRITERIA, float("inf"))
  fitting = _openturns().FittingTest
  # OpenTURNS reports its criteria per observation; rescaling to the usual total
  # puts a criterion on the usual scale for a whole sample.
  n = float(sample.getSize())
  return {
    key: n * float(getattr(fitting, method)(sample, distribution, k))
    for key, method in _FITTING_TEST.items()
  }


def _is_openturns_distribution(obj: Any) -> bool:
  """Whether ``obj`` is an ``openturns`` distribution.

  Parameters
  ----------
  obj : object
      Any object.

  Returns
  -------
  bool
      ``True`` for a concrete OpenTURNS distribution and for the
      ``Distribution`` interface object a factory returns, which share no base
      class beyond ``Object``.
  """
  return any(
    str(getattr(base, "__module__", "")).startswith("openturns")
    and base.__name__ in ("Distribution", "DistributionImplementation")
    for base in type(obj).__mro__
  )


def _adapt_openturns(obj: Any) -> MarginLike[Any]:
  """Adapt an ``openturns`` distribution.

  Parameters
  ----------
  obj : openturns.Distribution
      The distribution to wrap.

  Returns
  -------
  MarginLike
      An :class:`OpenTURNSMargin` around it, already fitted.
  """
  return OpenTURNSMargin.from_distribution(obj)


register_margin_adapter(_is_openturns_distribution, _adapt_openturns)
