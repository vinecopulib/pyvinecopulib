"""OpenTURNS distributions as univariate margins.

OpenTURNS spells the univariate surface its own way -- ``computePDF`` /
``computeCDF`` / ``computeQuantile``, over ``Point`` and ``Sample`` rather than
NumPy arrays -- so it needs an adapter of its own rather than one of the SciPy
ones. Three things ship here: that adapter, registered on import so
:func:`as_margin` accepts an OpenTURNS distribution; :class:`OpenTURNSMargin`,
one family fitted by a ``DistributionFactory``; and
:class:`OpenTURNSSelector`, a choice among several by an information criterion.

OpenTURNS is an optional dependency: it is imported when one of these is
constructed, not when this module is imported, and the adapter's predicate
never imports it at all, so :mod:`pyvinecopulib.margins` keeps working without
it.
"""

from __future__ import annotations

import time
import warnings
from typing import Any, Optional

import numpy as np

from ..core import Kde1d, MarginBase, MarginLike
from ._adapters import register_margin_adapter
from .selection import _CRITERIA, _criteria

__all__ = ["OpenTURNSMargin", "OpenTURNSSelector"]

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
  pyvinecopulib.margins.OpenTURNSSelector : Choose among several of these.
  pyvinecopulib.margins.ParametricMargin : The same thing over SciPy.

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
    if (factory is None) == (distribution is None):
      raise ValueError(
        "OpenTURNSMargin takes either a factory to estimate a family with or "
        "a distribution that already carries its parameters, not both and not "
        "neither"
      )
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
    self._loglik: Optional[float] = None

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
        f"OpenTURNSMargin({self.family_name!r}) is not fitted; call fit(x) or "
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
    """
    if self._distribution is not None:
      return str(self._distribution.getName())
    return _factory_name(self._factory)

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
  def fitted_loglik(self) -> float:
    """Log-likelihood attained at the fit.

    Returns
    -------
    float
        The maximized log-likelihood.

    Raises
    ------
    RuntimeError
        If the margin was never fitted. Use ``logpdf(x).sum()`` for the
        log-likelihood on an arbitrary sample.
    """
    if self._loglik is None:
      raise RuntimeError(
        f"OpenTURNSMargin({self.family_name!r}).fitted_loglik is only defined "
        "after fit(x); use logpdf(x).sum() for another sample"
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

  def fit(self, x: Any, *, weights: Optional[Any] = None) -> "OpenTURNSMargin":
    """Estimate the family's parameters with its OpenTURNS factory.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations; NaNs are dropped.
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
    if weights is not None:
      raise TypeError(
        f"OpenTURNSMargin({self.family_name!r}) cannot use observation "
        "weights: an OpenTURNS Sample carries none. Pass margins='kde' or "
        "margins='empirical' for a weighted fit, or drop weights="
      )
    if self._factory is None:
      raise TypeError(
        f"OpenTURNSMargin({self.family_name!r}) was built from a distribution "
        "that already carries its parameters; construct it with a factory to "
        "estimate one from data"
      )
    openturns = _openturns()
    data = np.asarray(x, dtype=float).ravel()
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
    return self

  # --- evaluation ---------------------------------------------------------- #

  def pdf(self, x: Any) -> np.ndarray:
    return _at_points(self.distribution.computePDF, x)

  def logpdf(self, x: Any) -> np.ndarray:
    """Log of :meth:`pdf`, from OpenTURNS' own log-density.

    Overrides the inherited ``log(pdf(x))``, which loses the tails to underflow
    exactly where a log-likelihood is most sensitive to them.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Log-density, ``-inf`` outside the support.
    """
    return _at_points(self.distribution.computeLogPDF, x)

  def cdf(self, x: Any) -> np.ndarray:
    return _at_points(self.distribution.computeCDF, x)

  def icdf(self, p: Any) -> np.ndarray:
    return _at_probabilities(self.distribution, p)

  def cdf_left(self, x: Any) -> np.ndarray:
    if self._var_type == "c":
      return self.cdf(x)
    return _left_limit(self.distribution, x)

  def _sample_uniform(self, n: int, seeds: list[int]) -> np.ndarray:
    return np.random.default_rng(seeds or None).random(n)

  def __repr__(self) -> str:
    if self._distribution is None:
      return f"OpenTURNSMargin({self.family_name!r}, unfitted)"
    shown = ", ".join(
      f"{name}={value:.4g}"
      for name, value in zip(self.parameter_names, self.parameters)
    )
    return f"OpenTURNSMargin({self.family_name!r}, {shown})"


# --- family selection ------------------------------------------------------ #


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
    return dict.fromkeys(_CRITERIA, float("inf"))
  fitting = _openturns().FittingTest
  # OpenTURNS reports its criteria per observation; rescaling to the usual total
  # puts a row on the same scale as `MarginSelector`'s.
  n = float(sample.getSize())
  return {
    key: n * float(getattr(fitting, method)(sample, distribution, k))
    for key, method in _FITTING_TEST.items()
  }


class OpenTURNSSelector(MarginBase[np.ndarray]):
  """Pick an OpenTURNS family by an information criterion.

  A selector is itself a margin: :meth:`fit` estimates every admissible
  candidate, keeps the best one on :attr:`selected_`, and forwards ``pdf`` /
  ``cdf`` / ``icdf`` to it, so a selected margin drops into a
  :class:`~pyvinecopulib.core.Vinedist` exactly like a fixed one. The
  per-candidate table is on :attr:`report_`.

  Parameters
  ----------
  candidates : str or sequence, optional
      ``"auto"`` (the default) takes OpenTURNS' own univariate factory list for
      the variable's type. A sequence selects explicitly, and may mix family
      names with ``DistributionFactory`` objects.
  criterion : {"aic", "bic", "aicc"}, optional
      What to minimize.
  var_type : {"c", "d"} or None, optional
      Force the data to be read as continuous or as discrete. Inferred when
      ``None``: integer-valued data are discrete.
  name : str or None, optional
      Name of the variable being fitted, reported in :attr:`report_` and in the
      fallback warning.

  Raises
  ------
  ImportError
      If OpenTURNS is not installed.
  ValueError
      If ``criterion`` is not one of the three, or ``var_type`` is not ``"c"``
      or ``"d"``.

  See Also
  --------
  pyvinecopulib.margins.OpenTURNSMargin : The candidates.
  pyvinecopulib.margins.MarginSelector : The same thing over SciPy.
  pyvinecopulib.core.Kde1d : The fallback, and the library default.

  Notes
  -----
  **Candidates are scored one at a time, not through
  ``FittingTest.BestModelBIC``.** That helper takes a single flat factory list,
  and drops without raising the candidates it cannot build -- including ones
  that would have won -- while happily comparing a continuous density's
  likelihood against discrete masses if the list mixes the two. Here the
  candidates are partitioned by ``isDiscrete()`` so that only the group
  matching the variable's type competes, and each one is scored on its own.

  **The criteria are OpenTURNS' own**, from ``FittingTest.AIC`` / ``BIC`` /
  ``AICC``, rescaled from per-observation to the usual total so a row sits
  alongside :attr:`MarginSelector.report_`'s.

  **Nothing is skipped silently.** A candidate that fails to build, that lands
  on the wrong side of the discrete / continuous split, or that reaches a
  non-finite log-likelihood still gets a row in :attr:`report_` with the
  reason. When every candidate fails, the selector falls back to
  :class:`~pyvinecopulib.core.Kde1d` and warns once, naming the variable and
  the reasons -- never to a normal, which would misspecify silently.

  **Observation weights are not supported**, because the candidates cannot use
  them; see the notes on :class:`OpenTURNSMargin`.

  Examples
  --------
  ::

      import numpy as np
      from pyvinecopulib.margins import OpenTURNSSelector

      x = np.random.default_rng(0).gamma(2.5, 1.5, size=500)
      sel = OpenTURNSSelector(criterion="bic").fit(x)
      sel.selected_.family_name
      [(r["family"], r["status"]) for r in sel.report_]
  """

  supports_weights: bool = False

  def __init__(
    self,
    candidates: Any = "auto",
    *,
    criterion: str = "aic",
    var_type: Optional[str] = None,
    name: Optional[str] = None,
  ) -> None:
    _openturns()
    if criterion not in _CRITERIA:
      raise ValueError(
        f"unknown criterion={criterion!r}; expected one of {list(_CRITERIA)}"
      )
    if var_type not in (None, "c", "d"):
      raise ValueError(
        f"unknown var_type={var_type!r}; expected 'c', 'd' or None"
      )
    self.candidates = candidates
    self.criterion = criterion
    self._forced_var_type = var_type
    self.name = name
    self._selected: Optional[Any] = None
    self._report: list[dict[str, Any]] = []

  # --- fitted state -------------------------------------------------------- #

  @property
  def selected_(self) -> Any:
    """The winning margin.

    Returns
    -------
    MarginLike
        The candidate with the smallest ``criterion``, or the
        :class:`~pyvinecopulib.core.Kde1d` fallback when every candidate was
        rejected.

    Raises
    ------
    RuntimeError
        If the selector has not been fitted.
    """
    if self._selected is None:
      raise RuntimeError("OpenTURNSSelector is not fitted; call fit(x)")
    return self._selected

  @property
  def report_(self) -> list[dict[str, Any]]:
    """One row per candidate, in the order they were tried.

    Rows carry the same keys as :attr:`MarginSelector.report_`'s, so tables
    from the two selectors concatenate.

    Returns
    -------
    list of dict
        The table; empty before :meth:`fit`.
    """
    return list(self._report)

  @property
  def is_fitted(self) -> bool:
    return self._selected is not None

  # --- estimation ---------------------------------------------------------- #

  def fit(
    self, x: Any, *, weights: Optional[Any] = None
  ) -> "OpenTURNSSelector":
    """Fit every candidate and keep the best.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations; NaNs are dropped.
    weights : array, shape (n,), or None, optional
        Not supported; the candidates cannot use them.

    Returns
    -------
    OpenTURNSSelector
        ``self``.

    Raises
    ------
    TypeError
        If ``weights`` is given.
    ValueError
        If no observation survives.

    Warns
    -----
    UserWarning
        Once, when every candidate is rejected and the kernel-density fallback
        is used instead.
    """
    if weights is not None:
      raise TypeError(
        "OpenTURNSSelector cannot use observation weights: an OpenTURNS "
        "Sample carries none. Pass margins='kde' or margins='empirical' for a "
        "weighted fit, or drop weights="
      )
    data = np.asarray(x, dtype=float).ravel()
    data = data[~np.isnan(data)]
    if data.size == 0:
      raise ValueError("OpenTURNSSelector.fit got no usable observation")

    discrete = self._reads_as_discrete(data)
    sample = _openturns().Sample(data.reshape(-1, 1))
    rows: list[dict[str, Any]] = []
    admissible: list[tuple[float, int, Any]] = []
    for index, factory in enumerate(self._candidates(discrete=discrete)):
      row, candidate = self._try(factory, data, sample, discrete=discrete)
      rows.append(row)
      if candidate is not None:
        # The index breaks ties toward the earlier candidate, which is what
        # keeps the outcome independent of OpenTURNS' list order changing.
        admissible.append((row[self.criterion], index, candidate))

    if admissible:
      _, index, winner = min(admissible, key=lambda item: item[:2])
      rows[index]["status"] = "selected"
      rows[index]["selected"] = True
    else:
      winner = self._fallback(data, rows, discrete=discrete)
    self._selected = winner
    self._report = rows
    return self

  def _reads_as_discrete(self, data: np.ndarray) -> bool:
    """Whether the data should be fitted with the discrete families.

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
    if self._forced_var_type is not None:
      return self._forced_var_type == "d"
    return bool(np.all(data == np.round(data)))

  def _candidates(self, *, discrete: bool) -> list[Any]:
    """Resolve the factories to try.

    Parameters
    ----------
    discrete : bool
        Whether the discrete group applies.

    Returns
    -------
    list
        OpenTURNS factories.

    Raises
    ------
    ValueError
        If ``candidates`` is a string other than ``"auto"``.
    """
    if isinstance(self.candidates, str):
      if self.candidates != "auto":
        raise ValueError(
          f"unknown candidates={self.candidates!r}; expected 'auto' or a "
          "sequence of family names and OpenTURNS factories"
        )
      registry = _openturns().DistributionFactory
      return list(
        registry.GetDiscreteUniVariateFactories()
        if discrete
        else registry.GetContinuousUniVariateFactories()
      )
    return [_resolve_factory(entry) for entry in self.candidates]

  def _try(
    self, factory: Any, data: np.ndarray, sample: Any, *, discrete: bool
  ) -> tuple[dict[str, Any], Optional[OpenTURNSMargin]]:
    """Fit one candidate and score it.

    Parameters
    ----------
    factory : openturns.DistributionFactory
        The candidate's factory.
    data : array, shape (n,), dtype float
        The observations.
    sample : openturns.Sample
        The same observations, as the ``(n, 1)`` sample the criteria read.
    discrete : bool
        Whether the data were read as discrete.

    Returns
    -------
    tuple
        The report row, and the fitted candidate when it is admissible.
    """
    family = _factory_name(factory)
    candidate: Optional[OpenTURNSMargin] = None
    status: Optional[str] = None
    started = time.perf_counter()
    # A factory's optimizer chatter says nothing a caller can act on -- the
    # outcome is the row below either way -- so collect it into the row rather
    # than emitting one warning per family tried.
    with warnings.catch_warnings(record=True) as caught:
      warnings.simplefilter("always")
      try:
        candidate = OpenTURNSMargin(factory).fit(data)
      except Exception as e:  # noqa: BLE001 - any build failure is a rejection
        status = f"{type(e).__name__}: {e}"
    elapsed = time.perf_counter() - started

    loglik, k, support = float("nan"), float("nan"), None
    criteria = dict.fromkeys(_CRITERIA, float("inf"))
    if candidate is not None:
      family = candidate.family_name
      loglik = candidate.fitted_loglik
      k = candidate.n_parameters
      support = candidate.support
      criteria = _openturns_criteria(
        sample, candidate.distribution, int(k), loglik
      )
      if (candidate.var_type == "d") != discrete:
        status = (
          f"a {candidate.var_type!r} family on data read as "
          f"{'d' if discrete else 'c'}; masses and densities are not comparable"
        )
      elif not np.isfinite(loglik):
        status = f"non-finite log-likelihood ({loglik})"
      if status is not None:
        candidate = None

    row: dict[str, Any] = {
      "column": self.name,
      "family": family,
      "n_parameters": k,
      "loglik": loglik,
      "support": support,
      "seconds": elapsed,
      "warnings": [str(w.message) for w in caught],
      "status": "ok" if status is None else status,
      "selected": False,
    }
    row.update(criteria)
    return row, candidate

  def _fallback(
    self, data: np.ndarray, rows: list[dict[str, Any]], *, discrete: bool
  ) -> Any:
    """Fit the kernel-density margin used when no candidate survives.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations.
    rows : list of dict
        The report so far; a row for the fallback is appended.
    discrete : bool
        Whether the data were read as discrete.

    Returns
    -------
    Kde1d
        The fitted fallback.
    """
    reasons = "; ".join(f"{row['family']}: {row['status']}" for row in rows)
    where = f"variable {self.name!r}" if self.name else "this variable"
    warnings.warn(
      f"no OpenTURNS family was admissible for {where}; falling back to a "
      f"kernel-density margin. Rejected: {reasons or 'no candidates'}",
      UserWarning,
      stacklevel=3,
    )
    margin = Kde1d(type="discrete" if discrete else "continuous")
    margin.fit(data)
    loglik = float(np.sum(margin.logpdf(data)))
    row: dict[str, Any] = {
      "column": self.name,
      "family": margin.family_name,
      "n_parameters": margin.n_parameters,
      "loglik": loglik,
      "support": margin.support,
      "seconds": float("nan"),
      "warnings": [],
      "status": "fallback",
      "selected": True,
    }
    row.update(_criteria(loglik, margin.n_parameters, data.size))
    rows.append(row)
    return margin

  # --- forwarding ---------------------------------------------------------- #

  @property
  def family_name(self) -> str:
    """Name of the selected family.

    Returns
    -------
    str
        The winner's own ``family_name``.
    """
    return str(getattr(self.selected_, "family_name", "unknown"))

  @property
  def n_parameters(self) -> float:
    """Number of freely estimated parameters of the selected margin.

    The selection itself is not counted, so an information criterion computed
    from this understates the model's complexity -- the usual post-selection
    caveat, not a property of this implementation.

    Returns
    -------
    float
        The winner's count.
    """
    return float(getattr(self.selected_, "n_parameters", float("nan")))

  @property
  def fitted_loglik(self) -> float:
    """Log-likelihood attained by the selected margin.

    Returns
    -------
    float
        The winner's fitted log-likelihood.
    """
    return float(self.selected_.fitted_loglik)

  @property
  def var_type(self) -> str:
    return str(getattr(self.selected_, "var_type", "c"))

  @property
  def support(self) -> tuple[float, float]:
    lo, hi = getattr(self.selected_, "support", (float("-inf"), float("inf")))
    return (float(lo), float(hi))

  def pdf(self, x: Any) -> np.ndarray:
    return np.asarray(self.selected_.pdf(x), dtype=float)

  def cdf(self, x: Any) -> np.ndarray:
    return np.asarray(self.selected_.cdf(x), dtype=float)

  def icdf(self, p: Any) -> np.ndarray:
    return np.asarray(self.selected_.icdf(p), dtype=float)

  def logpdf(self, x: Any) -> np.ndarray:
    return np.asarray(self.selected_.logpdf(x), dtype=float)

  def cdf_left(self, x: Any) -> np.ndarray:
    return np.asarray(self.selected_.cdf_left(x), dtype=float)

  def sample(self, n: int, *, seeds: Optional[list[int]] = None) -> np.ndarray:
    """Draw ``n`` samples from the selected margin.

    Parameters
    ----------
    n : int
        Number of samples.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n,), dtype float
        Samples on the original scale.
    """
    return np.asarray(self.selected_.sample(n, seeds=seeds), dtype=float)

  def __repr__(self) -> str:
    if self._selected is None:
      return f"OpenTURNSSelector(criterion={self.criterion!r}, unfitted)"
    return (
      f"OpenTURNSSelector(criterion={self.criterion!r}, "
      f"selected={self.family_name!r})"
    )


# --- coercion -------------------------------------------------------------- #


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
