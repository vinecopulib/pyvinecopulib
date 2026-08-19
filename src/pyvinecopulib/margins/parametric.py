"""Parametric margin: one SciPy family, fitted by maximum likelihood.

SciPy is an optional dependency; it is imported when a
``ParametricMargin`` is constructed, not when this module is imported, so
:mod:`pyvinecopulib.margins` keeps working without it.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

import numpy as np

from ..core import MarginBase
from ..core.margin_base import _reject_covariates

__all__ = ["ParametricMargin"]

#: Curated candidate families, grouped by the support they can represent.
#: ``MarginSelector`` draws its ``"auto"`` candidates from the groups the data
#: are compatible with. Counts are never compared against continuous families:
#: the two use different reference measures, so their log-likelihoods are not on
#: the same scale.
_PARTITIONS: dict[str, tuple[str, ...]] = {
  "real": ("norm", "t", "logistic", "laplace", "skewnorm", "gumbel_r"),
  "positive": ("gamma", "lognorm", "weibull_min", "invgauss", "expon"),
  "unit": ("beta",),
  "bounded": ("uniform", "beta"),
  "count": ("poisson", "nbinom", "geom"),
}

#: Open interval each fitted parameter must lie strictly inside. Two things
#: are caught here at once: a parameter that ran to a bound of its own domain
#: (``p = 1`` for ``nbinom`` is a point mass at 0), and one that ran to a
#: degeneracy threshold well inside it (a Student ``t`` with ``df = 2e10`` is a
#: normal spelled with an extra parameter, and beats it on log-likelihood by
#: rounding error alone). Anything not listed only has to be finite.
_LIMITS: dict[str, dict[str, tuple[float, float]]] = {
  "t": {"df": (0.0, 200.0)},
  "skewnorm": {"a": (-1e3, 1e3)},
  "gamma": {"a": (0.0, 1e4)},
  "lognorm": {"s": (1e-4, 1e3)},
  "weibull_min": {"c": (0.0, 1e3)},
  "invgauss": {"mu": (0.0, 1e4)},
  "beta": {"a": (0.0, 1e4), "b": (0.0, 1e4)},
  "poisson": {"mu": (0.0, 1e4)},
  "nbinom": {"n": (0.0, 1e4), "p": (0.0, 1.0)},
  "geom": {"p": (0.0, 1.0)},
  "binom": {"n": (0.0, 1e4), "p": (0.0, 1.0)},
}

#: Search bounds for the discrete families. Legacy ``rv_discrete`` has no
#: ``fit`` method at all, so those go through ``scipy.stats.fit``, which needs
#: a bound for every parameter it estimates.
_FIT_BOUNDS: dict[str, dict[str, tuple[float, float]]] = {
  "poisson": {"mu": (1e-8, 1e4)},
  "nbinom": {"n": (1.0, 1e4), "p": (1e-8, 1.0)},
  "geom": {"p": (1e-8, 1.0)},
  "binom": {"n": (1.0, 1e4), "p": (1e-8, 1.0)},
}

#: Families deliberately kept out of ``_PARTITIONS``, and why. Each is still
#: reachable by name — ``ParametricMargin("vonmises")`` works — but never enters
#: an automatic candidate set, because a blind sweep over all of SciPy's
#: continuous families ranks ``vonmises`` first on clean gamma data, beating the
#: true family by 701 AIC units, purely because it advertises the whole real
#: line as its support while its density integrates to 63.75 there. The reasons
#: are rendered into ``MarginSelector``'s docstring rather than restated there,
#: so the documented list and the honored one cannot drift apart.
_EXCLUDED: dict[str, str] = {
  "vonmises": "periodic; advertises an unbounded support its density does not "
  "integrate to 1 over",
  "wrapcauchy": "periodic; same support misstatement as vonmises",
  "levy_stable": "fit does not terminate",
  "studentized_range": "fit does not terminate",
  "pareto": "shape escapes to ~1e8 rather than converging",
  "f": "boundary escapes, and slow",
  "loggamma": "boundary escapes",
  "burr": "boundary escapes",
  "ncf": "boundary escapes, and slow",
  "johnsonsu": "boundary escapes",
  "nct": "boundary escapes, and slow",
  "tukeylambda": "boundary escapes",
  "chi2": "an exact reparameterization of gamma; ties its AIC, which makes the "
  "winner arbitrary",
  "erlang": "gamma with an integer shape; same tie as chi2",
  "irwinhall": "no usable fit",
  "kstwo": "no usable fit",
  "binom": "the total count is not identified alongside the success "
  "probability: the likelihood drifts to the Poisson limit (n -> inf, "
  "p -> 0). Usable with the count pinned, as ParametricMargin('binom', "
  "fn=20.0)",
}


def _stats() -> Any:
  """Import ``scipy.stats``, or explain which extra provides it.

  Returns
  -------
  module
      The ``scipy.stats`` module.

  Raises
  ------
  ImportError
      If SciPy is not installed.
  """
  try:
    import scipy.stats as stats
  except ImportError as e:  # pragma: no cover - exercised in a subprocess
    raise ImportError(
      "pyvinecopulib.margins.ParametricMargin requires SciPy. "
      "Install it with `pip install pyvinecopulib[scipy]`."
    ) from e
  return stats


#: Seed pinned into the discrete fitter. ``scipy.stats.fit`` optimizes with
#: ``differential_evolution``, which is stochastic, so an unseeded call returns a
#: different parameter vector -- and a different ``report_`` row -- every time it
#: runs on the same data.
_DISCRETE_FIT_SEED = 5489


def _seeded_optimizer(objective: Any, **kwargs: Any) -> Any:
  """Optimize as ``scipy.stats.fit`` does, with the RNG pinned.

  Parameters
  ----------
  objective : callable
      The negative log-likelihood, as ``scipy.stats.fit`` builds it.
  **kwargs
      The remaining arguments ``scipy.stats.fit`` supplies, ``bounds`` and
      ``integrality``.

  Returns
  -------
  scipy.optimize.OptimizeResult
      The optimizer result, carrying ``x`` and ``fun``.
  """
  import importlib
  import inspect

  de = importlib.import_module("scipy.optimize").differential_evolution
  # `seed` was renamed to `rng` across the supported SciPy range.
  key = "rng" if "rng" in inspect.signature(de).parameters else "seed"
  return de(objective, **{key: _DISCRETE_FIT_SEED}, **kwargs)


def _curated_margin(
  family: str,
  partition: str,
  bounds: Optional[tuple[float, float]] = None,
) -> "ParametricMargin":
  """Build a candidate with its group's fixed-parameter policy applied.

  Parameters
  ----------
  family : str
      A family name.
  partition : str
      One of the support groups.
  bounds : tuple of float, or None, optional
      Known support ``(a, b)``; required by the ``"bounded"`` group.

  Returns
  -------
  ParametricMargin
      An unfitted candidate.
  """
  return ParametricMargin(family, **_curated_fixed(family, partition, bounds))


def _curated_fixed(
  family: str,
  partition: str,
  bounds: Optional[tuple[float, float]] = None,
) -> dict[str, Any]:
  """Return the fixed-parameter policy for a curated family.

  Whether ``loc`` and ``scale`` are estimated or pinned is a property of the
  group a family is used in, not of the family itself: ``gamma`` on positive
  data has ``loc = 0`` by assumption, whereas the same family on an interval of
  unknown origin does not. Pinning matters twice over — a free ``loc`` is the
  boundary escape that makes ``weibull_min`` return a plausible triple whose
  support excludes the smallest observation, and every pinned parameter is one
  fewer parameter in the information criterion.

  Parameters
  ----------
  family : str
      A family name.
  partition : str
      One of the support groups.
  bounds : tuple of float, or None, optional
      Known support ``(a, b)``; required by the ``"bounded"`` group.

  Returns
  -------
  dict
      Keyword arguments to pass to :class:`ParametricMargin`. Typed loosely
      because they are splatted into a signature whose ``**fixed`` sits
      alongside named parameters of other types.

  """
  del family  # every family in a group is anchored the same way
  if partition == "real":
    return {}
  if partition == "positive":
    return {"floc": 0.0}
  if partition == "unit":
    return {"floc": 0.0, "fscale": 1.0}
  if partition == "count":
    return {"floc": 0.0}
  assert bounds is not None, "the bounded group is only used with bounds"
  return {"floc": float(bounds[0]), "fscale": float(bounds[1] - bounds[0])}


class ParametricMargin(MarginBase[np.ndarray]):
  """A univariate margin from one ``scipy.stats`` family.

  Construct it with a family name and, optionally, parameters to hold fixed in
  SciPy's own ``f``-prefixed spelling — ``ParametricMargin("gamma",
  floc=0.0)`` estimates the shape and scale of a gamma anchored at the origin.
  Supplying the full parameter vector instead makes the margin a *fixed* one,
  already fitted: ``ParametricMargin("norm", (0.0, 1.0))`` is the standard
  normal and needs no data.

  Parameters
  ----------
  family : str
      Name of a ``scipy.stats`` distribution, e.g. ``"gamma"``.
  params : sequence of float, or None, optional
      The full parameter vector in SciPy's order (shape parameters, then
      ``loc``, then ``scale`` for a continuous family). Given here, the margin
      is already fitted and :attr:`n_parameters` is 0, since nothing was
      estimated from data.
  bounds : mapping of str to tuple of float, or None, optional
      Search bounds per parameter name, used only by the discrete families.
      Curated families come with their own; pass this to widen them or to fit
      a family the curated set does not cover.
  **fixed : float
      Parameters to hold fixed, named as SciPy's legacy fitter names them:
      ``floc``, ``fscale``, and ``f`` followed by a shape parameter's name
      (``fa`` for ``gamma``, ``fdf`` for ``t``). Positional aliases ``f0``,
      ``f1``, … are accepted too.

  Raises
  ------
  ImportError
      If SciPy is not installed.
  ValueError
      If ``family`` names no ``scipy.stats`` distribution, if a keyword does
      not name one of its parameters, or if ``params`` has the wrong length.

  See Also
  --------
  pyvinecopulib.margins.MarginSelector : Choose among several of these.
  pyvinecopulib.core.Kde1d : The nonparametric default.

  Notes
  -----
  **Continuous families use SciPy's legacy fitter**,
  ``scipy.stats.<family>.fit``, not ``scipy.stats.fit``. The modern
  function defaults to *fixing* ``loc = 0`` and ``scale = 1`` rather than
  estimating them — it reports ``success=True`` for a ``norm`` it never
  fitted — and its ``differential_evolution`` search is some 800 times slower
  than the closed-form solutions SciPy ships for 14 of the 15 curated
  families. Discrete families have no legacy fitter, so those do go through
  ``scipy.stats.fit``, with explicit search bounds; their ``loc``
  (the lattice offset) is pinned to 0 unless ``floc`` says otherwise, because
  it is not identified alongside the shape parameters.

  **Observation weights are not supported.** Nothing in SciPy's estimation
  surface accepts them — not the legacy ``fit``, not ``scipy.stats.fit``,
  not ``nnlf`` — so a weighted objective would mean our own optimizer and the
  loss of the analytic fast paths above. :meth:`fit` raises rather than
  silently returning an unweighted fit.

  Examples
  --------
  ::

      import numpy as np
      from pyvinecopulib.margins import ParametricMargin

      x = np.random.default_rng(0).gamma(2.5, 1.5, size=500)
      m = ParametricMargin("gamma", floc=0.0).fit(x)
      m.parameters       # -> (shape, 0.0, scale)
      m.n_parameters     # -> 2, not 3: `loc` was pinned
  """

  supports_weights: bool = False

  def __init__(
    self,
    family: str,
    params: Optional[Sequence[float]] = None,
    *,
    bounds: Optional[Mapping[str, tuple[float, float]]] = None,
    **fixed: float,
  ) -> None:
    stats = _stats()
    dist = getattr(stats, family, None)
    if not isinstance(dist, (stats.rv_continuous, stats.rv_discrete)):
      raise ValueError(
        f"unknown scipy.stats family {family!r}; expected the name of a "
        "distribution such as 'gamma' or 'poisson'"
      )
    self._family = family
    self._discrete = isinstance(dist, stats.rv_discrete)
    self._names = _parameter_names(dist, discrete=self._discrete)
    self._fixed = _normalize_fixed(fixed, self._names, family)
    if self._discrete and "loc" not in self._fixed:
      # `scipy.stats.fit` needs a bound for every free parameter, and a free
      # lattice offset is not identified alongside the shape parameters.
      self._fixed["loc"] = 0.0
    self._bounds = (
      dict(bounds) if bounds is not None else dict(_FIT_BOUNDS.get(family, {}))
    )

    self._params: Optional[tuple[float, ...]] = None
    self._n_free = 0
    self._loglik: Optional[float] = None
    if params is not None:
      values = tuple(float(v) for v in params)
      if len(values) != len(self._names):
        raise ValueError(
          f"{family!r} takes {len(self._names)} parameters "
          f"{self._names}, got {len(values)}"
        )
      self._params = values

  # --- identity ------------------------------------------------------------ #

  @property
  def family_name(self) -> str:
    """Name of the SciPy family.

    Returns
    -------
    str
        The name passed to the constructor.
    """
    return self._family

  @property
  def parameter_names(self) -> tuple[str, ...]:
    """Parameter names in SciPy's order.

    Returns
    -------
    tuple of str
        The shape parameters, then ``"loc"``, then ``"scale"`` for a
        continuous family.
    """
    return self._names

  @property
  def parameters(self) -> tuple[float, ...]:
    """Parameter values in SciPy's order.

    Returns
    -------
    tuple of float
        One value per entry of :attr:`parameter_names`.

    Raises
    ------
    RuntimeError
        If the margin has neither been fitted nor been given parameters.
    """
    if self._params is None:
      raise RuntimeError(
        f"ParametricMargin({self._family!r}) is not fitted; call fit(y) or "
        "pass params="
      )
    return self._params

  @property
  def fixed_parameters(self) -> dict[str, float]:
    """Parameters held fixed rather than estimated.

    Returns
    -------
    dict
        Parameter name to pinned value.
    """
    return dict(self._fixed)

  @property
  def n_parameters(self) -> float:
    """Number of freely estimated parameters.

    Only what :meth:`fit` actually estimated counts. SciPy's legacy fitter
    always returns ``loc`` and ``scale``, so the length of the parameter vector
    overstates this whenever one of them was pinned — by enough to reverse a
    close comparison, since it moves AIC by 2 per parameter.

    Returns
    -------
    float
        The count; 0 for a margin whose parameters were given at construction.
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
        f"ParametricMargin({self._family!r}).loglik() is only defined after "
        "fit(y); pass a sample to evaluate one instead"
      )
    return self._loglik

  @property
  def is_fitted(self) -> bool:
    return self._params is not None

  @property
  def var_type(self) -> str:
    return "d" if self._discrete else "c"

  @property
  def support(self) -> tuple[float, float]:
    """Closed bounds of the fitted support.

    Returns
    -------
    tuple of float
        ``(lo, hi)`` at the current parameters, and ``(-inf, inf)`` while the
        margin is unfitted, since the support of most families depends on them.
    """
    if self._params is None:
      return (float("-inf"), float("inf"))
    lo, hi = self._dist.support(*self._params)
    return (float(lo), float(hi))

  # --- estimation ---------------------------------------------------------- #

  @property
  def _dist(self) -> Any:
    """Return the SciPy distribution object.

    Looked up by name rather than stored, so the margin pickles as a name and
    a parameter tuple.

    Returns
    -------
    rv_continuous or rv_discrete
        The family.
    """
    return getattr(_stats(), self._family)

  def fit(
    self, y: Any, *, x: Optional[Any] = None, weights: Optional[Any] = None
  ) -> "ParametricMargin":
    """Estimate the free parameters by maximum likelihood.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations; NaNs are dropped.
    x : array, shape (n, k), or None, optional
        Not supported; passing covariates raises rather than silently
        fitting an unconditional margin.
    weights : array, shape (n,), or None, optional
        Not supported; see the Notes on the class.

    Returns
    -------
    ParametricMargin
        ``self``.

    Raises
    ------
    TypeError
        If ``weights`` is given.
    ValueError
        If no observation survives, or a free discrete parameter has no search
        bound.
    """
    _reject_covariates(self, x)
    if weights is not None:
      raise TypeError(
        f"ParametricMargin({self._family!r}) cannot use observation weights: "
        "SciPy's estimators do not accept them. Pass margins='kde' or "
        "margins='empirical' for a weighted fit, or drop weights="
      )
    data = np.asarray(y, dtype=float).ravel()
    data = data[~np.isnan(data)]
    if data.size == 0:
      raise ValueError(
        f"ParametricMargin({self._family!r}).fit got no usable observation"
      )

    free = [name for name in self._names if name not in self._fixed]
    if not free:
      # SciPy's fitter refuses a problem with nothing to optimize.
      self._params = tuple(self._fixed[name] for name in self._names)
    elif self._discrete:
      self._params = self._fit_discrete(data, free)
    else:
      kwargs = {f"f{name}": value for name, value in self._fixed.items()}
      self._params = tuple(float(v) for v in self._dist.fit(data, **kwargs))
    self._n_free = len(free)
    self._loglik = float(np.sum(self.logpdf(data)))
    return self

  def _fit_discrete(
    self, data: np.ndarray, free: Sequence[str]
  ) -> tuple[float, ...]:
    """Fit a discrete family through ``scipy.stats.fit``.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        Observations, NaNs already dropped.
    free : sequence of str
        Names of the parameters to estimate.

    Returns
    -------
    tuple of float
        Parameter values in SciPy's order.

    Raises
    ------
    ValueError
        If a free parameter has no search bound.
    """
    missing = [name for name in free if name not in self._bounds]
    if missing:
      raise ValueError(
        f"fitting {self._family!r} needs search bounds for {missing}; pass "
        "bounds={'name': (lo, hi)} or pin the parameter with f<name>="
      )
    search: dict[str, tuple[float, float]] = {}
    for name in self._names:
      if name in self._fixed:
        value = self._fixed[name]
        search[name] = (value, value)
      else:
        search[name] = self._bounds[name]
    result = _stats().fit(
      self._dist, data, bounds=search, optimizer=_seeded_optimizer
    )
    return tuple(float(v) for v in result.params)

  # --- evaluation ---------------------------------------------------------- #

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    params = self.parameters
    dist = self._dist
    values = np.asarray(y, dtype=float)
    if self._discrete:
      return np.asarray(dist.pmf(values, *params), dtype=float)
    return np.asarray(dist.pdf(values, *params), dtype=float)

  def logpdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    """Log of :meth:`pdf`, from SciPy's own log-density.

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
    params = self.parameters
    dist = self._dist
    values = np.asarray(y, dtype=float)
    if self._discrete:
      return np.asarray(dist.logpmf(values, *params), dtype=float)
    return np.asarray(dist.logpdf(values, *params), dtype=float)

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return np.asarray(
      self._dist.cdf(np.asarray(y, dtype=float), *self.parameters),
      dtype=float,
    )

  def icdf(self, p: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return np.asarray(
      self._dist.ppf(np.asarray(p, dtype=float), *self.parameters),
      dtype=float,
    )

  def sample(
    self, n: int, *, x: Optional[Any] = None, seeds: Optional[list[int]] = None
  ) -> np.ndarray:
    """Draw ``n`` samples from the family.

    Parameters
    ----------
    n : int
        Number of samples.
    x : array, shape (n, k), or None, optional
        Ignored; this margin is unconditional.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n,), dtype float
        Samples on the original scale.
    """
    rng = np.random.default_rng(list(seeds) if seeds else None)
    return np.asarray(
      self._dist.rvs(*self.parameters, size=n, random_state=rng), dtype=float
    )

  def __repr__(self) -> str:
    if self._params is None:
      pinned = ", ".join(f"{k}={v!r}" for k, v in self._fixed.items())
      return f"ParametricMargin({self._family!r}, unfitted, {pinned})"
    shown = ", ".join(
      f"{name}={value:.4g}" for name, value in zip(self._names, self._params)
    )
    return f"ParametricMargin({self._family!r}, {shown})"


def _parameter_names(dist: Any, *, discrete: bool) -> tuple[str, ...]:
  """Return a family's parameter names in SciPy's order.

  Parameters
  ----------
  dist : rv_continuous or rv_discrete
      The family.
  discrete : bool
      Whether it is discrete, in which case it has no ``scale``.

  Returns
  -------
  tuple of str
      Shape parameter names, then ``"loc"``, then ``"scale"`` when continuous.
  """
  shapes = tuple(
    part.strip() for part in (dist.shapes or "").split(",") if part.strip()
  )
  return shapes + ("loc",) if discrete else shapes + ("loc", "scale")


def _normalize_fixed(
  fixed: Mapping[str, float], names: Sequence[str], family: str
) -> dict[str, float]:
  """Map SciPy's ``f``-prefixed keywords onto parameter names.

  Parameters
  ----------
  fixed : mapping of str to float
      Keywords as the caller spelled them: ``floc``, ``fscale``, ``f<name>``,
      or the positional aliases ``f0``, ``f1``, ….
  names : sequence of str
      Parameter names in SciPy's order.
  family : str
      Family name, for the error message.

  Returns
  -------
  dict
      Parameter name to pinned value.

  Raises
  ------
  ValueError
      If a keyword names no parameter of the family.
  """
  aliases = {f"f{name}": name for name in names}
  aliases.update({f"f{i}": name for i, name in enumerate(names)})
  resolved: dict[str, float] = {}
  for key, value in fixed.items():
    if key not in aliases:
      raise ValueError(
        f"{key!r} does not name a parameter of {family!r}; its parameters are "
        f"{tuple(names)}, pinned as {tuple('f' + n for n in names)}"
      )
    resolved[aliases[key]] = float(value)
  return resolved
