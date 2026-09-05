"""Parametric margin: one SciPy family, fitted by maximum likelihood.

SciPy is an optional dependency; it is imported when a
``SciPyMargin`` is constructed, not when this module is imported, so
:mod:`pyvinecopulib.margins` keeps working without it.
"""

from __future__ import annotations

import re
import warnings
from typing import (
  Any,
  Iterable,
  Mapping,
  Optional,
  Sequence,
)

import numpy as np

from ..core import MarginBase
from ..core.margin_base import _reject_covariates
from ..core._validation import validate_univariate
from .controls import CRITERIA, FitControlsMargin

__all__ = ["SciPyMargin"]

#: Curated candidate families, grouped by the support they can represent.
#: `SciPyMargin.select` draws its candidates from the groups the data are
#: compatible with. Counts are never compared against continuous families:
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
#: reachable by name — ``SciPyMargin("vonmises")`` works — but never enters
#: an automatic candidate set, because a blind sweep over all of SciPy's
#: continuous families ranks ``vonmises`` first on clean gamma data, beating the
#: true family by 701 AIC units, purely because it advertises the whole real
#: line as its support while its density integrates to 63.75 there. The reasons
#: are rendered into ``SciPyMargin``'s docstring rather than restated
#: there, so the documented list and the honored one cannot drift apart.
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
  "p -> 0). Usable with the count pinned, as SciPyMargin('binom', "
  "fn=20.0)",
}


def _excluded_block(indent: str = "  ") -> str:
  """Render the exclusion table as an RST bullet list.

  Parameters
  ----------
  indent : str, optional
      Leading whitespace for each bullet.

  Returns
  -------
  str
      One bullet per excluded family.
  """
  return "\n".join(
    f"{indent}- ``{family}``: {reason}"
    for family, reason in sorted(_EXCLUDED.items())
  )


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
      "pyvinecopulib.margins.SciPyMargin requires SciPy. "
      "Install it with `pip install pyvinecopulib[scipy]`."
    ) from e
  return stats


#: Seed pinned into the discrete fitter. ``scipy.stats.fit`` optimizes with
#: ``differential_evolution``, which is stochastic, so an unseeded call returns a
#: different parameter vector every time it
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


def _anchoring_group(family: str, groups: Sequence[str]) -> str:
  """The group whose fixed-parameter policy a named family should get.

  A narrowed search has to anchor a family exactly as the curated search would
  have anchored it, or the two estimate a different number of parameters for
  the same family and their criteria stop being comparable -- and the boundary
  escape the policy exists to prevent comes back precisely when the caller
  asks for one family. So the most specific applicable group that offers the
  family wins, and a family outside the curated partition falls back to the
  least specific, where nothing is pinned.

  Parameters
  ----------
  family : str
      A family name.
  groups : sequence of str
      Applicable support groups, least specific first.

  Returns
  -------
  str
      One of ``groups``.
  """
  for group in reversed(list(groups)):
    if family in _PARTITIONS.get(group, ()):
      return group
  return groups[0]


def _curated_margin(
  family: str,
  partition: str,
  bounds: Optional[tuple[float, float]] = None,
) -> "SciPyMargin":
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
  SciPyMargin
      An unfitted candidate.
  """
  return SciPyMargin(family, **_curated_fixed(family, partition, bounds))


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
      Keyword arguments to pass to :class:`SciPyMargin`. Typed loosely
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


def _criteria(loglik: float, k: float, n: int) -> dict[str, float]:
  """Evaluate every criterion at one fit.

  Parameters
  ----------
  loglik : float
      Maximized log-likelihood.
  k : float
      Number of freely estimated parameters.
  n : int
      Number of observations.

  Returns
  -------
  dict
      One entry per criterion name; ``inf`` where a criterion is not
      defined, so a candidate can never win by being undefined.
  """
  if not np.isfinite(loglik):
    return dict.fromkeys(CRITERIA, float("inf"))
  aic = -2.0 * loglik + 2.0 * k
  bic = -2.0 * loglik + k * float(np.log(n))
  tail = n - k - 1.0
  aicc = aic + (2.0 * k * (k + 1.0) / tail if tail > 0 else float("inf"))
  return {"aic": aic, "bic": bic, "aicc": aicc}


def _reject(candidate: Any, y: np.ndarray) -> Optional[str]:
  """Return why a fitted candidate is inadmissible, or ``None`` if it is fine.

  Four things have to hold, and none of them implies the others. Every
  parameter is finite and non-degenerate — a Student ``t`` fitted to normal
  data comes back with ``df`` on the order of ``1e10``, and a ``scale`` of 0 is
  a point mass. The support contains the data, which is what catches the
  boundary escapes: ``weibull_min`` fitted with a free ``loc`` returns an
  entirely plausible triple whose ``loc`` sits above the smallest observation.
  And the log-likelihood is finite, which is the same escape seen from the
  other side.

  Parameters
  ----------
  candidate : MarginLike
      A fitted candidate.
  y : array, shape (n,), dtype float
      The observations it was fitted to.

  Returns
  -------
  str or None
      A one-line reason, or ``None`` when the candidate is admissible.
  """
  # `parameters` / `fixed_parameters` / `loglik` are not on `MarginLike`; a
  # candidate may be any margin, so read them defensively.
  margin: Any = candidate
  names: Sequence[str] = getattr(margin, "parameter_names", ())
  values: Sequence[float] = getattr(margin, "parameters", ())
  if len(names) == len(values):
    limits = _LIMITS.get(getattr(margin, "family_name", ""), {})
    for name, value in zip(names, values):
      if not np.isfinite(value):
        return f"non-finite parameter {name}={value}"
      lo, hi = limits.get(name, (float("-inf"), float("inf")))
      if name == "scale":
        lo = max(lo, 0.0)
      if not lo < value < hi:
        return f"degenerate parameter {name}={value:.6g}, outside ({lo}, {hi})"

  # A continuous family whose scale collapses onto an atom is the other boundary
  # escape, and the one an information criterion rewards hardest: a density
  # concentrating on a repeated value diverges, and so does its likelihood. On
  # data that is half exact zeros, a `t` fitted with scale ~1e-9 beats every
  # honest candidate by thousands of AIC units.
  if len(names) == len(values):
    spread = float(np.max(y) - np.min(y))
    floor = 1e-6 * spread if spread > 0.0 else 0.0
    for name, value in zip(names, values):
      if name == "scale" and floor > 0.0 and value < floor:
        return (
          f"degenerate parameter scale={value:.6g}, below {floor:.6g} "
          "(collapsing onto an atom rather than fitting the spread)"
        )

  lo, hi = getattr(margin, "support", (float("-inf"), float("inf")))
  # A bound the caller pinned is allowed to coincide with the data range: a
  # uniform on a known [a, b] is not an escape. An estimated one is not, since
  # touching it is the escape.
  pinned = getattr(margin, "fixed_parameters", {})
  lo_ok = lo <= y.min() if "loc" in pinned else lo < y.min()
  hi_ok = y.max() <= hi if "scale" in pinned else y.max() < hi
  if not (lo_ok and hi_ok):
    return (
      f"support ({lo:.6g}, {hi:.6g}) does not contain the data range "
      f"[{y.min():.6g}, {y.max():.6g}]"
    )

  loglik = float(np.sum(margin.logpdf(y)))
  if not np.isfinite(loglik):
    return f"non-finite log-likelihood ({loglik})"
  return None


def _fit_candidate(candidate: Any, y: np.ndarray) -> Optional[str]:
  """Fit one candidate and report why it is inadmissible, or ``None``.

  Parameters
  ----------
  candidate : SciPyMargin
      An unfitted candidate.
  y : array, shape (n,), dtype float
      The observations.

  Returns
  -------
  str or None
      A one-line reason, or ``None`` when the candidate is admissible.
  """
  # A candidate's optimizer chatter says nothing a caller can act on -- the
  # outcome is admissible or not either way -- so it is collected rather than
  # emitted once per family tried.
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
      if not getattr(candidate, "is_fitted", True):
        candidate.fit(y)
      return _reject(candidate, y)
    except Exception as e:  # noqa: BLE001 - any fitter failure is a rejection
      return f"{type(e).__name__}: {e}"


def _dedupe(candidates: Iterable[Any]) -> list[Any]:
  """Drop candidates that would tie with one already present.

  Two unfitted candidates of the same family with the same pinned parameters
  and search bounds fit the same model, so they tie exactly on every criterion:
  the winner becomes an artifact of iteration order, and the report carries the
  row twice. Ready-made fitted candidates additionally include their parameter
  vectors in their identity. Anything whose identity cannot be read this way is
  kept, since dropping it would be a guess.

  Parameters
  ----------
  candidates : iterable
      Unfitted candidate margins, in preference order.

  Returns
  -------
  list
      The first candidate of each distinct family, pins, and search bounds,
      order preserved.
  """
  seen: set[tuple[Any, ...]] = set()
  out: list[Any] = []
  for margin in candidates:
    family = getattr(margin, "family_name", None)
    if family is None:
      out.append(margin)
      continue
    fixed = getattr(margin, "fixed_parameters", None) or {}
    search_bounds = (
      tuple(sorted(margin._bounds.items()))
      if isinstance(margin, SciPyMargin)
      else ()
    )
    is_fitted = bool(getattr(margin, "is_fitted", False))
    if is_fitted:
      raw_parameters = getattr(margin, "parameters", None)
      if raw_parameters is None:
        # Fitted state without a readable identity may represent any model.
        # Keeping it is the only ownership-safe choice.
        out.append(margin)
        continue
      parameters = tuple(raw_parameters)
    else:
      parameters = ()
    key = (
      str(family),
      tuple(sorted(fixed.items())),
      search_bounds,
      is_fitted,
      parameters,
    )
    if key in seen:
      continue
    seen.add(key)
    out.append(margin)
  return out


class SciPyMargin(MarginBase[np.ndarray]):
  """A univariate margin from one ``scipy.stats`` family.

  Construct it with a family name and, optionally, parameters to hold fixed in
  SciPy's own ``f``-prefixed spelling — ``SciPyMargin("gamma",
  floc=0.0)`` estimates the shape and scale of a gamma anchored at the origin.
  Supplying the full parameter vector instead makes the margin a *fixed* one,
  already fitted: ``SciPyMargin("norm", (0.0, 1.0))`` is the standard
  normal and needs no data.

  Parameters
  ----------
  family : str, or None, optional
      Name of a ``scipy.stats`` distribution, e.g. ``"gamma"``. Omitted, the
      margin has no family until :meth:`select` chooses one from the candidate
      set -- which is what ``margins="parametric"`` means.
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
  pyvinecopulib.margins.OpenTURNSMargin : The OpenTURNS counterpart.
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

  **The candidate set** :meth:`select` searches is curated rather than
  everything SciPy ships, and that is a statistical decision. A
  log-likelihood is only comparable across families normalized over the same
  set, so a blind sweep promotes a family that misstates its own support above
  the truth: over all of SciPy's continuous families, ``vonmises`` ranks first
  on clean gamma data, beating the true family by 701 AIC units, purely
  because it advertises the whole real line while its density integrates to
  63.75 there. The set is grouped by support instead -- real line, positive,
  unit interval, a caller-supplied interval, counts -- with ``loc`` (and
  ``scale`` where the support is known) pinned. Pinning pays twice: it removes
  the boundary escapes, and it removes a parameter the criterion would
  otherwise charge for. Every family stays reachable by name;
  the ones kept out of the automatic set are:

  {excluded}

  **Observation weights are not supported.** Nothing in SciPy's estimation
  surface accepts them — not the legacy ``fit``, not ``scipy.stats.fit``,
  not ``nnlf`` — so a weighted objective would mean our own optimizer and the
  loss of the analytic fast paths above. :meth:`fit` raises rather than
  silently returning an unweighted fit.

  Examples
  --------
  ::

      import numpy as np
      from pyvinecopulib.margins import SciPyMargin

      x = np.random.default_rng(0).gamma(2.5, 1.5, size=500)
      m = SciPyMargin("gamma", floc=0.0).fit(x)
      m.parameters       # -> (shape, 0.0, scale)
      m.n_parameters     # -> 2, not 3: `loc` was pinned
  """

  supports_weights: bool = False

  def __init__(
    self,
    family: Optional[str] = None,
    params: Optional[Sequence[float]] = None,
    *,
    bounds: Optional[Mapping[str, tuple[float, float]]] = None,
    **fixed: float,
  ) -> None:
    self._declared_var_type: Optional[str] = None
    self._declared_support: Optional[tuple[float, float]] = None
    if family is None:
      if params is not None:
        raise ValueError(
          "params= was given without a family; name the family it belongs to, "
          "or leave both out and let select(y) choose one"
        )
      self._unnamed(bounds, fixed)
      return
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
    self._nobs: Optional[int] = None
    if params is not None:
      values = tuple(float(v) for v in params)
      if len(values) != len(self._names):
        raise ValueError(
          f"{family!r} takes {len(self._names)} parameters "
          f"{self._names}, got {len(values)}"
        )
      self._params = values

  def _unnamed(
    self,
    bounds: Optional[Mapping[str, tuple[float, float]]],
    fixed: Mapping[str, float],
  ) -> None:
    """Initialize a margin whose family is not chosen yet.

    Parameters
    ----------
    bounds : mapping, or None
        Search bounds to apply to whichever family is chosen.
    fixed : mapping
        Parameters to pin on whichever family is chosen.

    Returns
    -------
    None
    """
    self._family = None
    self._discrete = False
    self._names: tuple[str, ...] = ()
    self._fixed: dict[str, float] = dict(fixed)
    self._bounds = dict(bounds) if bounds is not None else {}
    self._params: Optional[tuple[float, ...]] = None
    self._n_free = 0
    self._loglik: Optional[float] = None
    self._nobs: Optional[int] = None

  def declare(
    self,
    *,
    var_type: Optional[str] = None,
    support: Optional[tuple[float, float]] = None,
  ) -> "SciPyMargin":
    """Accept what the caller knows, to be honored by :meth:`select`.

    A named family already fixes both, so this only steers a search: the type
    decides whether the count families or the continuous ones are candidates,
    and the support selects the bounded group and pins its endpoints.

    Parameters
    ----------
    var_type : str or None, optional
        ``"c"``, ``"d"`` or ``"zi"``.
    support : tuple of float, or None, optional
        Declared bounds as ``(lo, hi)``.

    Returns
    -------
    SciPyMargin
        ``self``, so the call chains into :meth:`select`.
    """
    if var_type is not None:
      self._declared_var_type = "d" if var_type == "zi" else var_type
    if support is not None:
      lo, hi = support
      if (
        lo is not None
        and hi is not None
        and np.isfinite(lo)
        and np.isfinite(hi)
      ):
        self._declared_support = (float(lo), float(hi))
    return self

  # --- selection ----------------------------------------------------------- #

  def select(
    self,
    y: Any,
    /,
    controls: Optional[Any] = None,
    *,
    x: Optional[Any] = None,
    weights: Optional[Any] = None,
  ) -> "SciPyMargin":
    """Choose a family from the candidate set, fit it, and become it.

    Every admissible candidate is fitted and scored, and the best on the
    criterion wins; this margin then *is* that family. The default candidate
    set is curated rather than everything SciPy ships, and that is a
    statistical decision: a log-likelihood is only comparable across families
    normalized over the same set, so a blind sweep promotes a family that
    misstates its own support above the truth. See the class Notes for the
    grouping and the exclusions.

    A margin that already names a family reduces this to :meth:`fit`: naming
    the family *is* the choice, so ``family_set`` is how a caller asks for the
    search back on one.

    Nothing is dropped silently. A candidate that fails to fit, or fits to
    something inadmissible, is refused with its reason, and a variable where
    every candidate fails raises and names each one.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations; NaNs are dropped.
    controls : FitControlsMargin, or None, optional
        Bounds the search: ``family_set`` names the candidates,
        ``selection_criterion`` scores them, ``var_type`` and ``support`` say
        what the caller knows, and ``on_failure`` decides what an
        all-candidates-failed variable does.
    x : array, shape (n, k), or None, optional
        Not supported; passing covariates raises rather than silently
        selecting an unconditional margin.
    weights : array, shape (n,), or None, optional
        Not supported; SciPy's estimators do not accept them.

    Returns
    -------
    SciPyMargin
        ``self``, now carrying the chosen family.

    Raises
    ------
    TypeError
        If ``weights`` is given.
    ValueError
        If no observation survives, or if every candidate is inadmissible and
        ``on_failure`` is ``"raise"``.

    See Also
    --------
    fit : Estimate a named family, leaving the family alone.
    aic : Score the chosen fit.
    """
    _reject_covariates(self, x)
    if weights is not None:
      raise TypeError(
        "SciPyMargin cannot use observation weights: SciPy's estimators "
        "do not accept them. Pass margins='kde' for a weighted fit, or drop "
        "weights="
      )
    settings = controls if controls is not None else FitControlsMargin()
    criterion = getattr(settings, "selection_criterion", "aic")
    family_set = getattr(settings, "family_set", None)
    self.declare(
      var_type=getattr(settings, "var_type", None),
      support=getattr(settings, "support", None),
    )
    if self._family is not None and family_set is None:
      # Naming a family *is* the choice, so there is nothing left to select and
      # the base contract applies: reduce to `fit`. Replacing it silently would
      # answer a specification with a different model. A caller who does want
      # the search back asks for it by name, with `family_set`.
      return self.fit(y, x=x, weights=weights)

    data = validate_univariate(np.asarray(y, dtype=float))
    data = data[~np.isnan(data)]
    if data.size == 0:
      raise ValueError("SciPyMargin.select got no usable observation")

    counts = self._reads_as_counts(data)
    candidates = self._candidates(data, counts=counts, family_set=family_set)
    scored: list[tuple[float, SciPyMargin]] = []
    refused: list[str] = []
    for candidate in candidates:
      reason = _fit_candidate(candidate, data)
      if reason is not None:
        refused.append(f"{candidate.family_name}: {reason}")
        continue
      loglik = candidate._loglik
      score = _criteria(
        float("-inf") if loglik is None else float(loglik),
        candidate.n_parameters,
        int(data.size),
      )[criterion]
      scored.append((score, candidate))

    if not scored:
      return self._no_candidate(data, refused, settings)
    scored.sort(key=lambda pair: pair[0])
    return self._adopt(scored[0][1])

  def _reads_as_counts(self, data: np.ndarray) -> bool:
    """Whether the count families are the candidate group.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations.

    Returns
    -------
    bool
        ``True`` for the count group.
    """
    if self._declared_var_type is not None:
      return self._declared_var_type == "d"
    return bool(np.all(data >= 0.0) and np.all(data == np.round(data)))

  def _candidates(
    self,
    data: np.ndarray,
    *,
    counts: bool,
    family_set: Optional[Sequence[str]],
  ) -> list["SciPyMargin"]:
    """Build the candidates to fit.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations, which decide the applicable support groups.
    counts : bool
        Whether the count group applies.
    family_set : sequence of str, or None
        Families to consider, or ``None`` for the curated set.

    Returns
    -------
    list of SciPyMargin
        Unfitted candidates, in preference order.

    Raises
    ------
    ValueError
        If every named family is of the other kind than the data read as.
    """
    groups = self._groups(data, counts=counts)
    if family_set is None:
      # `unit` and `bounded` both offer `beta`, differing only in its pins,
      # and a future partition may overlap outright. Deduplicating keeps that
      # from producing two candidates that tie exactly on every criterion,
      # which would make the winner an artifact of iteration order.
      return _dedupe(
        _curated_margin(family, group, self._declared_support)
        for group in groups
        for family in _PARTITIONS[group]
      )

    want = "d" if counts else "c"
    out: list[SciPyMargin] = []
    dropped: list[str] = []
    for family in family_set:
      candidate = _curated_margin(
        family, _anchoring_group(family, groups), self._declared_support
      )
      # The partition the curated path applies, applied here too: a count
      # family reports a probability mass and a continuous one a Lebesgue
      # density, so ranking both on one criterion compares numbers that are
      # not commensurable -- and the density usually wins.
      if candidate.var_type == want:
        out.append(candidate)
      else:
        dropped.append(family)
    out = _dedupe(out)
    reads = "counts" if counts else "continuous"
    if dropped and not out:
      raise ValueError(
        f"the data read as {reads}, and every family in family_set is of the "
        f"other kind: {', '.join(dropped)}. Set var_type to say how the data "
        "should be read, or name families that match it."
      )
    if dropped:
      warnings.warn(
        f"the data read as {reads}, so {len(dropped)} family/families of the "
        f"other kind were dropped: {', '.join(dropped)}. A probability mass "
        "and a Lebesgue density are not comparable on one information "
        "criterion.",
        UserWarning,
        stacklevel=3,
      )
    return out

  def _groups(self, data: np.ndarray, *, counts: bool) -> list[str]:
    """Support groups the data are compatible with.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations.
    counts : bool
        Whether the count group applies.

    Returns
    -------
    list of str
        Keys into the curated partition table, in the order they are tried.
    """
    if counts:
      return ["count"]
    if self._declared_support is not None:
      return ["bounded"]
    groups = ["real"]
    if data.min() > 0.0:
      groups.append("positive")
      if data.max() < 1.0:
        groups.append("unit")
    return groups

  def _adopt(self, winner: "SciPyMargin") -> "SciPyMargin":
    """Become the selected candidate.

    Parameters
    ----------
    winner : SciPyMargin
        The candidate that won.

    Returns
    -------
    SciPyMargin
        ``self``.
    """
    self._family = winner._family
    self._discrete = winner._discrete
    self._names = winner._names
    self._fixed = dict(winner._fixed)
    self._bounds = dict(winner._bounds)
    self._params = winner._params
    self._n_free = winner._n_free
    self._loglik = winner._loglik
    self._nobs = winner._nobs
    return self

  def _no_candidate(
    self,
    data: np.ndarray,
    refused: Sequence[str],
    settings: Any,
  ) -> "SciPyMargin":
    """Answer a variable where every candidate was refused.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations.
    refused : sequence of str
        One ``family: reason`` line per candidate.
    settings : FitControlsMargin
        The resolved settings, read for ``on_failure``.

    Returns
    -------
    SciPyMargin
        Never returns; the annotation matches :meth:`select`.

    Raises
    ------
    ValueError
        Always. Substituting a different kind of margin is not something this
        class can do -- it would have to stop being parametric -- so
        ``on_failure="fallback"`` is honored one level up, where the margin
        for a column is chosen.
    """
    del data, settings
    detail = "\n  ".join(refused) or "(no candidate was admissible)"
    raise ValueError(
      "no parametric family fits this variable; every candidate was "
      f"refused:\n  {detail}\nName a family directly, widen family_set, or "
      'pass on_failure="fallback" for a kernel-density margin.'
    )

  # --- scoring ------------------------------------------------------------- #
  # --- identity ------------------------------------------------------------ #

  @property
  def family_name(self) -> str:
    """Name of the SciPy family.

    Returns
    -------
    str
        The name passed to the constructor, or chosen by :meth:`select`.

    Raises
    ------
    RuntimeError
        If the family is neither named nor selected yet.
    """
    if self._family is None:
      raise RuntimeError(
        "SciPyMargin() has no family yet; name one at construction or "
        "call select(y) to choose from the candidate set"
      )
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
        f"SciPyMargin({self._family!r}) is not fitted; call fit(y) or "
        "pass params="
      )
    return self._params

  @property
  def nobs(self) -> Optional[int]:
    """Number of observations the fit used.

    Returns
    -------
    int or None
        The sample size, or ``None`` on a margin given its parameters rather
        than fitted. It is what penalizes ``bic`` and ``aicc``.
    """
    return self._nobs

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
        f"SciPyMargin({self._family!r}).loglik() is only defined after "
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
  def supported_var_types(self) -> tuple[str, ...]:
    """Variable types this margin can represent.

    A SciPy family is either continuous or discrete, never both, so once a
    family is named the set is a singleton and a declared type cannot widen
    it. A margin whose family is still to be selected serves both, since the
    candidate set spans continuous and count families.

    Returns
    -------
    tuple of str
        ``("d",)`` for a discrete family, ``("c",)`` for a continuous one, and
        both while the family is unchosen.
    """
    if self._family is None:
      return ("c", "d")
    return ("d",) if self._discrete else ("c",)

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
    return getattr(_stats(), self.family_name)

  def to_json(self) -> dict[str, Any]:
    """Return this margin's JSON payload.

    Returns
    -------
    dict
        A JSON-serializable mapping that
        :func:`~pyvinecopulib.core.margin_from_json` reads back.
    """
    payload: dict[str, Any] = {
      "kind": "SciPyMargin",
      "family": self._family,
      "fixed": dict(self._fixed),
      "bounds": {k: list(v) for k, v in (self._bounds or {}).items()},
    }
    if self._params is not None:
      payload["params"] = list(self._params)
    if self._loglik is not None:
      payload["loglik"] = self._loglik
    payload["n_free"] = self._n_free
    return payload

  @classmethod
  def _from_json_payload(cls, payload: dict[str, Any]) -> "SciPyMargin":
    """Rebuild a margin from the payload :meth:`to_json` produced."""
    bounds = {
      k: (float(v[0]), float(v[1]))
      for k, v in (payload.get("bounds") or {}).items()
    } or None
    fixed: dict[str, float] = {
      str(k): float(v) for k, v in (payload.get("fixed") or {}).items()
    }
    # `params` positionally, so the `**fixed` unpacking cannot reach it.
    margin = cls(str(payload["family"]), None, bounds=bounds, **fixed)
    if "params" in payload:
      margin._params = tuple(float(v) for v in payload["params"])
    if payload.get("loglik") is not None:
      margin._loglik = float(payload["loglik"])
    # `n_free` is what `n_parameters` reports, and it separates a fitted margin
    # from one constructed with `params=` pinned.
    margin._n_free = int(payload.get("n_free", 0))
    return margin

  def fit(
    self,
    y: Any,
    /,
    controls: Optional[Any] = None,
    *,
    x: Optional[Any] = None,
    weights: Optional[Any] = None,
  ) -> "SciPyMargin":
    """Estimate the free parameters by maximum likelihood.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations; NaNs are dropped.
    controls : object, or None, optional
        Unused; the family is fixed here, so there is nothing to configure.
        A search over families is :meth:`select`.
    x : array, shape (n, k), or None, optional
        Not supported; passing covariates raises rather than silently
        fitting an unconditional margin.
    weights : array, shape (n,), or None, optional
        Not supported; see the Notes on the class.

    Returns
    -------
    SciPyMargin
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
    if self._family is None:
      # `RuntimeError`, as the property readers use for the same missing
      # state -- and deliberately not `ValueError`, which `fit_margin`'s
      # `on_failure="fallback"` catches: a margin with no family is a misuse,
      # not a variable no family fits.
      raise RuntimeError(
        "SciPyMargin() has no family yet; name one at construction, or call "
        "select(y) to choose from the candidate set"
      )
    if weights is not None:
      raise TypeError(
        f"SciPyMargin({self._family!r}) cannot use observation weights: "
        "SciPy's estimators do not accept them. Pass margins='kde' for a "
        "weighted fit, or drop weights="
      )
    data = validate_univariate(np.asarray(y, dtype=float))
    data = data[~np.isnan(data)]
    if data.size == 0:
      raise ValueError(
        f"SciPyMargin({self._family!r}).fit got no usable observation"
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
    self._nobs = int(data.size)
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
      return f"SciPyMargin({self._family!r}, unfitted, {pinned})"
    shown = ", ".join(
      f"{name}={value:.4g}" for name, value in zip(self._names, self._params)
    )
    return f"SciPyMargin({self._family!r}, {shown})"


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


# The exclusion table is rendered into the class docstring rather than restated
# there, so the documented list and the honored one cannot drift apart.
SciPyMargin.__doc__ = re.sub(
  r"^([ \t]*)\{excluded\}$",
  lambda m: _excluded_block(m.group(1)),
  SciPyMargin.__doc__ or "",
  flags=re.MULTILINE,
)
