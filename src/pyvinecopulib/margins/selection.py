"""Family selection for univariate margins.

The estimator is the selection: :class:`MarginSelector` fits every admissible
candidate, keeps the winner and forwards evaluation to it, the same shape
``GridSearchCV`` gives to hyper-parameter search.
"""

from __future__ import annotations

import re
import textwrap
import time
import warnings
from typing import Any, Iterable, Optional, Sequence

import numpy as np

from ..core import MarginBase
from ..core.margin_base import _reject_covariates
from ..core import Kde1d
from .parametric import (
  _EXCLUDED,
  _LIMITS,
  _PARTITIONS,
  ParametricMargin,
  _curated_margin,
)

__all__ = ["MarginSelector"]

#: Information criteria, as functions of the log-likelihood, the number of
#: freely estimated parameters and the sample size. All three are minimized.
_CRITERIA: tuple[str, ...] = ("aic", "bic", "aicc")

#: What :meth:`MarginSelector.fit` does when no candidate is admissible.
#: ``"raise"`` reports every family and its cause; ``"fallback"`` substitutes
#: a kernel-density margin and warns once.
_ON_FAILURE = ("raise", "fallback")


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
    return dict.fromkeys(_CRITERIA, float("inf"))
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


def _dedupe(candidates: Iterable[Any]) -> list[Any]:
  """Drop candidates that would tie with one already present.

  Two candidates of the same family with the same pinned parameters fit to the
  same numbers, so they tie exactly on every criterion: the winner becomes an
  artifact of iteration order, and the report carries the row twice. Anything
  whose identity cannot be read this way is kept, since dropping it would be a
  guess.

  Parameters
  ----------
  candidates : iterable
      Unfitted candidate margins, in preference order.

  Returns
  -------
  list
      The first candidate of each distinct family-and-pins, order preserved.
  """
  seen: set[tuple[str, tuple[tuple[str, float], ...]]] = set()
  out: list[Any] = []
  for margin in candidates:
    family = getattr(margin, "family_name", None)
    if family is None:
      out.append(margin)
      continue
    fixed = getattr(margin, "fixed_parameters", None) or {}
    key = (str(family), tuple(sorted(fixed.items())))
    if key in seen:
      continue
    seen.add(key)
    out.append(margin)
  return out


class _SelectorBase(MarginBase[np.ndarray]):
  """State and forwarding shared by the family selectors.

  A selector is a margin that fits several candidates, keeps the winner on
  ``selected_`` and forwards evaluation to it. Everything in that sentence
  except *which* candidates and *how* one is fitted is the same whichever
  ecosystem supplies the families, so it lives here: the fitted state, the
  report, and the eleven members that forward to the winner. Subclasses own
  ``__init__``, ``fit`` and the candidate machinery.

  Names in messages come from ``type(self).__name__``, so a subclass needs no
  override to identify itself.
  """

  supports_weights = False

  #: Set by every subclass's ``__init__``; declared here because the shared
  #: members below read them.
  criterion: str
  _selected: Optional[Any]
  _report: list[dict[str, Any]]

  # --- fitted state -------------------------------------------------------- #

  @property
  def selected_(self) -> Any:
    """The winning margin.

    Returns
    -------
    MarginLike
        The candidate with the smallest ``criterion``, or the fallback margin
        when every candidate was rejected.

    Raises
    ------
    RuntimeError
        If the selector has not been fitted.
    """
    if self._selected is None:
      raise RuntimeError(f"{type(self).__name__} is not fitted; call fit(y)")
    return self._selected

  @property
  def report_(self) -> list[dict[str, Any]]:
    """One row per candidate, in the order they were tried.

    Each row carries ``column``, ``family``, ``n_parameters``, ``loglik``,
    ``aic``, ``bic``, ``aicc``, ``criterion``, ``support``, ``seconds``,
    ``warnings`` and ``status`` / ``selected``. ``status`` is ``"ok"`` for an
    admissible candidate, ``"selected"`` for the winner, ``"fallback"`` for the
    margin used when nothing else survived, and otherwise the reason it was
    rejected. Rows from several variables concatenate into one tidy table, e.g.
    ``pandas.DataFrame(rows)``.

    Returns
    -------
    list of dict
        The table; empty before :meth:`fit`.
    """
    return list(self._report)

  @property
  def is_fitted(self) -> bool:
    return self._selected is not None

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
  def _fitted_loglik(self) -> float:
    """Log-likelihood attained by the selected margin.

    Returns
    -------
    float
        The winner's fitted log-likelihood.
    """
    return float(self.selected_.loglik())

  @property
  def var_type(self) -> str:
    return str(getattr(self.selected_, "var_type", "c"))

  @property
  def support(self) -> tuple[float, float]:
    lo, hi = getattr(self.selected_, "support", (float("-inf"), float("inf")))
    return (float(lo), float(hi))

  def pdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return np.asarray(self.selected_.pdf(y), dtype=float)

  def cdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return np.asarray(self.selected_.cdf(y), dtype=float)

  def icdf(self, p: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return np.asarray(self.selected_.icdf(p), dtype=float)

  def logpdf(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return np.asarray(self.selected_.logpdf(y), dtype=float)

  def cdf_left(self, y: Any, *, x: Optional[Any] = None) -> np.ndarray:
    return np.asarray(self.selected_.cdf_left(y), dtype=float)

  def sample(
    self, n: int, *, x: Optional[Any] = None, seeds: Optional[list[int]] = None
  ) -> np.ndarray:
    """Draw ``n`` samples from the selected margin.

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
    return np.asarray(self.selected_.sample(n, seeds=seeds), dtype=float)

  def __repr__(self) -> str:
    name = type(self).__name__
    if self._selected is None:
      return f"{name}(criterion={self.criterion!r}, unfitted)"
    return (
      f"{name}(criterion={self.criterion!r}, selected={self.family_name!r})"
    )


class MarginSelector(_SelectorBase):
  """Pick a univariate family by an information criterion.

  A selector is itself a margin: :meth:`fit` estimates every admissible
  candidate, keeps the best one on :attr:`selected_`, and forwards ``pdf`` /
  ``cdf`` / ``icdf`` to it, so a selected margin drops into a
  :class:`~pyvinecopulib.core.Vinedist` exactly like a fixed one. The
  per-candidate table is on :attr:`report_`.

  Parameters
  ----------
  candidates : str or sequence, optional
      ``"auto"`` (the default) uses the curated families whose support the data
      are compatible with. A sequence selects explicitly, and may mix family
      names with ready-made margins, e.g. ``["gamma", "lognorm"]``, which is
      how a family kept out of the curated set is opted into.
  criterion : {"aic", "bic", "aicc"}, optional
      What to minimize.
  bounds : tuple of float, or None, optional
      A known support ``(a, b)``. Given, the candidates are the families that
      live on a bounded interval, anchored there, and the families with
      unbounded support are not considered — they would put mass outside a
      support the caller has stated.
  var_type : {"c", "d"} or None, optional
      Force the data to be read as continuous or as counts. Inferred when
      ``None``: non-negative integer-valued data are counts.
  on_failure : {"raise", "fallback"}, optional
      What to do when every candidate is rejected. ``"raise"`` (the default)
      reports each family and why it was rejected. ``"fallback"`` substitutes a
      kernel-density margin and warns once -- available, but not the default,
      because answering a parametric request with a nonparametric fit is a
      silent downgrade, and a warning is easy to lose in a loop over fifty
      columns.
  name : str or None, optional
      Name of the variable being fitted, reported in :attr:`report_` and in the
      fallback warning.

  Raises
  ------
  ValueError
      If ``criterion`` is not one of the three, ``var_type`` is not ``"c"`` or
      ``"d"``, or ``on_failure`` is neither.

  See Also
  --------
  pyvinecopulib.margins.ParametricMargin : The candidates.
  pyvinecopulib.core.Kde1d : The fallback, and the library default.

  Notes
  -----
  **The candidate set is ours, not SciPy's.** Sweeping all 110 continuous
  families is worse than useless: on clean gamma data it ranks ``vonmises``
  first, ahead of the true family by 701 AIC units, because ``vonmises``
  advertises the whole real line while its density integrates to 63.75 there.
  ``levy_stable`` does not return. So the curated set is small and grouped by
  support:

  * **real line** — ``norm``, ``t``, ``logistic``, ``laplace``, ``skewnorm``,
    ``gumbel_r``, with ``loc`` and ``scale`` free.
  * **positive** — ``gamma``, ``lognorm``, ``weibull_min``, ``invgauss``,
    ``expon``, anchored at ``loc = 0``.
  * **unit interval** — ``beta``, anchored at ``loc = 0``, ``scale = 1``.
  * **bounded** on a given ``(a, b)`` — ``uniform`` and ``beta``, anchored at
    ``loc = a``, ``scale = b - a``.
  * **counts** — ``poisson``, ``nbinom``, ``geom``, with ``loc = 0``.

  Every ``loc`` above is pinned for two reasons at once. It is the boundary
  escape a free ``loc`` invites, and it is one fewer parameter in the
  criterion: the same gamma scores 4289.5 with a free ``loc`` and 4287.8 with
  ``loc = 0``, which is enough to reverse a close comparison.

  Counts form their own group and are never compared against continuous
  families, whose log-likelihoods are on a different scale. ``geom`` is
  parameterized as a number of trials, so it is inadmissible — and reported as
  such — for data containing 0; ``nbinom`` covers the zero-inclusive
  geometric. Every other family SciPy offers is left out on purpose, and each
  exclusion carries its reason:

  {excluded}

  **Nothing is skipped silently.** A candidate that fails to fit, or that fits
  to something inadmissible, still gets a row in :attr:`report_` with the
  reason. When every candidate fails, the selector falls back to
  ``Kde1d`` and warns once, naming the variable and the reasons —
  never to ``norm``, which would misspecify silently. Marginal
  misspecification is the more damaging of the two errors here: it distorts the
  probability transform and biases the copula, where a conservative margin does
  not.

  **Observation weights are not supported**, because the parametric candidates
  cannot use them; see the notes on :class:`ParametricMargin`.

  Examples
  --------
  ::

      import numpy as np
      from pyvinecopulib.margins import MarginSelector

      y = np.random.default_rng(0).gamma(2.5, 1.5, size=500)
      sel = MarginSelector(criterion="bic").fit(y)
      sel.selected_.family_name              # -> 'gamma'
      [(r["family"], r["status"]) for r in sel.report_]
  """

  supports_weights: bool = False

  def __init__(
    self,
    candidates: Any = "auto",
    *,
    criterion: str = "aic",
    bounds: Optional[tuple[float, float]] = None,
    var_type: Optional[str] = None,
    on_failure: str = "raise",
    name: Optional[str] = None,
  ) -> None:
    if criterion not in _CRITERIA:
      raise ValueError(
        f"unknown criterion={criterion!r}; expected one of {list(_CRITERIA)}"
      )
    if on_failure not in _ON_FAILURE:
      raise ValueError(
        f"unknown on_failure={on_failure!r}; expected one of "
        f"{list(_ON_FAILURE)}"
      )
    if var_type not in (None, "c", "d"):
      raise ValueError(
        f"unknown var_type={var_type!r}; expected 'c', 'd' or None"
      )
    self.candidates = candidates
    self.criterion = criterion
    self.bounds = (
      None if bounds is None else (float(bounds[0]), float(bounds[1]))
    )
    self._var_type = var_type
    self.on_failure = on_failure
    # A constructor argument is authoritative; `declare` fills in only what the
    # caller left open.
    self._pinned = (var_type is not None, bounds is not None)
    self.name = name
    self._selected: Optional[Any] = None
    self._report: list[dict[str, Any]] = []

  #: The selector partitions its candidates by kind, so it serves either.
  supported_var_types: tuple[str, ...] = ("c", "d")

  def declare(
    self,
    *,
    var_type: Optional[str] = None,
    support: Optional[tuple[float, float]] = None,
  ) -> "MarginSelector":
    """Adopt the caller's variable type and support where none was pinned.

    This is the selector's only route to information it cannot recover from the
    sample. Without it, a column the caller knows to be a bounded count is
    retyped by ``_is_discrete``'s heuristic and scored against candidates
    with no bounds at all.

    Parameters
    ----------
    var_type : str or None, optional
        ``"c"``, ``"d"`` or ``"zi"``. ``"zi"`` is reduced to ``"d"``, which is
        what the candidate partition distinguishes.
    support : tuple of float, or None, optional
        Declared bounds. Only a pair finite at both ends is adopted, since
        ``bounds`` selects the bounded candidate group, whose families need a
        finite interval to be rescaled onto. A half-bounded support is still
        read off the data.

    Returns
    -------
    MarginSelector
        ``self``, so the call chains into :meth:`fit`.
    """
    pinned_type, pinned_bounds = self._pinned
    if var_type is not None and not pinned_type:
      self._var_type = "d" if var_type == "zi" else var_type
    if support is not None and not pinned_bounds:
      lo, hi = float(support[0]), float(support[1])
      if np.isfinite(lo) and np.isfinite(hi):
        self.bounds = (lo, hi)
    return self

  # --- fitted state -------------------------------------------------------- #

  # --- estimation ---------------------------------------------------------- #

  def fit(
    self, y: Any, *, x: Optional[Any] = None, weights: Optional[Any] = None
  ) -> "MarginSelector":
    """Fit every candidate and keep the best.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations; NaNs are dropped.
    x : array, shape (n, k), or None, optional
        Not supported; passing covariates raises rather than silently
        fitting an unconditional margin.
    weights : array, shape (n,), or None, optional
        Not supported; the parametric candidates cannot use them.

    Returns
    -------
    MarginSelector
        ``self``.

    Raises
    ------
    TypeError
        If ``weights`` is given.
    ValueError
        If no observation survives, or -- with ``on_failure="raise"`` -- if every
        candidate was rejected.

    Warns
    -----
    UserWarning
        Once, when every candidate is rejected and ``on_failure="fallback"``
        substitutes a kernel-density margin.
    """
    _reject_covariates(self, x)
    if weights is not None:
      raise TypeError(
        "MarginSelector cannot use observation weights: its parametric "
        "candidates are fitted by SciPy, which does not accept them. Pass "
        "margins='kde' for a weighted fit, or drop weights="
      )
    data = np.asarray(y, dtype=float).ravel()
    data = data[~np.isnan(data)]
    if data.size == 0:
      raise ValueError("MarginSelector.fit got no usable observation")

    counts = self._reads_as_counts(data)
    rows: list[dict[str, Any]] = []
    admissible: list[tuple[float, int, Any]] = []
    for index, candidate in enumerate(
      self._candidate_margins(data, counts=counts)
    ):
      row, ok = self._try(candidate, data)
      rows.append(row)
      if ok:
        # The index breaks ties toward the earlier candidate, which is what
        # keeps the outcome independent of dictionary iteration order.
        admissible.append((row[self.criterion], index, candidate))

    if admissible:
      _, index, winner = min(admissible, key=lambda item: item[:2])
      rows[index]["status"] = "selected"
      rows[index]["selected"] = True
    else:
      winner = self._fallback(data, rows, counts=counts)
    self._selected = winner
    self._report = rows
    return self

  def _reads_as_counts(self, data: np.ndarray) -> bool:
    """Whether the data should be fitted with the count families.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations.

    Returns
    -------
    bool
        ``True`` for the count group.
    """
    if self._var_type is not None:
      return self._var_type == "d"
    return bool(np.all(data >= 0.0) and np.all(data == np.round(data)))

  def _candidate_margins(self, data: np.ndarray, *, counts: bool) -> list[Any]:
    """Build the candidates to fit.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations, used to decide which support groups apply.
    counts : bool
        Whether the count group applies.

    Returns
    -------
    list
        Unfitted candidate margins.

    Raises
    ------
    ValueError
        If ``candidates`` is neither ``"auto"`` nor a sequence.
    """
    if self.candidates == "auto":
      # `unit` and `bounded` both offer `beta`, differing only in its pins, and
      # a future partition may overlap outright. Deduplicating here keeps that
      # from producing two candidates that tie exactly on every criterion,
      # which would make the winner an artifact of iteration order.
      return _dedupe(
        _curated_margin(family, group, self.bounds)
        for group in self._groups(data, counts=counts)
        for family in _PARTITIONS[group]
      )
    if isinstance(self.candidates, str):
      raise ValueError(
        f"unknown candidates={self.candidates!r}; expected 'auto' or a "
        "sequence of family names and margins"
      )
    out: list[Any] = []
    dropped: list[str] = []
    want = "d" if counts else "c"
    for entry in self.candidates:
      margin = ParametricMargin(entry) if isinstance(entry, str) else entry
      # The partition the "auto" path applies, applied here too: a count family
      # reports a probability mass and a continuous one a Lebesgue density, so
      # ranking both on one criterion compares numbers that are not
      # commensurable -- and the density usually wins.
      #
      # A family that serves both says so through `supported_var_types`, which
      # is a property of the family; `var_type` only reports what one instance
      # currently is. Such a candidate is given the wanted type and then
      # re-checked, so a family that declares the set but ignores `declare` is
      # dropped rather than fitted as the wrong kind.
      served = getattr(
        margin, "supported_var_types", (getattr(margin, "var_type", "c"),)
      )
      if want in served:
        configure = getattr(margin, "declare", None)
        if configure is not None:
          configure(var_type=want)
      if getattr(margin, "var_type", "c") == want:
        out.append(margin)
      else:
        dropped.append(str(getattr(margin, "family_name", margin)))
    out = _dedupe(out)
    reads = "counts" if counts else "continuous"
    if dropped and not out:
      raise ValueError(
        f"the data read as {reads}, and every candidate is of the other kind: "
        f"{', '.join(dropped)}. Pass var_type to say how the data should be "
        "read, or candidates that match it."
      )
    if dropped:
      warnings.warn(
        f"the data read as {reads}, so {len(dropped)} candidate(s) of the "
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
    if self.bounds is not None:
      return ["bounded"]
    groups = ["real"]
    if data.min() > 0.0:
      groups.append("positive")
      if data.max() < 1.0:
        groups.append("unit")
    return groups

  def _try(
    self, candidate: Any, data: np.ndarray
  ) -> tuple[dict[str, Any], bool]:
    """Fit one candidate and score it.

    Parameters
    ----------
    candidate : MarginLike
        An unfitted candidate.
    data : array, shape (n,), dtype float
        The observations.

    Returns
    -------
    tuple
        The report row, and whether the candidate is admissible.
    """
    family = str(getattr(candidate, "family_name", type(candidate).__name__))
    started = time.perf_counter()
    # A candidate's optimizer chatter says nothing a caller can act on -- the
    # outcome is the row below either way -- so collect it into the row rather
    # than emitting one warning per family tried.
    fitted = True
    with warnings.catch_warnings(record=True) as caught:
      warnings.simplefilter("always")
      try:
        candidate.fit(data)
        status = _reject(candidate, data)
      except Exception as e:  # noqa: BLE001 - any fitter failure is a rejection
        status, fitted = f"{type(e).__name__}: {e}", False
    elapsed = time.perf_counter() - started

    ok = status is None
    # A rejected candidate still reports the log-likelihood and criteria it
    # reached, so a table row explains why it lost rather than only that it did.
    # A candidate that never fitted has neither; most margins raise rather than
    # report a parameter count before they have parameters.
    loglik, k, support = float("nan"), float("nan"), None
    if fitted:
      loglik = float(np.sum(candidate.logpdf(data)))
      k = float(getattr(candidate, "n_parameters", float("nan")))
      support = getattr(candidate, "support", None)
    row: dict[str, Any] = {
      "column": self.name,
      "family": family,
      "n_parameters": k,
      "loglik": loglik,
      "support": support,
      "seconds": elapsed,
      "warnings": [str(w.message) for w in caught],
      "status": "ok" if ok else status,
      "selected": False,
    }
    row.update(_criteria(loglik, k, data.size))
    return row, ok

  def _fallback(
    self, data: np.ndarray, rows: list[dict[str, Any]], *, counts: bool
  ) -> Any:
    """Fit the kernel-density margin used when no candidate survives.

    Parameters
    ----------
    data : array, shape (n,), dtype float
        The observations.
    rows : list of dict
        The report so far; a row for the fallback is appended.
    counts : bool
        Whether the data were read as counts.

    Returns
    -------
    Kde1d
        The fitted fallback.

    Raises
    ------
    ValueError
        If ``on_failure`` is ``"raise"``, naming every rejected family and its
        cause.
    """
    reasons = "; ".join(f"{row['family']}: {row['status']}" for row in rows)
    where = f"variable {self.name!r}" if self.name else "this variable"
    if self.on_failure == "raise":
      raise ValueError(
        f"no parametric family was admissible for {where}. Rejected: "
        f"{reasons or 'no candidates'}. Pass candidates= to widen the search, "
        "var_type= or bounds= to say how the data should be read, or "
        'on_failure="fallback" to accept a kernel-density margin instead.'
      )
    warnings.warn(
      f"no parametric family was admissible for {where}; falling back to a "
      f"kernel-density margin. Rejected: {reasons or 'no candidates'}",
      UserWarning,
      stacklevel=3,
    )
    margin = Kde1d(type="discrete" if counts else "continuous")
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


def _excluded_block(indent: str) -> str:
  """Render the excluded-family table as an RST bullet list.

  The exclusions are data, not prose: rendering the docstring from them keeps
  the documented list and the list the selector actually honors identical.

  Parameters
  ----------
  indent : str
      Leading whitespace of the line the list replaces. Python 3.13 dedents
      docstrings at compile time and earlier versions do not, so the indent is
      read off the placeholder rather than assumed.

  Returns
  -------
  str
      A bullet list, one entry per excluded family.
  """
  return "\n".join(
    textwrap.fill(
      f"* ``{family}`` \N{EM DASH} {reason}",
      width=78,
      initial_indent=indent,
      subsequent_indent=indent + "  ",
    )
    for family, reason in _EXCLUDED.items()
  )


MarginSelector.__doc__ = re.sub(
  r"^([ \t]*)\{excluded\}$",
  lambda m: _excluded_block(m.group(1)),
  MarginSelector.__doc__ or "",
  flags=re.MULTILINE,
)
