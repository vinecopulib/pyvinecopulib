"""Empirical margin — the rank transform, as a first-class distribution."""

from __future__ import annotations

from typing import Any, Optional

import numpy as np

from ..core import MarginBase

__all__ = ["EmpiricalMargin"]

# Raised by `pdf` / `logpdf`, and so by everything downstream of them. Spelled
# out because the mistake it catches is silent otherwise: the natural use of
# this margin is on continuous data, where a mass is not a density, and a wrong
# number would survive all the way into a log-likelihood.
_NO_DENSITY = (
  "EmpiricalMargin has no density: the empirical distribution is atomic, so "
  "there is no density with respect to Lebesgue measure and no meaningful "
  "log-likelihood on the original scale. Use Kde1d for a density — it is the "
  'default, and `margins="kde"` selects it explicitly. EmpiricalMargin remains '
  "valid for cdf, icdf and sample, and as the rank transform feeding a copula."
)


class EmpiricalMargin(MarginBase[np.ndarray]):
  r"""Rescaled empirical distribution function of the observed sample.

  Fitting a vine copula to :func:`pyvinecopulib.utils.to_pseudo_obs` output is
  the same model as fitting a vine *distribution* with these margins: that
  helper is exactly ``n / (n + 1)`` times the empirical cdf, and the ``n + 1``
  denominator is the standard device that keeps the pseudo-observations off the
  boundary, where pair-copula densities diverge. Making it a margin turns the
  long-standing copula-only workflow into one configuration of the general one.

  Parameters
  ----------
  weights : array, shape (n,), or None, optional
      Ignored here; pass weights to :meth:`fit` instead.

  See Also
  --------
  pyvinecopulib.utils.to_pseudo_obs : The equivalent free function.

  Notes
  -----
  Honest limitations. **There is no density**: the empirical distribution is
  atomic, so :meth:`pdf` and :meth:`logpdf` raise, and so does anything built
  on them — a :class:`~pyvinecopulib.core.Vinedist` on these margins has a
  ``cdf`` but no ``pdf``, ``logpdf`` or ``loglik``. Use ``Kde1d`` — the default
  — when a density is wanted. :meth:`icdf` is a step function, so sampling
  through this margin resamples observed values and never produces a new one.
  The whole sample is retained, so memory is ``O(n)``. Out of sample, ``cdf``
  saturates at ``1 / (n + 1)`` and ``n / (n + 1)`` rather than at 0 and 1 —
  which is the point.

  Ties are broken upward: ``cdf`` counts observations ``<= x``, matching
  :func:`~pyvinecopulib.utils.to_pseudo_obs` with ``ties_method="max"``. On a
  sample without ties — the case for continuous data — it agrees with that
  helper's ``"average"`` default exactly.

  What this margin is *for* is the copula: ``cdf`` on the observed sample is
  exactly :func:`~pyvinecopulib.utils.to_pseudo_obs`, so a
  :class:`~pyvinecopulib.core.Vinedist` on these margins fits and evaluates the
  same copula as the long-standing copula-only workflow, and
  :meth:`~pyvinecopulib.core.Vinedist.copula_data` reproduces its input.
  """

  supports_weights: bool = True

  def __init__(self, weights: Optional[Any] = None) -> None:
    del weights  # accepted for symmetry with the other margins; use fit()
    self._sorted: Optional[np.ndarray] = None
    self._cum: Optional[np.ndarray] = None
    self._total: float = 0.0

  @property
  def is_fitted(self) -> bool:
    return self._sorted is not None

  @property
  def var_type(self) -> str:
    """Variable type, reported as continuous.

    The empirical distribution is atomic, but the ``n / (n + 1)`` scaling is
    precisely the device for treating its output as continuous, and that is what
    fitting a vine copula to pseudo-observations means. Declaring ``"d"`` would
    instead fit a discrete copula against a left-limit column, which is a
    different model.

    Returns
    -------
    str
        Always ``"c"``.
    """
    return "c"

  @property
  def support(self) -> tuple[float, float]:
    if self._sorted is None:
      return (float("-inf"), float("inf"))
    return (float(self._sorted[0]), float(self._sorted[-1]))

  @property
  def n_parameters(self) -> float:
    """Effective parameter count, taken to be the sample size.

    Returns
    -------
    float
        ``n``; an empirical margin spends one parameter per observation, so it
        is never competitive on an information criterion and should not be
        compared with parametric families on one.
    """
    return float(self._values().size)

  @property
  def family_name(self) -> str:
    """Name used for this margin in selection reports.

    Returns
    -------
    str
        Always ``"empirical"``.
    """
    return "empirical"

  def _values(self) -> np.ndarray:
    """Return the sorted sample, or raise if unfitted.

    Returns
    -------
    array, shape (m,), dtype float
        The distinct observed values, ascending.

    Raises
    ------
    RuntimeError
        If the margin has not been fitted.
    """
    if self._sorted is None:
      raise RuntimeError("EmpiricalMargin is not fitted; call fit(x)")
    return self._sorted

  def fit(self, x: Any, *, weights: Optional[Any] = None) -> "EmpiricalMargin":
    """Record the sample.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations; NaNs are dropped, as
        :func:`~pyvinecopulib.utils.to_pseudo_obs` drops them.
    weights : array, shape (n,), or None, optional
        Observation weights.

    Returns
    -------
    EmpiricalMargin
        ``self``.
    """
    x_arr = np.asarray(x, dtype=float).ravel()
    w = (
      np.ones_like(x_arr)
      if weights is None
      else np.asarray(weights, dtype=float).ravel()
    )
    keep = ~np.isnan(x_arr)
    x_arr, w = x_arr[keep], w[keep]

    order = np.argsort(x_arr, kind="stable")
    x_sorted, w_sorted = x_arr[order], w[order]
    # Collapse ties so `cdf` is a lookup on distinct values.
    self._sorted, first = np.unique(x_sorted, return_index=True)
    self._cum = np.add.reduceat(w_sorted, first).cumsum()
    self._total = float(w.sum())
    return self

  def cdf(self, x: Any) -> np.ndarray:
    values = self._values()
    assert self._cum is not None  # guaranteed once fitted
    idx = np.searchsorted(values, np.asarray(x, dtype=float), side="right")
    below = np.concatenate(([0.0], self._cum))
    return below[idx] / (self._total + 1.0)

  def cdf_left(self, x: Any) -> np.ndarray:
    """Left limit of the empirical cdf.

    Counts strictly below ``x``, so the jump at an observed value is exactly
    its mass. Overrides the base class rather than stepping back a lattice
    point, since the observed values need not lie on one.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations.

    Returns
    -------
    array, shape (n,), dtype float
        Left limits.
    """
    values = self._values()
    assert self._cum is not None
    idx = np.searchsorted(values, np.asarray(x, dtype=float), side="left")
    below = np.concatenate(([0.0], self._cum))
    return below[idx] / (self._total + 1.0)

  def pdf(self, x: Any) -> np.ndarray:
    """Raise: the empirical distribution has no density.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Unused.

    Returns
    -------
    array, shape (n,), dtype float
        Never returns.

    Raises
    ------
    NotImplementedError
        Always.
    """
    del x
    raise NotImplementedError(_NO_DENSITY)

  def logpdf(self, x: Any) -> np.ndarray:
    """Raise: the empirical distribution has no density.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Unused.

    Returns
    -------
    array, shape (n,), dtype float
        Never returns.

    Raises
    ------
    NotImplementedError
        Always.
    """
    del x
    raise NotImplementedError(_NO_DENSITY)

  def icdf(self, p: Any) -> np.ndarray:
    values = self._values()
    assert self._cum is not None
    probs = np.asarray(p, dtype=float)
    targets = np.clip(probs, 0.0, 1.0) * (self._total + 1.0)
    idx = np.searchsorted(self._cum, targets, side="left")
    return values[np.clip(idx, 0, values.size - 1)]

  def _sample_uniform(self, n: int, seeds: list[int]) -> np.ndarray:
    """Draw ``n`` uniforms for the inherited ``sample``.

    Parameters
    ----------
    n : int
        Number of draws.
    seeds : list of int
        RNG seeds.

    Returns
    -------
    array, shape (n,), dtype float
        Uniforms in ``(0, 1)``.
    """
    rng = np.random.default_rng(seeds if seeds else None)
    return rng.uniform(size=n)

  def __repr__(self) -> str:
    state = f"n={int(self._total)}" if self.is_fitted else "unfitted"
    return f"EmpiricalMargin({state})"
