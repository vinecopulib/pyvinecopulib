"""Kernel-density margin, wrapping :class:`pyvinecopulib.utils.Kde1d`."""

from __future__ import annotations

import math
from typing import Any, Optional

import numpy as np

from ..core import MarginBase
from ..utils import Kde1d

__all__ = ["Kde1dMargin"]

#: kde1d's spellings of its variable types, mapped onto the margin vocabulary.
_TYPE_TO_VAR_TYPE = {
  "continuous": "c",
  "discrete": "d",
  "zero-inflated": "zi",
  "zero_inflated": "zi",
}


class Kde1dMargin(MarginBase[np.ndarray]):
  """Boundary-corrected kernel-density margin, the library's default.

  Wraps :class:`pyvinecopulib.utils.Kde1d`, whose ``pdf`` already returns the
  density with respect to the right reference measure in each of its three
  modes — a density when continuous, a probability mass when discrete, and the
  point mass at zero when zero-inflated — which is exactly what
  :class:`~pyvinecopulib.core.MarginLike` asks for. The inherited ``cdf_left``
  therefore needs no override either: stepping back a lattice point is correct
  for ``"discrete"``, and removing ``pdf(0)`` from ``cdf(0)`` is correct for
  ``"zero-inflated"``.

  Parameters
  ----------
  xmin, xmax : float or None, optional
      Bounds of the support; ``None`` means unbounded. Worth setting whenever
      they are known: a bound switches on the boundary-corrected estimator, and
      without one the fitted grid is padded past the data, which for count data
      puts mass on values that cannot occur.
  type : {"continuous", "discrete", "zero-inflated"}, optional
      Variable type.
  multiplier : float, optional
      Bandwidth multiplier.
  bandwidth : float or None, optional
      Bandwidth; ``None`` selects it automatically at each fit.
  degree : int, optional
      Degree of the local polynomial (0, 1 or 2).
  grid_size : int, optional
      Number of interpolation grid points.

  See Also
  --------
  pyvinecopulib.utils.Kde1d : The estimator this wraps.
  pyvinecopulib.core.MarginBase : The contract it fills in.

  Notes
  -----
  Each :meth:`fit` builds a fresh ``Kde1d``, so refitting one margin on new
  data re-selects the bandwidth. Refitting a ``Kde1d`` in place does not: it
  feeds the previous fit's bandwidth back into the selector and silently reuses
  it (fixed upstream in kde1d-cpp#28). That is also why broadcasting one margin
  across several columns must copy it rather than share it.
  """

  supports_weights: bool = True

  def __init__(
    self,
    *,
    xmin: Optional[float] = None,
    xmax: Optional[float] = None,
    type: str = "continuous",
    multiplier: float = 1.0,
    bandwidth: Optional[float] = None,
    degree: int = 2,
    grid_size: int = 400,
  ) -> None:
    if type not in _TYPE_TO_VAR_TYPE:
      raise ValueError(
        f"unknown type={type!r}; expected one of "
        f"{sorted(set(_TYPE_TO_VAR_TYPE))}"
      )
    self.xmin = xmin
    self.xmax = xmax
    self.type = type
    self.multiplier = multiplier
    self.bandwidth = bandwidth
    self.degree = degree
    self.grid_size = grid_size
    self._kde: Optional[Kde1d] = None

  @classmethod
  def from_kde1d(cls, kde: Kde1d) -> "Kde1dMargin":
    """Adopt an already-fitted ``Kde1d``.

    Parameters
    ----------
    kde : Kde1d
        A fitted estimator.

    Returns
    -------
    Kde1dMargin
        A margin evaluating through ``kde``.

    Raises
    ------
    ValueError
        If ``kde`` has not been fitted.
    """
    if not kde.is_fitted:
      raise ValueError("Kde1dMargin.from_kde1d needs a fitted Kde1d")
    margin = cls(
      xmin=None if math.isnan(kde.xmin) else kde.xmin,
      xmax=None if math.isnan(kde.xmax) else kde.xmax,
      type=kde.type,
    )
    margin._kde = kde
    return margin

  @property
  def kde1d(self) -> Kde1d:
    """The underlying fitted estimator.

    Returns
    -------
    Kde1d
        The estimator built by :meth:`fit`.

    Raises
    ------
    RuntimeError
        If the margin has not been fitted.
    """
    if self._kde is None:
      raise RuntimeError(
        "Kde1dMargin is not fitted; call fit(x) or use from_kde1d(...)"
      )
    return self._kde

  @property
  def is_fitted(self) -> bool:
    return self._kde is not None and self._kde.is_fitted

  @property
  def var_type(self) -> str:
    return _TYPE_TO_VAR_TYPE[self.type]

  @property
  def support(self) -> tuple[float, float]:
    lo = float("-inf") if self.xmin is None else float(self.xmin)
    hi = float("inf") if self.xmax is None else float(self.xmax)
    return (lo, hi)

  def fit(self, x: Any, *, weights: Optional[Any] = None) -> "Kde1dMargin":
    """Fit the kernel density estimate to ``x``.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations.
    weights : array, shape (n,), or None, optional
        Observation weights.

    Returns
    -------
    Kde1dMargin
        ``self``.
    """
    kde = Kde1d(
      xmin=self.xmin,
      xmax=self.xmax,
      type=self.type,
      multiplier=self.multiplier,
      bandwidth=self.bandwidth,
      degree=self.degree,
      grid_size=self.grid_size,
    )
    x_arr = np.asarray(x, dtype=float).ravel()
    if weights is None:
      kde.fit(x_arr)
    else:
      kde.fit(x_arr, np.asarray(weights, dtype=float).ravel())
    self._kde = kde
    return self

  def pdf(self, x: Any) -> np.ndarray:
    return np.asarray(self.kde1d.pdf(np.asarray(x, dtype=float)))

  def cdf(self, x: Any) -> np.ndarray:
    return np.asarray(self.kde1d.cdf(np.asarray(x, dtype=float)))

  def icdf(self, p: Any) -> np.ndarray:
    return np.asarray(self.kde1d.icdf(np.asarray(p, dtype=float)))

  def simulate(
    self, n: int, *, seeds: Optional[list[int]] = None
  ) -> np.ndarray:
    """Draw ``n`` samples from the fitted estimate.

    Parameters
    ----------
    n : int
        Number of samples.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n,), dtype float
        Samples.
    """
    return np.asarray(
      self.kde1d.simulate(n, seeds=list(seeds) if seeds else [])
    )

  @property
  def fitted_loglik(self) -> float:
    """Log-likelihood attained at the fit.

    Named apart from :meth:`~pyvinecopulib.core.MarginBase.loglik`, which is a
    *method* evaluating the log-likelihood of data: one name cannot be both. The
    value is the estimator's own, including the boundary-transform and
    zero-inflation corrections, so it is not reproducible by calling ``loglik``
    on the training sample.

    Returns
    -------
    float
        The fitted log-likelihood.
    """
    return float(self.kde1d.loglik)

  @property
  def n_parameters(self) -> float:
    """Effective degrees of freedom of the fitted estimate.

    Returns
    -------
    float
        The trace-of-smoother effective d.f., comparable with a parametric
        parameter count only as a heuristic.
    """
    return float(self.kde1d.edf)

  @property
  def family_name(self) -> str:
    """Name used for this margin in selection reports.

    Returns
    -------
    str
        Always ``"kde1d"``.
    """
    return "kde1d"

  def __repr__(self) -> str:
    state = "fitted" if self.is_fitted else "unfitted"
    return f"Kde1dMargin(type={self.type!r}, {state})"
