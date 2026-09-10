"""What a margin plot draws, for every variable type it can describe.

``margin_plot`` is bound to ``Kde1d.plot`` by the binding, which looks it up by
module path, and is what ``MarginBase.plot`` forwards to -- so a compiled
``Kde1d``, a ``TorchKde1d``, a ``SciPyMargin`` and a hand-written subclass all
draw the same picture. It reads only ``var_type``, ``support``, ``is_fitted``
and the requested one of ``pdf`` / ``cdf``, each of which ``MarginBase``
supplies with a continuous-correct default, and prefers an optional
``grid_points`` for the x-range where the margin has one.
"""

from typing import Any, Callable, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from ..pyvinecopulib_ext import Kde1d
from ._covariates import covariate_row
from ._placement import to_numpy
from .protocols import ArrayT, MarginLike


#: `Kde1d` is named outright because it satisfies `MarginLike` *nominally*
#: only: its `pdf` / `cdf` take no covariate `x`, so a static check rejects it.
_Margin = Union[MarginLike[ArrayT], Kde1d]


#: Shared with `MarginBase.plot`, which adds `x` and a `Raises`.
MARGIN_PLOT_PARAMS = """    xlim : tuple of float, or None, optional
        Limits for the x axis. ``None`` uses the declared support where it is
        finite and pads the drawn range otherwise.
    ylim : tuple of float, or None, optional
        Limits for the y axis. ``None`` fits them to what was drawn.
    grid_size : int, default=200
        Number of grid points, for a continuous or zero-inflated variable. A
        discrete one is drawn on its lattice, whose length it does not choose.
    show_zero_mass : bool, default=True
        Whether a zero-inflated variable's point mass at zero is marked.
    kind : str, default="density"
        What to draw: ``"density"`` or ``"cdf"``.
"""

#: Shared with `MarginBase.plot`.
MARGIN_PLOT_SUMMARY = """
    Plot the density or the distribution function of this margin.

    Each variable type is drawn as what it is: a continuous variable as a
    curve, a discrete one as marks on its integer support, and a
    zero-inflated one as the curve with zero excised plus the one emphasized
    point mass.
"""

#: `Kde1d.plot`'s docstring, which the binding reads from here by name.
MARGIN_PLOT_DOC = (
  MARGIN_PLOT_SUMMARY
  + """
    Parameters
    ----------
"""
  + MARGIN_PLOT_PARAMS
  + """
    Returns
    -------
    None
        The figure is drawn with matplotlib.

    Examples
    --------
    >>> import numpy as np
    >>> import pyvinecopulib as pv
    >>> rng = np.random.default_rng(123)
    >>> kde = pv.core.Kde1d().fit(rng.beta(0.5, 2.0, 100))
    >>> kde.plot()
    >>> kde.plot(kind="cdf")
    >>> counts = pv.core.Kde1d(type="discrete")
    >>> counts.fit(rng.poisson(3, 100).astype(float))
    >>> counts.plot()
    >>> y = rng.exponential(2, 100)
    >>> y[rng.choice(100, 30, replace=False)] = 0.0
    >>> pv.core.Kde1d(xmin=0, type="zero-inflated").fit(y).plot()
"""
)

#: The two quantiles a margin with an unbounded support is drawn between,
#: reached only where there is no fitted grid to read the range off.
_TAIL = (1e-3, 1 - 1e-3)


def _support(margin: _Margin[Any]) -> tuple[float, float]:
  """The margin's declared support, as two floats."""
  lo, hi = getattr(margin, "support", (-np.inf, np.inf))
  return (float(lo), float(hi))


def _grid_points(margin: _Margin[Any]) -> Optional[np.ndarray]:
  """The fitted evaluation grid, where the margin publishes one."""
  gp = getattr(margin, "grid_points", None)
  if gp is None:
    return None
  points = to_numpy(gp)
  return points if points.size else None


def _draw_range(
  margin: _Margin[ArrayT],
  place: Optional[Callable[[np.ndarray], Any]],
) -> tuple[float, float]:
  """The interval to evaluate over.

  A declared bound wins wherever there is one. What fills an undeclared bound
  is the fitted grid where the margin publishes one -- that is the range the
  density was estimated on -- and otherwise the quantile at ``_TAIL``, since
  an unbounded support has no last point to draw.
  """
  lo, hi = _support(margin)
  points = _grid_points(margin)
  if points is not None:
    return (
      lo if np.isfinite(lo) else float(points.min()),
      hi if np.isfinite(hi) else float(points.max()),
    )
  if np.isfinite(lo) and np.isfinite(hi):
    return (lo, hi)
  icdf = getattr(margin, "icdf", None)
  if icdf is None:
    raise ValueError(
      "cannot choose an x-range: the margin's support is unbounded and it "
      "has no `icdf` to read a quantile from; pass `xlim=`"
    )
  p = np.asarray(_TAIL, dtype=float)
  tails = to_numpy(icdf(p if place is None else place(p)))
  return (
    lo if np.isfinite(lo) else float(tails[0]),
    hi if np.isfinite(hi) else float(tails[1]),
  )


def make_plotting_grid(
  margin: _Margin[ArrayT],
  grid_size: int = 200,
  *,
  place: Optional[Callable[[np.ndarray], Any]] = None,
) -> np.ndarray:
  """The points a margin's plot is evaluated on, by variable type."""
  var_type = getattr(margin, "var_type", "c")
  lo, hi = _draw_range(margin, place)

  if var_type == "d":
    # The integer support, which is what carries mass. A fitted grid runs half
    # a unit wider than that support, because that is where the jitter cells
    # end, so it is rounded to the nearest level rather than outwards -- which
    # would plot one level below and one above that cannot occur.
    first, last = int(np.round(lo)), int(np.round(hi))
    return np.arange(first, max(first, last) + 1, dtype=float)

  ev = np.linspace(lo, hi, grid_size)
  if var_type == "zi":
    # The atom at zero is drawn on its own, so the continuous part excludes it.
    ev = ev[ev != 0]
  return np.asarray(ev)


def margin_plot(
  margin: _Margin[ArrayT],
  xlim: Optional[tuple[float, float]] = None,
  ylim: Optional[tuple[float, float]] = None,
  grid_size: int = 200,
  show_zero_mass: bool = True,
  *,
  kind: str = "density",
  x: Optional[ArrayT] = None,
  place: Optional[Callable[[np.ndarray], Any]] = None,
) -> None:
  """{}""".format(MARGIN_PLOT_DOC)

  if kind not in ("density", "cdf"):
    raise ValueError(f"kind must be 'density' or 'cdf'; got {kind!r}")

  ## A margin that publishes a grid says it is unfitted by leaving it empty,
  ## which is the only signal one carrying no `is_fitted` gives.
  empty_grid = (
    getattr(margin, "grid_points", None) is not None
    and _grid_points(margin) is None
  )
  if not getattr(margin, "is_fitted", True) or empty_grid:
    raise ValueError("the margin must be fitted before plotting")

  ## A conditional margin's density is a different curve at every covariate
  ## value, so a 2-d plot shows one slice: a single row, repeated across the
  ## grid. Placed but not clamped -- covariates are reals.
  row: Optional[np.ndarray] = None
  if x is not None:
    if not getattr(margin, "supports_covariates", False):
      raise ValueError(
        "this margin declares no `supports_covariates`, so it models f(y) "
        "and there is no covariate value to draw at; drop `x=`"
      )
    row = covariate_row(np.asarray(x, dtype=float))

  var_type = getattr(margin, "var_type", "c")
  ev = make_plotting_grid(margin, grid_size, place=place)

  ## The grid is manufactured here, so this is the one place a margin is
  ## handed an array it did not supply and cannot infer a namespace from.
  ## `place` brings it onto the margin's own namespace, dtype and device; a
  ## compiled `Kde1d` passes none and keeps the NumPy grid it always had.
  def evaluate(y: np.ndarray) -> np.ndarray:
    ya: Any = y if place is None else place(y)
    kwargs: dict[str, Any] = {}
    if row is not None:
      tiled = np.repeat(row, y.shape[0], axis=0)
      kwargs["x"] = tiled if place is None else place(tiled)
    fn = margin.pdf if kind == "density" else margin.cdf
    return to_numpy(fn(ya, **kwargs))

  vals = evaluate(ev)
  zero = evaluate(np.zeros(1)) if var_type == "zi" and show_zero_mass else None

  ## A discrete margin carries mass only on its lattice, so it is drawn as
  ## marks rather than a curve; a zero-inflated one is the continuous curve
  ## plus the one emphasized atom.
  if var_type == "d":
    plt.plot(ev, vals, marker="o", linestyle="None", markersize=6)
  else:
    plt.plot(ev, vals, linestyle="-", linewidth=2)
  if zero is not None:
    plt.plot(0, zero[0], "o", markersize=8, color="C0")

  if xlim is not None:
    plt.xlim(xlim)
  else:
    lo, hi = _support(margin)
    if np.isfinite(lo) and np.isfinite(hi):
      plt.xlim(lo, hi)
    else:
      pad = 0.05 * (ev.max() - ev.min())
      plt.xlim(ev.min() - pad, ev.max() + pad)

  if ylim is not None:
    plt.ylim(ylim)
  elif kind == "cdf":
    plt.ylim(0, 1.05)
  else:
    top = vals.max() if vals.size else 1.0
    if zero is not None:
      top = max(top, zero[0])
    plt.ylim(0, 1.1 * top)

  plt.xlabel("x")
  plt.ylabel("density" if kind == "density" else "probability")
  plt.grid(True, alpha=0.3)
  plt.show()
