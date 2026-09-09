"""What a margin plot draws, for every variable type it can describe.

``margin_plot`` is bound to ``Kde1d.plot`` by the binding, which looks it up by
module path, and is what ``MarginBase.plot`` forwards to -- so a compiled
``Kde1d``, a ``TorchKde1d``, a ``SciPyMargin`` and a hand-written subclass all
draw the same picture. It reads only ``var_type``, ``support``, ``is_fitted``
and the requested one of ``pdf`` / ``cdf``, each of which ``MarginBase``
supplies with a continuous-correct default, and prefers an optional
``grid_points`` for the x-range where the margin has one.
"""

from typing import Any, Callable, Optional

import matplotlib.pyplot as plt
import numpy as np

MARGIN_PLOT_DOC = """
    Generates a plot for the Kde1d object.

    This method creates a line plot for continuous data, a point plot for discrete data,
    and handles zero-inflated data with special point marking at zero.

    Parameters
    ----------
    xlim : tuple (default=None)
        The limits for the x axis. Automatically set if None.
    ylim : tuple (default=None)
        The limits for the y axis. Automatically set if None.
    grid_size : int (default=200)
        The number of grid points to use for continuous data.
    show_zero_mass : bool (default=True)
        Whether to show the point mass at zero for zero-inflated data.
    kind : str (default="density")
        What to draw: `"density"` or `"cdf"`.

    Returns
    -------
    Nothing, the function generates a plot and shows it using matplotlib.

    Examples
    --------
    >>> import pyvinecopulib as pv
    >>> import numpy as np
    >>> # Continuous data
    >>> np.random.seed(123)
    >>> x = np.random.beta(0.5, 2.0, 100)
    >>> kde = pv.core.Kde1d()
    >>> kde.fit(x)
    >>> kde.plot()
    >>> kde.plot(kind="cdf")
    >>> # Discrete data
    >>> x_discrete = np.random.poisson(3, 100)
    >>> kde_discrete = pv.core.Kde1d(type="discrete")
    >>> kde_discrete.fit(x_discrete)
    >>> kde_discrete.plot()
    >>> # Zero-inflated data
    >>> x_zi = np.random.exponential(2, 100)
    >>> x_zi[np.random.choice(100, 30, replace=False)] = 0
    >>> kde_zi = pv.core.Kde1d(xmin=0, type="zero-inflated")
    >>> kde_zi.fit(x_zi)
    >>> kde_zi.plot()
"""

#: The two quantiles a margin with an unbounded support is drawn between,
#: reached only where there is no fitted grid to read the range off.
_TAIL = (1e-3, 1 - 1e-3)


def _to_numpy(a: Any) -> np.ndarray:  # noqa: ANN401
  """Bring one array of any namespace onto NumPy, host memory included."""
  detach = getattr(a, "detach", None)
  if detach is not None:
    a = detach()
  to_cpu = getattr(a, "cpu", None)
  if to_cpu is not None:
    a = to_cpu()
  return np.asarray(a, dtype=float)


def _support(margin: Any) -> tuple[float, float]:  # noqa: ANN401
  """The margin's declared support, as two floats."""
  lo, hi = getattr(margin, "support", (-np.inf, np.inf))
  return (float(lo), float(hi))


def _grid_points(margin: Any) -> Optional[np.ndarray]:  # noqa: ANN401
  """The fitted evaluation grid, where the margin publishes one."""
  gp = getattr(margin, "grid_points", None)
  if gp is None:
    return None
  points = _to_numpy(gp)
  return points if points.size else None


def _draw_range(
  margin: Any,  # noqa: ANN401
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
  tails = _to_numpy(icdf(p if place is None else place(p)))
  return (
    lo if np.isfinite(lo) else float(tails[0]),
    hi if np.isfinite(hi) else float(tails[1]),
  )


def make_plotting_grid(
  margin: Any,  # noqa: ANN401
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
  # A margin -- a `MarginLike`, or the `Kde1d` the binding hands here, which
  # it looks this function up by name to reach. Importing the extension to
  # name that type would invert the layering.
  margin: Any,  # noqa: ANN401
  xlim: Optional[tuple[float, float]] = None,
  ylim: Optional[tuple[float, float]] = None,
  grid_size: int = 200,
  show_zero_mass: bool = True,
  *,
  kind: str = "density",
  x: Optional[Any] = None,  # noqa: ANN401
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
    row = np.asarray(x, dtype=float)
    if row.ndim == 1:
      row = row.reshape(1, -1)
    if row.ndim != 2 or row.shape[0] != 1:
      raise ValueError(
        "x must be a single covariate row, shape (p,) or (1, p): a 2-d plot "
        f"shows the margin at one covariate value; got {tuple(row.shape)}"
      )

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
    return _to_numpy(fn(ya, **kwargs))

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
