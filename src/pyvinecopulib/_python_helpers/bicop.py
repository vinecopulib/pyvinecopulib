from typing import Any, Optional, cast

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d.axis3d import XAxis as XAxis3D, YAxis as YAxis3D

from ..core._covariates import pair_eval
from .stats import expon_cdf, expon_pdf, expon_ppf, norm_cdf, norm_pdf, norm_ppf

BICOP_PLOT_DOC = """
    Generates a plot for the Bicop object.

    This method generates a contour or surface plot of the copula density. It can be used to visualize the copula density with different types of margins.

    Parameters
    ----------
    plot_type : str (default="contour")
        The type of plot to generate. Either `"contour"` or `"surface"`.
    margin_type : str (default="unif")
        The type of margins to use. Either `"unif"`, `"norm"`, or `"exp"`.
    xylim : tuple (default=None)
        The limits for the x and y axes. Automatically set if None.
    grid_size : int (default=None)
        The number of grid points to use. Automatically set if None.

    Returns
    -------
    Nothing, the function generates a plot and shows it using matplotlib.

    Examples
    --------
    >>> import pyvinecopulib as pv
    >>> import numpy as np
    >>> cop = pv.Bicop(
    ...     family=pv.BicopFamily.gaussian, parameters=np.array([[0.5]]),
    ... )
    >>> cop.plot()  # surface plot of copula density
    >>> cop.plot(plot_type="contour", margin_type="norm")
    >>> cop.plot(plot_type="contour", margin_type="unif")
"""


def get_default_xylim(margin_type: str) -> tuple[float, float]:
  if margin_type == "unif":
    return (1e-2, 1 - 1e-2)
  elif margin_type == "norm":
    return (-3, 3)
  elif margin_type == "exp":
    return (0, 6)
  else:
    raise ValueError("Unknown margin type")


def get_default_grid_size(plot_type: str) -> int:
  if plot_type == "contour":
    return 100
  elif plot_type == "surface":
    return 40
  else:
    raise ValueError("Unknown plot type")


def bicop_plot(
  cop: Any,
  plot_type: str = "surface",
  margin_type: str = "unif",
  xylim: Optional[tuple[float, float]] = None,
  grid_size: Optional[int] = None,
  *,
  x: Optional[Any] = None,
  place: Optional[Any] = None,
) -> None:
  """{}""".format(BICOP_PLOT_DOC)

  if plot_type not in ["contour", "surface"]:
    raise ValueError("Unknown type")

  if margin_type not in ["unif", "norm", "exp"]:
    raise ValueError("Unknown margin type")

  if xylim is None:
    xylim = get_default_xylim(margin_type)

  if grid_size is None:
    grid_size = get_default_grid_size(plot_type)

  if margin_type == "unif":
    if plot_type == "contour":
      points = np.linspace(1e-5, 1 - 1e-5, grid_size)
    else:
      points = np.linspace(1, grid_size, grid_size) / (grid_size + 1)

    g = np.meshgrid(points, points)
    points = g[0][0]
    adj = 1
    levels = [0.2, 0.6, 1, 1.5, 2, 3, 5, 10, 20]
    xlabel = "u1"
    ylabel = "u2"
  elif margin_type == "norm":
    points = norm_cdf(np.linspace(xylim[0], xylim[1], grid_size))
    g = np.meshgrid(points, points)
    points = norm_ppf(g[0][0])
    adj = np.outer(norm_pdf(points), norm_pdf(points))
    levels = [0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    xlabel = "z1"
    ylabel = "z2"
  elif margin_type == "exp":
    ll = 1e-2 if plot_type == "contour" else 1e-1
    points = expon_cdf(np.linspace(ll, xylim[1], grid_size))
    g = np.meshgrid(points, points)
    points = expon_ppf(g[0][0])
    adj = np.outer(expon_pdf(points), expon_pdf(points))
    levels = [0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5]
    xlabel = "e1"
    ylabel = "e2"
  else:
    raise ValueError("Unknown margin type")

  ## The grid is manufactured here, so this is the one place a pair copula is
  ## handed an array it did not supply and cannot infer a namespace from.
  ## `place` brings it onto the pair's own namespace, dtype and device; a
  ## compiled `Bicop` passes none and keeps the NumPy grid it always had.
  grid = np.stack(g, axis=-1).reshape(-1, 2)
  u_grid: Any = grid if place is None else place(grid)
  ## A conditional pair copula's density is a different surface for every
  ## covariate value, so a 2-d plot shows one slice: a single row, repeated
  ## across the grid. Placed but not clamped -- covariates are reals.
  x_grid: Optional[Any] = None
  if x is not None:
    row = np.asarray(x, dtype=float)
    if row.ndim == 1:
      row = row.reshape(1, -1)
    if row.ndim != 2 or row.shape[0] != 1:
      raise ValueError(
        "x must be a single covariate row, shape (p,) or (1, p): a 2-d plot "
        f"shows the density at one covariate value; got {tuple(row.shape)}"
      )
    tiled = np.repeat(row, grid.shape[0], axis=0)
    x_grid = tiled if place is None else place(tiled)

  ## evaluate on grid. Use a continuous copy when the pair stores discrete
  ## variable types. A third-party pair without that capability is restored in
  ## a finally block, so plotting can never leave caller-owned model state
  ## changed when density evaluation or plotting raises.
  vt = getattr(cop, "var_types", None)
  if vt is not None:
    # Read the capability from the type: permissive proxy objects such as
    # mocks synthesize arbitrary instance attributes on demand.
    as_continuous = getattr(type(cop), "as_continuous", None)
    if callable(as_continuous):
      eval_cop = as_continuous(cop)
      vals = pair_eval(eval_cop.pdf, u_grid, x_grid)
    else:
      cop.var_types = ["c", "c"]
      try:
        vals = pair_eval(cop.pdf, u_grid, x_grid)
      finally:
        cop.var_types = vt
  else:
    vals = pair_eval(cop.pdf, u_grid, x_grid)
  # Coerce the density to a NumPy array so a torch-tensor return reshapes
  # cleanly -- via the host, since ``np.asarray`` raises on a device tensor.
  detach = getattr(vals, "detach", None)
  if detach is not None:
    vals = detach()
  to_cpu = getattr(vals, "cpu", None)
  if to_cpu is not None:
    vals = to_cpu()
  vals = np.asarray(vals)
  cop = np.reshape(vals, (grid_size, grid_size))

  ## adjust for margins
  dens = cop * adj
  if len(np.unique(dens)) == 1:
    dens[0] = 1.000001 * dens[0]

  if margin_type == "unif":
    zlim = (0, max(3, 1.1 * max(dens.flatten())))
  elif margin_type == "norm":
    zlim = (0, max(0.4, 1.1 * max(dens.flatten())))
  elif margin_type == "exp":
    dens = np.minimum(dens, 6)
    zlim = (0, max(1, 1.1 * max(dens.flatten())))

  # Define the colors as in the R code
  colors = [
    "#00007F",
    "blue",
    "#007FFF",
    "cyan",
    "#7FFF7F",
    "yellow",
    "#FF7F00",
    "red",
    "#7F0000",
  ]

  # Create the custom colormap
  jet_colors = LinearSegmentedColormap.from_list("jet_colors", colors, N=100)

  ## plot
  if plot_type == "contour":
    contour = plt.contour(points, points, dens, levels=levels, cmap="gray")
    plt.clabel(contour, inline=True, fontsize=8, fmt="%1.2f")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
  elif plot_type == "surface":
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(elev=30, azim=-110)
    X, Y = np.meshgrid(points, points)
    ax.plot_surface(X, Y, dens, cmap=jet_colors, edgecolor="none", shade=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xylim)
    ax.set_ylim(xylim)
    ax.set_zlim(zlim)
    ax.set_box_aspect([1, 1, 1])
    # ax is Axes3D, but matplotlib's stubs declare ax.xaxis / ax.yaxis as
    # the 2D XAxis / YAxis (which lack `pane`). Cast to the 3D variants.
    # ax.zaxis is already typed correctly as the 3D ZAxis.
    cast(XAxis3D, ax.xaxis).pane.fill = False
    cast(YAxis3D, ax.yaxis).pane.fill = False
    ax.zaxis.pane.fill = False
    ax.grid(False)
    plt.draw()
    plt.show()
  else:
    raise ValueError("Unknown plot type")
