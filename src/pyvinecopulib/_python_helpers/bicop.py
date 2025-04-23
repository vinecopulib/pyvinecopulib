from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

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

    Usage
    -----
    .. code-block:: python

        import pyvinecopulib as pv
        import numpy as np
        cop = pv.Bicop(family=pv.BicopFamily.gaussian, parameters=np.array([[0.5]]))
        cop.plot() # surface plot of copula density
        cop.plot(plot_type="contour", margin_type="norm") # contour plot with normal margins
        cop.plot(plot_type="contour", margin_type="unif") # contour plot of copula density
"""


def get_default_xylim(margin_type: str) -> tuple:
  if margin_type == "unif":
    return (1e-2, 1 - 1e-2)
  elif margin_type == "norm":
    return (-3, 3)
  elif margin_type == "exp":
    return (0, 6)
  else:
    raise ValueError("Unknown margin type")


def get_default_grid_size(plot_type: str) -> tuple:
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
  xylim: tuple = None,
  grid_size: int = None,
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

  ## evaluate on grid
  vt = cop.var_types
  cop.var_types = ["c", "c"]
  vals = cop.pdf(np.stack(g, axis=-1).reshape(-1, 2))
  cop.var_types = vt
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
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.grid(False)
    plt.draw()
    plt.show()
  else:
    raise ValueError("Unknown plot type")
