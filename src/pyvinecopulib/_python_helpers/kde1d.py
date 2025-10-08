from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np

KDE1D_PLOT_DOC = """
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
    **kwargs
        Additional keyword arguments passed to matplotlib plotting functions.

    Returns
    -------
    Nothing, the function generates a plot and shows it using matplotlib.

    Usage
    -----
    .. code-block:: python

        import pyvinecopulib as pv
        import numpy as np
        import matplotlib.pyplot as plt

        # Continuous data
        np.random.seed(123)
        x = np.random.beta(0.5, 2.0, 100)
        kde = pv.Kde1d()
        kde.fit(x)

        plt.figure(figsize=(10, 6))
        kde.plot()

        # Discrete data
        x_discrete = np.random.poisson(3, 100)
        kde_discrete = pv.Kde1d(type="discrete")
        kde_discrete.fit(x_discrete)
        kde_discrete.plot()

        # Zero-inflated data
        x_zi = np.random.exponential(2, 100)
        x_zi[np.random.choice(100, 30, replace=False)] = 0
        kde_zi = pv.Kde1d(xmin=0, type="zero-inflated")
        kde_zi.fit(x_zi)
        kde_zi.plot()
"""


def make_plotting_grid(kde: Any, grid_size: int = 200) -> np.ndarray:
  """Create appropriate plotting grid based on kde type and data."""

  if kde.type == "discrete":
    # For discrete data, use integer grid points
    grid_points = kde.grid_points
    return np.arange(
      int(np.floor(grid_points.min())),
      int(np.ceil(grid_points.max())) + 1,
      dtype=float,
    )
  else:
    # For continuous data, create smooth grid
    grid_points = kde.grid_points
    ev = np.linspace(grid_points.min(), grid_points.max(), grid_size)

    # Adjust boundaries if specified
    try:
      if not np.isnan(kde.xmin):
        ev[0] = kde.xmin
      if not np.isnan(kde.xmax):
        ev[-1] = kde.xmax
    except (ValueError, TypeError):
      # Handle case where xmin/xmax might be arrays or not available
      pass

    # For zero-inflated, exclude zero from the main grid
    if kde.type == "zero-inflated":
      ev = ev[ev != 0]

    return np.asarray(ev)


def kde1d_plot(
  kde: Any,
  xlim: Optional[tuple[float, float]] = None,
  ylim: Optional[tuple[float, float]] = None,
  grid_size: int = 200,
  show_zero_mass: bool = True,
) -> None:
  """{}""".format(KDE1D_PLOT_DOC)

  # Check if kde is fitted
  if kde.grid_points.size == 0:
    raise ValueError("Kde1d object must be fitted before plotting")

  # Create plotting grid
  ev = make_plotting_grid(kde, grid_size)

  # Evaluate density
  vals = kde.pdf(ev)

  # Create the main plot based on type
  if kde.type == "discrete":
    plt.plot(ev, vals, marker="o", linestyle="None", markersize=6)
  else:
    plt.plot(ev, vals, linestyle="-", linewidth=2)

  # Handle zero-inflated case
  if kde.type == "zero-inflated" and show_zero_mass:
    zero_density = kde.pdf(np.array([0]))
    plt.plot(0, zero_density[0], "o", markersize=8, color="C0")

  # Set axis limits
  if xlim is not None:
    plt.xlim(xlim)
  else:
    # Auto-set x limits with some padding
    try:
      # Safely check if xmin and xmax are valid scalar values
      xmin_val = getattr(kde, "xmin", np.nan)
      xmax_val = getattr(kde, "xmax", np.nan)

      # Convert to scalar if needed and check if they're finite
      if hasattr(xmin_val, "size") and hasattr(xmin_val, "item"):
        xmin_val = xmin_val.item() if xmin_val.size > 0 else np.nan
      if hasattr(xmax_val, "size") and hasattr(xmax_val, "item"):
        xmax_val = xmax_val.item() if xmax_val.size > 0 else np.nan

      if (
        np.isscalar(xmin_val)
        and np.isscalar(xmax_val)
        and np.isfinite(xmin_val)
        and np.isfinite(xmax_val)
      ):
        plt.xlim(xmin_val, xmax_val)
      else:
        x_range = ev.max() - ev.min()
        plt.xlim(ev.min() - 0.05 * x_range, ev.max() + 0.05 * x_range)
    except (ValueError, TypeError, AttributeError):
      # Handle case where xmin/xmax might be arrays or other types
      x_range = ev.max() - ev.min()
      plt.xlim(ev.min() - 0.05 * x_range, ev.max() + 0.05 * x_range)

  if ylim is not None:
    plt.ylim(ylim)
  else:
    # Auto-set y limits
    max_val = vals.max()
    if kde.type == "zero-inflated" and show_zero_mass:
      zero_density = kde.pdf(np.array([0]))
      max_val = max(max_val, zero_density[0])
    plt.ylim(0, 1.1 * max_val)

  # Set labels
  plt.xlabel("x")
  plt.ylabel("density")

  # Add grid for better readability
  plt.grid(True, alpha=0.3)

  # Show the plot
  plt.show()
