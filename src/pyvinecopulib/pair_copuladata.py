from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import ArrayLike

from ._python_helpers.stats import norm_cdf, norm_pdf
from .pyvinecopulib_ext import Bicop, BicopFamily, FitControlsBicop, wdm


def pairs_copula_data(
  data: ArrayLike,
  main: str = "",
  cols: Optional[list[str]] = None,
  grid_size: int = 50,
  bins: int = 20,
  scatter_size: float = 6.0,
) -> tuple[Figure, Axes]:
  """
  Pair plot for copula data U in (0,1)^d using pure Matplotlib.
  - Lower: bivariate copula density contours (fitted with pyvinecopulib), drawn in z-space.
  - Upper: scatter with Kendall's tau annotation (copula space).
  - Diagonal: histograms (copula space).

  Parameters
  ----------
  data : (n,d) array-like
      Copula data with entries strictly in (0,1).
  main : str
      Figure title.
  grid_size : int
      Resolution of the contour grid per dimension (lower panels). Must be positive.
  bins : int
      Number of histogram bins (diagonal). Must be positive.
  scatter_size : float
      Marker size for upper-panel scatter. Must be positive.

  Returns
  -------
  fig, axes : matplotlib Figure and Axes array of shape (d, d)
  """
  # Input validation
  if data is None:
    raise ValueError("`data` cannot be None.")

  try:
    U = np.asarray(data, dtype=float)
  except (ValueError, TypeError) as e:
    raise ValueError(f"Could not convert `data` to numeric array: {e}")

  if U.ndim != 2:
    raise ValueError("`data` must be a 2D array-like (n,d).")

  if U.size == 0:
    raise ValueError("`data` cannot be empty.")

  if not (np.all(U > 0.0) and np.all(U < 1.0)):
    raise ValueError("All values must lie strictly in (0,1).")

  # Parameter validation
  if not isinstance(grid_size, int) or grid_size <= 0:
    raise ValueError("`grid_size` must be a positive integer.")

  if not isinstance(bins, int) or bins <= 0:
    raise ValueError("`bins` must be a positive integer.")

  if not isinstance(scatter_size, (int, float)) or scatter_size <= 0:
    raise ValueError("`scatter_size` must be a positive number.")

  n, d = U.shape

  if n < 2:
    raise ValueError(f"Need at least 2 observations, got {n}.")

  if cols is not None and (
    not isinstance(cols, list)
    or len(cols) != d
    or not all(isinstance(c, str) for c in cols)
  ):
    raise ValueError(f"`cols` must be a list of {d} strings or None.")
  if cols is None:
    cols = [f"var{i + 1}" for i in range(d)]

  # Prepare z-grid once
  z = np.linspace(-3.0, 3.0, grid_size)
  Z1, Z2 = np.meshgrid(z, z, indexing="xy")
  U1 = norm_cdf(Z1)
  U2 = norm_cdf(Z2)
  grid_u = np.column_stack([U1.ravel(), U2.ravel()])

  fig, axes = plt.subplots(
    d, d, figsize=(2.8 * d, 2.8 * d), sharex=False, sharey=False
  )
  if d == 1:
    axes = np.array([[axes]])

  # Helpers for consistent styling
  def set_zspace(ax: Axes) -> None:
    ax.set_xlim(-3.0, 3.0)
    ax.set_ylim(-3.0, 3.0)
    ax.set_xticks([-3.0, 0.0, 3.0])
    ax.set_yticks([-3.0, 0.0, 3.0])

  def set_copula_space(ax: Axes) -> None:
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([0.0, 0.5, 1.0])

  # Main loop
  for i in range(d):
    for j in range(d):
      ax = axes[i, j]
      if i == j:
        # Diagonal: histogram in copula space
        x = U[:, j].flatten()
        ax.hist(x, bins=bins, range=(0.0, 1.0), density=True, edgecolor="white")
        ax.hlines(1.0, 0.0, 1.0, linestyles="dashed", linewidth=0.8)
        set_copula_space(ax)
        if j == 0:
          ax.set_ylabel(cols[i])
        if i == d - 1:
          ax.set_xlabel(cols[j])

      elif i < j:
        # Upper: scatter with Kendall's tau (copula space)
        x = U[:, j].flatten()
        y = U[:, i].flatten()
        ax.scatter(x, y, s=scatter_size, alpha=0.6)

        try:
          tau = wdm(x, y, "kendall")
          tau_text = f"τ = {tau:.2f}"
          fontsize = 10 + 8 * abs(tau)
        except Exception:
          tau_text = "τ = N/A"
          fontsize = 10

        ax.text(
          0.5,
          0.5,
          tau_text,
          transform=ax.transAxes,
          ha="center",
          va="center",
          fontsize=fontsize,
          weight="bold",
        )
        set_copula_space(ax)
        if i == d - 1:
          ax.set_xlabel(cols[j])
        if j == 0:
          ax.set_ylabel(cols[i])

      else:
        # Lower: bicop contours in z-space
        x = U[:, j].flatten()
        y = U[:, i].flatten()
        uv = np.column_stack([x, y])

        try:
          controls = FitControlsBicop(family_set=[BicopFamily.tll])
          cop = Bicop.from_data(uv, controls=controls)

          # Temporarily enforce continuous var-types for pdf evaluation
          vt = cop.var_types
          cop.var_types = ["c", "c"]
          cvals = cop.pdf(grid_u).reshape(grid_size, grid_size)
          cop.var_types = vt

          dens = (
            cvals * norm_pdf(Z1) * norm_pdf(Z2)
          )  # z-space density via Jacobian

          # Safeguard for flat fields (to avoid contour errors)
          if np.allclose(dens.min(), dens.max()):
            dens = dens.copy()
            dens.flat[0] *= 1.000001

          # Contours in z-space
          ax.contour(
            Z1,
            Z2,
            dens,
            levels=[0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5],
            linewidths=0.8,
          )
        except Exception as e:
          # If copula fitting or plotting fails, show a simple message
          ax.text(
            0.5,
            0.5,
            f"Fit failed:\n{type(e).__name__}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
            weight="bold",
          )

        set_zspace(ax)

        if i == d - 1:
          ax.set_xlabel(cols[j])
        if j == 0:
          ax.set_ylabel(cols[i])

      # Keep all ticks/labels visible (no axis sharing side-effects)
      ax.tick_params(labelbottom=True, labelleft=True)

  if main:
    fig.suptitle(main, y=1.02)

  plt.tight_layout(rect=(0, 0, 1, 0.97))
  return fig, axes
