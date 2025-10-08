import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ._python_helpers.stats import norm_cdf, norm_pdf
from .pyvinecopulib_ext import Bicop, BicopFamily, FitControlsBicop


def pairs_copula_data(
  data: np.ndarray,
  main: str = "",
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
      Resolution of the contour grid per dimension (lower panels).
  bins : int
      Number of histogram bins (diagonal).
  scatter_size : float
      Marker size for upper-panel scatter.

  Returns
  -------
  fig, axes : matplotlib Figure and Axes array of shape (d, d)
  """
  U = np.asarray(data, dtype=float)
  if U.ndim != 2:
    raise ValueError("`data` must be a 2D array-like (n,d).")
  if not (np.all(U > 0.0) and np.all(U < 1.0)):
    raise ValueError("All values must lie strictly in (0,1).")

  n, d = U.shape
  cols = [f"var{i + 1}" for i in range(d)]
  # df = pd.DataFrame(U, columns=cols)

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
        x = U[:, j]
        ax.hist(x, bins=bins, range=(0.0, 1.0), density=True, edgecolor="white")
        ax.hlines(1.0, 0.0, 1.0, linestyles="dashed", linewidth=0.8)
        set_copula_space(ax)
        if j == 0:
          ax.set_ylabel(cols[i])
        if i == d - 1:
          ax.set_xlabel(cols[j])

      elif i < j:
        # Upper: scatter with Kendall's tau (copula space)
        x = U[:, j]
        y = U[:, i]
        ax.scatter(x, y, s=scatter_size, alpha=0.6)
        # tau, _ = kendalltau(x, y)
        tau = 0.0  # TODO
        ax.text(
          0.5,
          0.5,
          f"τ = {tau:.2f}",
          transform=ax.transAxes,
          ha="center",
          va="center",
          fontsize=10 + 8 * abs(tau),
          weight="bold",
        )
        set_copula_space(ax)
        if i == d - 1:
          ax.set_xlabel(cols[j])
        if j == 0:
          ax.set_ylabel(cols[i])

      else:
        # Lower: bicop contours in z-space
        x = U[:, j]
        y = U[:, i]
        uv = np.column_stack([x, y])

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
        set_zspace(ax)

        if i == d - 1:
          ax.set_xlabel(cols[j] + " (z)")
        if j == 0:
          ax.set_ylabel(cols[i] + " (z)")

      # Keep all ticks/labels visible (no axis sharing side-effects)
      ax.tick_params(labelbottom=True, labelleft=True)

  if main:
    fig.suptitle(main, y=1.02)

  plt.tight_layout(rect=(0, 0, 1, 0.97))
  return fig, axes
