import math
from typing import Literal, Optional, Sequence


def grid_size_2d(n: int, h: Optional[float] = None, c: float = 0.7) -> int:
  """
  Compute the number of grid points per axis for a 2D KDE table
  with bilinear interpolation.

  Parameters
  ----------
  n : int
      Sample size.
  h : float or None, optional
      Bandwidth. If None, uses rule-of-thumb h ~ n^(-1/6).
  c : float, optional
      Safety constantSequence (< 1). Default 0.7.

  Returns
  -------
  m : int
      Grid size per axis.
  """
  if h is None:
    # rule-of-thumb bandwidth for 2D
    h = n ** (-1 / 6)

  # grid spacing Δ = c * h * (n h^2)^(-1/4)
  delta = c * h * (n * h**2) ** (-0.25)

  # number of intervals on [0,1]
  m = math.ceil(1.0 / delta)
  return m


def grid_size_1d(
  n: int, h: Optional[float] = None, c: float = 0.7, length: float = 1.0
) -> int:
  """
  Number of grid points for a 1D KDE lookup table over an interval of length `length`,
  using linear interpolation.

  Rule: Δ = c * h * (n*h)^(-1/4)  (variance-matched),
        or Δ = min(Δ, c*h*h)      (cap by bias-matched spacing).

  Parameters
  ----------
  n : int
      Sample size.
  h : float or None
      Bandwidth. If None, uses 1D rate-optimal h ~ n^(-1/5).
  c : float
      Safety factor in (0,1]; smaller -> finer grid. Default 0.7.
  length : float
      Interval length (domain width). Default 1.0.

  Returns
  -------
  m : int
      Grid size (number of knots) across the interval.
  """
  if h is None:
    h = n ** (-1 / 5)  # rate-optimal scaling in 1D

  delta_var = c * h * (n * h) ** (-0.25)  # match stochastic error
  delta_bias = c * h * h  # match bias (if bias-dominated)
  delta = min(delta_var, delta_bias)

  m = max(2, math.ceil(length / delta))
  return m


def grid_size_1d_many(
  n: int,
  lengths: Sequence[float],
  h: Optional[float] = None,
  mode: Literal["rate", "bandwidth_steps"] = "bandwidth_steps",
  # parameters for "rate" mode
  c: float = 0.7,
  # parameters for "bandwidth_steps" mode
  K_min: int = 8,
  K_max: int = 16,
  grow_with_n: bool = True,
  # global caps
  m_cap: Optional[int] = 2048,
) -> list[int]:
  """
  Compute grid sizes per axis for 1D KDE lookup tables over given axis lengths.

  Parameters
  ----------
  n : int
      Sample size.
  lengths : Sequence[float]
      Interval lengths per feature, e.g. max - min for each column.
  h : float or None
      Bandwidth. If None, uses 1D rate-optimal scaling h ~ n^(-1/5).
      (You can pass your plug-in/Silverman bandwidth if available.)
  mode : {"rate","bandwidth_steps"}
      "rate":   Δ = c * h * (n*h)^(-1/4)   (theoretical, can be large for big n)
      "bandwidth_steps": choose Δ = h / K with K in [K_min, K_max],
                         optionally slowly growing with n.
  c : float
      Safety factor for "rate" mode (smaller -> finer grid).
  K_min, K_max : int
      Min/max steps per bandwidth for "bandwidth_steps" mode.
  grow_with_n : bool
      If True, use K = min(K_max, K_min * n**0.1) (very slow growth in n).
  m_cap : int or None
      Hard cap on grid size per axis (None to disable).

  Returns
  -------
  list[int]
      Grid size per axis.
  """
  if h is None:
    h = n ** (-1 / 5)  # rate-optimal scaling in 1D

  if mode == "rate":
    delta = c * h * (n * h) ** (-0.25)  # match stochastic term
    # (Bias-matched spacing equals this at h ~ n^(-1/5); you could do: delta = min(delta, c*h*h))
    m_list = [max(2, math.ceil(L / delta)) for L in lengths]
  else:
    # steps per bandwidth; grow very slowly with n (n**0.1) but clamp to K_max
    K = K_min * (n**0.1) if grow_with_n else K_min
    K = max(K_min, min(int(round(K)), K_max))
    delta = h / K
    m_list = [max(2, math.ceil(L / delta)) for L in lengths]

  if m_cap is not None:
    m_list = [min(m, m_cap) for m in m_list]
  return m_list
