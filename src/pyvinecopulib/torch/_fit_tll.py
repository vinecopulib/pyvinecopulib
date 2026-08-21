"""Pure-PyTorch port of the C++ ``TllBicop::fit`` for the ``constant`` method.

Mirrors the algorithm in ``lib/vinecopulib/.../bicop/implementation/tll.ipp``
end-to-end so :meth:`TorchBicop.from_data` produces an ``(m, m)`` density
grid that matches the C++ ``pv.Bicop.from_data`` output to machine
precision after :meth:`InterpolationGrid2D.normalize_margins(25)`.

Three pieces:

* :func:`_to_pseudo_obs_continuous` — empirical CDF on continuous data
  (rank/(n+1); no tie handling, mirrors C++ ``wdm`` ranks for
  jitter-free input).
* :func:`_ace` — alternating conditional expectations for the maximal-
  correlation coefficient. Outer/inner convergence loop matches the C++
  tolerances (``2e-15`` / ``1e-4``) and uses an FFT-free moving-average
  window smoother (``F.conv1d``) for portability.
* :func:`fit_tll_constant` — the user-facing entry point. Wraps the
  bandwidth selection (``select_bandwidth_constant``) and the local-
  constant kernel density evaluation (``fit_local_likelihood_constant``)
  into a single call returning the unnormalized ``(m, m)`` values grid.

Only the ``constant`` method is supported here; the ``linear`` and
``quadratic`` variants need more elaborate per-grid-point bookkeeping
(see C++ ``TllBicop::fit_local_likelihood`` and ``calculate_infl``).
"""

from __future__ import annotations

import math

import torch
import torch.nn.functional as F
from torch import Tensor

from ._interp import InterpolationGrid2D

_SQRT_2 = math.sqrt(2.0)
_SQRT_2PI_INV = 1.0 / math.sqrt(2.0 * math.pi)


# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #


def _qnorm(p: Tensor) -> Tensor:
  """Inverse standard normal CDF; matches C++ ``tools_stats::qnorm``."""
  return torch.special.ndtri(p)


def _dnorm(z: Tensor) -> Tensor:
  """Standard normal PDF."""
  return torch.exp(-0.5 * z * z) * _SQRT_2PI_INV


def _to_pseudo_obs_continuous(x: Tensor) -> Tensor:
  """Empirical CDF on continuous data: ``rank/(n+1)``.

  Mirrors ``tools_stats::to_pseudo_obs`` with no ties (continuous input,
  no jittering needed). C++ uses ``wdm`` ranks; for unique data those are
  the same as :func:`torch.argsort` ranks.
  """
  n = x.shape[0]
  ranks = x.argsort(dim=0).argsort(dim=0)
  return (ranks + 1).to(x.dtype) / (n + 1)


def _win_smoother(x: Tensor, wl: int) -> Tensor:
  """Centered moving average of half-window ``wl``, with edge clamps.

  Mirrors ``tools_stats::win`` semantically: pad ``x`` with ``wl`` zeros
  on each side, convolve with a uniform kernel of length ``2*wl + 1``,
  then clamp the leading / trailing ``wl`` entries to ``out[wl]`` and
  ``out[n - wl - 1]`` so the edges are flat. C++ uses FFT for the
  convolution; we use ``F.conv1d`` here for portability — same
  mathematical answer, possibly different rounding order at ~ulp level.
  """
  n = x.shape[0]
  weight = torch.full(
    (1, 1, 2 * wl + 1), 1.0 / (2 * wl + 1), dtype=x.dtype, device=x.device
  )
  x_padded = F.pad(x.view(1, 1, n), (wl, wl))
  out = F.conv1d(x_padded, weight).view(n).clone()
  if wl > 0:
    out[:wl] = out[wl]
    out[-wl:] = out[n - wl - 1]
  return out


def _cef(x: Tensor, ind: Tensor, ranks: Tensor, wl: int) -> Tensor:
  """``cef`` helper: ``win(x[ind], wl)[ranks]``.

  Smooths ``x`` in sorted order, then maps back to the original order.
  Mirrors ``tools_stats::cef``.
  """
  return _win_smoother(x[ind], wl)[ranks]


def _ace(
  data: Tensor,
  *,
  outer_iter_max: int = 100,
  inner_iter_max: int = 10,
  # Iteration tolerances are hardcoded to match ``tools_stats::ace`` in
  # vinecopulib; if upstream changes them, the C++↔torch parity guard
  # (tests/test_torch_bicop.py::test_from_data_matches_cpp) fails loudly.
  outer_abs_tol: float = 2e-15,
  inner_abs_tol: float = 1e-4,
) -> Tensor:
  """Alternating conditional expectations.

  Mirrors ``tools_stats::ace`` for the unweighted bivariate case (the
  only one that the TLL ``constant`` bandwidth path needs).

  Args:
    data: shape ``(n, 2)``.

  Returns:
    ``(n, 2)`` tensor of the ACE-transformed scores ``phi``.
  """
  n = data.shape[0]
  dtype, device = data.dtype, data.device
  wl = int(math.ceil(n / 5.0))

  ind = torch.empty(n, 2, dtype=torch.long, device=device)
  ranks = torch.empty(n, 2, dtype=torch.long, device=device)
  for i in range(2):
    order = data[:, i].argsort(stable=True)
    ind[:, i] = order
    ranks[order, i] = torch.arange(n, device=device)

  phi = ranks.to(dtype).clone()
  phi -= (n - 1) / 2.0 - 1.0
  phi /= math.sqrt(n * (n - 1) / 12.0)

  outer_iter, outer_eps, outer_abs_err = 1, 1.0, 1.0
  while outer_iter <= outer_iter_max and outer_abs_err > outer_abs_tol:
    inner_iter, inner_eps, inner_abs_err = 1, 1.0, 1.0
    while inner_iter <= inner_iter_max and inner_abs_err > inner_abs_tol:
      phi[:, 1] = _cef(phi[:, 0], ind[:, 1], ranks[:, 1], wl)
      phi[:, 1] = phi[:, 1] - phi[:, 1].sum() / n
      phi[:, 1] = phi[:, 1] / ((phi[:, 1] ** 2).sum() / (n - 1)).sqrt()
      prev = inner_eps
      inner_eps = ((phi[:, 1] - phi[:, 0]) ** 2).sum().item() / n
      inner_abs_err = abs(prev - inner_eps)
      inner_iter += 1
    phi[:, 0] = _cef(phi[:, 1], ind[:, 0], ranks[:, 0], wl)
    phi[:, 0] = phi[:, 0] - phi[:, 0].sum() / n
    phi[:, 0] = phi[:, 0] / ((phi[:, 0] ** 2).sum() / (n - 1)).sqrt()
    prev = outer_eps
    outer_eps = ((phi[:, 1] - phi[:, 0]) ** 2).sum().item() / n
    outer_abs_err = abs(prev - outer_eps)
    outer_iter += 1

  return phi


def _pearson_cor(x: Tensor) -> Tensor:
  """Pearson correlation of two columns of ``x: (n, 2)``. Returns a 0-D tensor."""
  x0 = x[:, 0] - x[:, 0].mean()
  x1 = x[:, 1] - x[:, 1].mean()
  return (x0 * x1).sum() / ((x0**2).sum().sqrt() * (x1**2).sum().sqrt())


def _pairwise_mcor(x: Tensor) -> float:
  """Maximal correlation via ACE + Pearson. Returns a Python float."""
  phi = _ace(x)
  return _pearson_cor(phi).item()


def _chol22(B: Tensor) -> Tensor:
  """Cholesky factor of a 2x2 SPD matrix; lower-triangular."""
  rB = torch.zeros_like(B)
  rB[0, 0] = B[0, 0].sqrt()
  rB[1, 0] = B[1, 0] / rB[0, 0]
  rB[1, 1] = (B[1, 1] - rB[1, 0] ** 2).sqrt()
  return rB


def _select_bandwidth_constant(z: Tensor) -> Tensor:
  """Bandwidth matrix for the constant-method local-likelihood KDE.

  Mirrors ``TllBicop::select_bandwidth`` for ``method == "constant"``.
  """
  n = z.shape[0]
  cor = _pearson_cor(z).clamp(-0.95, 0.95).item()
  cov = torch.tensor([[1.0, cor], [cor, 1.0]], dtype=z.dtype, device=z.device)
  mult = n ** (-1.0 / 3.0)
  mcor = _pairwise_mcor(z)
  scale = abs(cor / mcor) ** (0.5 * mcor)
  return mult * cov * scale


def _fit_local_likelihood_constant(
  z: Tensor, z_data: Tensor, B: Tensor
) -> Tensor:
  """Local-constant kernel density estimate at ``z`` from ``z_data``.

  Mirrors ``TllBicop::fit_local_likelihood`` for ``method == "constant"``,
  vectorized over all grid points at once (the per-grid-point ``for`` loop
  in C++ becomes a single broadcast).
  """
  irB = torch.linalg.inv(_chol22(B))
  det_irB = torch.linalg.det(irB)
  z_dec = (irB @ z.T).T  # (m², 2)
  z_data_dec = (irB @ z_data.T).T  # (n, 2)
  zz = z_data_dec.unsqueeze(0) - z_dec.unsqueeze(1)  # (m², n, 2)
  kernels = (
    torch.exp(-0.5 * (zz[..., 0] ** 2 + zz[..., 1] ** 2))
    * (_SQRT_2PI_INV * _SQRT_2PI_INV)
    * det_irB
  )
  return kernels.mean(dim=1)  # (m²,)


# --------------------------------------------------------------------------- #
# Public entry point                                                           #
# --------------------------------------------------------------------------- #


def fit_tll_constant(
  u: Tensor,
  grid_size: int = 30,
  mult: float = 1.0,
  grid_type: str = "normal",
) -> tuple[Tensor, Tensor]:
  """Fit a TLL pair-copula via local-constant kernel density estimation.

  Pure-PyTorch port of ``TllBicop::fit`` for ``method == "constant"``.

  Args:
    u: ``(n, 2)`` tensor of pseudo-observations in ``(0, 1)``.
    grid_size: number of grid points per axis (default 30; matches C++).
    mult: bandwidth multiplier passed through to ``select_bandwidth``;
      the C++ default is 1.
    grid_type: either ``"normal"`` (default, matches C++ — grid_points are
      ``pnorm(linspace(-3.25, 3.25, m))``) or ``"linear"`` (grid_points
      are uniformly spaced on ``[Phi(-3.25), 1-Phi(-3.25)]``, with the
      InterpolationGrid2D constructor forcing the endpoints to exactly
      0 / 1). The linear option matches the z-range of the normal grid so
      bandwidth selection is unaffected; only the *storage* grid changes,
      enabling O(1) cell-finding at eval time.

  Returns:
    A ``(grid_points, values)`` pair. ``values`` is the unnormalized
    ``(m, m)`` density; callers should pass it through
    ``InterpolationGrid2D(grid_points, values, norm_maxiter=25,
    is_linear=(grid_type == "linear"))`` to match the C++
    ``Bicop.parameters`` output to machine precision for ``"normal"``.
  """
  if u.ndim != 2 or u.shape[1] != 2:
    raise ValueError(f"u must have shape (n, 2); got {tuple(u.shape)}")
  if grid_size < 2:
    raise ValueError(f"grid_size must be >= 2; got {grid_size}")
  if mult <= 0:
    raise ValueError(f"mult must be > 0; got {mult}")
  dtype, device = u.dtype, u.device

  # Pseudo-observations + qnorm to z-space.
  psobs = _to_pseudo_obs_continuous(u)
  z_data = _qnorm(psobs)

  # Bandwidth selection.
  B = _select_bandwidth_constant(z_data) * mult

  # Storage and KDE-eval grids come from the centralized factory on
  # InterpolationGrid2D so any future grid type lives in one place.
  # For "normal" the two coincide and match C++; for "linear" the KDE
  # eval grid is uniform-then-clamped so qnorm stays finite while the
  # storage grid is the unclipped linspace(0, 1, m) — the resulting
  # densities at u=0/1 are stored at z=±3.25 (same trick as C++).
  grid_points = InterpolationGrid2D.make_grid_points(
    grid_type, grid_size, dtype=dtype
  ).to(device=device)
  kde_points = InterpolationGrid2D.make_kde_eval_points(
    grid_type, grid_size, dtype=dtype
  ).to(device=device)

  # Expand grid (row-major matches the C++ Eigen::Map<...>().transpose()
  # reshape in TllBicop::fit at the values-assembly step).
  g_i = kde_points.repeat_interleave(grid_size)
  g_j = kde_points.repeat(grid_size)
  grid_2d = torch.stack([g_i, g_j], dim=-1)
  z = _qnorm(grid_2d)

  # Local-constant KDE on the z-scale grid.
  f0 = _fit_local_likelihood_constant(z, z_data, B)

  # Transform z-space density to copula-scale density.
  phi_z = (
    torch.exp(-0.5 * (z[:, 0] ** 2 + z[:, 1] ** 2))
    * _SQRT_2PI_INV
    * _SQRT_2PI_INV
  )
  c = f0 / phi_z
  values = c.reshape(grid_size, grid_size)

  # The canonical TorchBicop builds an InterpolationGrid2D from
  # (grid_points, values) with norm_maxiter=25 — that's what matches C++
  # to machine precision. The grid_points returned here are the
  # un-forced ones used for the fit positions; InterpolationGrid2D's
  # constructor will clamp the endpoints to 0/1 internally.
  return grid_points, values
