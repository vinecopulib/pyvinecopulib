"""PyTorch ``BiCop`` module — port of the kernel/TLL bivariate copula.

Consumes the fitted grid from :class:`pyvinecopulib.Bicop` (with the ``tll``
family). The evaluation chain (``pdf`` / ``cdf`` / ``hfunc`` / ``hinv`` /
``sample``) lives entirely in PyTorch, so it can be moved to GPU and
combined with autograd-aware downstream code.

Fitting is intentionally not provided here: the underlying TLL kernel
estimator is computed in the C++ library, and the user is expected to
construct an instance via :meth:`TorchBicop.from_bicop` after fitting with
``pv.Bicop(family=tll)``.
"""

from __future__ import annotations

from typing import Optional

import torch
from torch import Tensor

from ..pyvinecopulib_ext import Bicop, tll as _TLL_FAMILY
from ._interp import InterpolationGrid2D, _TRIM_LO, _TRIM_HI
from ._util import make_normal_grid, solve_ITP

_LOG_FLOOR: float = -13.815510557964274  # log(1e-6); same as torchvinecopulib

_VALID_ROTATIONS = (0, 90, 180, 270)


class TorchBicop(torch.nn.Module):
  """PyTorch bivariate copula on an interpolation grid.

  This is the PyTorch analogue of ``KernelBicop`` (and the TLL family): it
  stores the fitted density on an ``m × m`` grid in ``[0, 1]^2`` and exposes
  the standard copula evaluation API.

  Parameters
  ----------
  grid_points:
    1-D tensor of strictly increasing grid points on ``[0, 1]`` (the same
    grid is used along both axes). Endpoints will be clipped to exactly
    ``0`` and ``1`` to avoid extrapolation.
  values:
    Square ``(m, m)`` tensor of density values on the tensor-product grid.
  rotation:
    Counter-clockwise copula rotation in degrees; one of
    ``{0, 90, 180, 270}``.
  cache_integrals:
    If ``True``, precompute ``cdf`` / ``hfunc1`` / ``hfunc2`` at every grid
    node (three extra ``m × m`` buffers, ≈21 KiB at ``m=30``). ``cdf``/
    ``hfunc`` calls then use one bilinear lookup on the cached grid —
    much faster for large batches, with a small interpolation gap between
    grid nodes vs. the on-the-fly trapezoidal integration. Default
    ``False`` matches the C++ implementation exactly.
  norm_times:
    Number of margin-normalization rounds; passed through to
    :class:`InterpolationGrid2D`. The C++ default is 3; pass 0 to skip
    when the grid already integrates to uniform margins.
  device, dtype:
    Standard module placement / precision controls. ``dtype`` defaults to
    ``torch.float64`` for parity with the C++ evaluation.
  """

  rotation: int
  is_indep: bool

  def __init__(
    self,
    grid_points: Optional[Tensor] = None,
    values: Optional[Tensor] = None,
    rotation: int = 0,
    cache_integrals: bool = False,
    norm_times: int = 3,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> None:
    super().__init__()
    if rotation not in _VALID_ROTATIONS:
      raise ValueError(
        f"rotation must be one of {_VALID_ROTATIONS}; got {rotation}"
      )
    self.rotation = int(rotation)

    if grid_points is None and values is None:
      # Independent-copula short-circuit: pdf=1 everywhere, cdf=u1*u2, etc.
      self.is_indep = True
      m = 2
      grid_points = torch.tensor([0.0, 1.0], dtype=dtype, device=device)
      values = torch.ones(m, m, dtype=dtype, device=device)
      norm_times_eff = 0
    else:
      if grid_points is None or values is None:
        raise ValueError(
          "Provide either both grid_points and values, or neither (for the "
          "independence copula)."
        )
      self.is_indep = False
      grid_points = torch.as_tensor(grid_points, dtype=dtype, device=device)
      values = torch.as_tensor(values, dtype=dtype, device=device)
      norm_times_eff = norm_times

    self.interp_grid = InterpolationGrid2D(
      grid_points=grid_points,
      values=values,
      norm_times=norm_times_eff,
    )

    self._cache_integrals = bool(cache_integrals)
    if self._cache_integrals and not self.is_indep:
      cdf_vals, h1_vals, h2_vals = self.interp_grid.build_caches()
      self.register_buffer("_cdf_cache", cdf_vals)
      self.register_buffer("_hfunc1_cache", h1_vals)
      self.register_buffer("_hfunc2_cache", h2_vals)
    else:
      self._cdf_cache = None
      self._hfunc1_cache = None
      self._hfunc2_cache = None

  # --------------------------------------------------------------------- #
  # Constructors                                                           #
  # --------------------------------------------------------------------- #

  @classmethod
  def from_bicop(
    cls,
    cop: Bicop,
    cache_integrals: bool = False,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchBicop":
    """Build a ``TorchBicop`` from a fitted :class:`pyvinecopulib.Bicop`.

    The source ``cop`` is expected to be a kernel-based family (``tll``).
    The interpolation grid is reconstructed from
    ``cop.parameters`` (the fitted ``m × m`` density grid) and the
    canonical normal-scale grid points used by the C++ library.
    """
    if cop.family != _TLL_FAMILY:
      raise ValueError(
        f"from_bicop expects a kernel-family bicop (e.g. pv.tll); "
        f"got family={cop.family!r}"
      )

    values_np = cop.parameters
    if values_np.ndim != 2 or values_np.shape[0] != values_np.shape[1]:
      raise ValueError(
        "cop.parameters must be a square 2D array for a kernel bicop; "
        f"got shape {values_np.shape}"
      )
    m = values_np.shape[0]
    grid_points = make_normal_grid(m, dtype=dtype).to(device=device)
    values = torch.as_tensor(values_np, dtype=dtype, device=device)
    # The grid stored on cop is already normalized; skip renormalization
    # to avoid drifting away from the C++ values.
    return cls(
      grid_points=grid_points,
      values=values,
      rotation=int(cop.rotation),
      cache_integrals=cache_integrals,
      norm_times=0,
      device=device,
      dtype=dtype,
    )

  # --------------------------------------------------------------------- #
  # Rotation plumbing                                                      #
  # --------------------------------------------------------------------- #

  def _rotate_input(self, u: Tensor) -> Tensor:
    """Counter-clockwise rotation of ``u`` to the 0deg frame.

    Mirrors ``Bicop::rotate_data`` (class.ipp). The underlying interpolation
    grid is always queried at the rotated coordinates.
    """
    r = self.rotation
    if r == 0:
      return u
    if r == 90:
      # swap columns then 1 - new col 1
      out = u[:, [1, 0]].clone()
      out[:, 1] = 1.0 - out[:, 1]
      return out
    if r == 180:
      return 1.0 - u
    # r == 270
    out = u[:, [1, 0]].clone()
    out[:, 0] = 1.0 - out[:, 0]
    return out

  # --------------------------------------------------------------------- #
  # Densities                                                              #
  # --------------------------------------------------------------------- #

  def _prep(self, u: Tensor) -> Tensor:
    u = torch.as_tensor(
      u,
      dtype=self.interp_grid.values.dtype,
      device=self.interp_grid.values.device,
    )
    if u.ndim != 2 or u.shape[1] != 2:
      raise ValueError(f"u must have shape (n, 2); got {tuple(u.shape)}")
    # Trim to (1e-10, 1 - 1e-10) before rotation, mirroring Bicop::prep_for_abstract.
    return u.clamp(_TRIM_LO, _TRIM_HI)

  def pdf(self, u: Tensor) -> Tensor:
    """Bivariate copula density at ``u``."""
    u = self._prep(u)
    if self.is_indep:
      return torch.ones(u.shape[0], dtype=u.dtype, device=u.device)
    u_rot = self._rotate_input(u)
    return self.interp_grid.interpolate(u_rot).clamp_min(1e-20)

  def log_pdf(self, u: Tensor) -> Tensor:
    """``log pdf`` with safe handling of ``-inf`` / ``nan``."""
    return self.pdf(u).log().nan_to_num(neginf=_LOG_FLOOR, posinf=0.0)

  def cdf(self, u: Tensor) -> Tensor:
    """Bivariate copula CDF at ``u``."""
    u = self._prep(u)
    if self.is_indep:
      return (u[:, 0] * u[:, 1]).clamp(_TRIM_LO, _TRIM_HI)
    u_rot = self._rotate_input(u)
    if self._cdf_cache is not None:
      p = self.interp_grid.interp_at(self._cdf_cache, u_rot)
    else:
      p = self.interp_grid.integrate_2d(u_rot)
    r = self.rotation
    if r == 90:
      return u[:, 1] - p
    if r == 180:
      return p - 1.0 + u[:, 0] + u[:, 1]
    if r == 270:
      return u[:, 0] - p
    return p

  # --------------------------------------------------------------------- #
  # h-functions                                                            #
  # --------------------------------------------------------------------- #

  def _hfunc_raw(self, u_rot: Tensor, cond_var: int) -> Tensor:
    cache = self._hfunc1_cache if cond_var == 1 else self._hfunc2_cache
    if cache is not None:
      return self.interp_grid.interp_at(cache, u_rot).clamp(_TRIM_LO, _TRIM_HI)
    return self.interp_grid.integrate_1d(u_rot, cond_var=cond_var)

  def hfunc1(self, u: Tensor) -> Tensor:
    """``H1(u1, u2) = P(U2 <= u2 | U1 = u1)``."""
    u = self._prep(u)
    if self.is_indep:
      return u[:, 1].clamp(_TRIM_LO, _TRIM_HI)
    u_rot = self._rotate_input(u)
    r = self.rotation
    if r == 0:
      h = self._hfunc_raw(u_rot, 1)
    elif r == 90:
      h = self._hfunc_raw(u_rot, 2)
    elif r == 180:
      h = 1.0 - self._hfunc_raw(u_rot, 1)
    else:  # r == 270
      h = 1.0 - self._hfunc_raw(u_rot, 2)
    return h.clamp(0.0, 1.0)

  def hfunc2(self, u: Tensor) -> Tensor:
    """``H2(u1, u2) = P(U1 <= u1 | U2 = u2)``."""
    u = self._prep(u)
    if self.is_indep:
      return u[:, 0].clamp(_TRIM_LO, _TRIM_HI)
    u_rot = self._rotate_input(u)
    r = self.rotation
    if r == 0:
      h = self._hfunc_raw(u_rot, 2)
    elif r == 90:
      h = 1.0 - self._hfunc_raw(u_rot, 1)
    elif r == 180:
      h = 1.0 - self._hfunc_raw(u_rot, 2)
    else:  # r == 270
      h = self._hfunc_raw(u_rot, 1)
    return h.clamp(0.0, 1.0)

  # --------------------------------------------------------------------- #
  # Inverse h-functions (numerical via ITP).                               #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def _hinv1_raw(self, u_rot: Tensor) -> Tensor:
    """At rotation 0: solve ``hfunc1((u_rot[:, 0], x)) = u_rot[:, 1]`` for ``x``."""
    a = u_rot[:, 0:1]
    p = u_rot[:, 1:2]

    def fun(x: Tensor) -> Tensor:
      u_eval = torch.cat([a, x], dim=-1)
      return self._hfunc_raw(u_eval, 1).unsqueeze(-1) - p

    x = solve_ITP(fun, x_a=torch.zeros_like(p), x_b=torch.ones_like(p))
    return x.squeeze(-1).clamp(0.0, 1.0)

  @torch.no_grad()
  def _hinv2_raw(self, u_rot: Tensor) -> Tensor:
    """At rotation 0: solve ``hfunc2((x, u_rot[:, 1])) = u_rot[:, 0]`` for ``x``."""
    p = u_rot[:, 0:1]
    b = u_rot[:, 1:2]

    def fun(x: Tensor) -> Tensor:
      u_eval = torch.cat([x, b], dim=-1)
      return self._hfunc_raw(u_eval, 2).unsqueeze(-1) - p

    x = solve_ITP(fun, x_a=torch.zeros_like(p), x_b=torch.ones_like(p))
    return x.squeeze(-1).clamp(0.0, 1.0)

  @torch.no_grad()
  def hinv1(self, u: Tensor) -> Tensor:
    """Inverse of ``hfunc1`` w.r.t. the second argument.

    Given ``u = [u1, p]`` with shape ``(n, 2)``, returns ``u2`` such that
    ``hfunc1((u1, u2)) = p``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 1].clamp(_TRIM_LO, _TRIM_HI)
    u_rot = self._rotate_input(u)
    r = self.rotation
    if r == 0:
      return self._hinv1_raw(u_rot)
    if r == 90:
      return self._hinv2_raw(u_rot)
    if r == 180:
      return (1.0 - self._hinv1_raw(u_rot)).clamp(0.0, 1.0)
    # r == 270
    return (1.0 - self._hinv2_raw(u_rot)).clamp(0.0, 1.0)

  @torch.no_grad()
  def hinv2(self, u: Tensor) -> Tensor:
    """Inverse of ``hfunc2`` w.r.t. the first argument.

    Given ``u = [p, u2]`` with shape ``(n, 2)``, returns ``u1`` such that
    ``hfunc2((u1, u2)) = p``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 0].clamp(_TRIM_LO, _TRIM_HI)
    u_rot = self._rotate_input(u)
    r = self.rotation
    if r == 0:
      return self._hinv2_raw(u_rot)
    if r == 90:
      return (1.0 - self._hinv1_raw(u_rot)).clamp(0.0, 1.0)
    if r == 180:
      return (1.0 - self._hinv2_raw(u_rot)).clamp(0.0, 1.0)
    # r == 270
    return self._hinv1_raw(u_rot)

  # --------------------------------------------------------------------- #
  # Sampling                                                               #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def sample(
    self,
    num_sample: int = 100,
    seed: Optional[int] = 42,
    is_sobol: bool = False,
  ) -> Tensor:
    """Draw ``num_sample`` samples via inverse Rosenblatt (``hinv1``)."""
    device = self.interp_grid.values.device
    dtype = self.interp_grid.values.dtype
    if is_sobol:
      kwargs = {"dimension": 2, "scramble": True}
      if seed is not None:
        kwargs["seed"] = seed
      u = (
        torch.quasirandom.SobolEngine(**kwargs)
        .draw(n=num_sample, dtype=dtype)
        .to(device=device)
      )
    else:
      if seed is not None:
        torch.manual_seed(seed)
      u = torch.rand(num_sample, 2, dtype=dtype, device=device)
    if self.is_indep:
      return u
    u2 = self.hinv1(u).unsqueeze(-1)
    return torch.cat([u[:, 0:1], u2], dim=-1)
