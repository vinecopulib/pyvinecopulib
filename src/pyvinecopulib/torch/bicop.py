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
from ._util import solve_itp

_LOG_FLOOR: float = -13.815510557964274  # log(1e-6); same as torchvinecopulib


class TorchBicop(torch.nn.Module):
  """PyTorch bivariate copula on an interpolation grid.

  This is the PyTorch analogue of ``KernelBicop`` (and the TLL family): it
  stores the fitted density on an ``m × m`` grid in ``[0, 1]^2`` and exposes
  the standard copula evaluation API. Non-zero copula rotations are not
  supported — TLL/kernel pair-copulas in pyvinecopulib always have
  ``rotation=0`` (the family is non-parametric).

  Parameters
  ----------
  grid_points:
    1-D tensor of strictly increasing grid points on ``[0, 1]`` (the same
    grid is used along both axes). Endpoints will be clipped to exactly
    ``0`` and ``1`` to avoid extrapolation.
  values:
    Square ``(m, m)`` tensor of density values on the tensor-product grid.
  cache_integrals:
    If ``True``, precompute ``cdf`` / ``hfunc1`` / ``hfunc2`` / ``hinv1`` /
    ``hinv2`` at every grid node (five extra ``m × m`` buffers, ≈36 KiB
    at ``m=30``). ``cdf`` / ``hfunc`` / ``hinv`` calls then use one
    bilinear lookup on the cached grid — much faster for large batches
    (``hinv`` in particular drops from 35 bisection iterations to a
    single interp), with a small interpolation gap between grid nodes
    (~1e-3 mean / ~1e-2 max relative to the on-the-fly trapezoidal +
    bisection path). Default ``False`` matches the C++ implementation
    exactly.
  norm_times:
    Number of margin-normalization rounds; passed through to
    :class:`InterpolationGrid2D`. The C++ default is 3; pass 0 to skip
    when the grid already integrates to uniform margins.
  device, dtype:
    Standard module placement / precision controls. ``dtype`` defaults to
    ``torch.float64`` for parity with the C++ evaluation.
  """

  is_indep: bool

  def __init__(
    self,
    grid_points: Optional[Tensor] = None,
    values: Optional[Tensor] = None,
    cache_integrals: bool = False,
    norm_times: int = 3,
    is_linear: bool = False,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> None:
    super().__init__()
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
      is_linear=is_linear,
    )

    self._cache_integrals = bool(cache_integrals)
    if self._cache_integrals and not self.is_indep:
      cdf_vals, h1_vals, h2_vals, hinv1_vals, hinv2_vals = (
        self.interp_grid.build_caches()
      )
      self.register_buffer("_cdf_cache", cdf_vals)
      self.register_buffer("_hfunc1_cache", h1_vals)
      self.register_buffer("_hfunc2_cache", h2_vals)
      self.register_buffer("_hinv1_cache", hinv1_vals)
      self.register_buffer("_hinv2_cache", hinv2_vals)
    else:
      self._cdf_cache = None
      self._hfunc1_cache = None
      self._hfunc2_cache = None
      self._hinv1_cache = None
      self._hinv2_cache = None

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
    if int(cop.rotation) != 0:
      raise ValueError(
        "TorchBicop is rotation-less. TLL pair-copulas in pyvinecopulib "
        f"always have rotation=0; got rotation={int(cop.rotation)!r}."
      )

    values_np = cop.parameters
    if values_np.ndim != 2 or values_np.shape[0] != values_np.shape[1]:
      raise ValueError(
        "cop.parameters must be a square 2D array for a kernel bicop; "
        f"got shape {values_np.shape}"
      )
    m = values_np.shape[0]
    grid_points = InterpolationGrid2D.make_grid_points(
      "normal", m, dtype=dtype
    ).to(device=device)
    values = torch.as_tensor(values_np, dtype=dtype, device=device)
    # The grid stored on cop is already normalized; skip renormalization
    # to avoid drifting away from the C++ values.
    return cls(
      grid_points=grid_points,
      values=values,
      cache_integrals=cache_integrals,
      norm_times=0,
      device=device,
      dtype=dtype,
    )

  @classmethod
  def from_data(
    cls,
    u,
    grid_size: int = 30,
    mult: float = 1.0,
    cache_integrals: bool = False,
    grid_type: str = "normal",
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchBicop":
    """Fit a TLL bicop on pseudo-observations and wrap in a ``TorchBicop``.

    Pure-PyTorch port of the C++ ``TllBicop::fit`` for the ``constant``
    method (the C++ default). With ``grid_type="normal"`` the output
    matches ``pv.Bicop.from_data(u,
    controls=FitControlsBicop(family_set=[tll]))`` to machine precision
    (worst case observed: ~1e-12 across (n, ρ) in
    ``{500, 2000} × {0.3, 0.6, 0.9}``).

    Args:
      u: ``(n, 2)`` pseudo-observations; np.ndarray or Tensor.
      grid_size: density grid size per axis (default 30; matches C++).
      mult: bandwidth multiplier (default 1; matches C++).
      grid_type: ``"normal"`` (default, matches C++ — Phi-spaced grid) or
        ``"linear"`` (uniform grid on [0, 1] with O(1) cell-finding at
        eval time). Both build the KDE on the same z-range so bandwidth
        selection is unaffected; only the storage grid changes.
      cache_integrals, device, dtype: same semantics as :meth:`__init__`.
    """
    from ._fit_tll import fit_tll_constant

    u_t = torch.as_tensor(u, dtype=dtype, device=device)
    grid_points, values = fit_tll_constant(
      u_t, grid_size=grid_size, mult=mult, grid_type=grid_type
    )
    return cls(
      grid_points=grid_points,
      values=values,
      cache_integrals=cache_integrals,
      norm_times=3,
      is_linear=(grid_type == "linear"),
      device=device,
      dtype=dtype,
    )

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
    # Trim to (1e-10, 1 - 1e-10), mirroring Bicop::prep_for_abstract.
    return u.clamp(_TRIM_LO, _TRIM_HI)

  def pdf(self, u: Tensor) -> Tensor:
    """Bivariate copula density ``c(u1, u2)`` at the query points ``u``.

    Args:
      u: ``(n, 2)`` tensor of pseudo-observations in ``[0, 1]^2``. Inputs
        outside the unit square are clamped to ``[1e-10, 1 - 1e-10]``;
        ``NaN`` propagates through the bilinear interpolation.

    Returns:
      ``(n,)`` tensor of density values; clamped to a strictly-positive
      floor of ``1e-20`` so subsequent ``log`` calls stay finite.
    """
    u = self._prep(u)
    if self.is_indep:
      return torch.ones(u.shape[0], dtype=u.dtype, device=u.device)
    return self.interp_grid.interpolate(u).clamp_min(1e-20)

  def log_pdf(self, u: Tensor) -> Tensor:
    """``log c(u1, u2)`` with safe handling of ``-inf`` / ``NaN``.

    Equivalent to ``pdf(u).log()`` but replaces ``-inf`` (from the
    density floor) with a fixed lower bound and ``+inf`` / ``NaN`` with
    finite sentinels — useful when the result feeds into an autograd
    loss.
    """
    return self.pdf(u).log().nan_to_num(neginf=_LOG_FLOOR, posinf=0.0)

  def cdf(self, u: Tensor) -> Tensor:
    """Bivariate copula CDF ``C(u1, u2) = ∫_0^{u1} ∫_0^{u2} c(s, t) ds dt``.

    Computed via :meth:`InterpolationGrid2D.integrate_2d` (nested
    trapezoidal integration) when ``cache_integrals=False``, or via a
    single bilinear interp on the precomputed cache when
    ``cache_integrals=True``.

    Args:
      u: ``(n, 2)`` tensor of pseudo-observations in ``[0, 1]^2``.

    Returns:
      ``(n,)`` tensor of CDF values in ``[1e-10, 1 - 1e-10]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return (u[:, 0] * u[:, 1]).clamp(_TRIM_LO, _TRIM_HI)
    if self._cdf_cache is not None:
      return self.interp_grid.interp_at(self._cdf_cache, u)
    return self.interp_grid.integrate_2d(u)

  # --------------------------------------------------------------------- #
  # h-functions                                                            #
  # --------------------------------------------------------------------- #

  def _hfunc_raw(self, u: Tensor, cond_var: int) -> Tensor:
    cache = self._hfunc1_cache if cond_var == 1 else self._hfunc2_cache
    if cache is not None:
      return self.interp_grid.interp_at(cache, u).clamp(_TRIM_LO, _TRIM_HI)
    return self.interp_grid.integrate_1d(u, cond_var=cond_var)

  def hfunc1(self, u: Tensor) -> Tensor:
    """First h-function: ``H1(u1, u2) = P(U2 ≤ u2 | U1 = u1)``.

    Computed via :meth:`InterpolationGrid2D.integrate_1d` with
    ``cond_var=1`` when ``cache_integrals=False``, or via a single
    bilinear interp on the precomputed ``hfunc1`` cache when
    ``cache_integrals=True``.

    Args:
      u: ``(n, 2)`` tensor of pseudo-observations in ``[0, 1]^2``.

    Returns:
      ``(n,)`` tensor of conditional CDF values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 1].clamp(_TRIM_LO, _TRIM_HI)
    return self._hfunc_raw(u, 1).clamp(0.0, 1.0)

  def hfunc2(self, u: Tensor) -> Tensor:
    """Second h-function: ``H2(u1, u2) = P(U1 ≤ u1 | U2 = u2)``.

    Symmetric to :meth:`hfunc1` with the conditioning variable swapped.

    Args:
      u: ``(n, 2)`` tensor of pseudo-observations in ``[0, 1]^2``.

    Returns:
      ``(n,)`` tensor of conditional CDF values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 0].clamp(_TRIM_LO, _TRIM_HI)
    return self._hfunc_raw(u, 2).clamp(0.0, 1.0)

  # --------------------------------------------------------------------- #
  # Inverse h-functions (numerical via bisection).                         #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def hinv1(self, u: Tensor) -> Tensor:
    """Inverse of :meth:`hfunc1` w.r.t. the second argument.

    Given ``u = [u1, p]`` of shape ``(n, 2)``, returns ``u2`` such that
    ``H1(u1, u2) = p``. With ``cache_integrals=True`` this is a single
    bilinear interp on the precomputed ``hinv1`` cache; otherwise each
    call runs the fixed-iter vectorized ITP root-finder
    (:func:`._util.solve_itp`) over the on-the-fly h-function.

    Args:
      u: ``(n, 2)`` tensor where column 0 is ``u1`` and column 1 is the
        target probability ``p``.

    Returns:
      ``(n,)`` tensor of ``u2`` values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 1].clamp(_TRIM_LO, _TRIM_HI)
    if self._hinv1_cache is not None:
      return self.interp_grid.interp_at(self._hinv1_cache, u).clamp(0.0, 1.0)
    a = u[:, 0:1]
    p = u[:, 1:2]

    def fun(x: Tensor) -> Tensor:
      u_eval = torch.cat([a, x], dim=-1)
      return self._hfunc_raw(u_eval, 1).unsqueeze(-1) - p

    x = solve_itp(fun, x_a=torch.zeros_like(p), x_b=torch.ones_like(p))
    return x.squeeze(-1).clamp(0.0, 1.0)

  @torch.no_grad()
  def hinv2(self, u: Tensor) -> Tensor:
    """Inverse of :meth:`hfunc2` w.r.t. the first argument.

    Given ``u = [p, u2]`` of shape ``(n, 2)``, returns ``u1`` such that
    ``H2(u1, u2) = p``. See :meth:`hinv1` for the cache vs. ITP-bisection
    semantics.

    Args:
      u: ``(n, 2)`` tensor where column 0 is the target probability ``p``
        and column 1 is ``u2``.

    Returns:
      ``(n,)`` tensor of ``u1`` values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 0].clamp(_TRIM_LO, _TRIM_HI)
    if self._hinv2_cache is not None:
      return self.interp_grid.interp_at(self._hinv2_cache, u).clamp(0.0, 1.0)
    p = u[:, 0:1]
    b = u[:, 1:2]

    def fun(x: Tensor) -> Tensor:
      u_eval = torch.cat([x, b], dim=-1)
      return self._hfunc_raw(u_eval, 2).unsqueeze(-1) - p

    x = solve_itp(fun, x_a=torch.zeros_like(p), x_b=torch.ones_like(p))
    return x.squeeze(-1).clamp(0.0, 1.0)

  # --------------------------------------------------------------------- #
  # Sampling                                                               #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def simulate(
    self,
    n: int = 100,
    qrng: bool = False,
    seeds: list[int] = [],
  ) -> Tensor:
    """Draw ``n`` joint samples from the fitted copula.

    Mirror of :meth:`pyvinecopulib.Bicop.simulate`. Uses the inverse
    Rosenblatt scheme: sample two independent uniforms ``(U1, P)``,
    then set ``U2 = hinv1((U1, P))`` so that ``(U1, U2)`` has the
    fitted joint distribution.

    Args:
      n: Number of samples to draw (must be ``> 0``).
      qrng: If ``True``, draw the base uniforms from a scrambled Sobol
        sequence (better low-discrepancy in 2-D) instead of
        pseudo-random uniforms.
      seeds: When ``qrng=True``, the first entry seeds the
        :class:`torch.quasirandom.SobolEngine` scramble. When
        ``qrng=False``, the first entry seeds the global torch RNG
        before the ``torch.rand`` call. Empty list keeps the existing
        global state.

    Returns:
      ``(n, 2)`` tensor of samples in ``(0, 1)^2``.
    """
    if n <= 0:
      raise ValueError(f"n must be > 0; got {n}")
    device = self.interp_grid.values.device
    dtype = self.interp_grid.values.dtype
    seed = int(seeds[0]) if seeds else None
    if qrng:
      kwargs = {"dimension": 2, "scramble": True}
      if seed is not None:
        kwargs["seed"] = seed
      u = (
        torch.quasirandom.SobolEngine(**kwargs)
        .draw(n=n, dtype=dtype)
        .to(device=device)
      )
    else:
      if seed is not None:
        torch.manual_seed(seed)
      u = torch.rand(n, 2, dtype=dtype, device=device)
    if self.is_indep:
      return u
    u2 = self.hinv1(u).unsqueeze(-1)
    return torch.cat([u[:, 0:1], u2], dim=-1)

  @torch.no_grad()
  def sample(
    self,
    num_sample: int = 100,
    seed: Optional[int] = 42,
    is_sobol: bool = False,
  ) -> Tensor:
    """Deprecated alias for :meth:`simulate`.

    .. deprecated::
       Use :meth:`simulate` with ``n``, ``qrng``, ``seeds``. The old
       parameter names are kept as a thin pass-through within the
       ``feature/torch-bicop`` branch and will be removed before the
       next stable release.
    """
    return self.simulate(
      n=num_sample,
      qrng=is_sobol,
      seeds=[seed] if seed is not None else [],
    )
