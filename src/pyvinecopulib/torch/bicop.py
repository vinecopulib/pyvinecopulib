"""PyTorch evaluator for a bivariate pair copula on a density grid.

The evaluation chain (``pdf`` / ``cdf`` / ``hfunc`` / ``hinv`` /
``simulate``) lives entirely in PyTorch, so the pair copula can be
moved to GPU with ``.to("cuda")`` and composed with autograd-aware
downstream code.

Three constructors are provided:

* :meth:`TorchBicop.from_data` — fit on pseudo-observations directly
  in PyTorch. Dispatches on the ``method`` field of a
  :class:`FitControlsTorchBicop`. ``"tll"`` (default) is the
  *Transformed Local Likelihood* kernel density estimator (Geenens
  2014; Nagler 2018), the non-parametric family exposed as
  ``tll``; this path matches the
  ``Bicop.from_data()`` TLL fit to machine precision.
* :meth:`TorchBicop.from_bicop` — lift a fitted
  ``Bicop`` (TLL family) into a ``TorchBicop``.
  Useful when you already fit with ``Bicop`` and
  want GPU / autograd evaluation downstream.
* ``TorchBicop(grid_points=..., values=...)`` — construct from an
  externally-prepared ``(m, m)`` density grid.

See Also
--------
pyvinecopulib.core.Bicop : Reference pair copula.
FitControlsTorchBicop : Fit-time controls.
"""

from __future__ import annotations

from typing import Optional

import torch
from torch import Tensor

from ..core import BicopBase
from ..pyvinecopulib_ext import Bicop, tll as _TLL_FAMILY
from ._interp import InterpolationGrid2D, _TRIM_LO, _TRIM_HI
from .controls import FitControlsTorchBicop


class TorchBicop(BicopBase[torch.Tensor], torch.nn.Module):
  """PyTorch evaluator for a bivariate copula stored as a density grid.

  Torch counterpart of the non-parametric pair-copula path in
  ``Bicop`` (the *Transformed Local Likelihood* /
  ``tll`` family). The fitted density lives on
  an ``m x m`` grid in ``[0, 1]^2`` and is evaluated by bilinear
  interpolation. Non-zero copula rotations are not supported — TLL
  pair-copulas always have ``rotation=0``.

  Parameters
  ----------
  grid_points : Tensor, shape (m,), dtype float
      Strictly increasing 1-D tensor of grid points on ``[0, 1]``
      (the same grid is used along both axes). Endpoints are
      clipped to exactly ``0`` and ``1`` to avoid extrapolation.
  values : Tensor, shape (m, m), dtype float
      Density values on the tensor-product grid.
  cache_integrals : bool, default=True
      If ``True``, precompute ``cdf`` / ``hfunc1`` / ``hfunc2`` /
      ``hinv1`` / ``hinv2`` at every grid node so subsequent calls
      are a single bilinear lookup (~80–300x faster with a ~1e-3
      mean IAE cost relative to the on-the-fly trapezoidal +
      bisection path). The caches are built once at construction
      from ``values`` and are not refreshed afterwards; mutating
      ``interp_grid.values`` in place on a cached instance is
      unsupported (rebuild a new ``TorchBicop`` instead).
  norm_times : int, default=3
      Number of margin-normalization rounds. Matches the ``Bicop``
      TLL default. Pass ``0`` to skip when the grid already integrates
      to uniform margins.
  is_linear : bool, default=False
      Internal flag selecting the linear-grid fast-path in the
      underlying :class:`InterpolationGrid2D`. Set by
      :meth:`from_data` when ``grid_type="linear"``; users
      normally do not pass it directly.
  device : torch.device or None, default=None
      Placement of the underlying tensors. ``None`` keeps them on
      the input's device.
  dtype : torch.dtype, default=torch.float64
      Precision of the underlying tensors. ``torch.float64`` mirrors
      the ``Bicop`` evaluation.
  """

  is_indep: bool
  #: TorchBicop exposes the grid/cache internals the batched vine path needs.
  supports_batched: bool = True

  def __init__(
    self,
    grid_points: Optional[Tensor] = None,
    values: Optional[Tensor] = None,
    cache_integrals: bool = True,
    norm_times: int = 3,
    is_linear: bool = False,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> None:
    # Initialise nn.Module explicitly: TorchBicop also subclasses BicopBase
    # (a Protocol-derived ABC), whose __init__ chain would otherwise shadow
    # nn.Module's under super().
    torch.nn.Module.__init__(self)
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
    cache_integrals: bool = True,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchBicop":
    """Lifts a fitted ``Bicop`` into a ``TorchBicop``.

    The resulting ``TorchBicop`` is a ``torch.nn.Module`` (so
    ``.to("cuda")`` moves the density grid in one line) built from a
    differentiable bilinear interpolator, with autograd flowing
    through every ``pdf`` / ``cdf`` / ``hfunc`` / ``hinv`` call.

    Parameters
    ----------
    cop : Bicop
        A fitted ``Bicop`` of the TLL family at
        ``rotation=0`` (the only shape ``TorchBicop`` represents
        directly). The density values are taken straight from
        ``cop.parameters``; the grid coordinates come from
        ``InterpolationGrid2D.make_grid_points`` with the canonical
        Phi-spaced normal grid ``Bicop`` uses. The grid is
        already normalised, so renormalisation is skipped (parity
        is typically ``< 1e-12`` per cell on fresh fits).
    cache_integrals : bool, default=True
        See ``TorchBicop.__init__``.
    device : torch.device or None, default=None
        See ``TorchBicop.__init__``.
    dtype : torch.dtype, default=torch.float64
        See ``TorchBicop.__init__``.

    Returns
    -------
    TorchBicop
        A ``TorchBicop`` mirroring ``cop``.
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
    # to avoid drifting away from the reference ``Bicop`` values.
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
    controls: Optional[FitControlsTorchBicop] = None,
    *,
    cache_integrals: bool = True,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchBicop":
    """Fits a bicop on pseudo-observations and wraps in a ``TorchBicop``.

    Dispatches on ``controls.method`` (``"tll"``: pure-torch
    Transformed Local Likelihood, matching the ``Bicop``
    TLL fit to machine precision).

    Parameters
    ----------
    u : ndarray or Tensor, shape (n, 2), dtype float
        Pseudo-observations.
    controls : FitControlsTorchBicop or None, default=None
        Fit-time controls. `None` defaults to TLL with
        ``grid_size=30`` on the normal-spaced grid.
    cache_integrals : bool, default=True
        See ``TorchBicop.__init__``.
    device : torch.device or None, default=None
        See ``TorchBicop.__init__``.
    dtype : torch.dtype, default=torch.float64
        See ``TorchBicop.__init__``.

    Returns
    -------
    TorchBicop
        A fitted ``TorchBicop``.
    """
    if controls is None:
      controls = FitControlsTorchBicop()

    u_t = torch.as_tensor(u, dtype=dtype, device=device)
    if u_t.ndim != 2 or u_t.shape[1] != 2:
      raise ValueError(f"u must have shape (n, 2); got {tuple(u_t.shape)}")

    # ``method`` is validated to be "tll" by FitControlsTorchBicop; it is
    # kept as the dispatch seam for future torch fitters.
    from ._fit_tll import fit_tll_constant

    grid_points, values = fit_tll_constant(
      u_t,
      grid_size=controls.grid_size,
      mult=controls.mult,
      grid_type=controls.grid_type,
    )
    return cls(
      grid_points=grid_points,
      values=values,
      cache_integrals=cache_integrals,
      norm_times=3,
      is_linear=(controls.grid_type == "linear"),
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

  def pdf(self, u: Tensor, x: Optional[Tensor] = None) -> Tensor:
    """Evaluates the bivariate copula density ``c(u1, u2)``.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Pseudo-observations in ``[0, 1]^2``. Inputs outside the
        unit square are clamped to ``[1e-10, 1 - 1e-10]``; ``NaN``
        propagates through the interpolation.
    x : Tensor or None, optional
        Conditioning variables, shape ``(n, k)``. Ignored by ``TorchBicop``
        (an unconditional pair copula); accepted so the class satisfies the
        :class:`~pyvinecopulib.core.BicopLike` contract.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Density values, clamped to a strictly-positive floor of
        ``1e-20``.
    """
    u = self._prep(u)
    if self.is_indep:
      return torch.ones(u.shape[0], dtype=u.dtype, device=u.device)
    return self.interp_grid.interpolate(u).clamp_min(1e-20)

  def cdf(self, u: Tensor, x: Optional[Tensor] = None) -> Tensor:
    """Evaluates the bivariate copula CDF.

    .. math::

       C(u_1, u_2) = \\int_0^{u_1} \\int_0^{u_2} c(s, t)\\, ds\\, dt.

    Computed via nested trapezoidal integration when
    ``cache_integrals=False``, or via a single bilinear
    interpolation on the precomputed cache when
    ``cache_integrals=True``.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Pseudo-observations in ``[0, 1]^2``.
    x : Tensor or None, optional
        Conditioning variables, shape ``(n, k)``. Ignored by ``TorchBicop``
        (an unconditional pair copula); accepted so the class satisfies the
        :class:`~pyvinecopulib.core.BicopLike` contract.

    Returns
    -------
    Tensor, shape (n,), dtype float
        CDF values in ``[1e-10, 1 - 1e-10]``.
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

  def hfunc1(self, u: Tensor, x: Optional[Tensor] = None) -> Tensor:
    """Evaluates the first h-function.

    .. math::

       h_1(u_1, u_2) = \\mathbb{P}(U_2 \\le u_2 \\mid U_1 = u_1).

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Pseudo-observations in ``[0, 1]^2``.
    x : Tensor or None, optional
        Conditioning variables, shape ``(n, k)``. Ignored by ``TorchBicop``
        (an unconditional pair copula); accepted so the class satisfies the
        :class:`~pyvinecopulib.core.BicopLike` contract.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Conditional CDF values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 1].clamp(_TRIM_LO, _TRIM_HI)
    return self._hfunc_raw(u, 1).clamp(0.0, 1.0)

  def hfunc2(self, u: Tensor, x: Optional[Tensor] = None) -> Tensor:
    """Evaluates the second h-function.

    .. math::

       h_2(u_1, u_2) = \\mathbb{P}(U_1 \\le u_1 \\mid U_2 = u_2).

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Pseudo-observations in ``[0, 1]^2``.
    x : Tensor or None, optional
        Conditioning variables, shape ``(n, k)``. Ignored by ``TorchBicop``
        (an unconditional pair copula); accepted so the class satisfies the
        :class:`~pyvinecopulib.core.BicopLike` contract.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Conditional CDF values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 0].clamp(_TRIM_LO, _TRIM_HI)
    return self._hfunc_raw(u, 2).clamp(0.0, 1.0)

  # --------------------------------------------------------------------- #
  # Inverse h-functions (closed-form conditional quantiles).               #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def _hinv_raw(self, u: Tensor, cond_var: int) -> Tensor:
    """Shared inverse-h-function body for `hinv1` / `hinv2`.

    ``cond_var=1`` solves ``H1(u1, x) = p`` for ``x`` (free column 1,
    ``u = [u1, p]``); ``cond_var=2`` solves ``H2(x, u2) = p`` (free
    column 0, ``u = [p, u2]``). With a precomputed cache this is one
    bilinear interpolation; otherwise the exact closed-form inversion of
    the piecewise-quadratic conditional cdf
    (:meth:`InterpolationGrid2D.inverse_integrate_1d`, mirroring
    vinecopulib#691).
    """
    cache = self._hinv1_cache if cond_var == 1 else self._hinv2_cache
    if cache is not None:
      return self.interp_grid.interp_at(cache, u).clamp(0.0, 1.0)
    return self.interp_grid.inverse_integrate_1d(u, cond_var).clamp(0.0, 1.0)

  @torch.no_grad()
  def hinv1(self, u: Tensor, x: Optional[Tensor] = None) -> Tensor:
    """Inverts `hfunc1` w.r.t. the second argument.

    Given ``u = [u1, p]``, returns ``u2`` such that
    ``H1(u1, u2) = p``. With ``cache_integrals=True`` this is a
    single bilinear interpolation on the precomputed cache;
    otherwise each call inverts the piecewise-quadratic conditional
    cdf in closed form (exact inverse of the on-the-fly h-function).

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Column 0 is ``u1``; column 1 is the target probability
        ``p``.
    x : Tensor or None, optional
        Conditioning variables, shape ``(n, k)``. Ignored by ``TorchBicop``
        (an unconditional pair copula); accepted so the class satisfies the
        :class:`~pyvinecopulib.core.BicopLike` contract.

    Returns
    -------
    Tensor, shape (n,), dtype float
        ``u2`` values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 1].clamp(_TRIM_LO, _TRIM_HI)
    return self._hinv_raw(u, 1)

  @torch.no_grad()
  def hinv2(self, u: Tensor, x: Optional[Tensor] = None) -> Tensor:
    """Inverts `hfunc2` w.r.t. the first argument.

    Given ``u = [p, u2]``, returns ``u1`` such that
    ``H2(u1, u2) = p``. See `hinv1` for the cache vs. closed-form
    semantics.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Column 0 is the target probability ``p``; column 1 is
        ``u2``.
    x : Tensor or None, optional
        Conditioning variables, shape ``(n, k)``. Ignored by ``TorchBicop``
        (an unconditional pair copula); accepted so the class satisfies the
        :class:`~pyvinecopulib.core.BicopLike` contract.

    Returns
    -------
    Tensor, shape (n,), dtype float
        ``u1`` values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return u[:, 0].clamp(_TRIM_LO, _TRIM_HI)
    return self._hinv_raw(u, 2)

  # --------------------------------------------------------------------- #
  # Sampling                                                               #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def simulate(
    self,
    n: int = 100,
    *,
    x: Optional[Tensor] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> Tensor:
    """Draws ``n`` joint samples from the fitted copula.

    Uses the inverse Rosenblatt scheme: sample two independent
    uniforms ``(U1, P)`` and set ``U2 = hinv1((U1, P))`` so that
    ``(U1, U2)`` has the fitted joint distribution.

    Parameters
    ----------
    n : int, default=100
        Number of samples to draw (must be ``> 0``).
    x : Tensor or None, optional
        Conditioning variables. Ignored by ``TorchBicop`` (an unconditional
        pair copula); accepted for signature parity with the ``BicopLike``
        contract.
    qrng : bool, default=False
        If ``True``, draw the base uniforms from a scrambled Sobol
        sequence instead of pseudo-random uniforms.
    seeds : list of int or None, optional
        When ``qrng=True`` the first entry seeds the
        ``torch.quasirandom.SobolEngine`` scramble; when
        ``qrng=False`` it seeds the global torch RNG before the
        ``torch.rand`` call. ``None`` keeps the existing global state.

    Returns
    -------
    Tensor, shape (n, 2), dtype float
        Samples in ``(0, 1)^2``.
    """
    del x
    seeds = list(seeds) if seeds else []
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
