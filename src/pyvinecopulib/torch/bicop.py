"""PyTorch evaluator for a bivariate pair copula on a density grid.

The evaluation chain (``pdf`` / ``cdf`` / ``hfunc`` / ``hinv`` /
``sample``) lives entirely in PyTorch, so the pair copula can be
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
from ._interp import InterpolationGrid2D, _trim
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
      unsupported (rebuild a new ``TorchBicop`` instead). For the same
      reason the cached members are exact, differentiable functions of
      ``u`` but carry **no** gradient with respect to ``values`` -- the
      tables are constants. ``pdf`` never uses a cache and always does.
      Pass ``cache_integrals=False`` to differentiate the integrals with
      respect to the density grid.
  norm_maxiter : int, default=25
      Maximum number of margin-rescaling passes; rescaling also stops
      as soon as both margins integrate to 1 within ``1e-10``. Matches
      the ``Bicop`` TLL default. Pass ``0`` to skip when the grid
      already integrates to uniform margins.
  is_linear : bool, default=False
      Internal flag selecting the linear-grid fast-path in the
      underlying ``InterpolationGrid2D``. Set by
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
  # Class-level hints so the prefix tables registered in __init__ are
  # statically typed as Tensors instead of nn.Module. They are set together:
  # all three are present or all three are None.
  _sy: Tensor | None
  _sx: Tensor | None
  _prefix: Tensor | None
  #: TorchBicop exposes the grid/cache internals the batched vine path needs.
  supports_batched: bool = True

  def __init__(
    self,
    grid_points: Optional[Tensor] = None,
    values: Optional[Tensor] = None,
    cache_integrals: bool = True,
    norm_maxiter: int = 25,
    is_linear: bool = False,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> None:
    # Initialize nn.Module explicitly: TorchBicop also subclasses BicopBase
    # (a Protocol-derived ABC), whose __init__ chain would otherwise shadow
    # nn.Module's under super().
    torch.nn.Module.__init__(self)
    if grid_points is None and values is None:
      # Independent-copula short-circuit: pdf=1 everywhere, cdf=u1*u2, etc.
      self.is_indep = True
      m = 2
      grid_points = torch.tensor([0.0, 1.0], dtype=dtype, device=device)
      values = torch.ones(m, m, dtype=dtype, device=device)
      norm_maxiter_eff = 0
    else:
      if grid_points is None or values is None:
        raise ValueError(
          "Provide either both grid_points and values, or neither (for the "
          "independence copula)."
        )
      self.is_indep = False
      grid_points = torch.as_tensor(grid_points, dtype=dtype, device=device)
      values = torch.as_tensor(values, dtype=dtype, device=device)
      norm_maxiter_eff = norm_maxiter

    self.interp_grid = InterpolationGrid2D(
      grid_points=grid_points,
      values=values,
      norm_maxiter=norm_maxiter_eff,
      is_linear=is_linear,
    )

    self._cache_integrals = bool(cache_integrals)
    if self._cache_integrals and not self.is_indep:
      # Detached, because a buffer is a cache and holding a graph in one would
      # keep it alive for the module's lifetime. `_tables` rebuilds them inside
      # the graph whenever a gradient is actually being taken.
      sy, sx, pref = (t.detach() for t in self.interp_grid.build_caches())
      self.register_buffer("_sy", sy)
      self.register_buffer("_sx", sx)
      self.register_buffer("_prefix", pref)
    else:
      self._sy = None
      self._sx = None
      self._prefix = None

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
    differentiable bilinear interpolator: autograd flows through every
    ``pdf`` / ``cdf`` / ``hfunc`` / ``hinv`` call in ``u``, and through
    ``pdf`` in the density grid. With ``cache_integrals=True`` the other
    four read a frozen table, so they carry no gradient in the grid; see
    that parameter on ``TorchBicop.__init__``.

    Parameters
    ----------
    cop : Bicop
        A fitted ``Bicop`` of the TLL family at
        ``rotation=0`` (the only shape ``TorchBicop`` represents
        directly). The density values are taken straight from
        ``cop.parameters``; the grid coordinates come from
        ``InterpolationGrid2D.make_grid_points`` with the canonical
        Phi-spaced normal grid ``Bicop`` uses. The grid is
        already normalized, so renormalization is skipped (parity
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
      "normal", m, dtype=dtype, device=device
    )
    values = torch.as_tensor(values_np, dtype=dtype, device=device)
    # The grid stored on cop is already normalized; skip renormalization
    # to avoid drifting away from the reference ``Bicop`` values.
    return cls(
      grid_points=grid_points,
      values=values,
      cache_integrals=cache_integrals,
      norm_maxiter=0,
      device=device,
      dtype=dtype,
    )

  @classmethod
  def from_data_batched(
    cls,
    u: Tensor,
    controls: Optional[FitControlsTorchBicop] = None,
    *,
    cache_integrals: bool = True,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "list[TorchBicop]":
    """Fit ``P`` pair copulas from one stacked sample, in one call.

    The pairs share the sample size, the interpolation grid and the
    evaluation points, so everything that does not depend on the data is
    computed once and the two convergence loops of the bandwidth search
    advance all ``P`` lanes together, each freezing as it converges. What
    comes back is ``P`` ordinary :class:`TorchBicop` objects; the batching
    is confined to the fit.

    Nothing here knows what the pairs are *for*: the leading axis is just
    ``P`` independent pairs on shared rows. One vine's tree level is the
    obvious caller, but several vines' levels concatenate into the same
    axis just as well.

    Parameters
    ----------
    u : Tensor, shape (P, n, 2), dtype float
        Copula-scale observations, one ``(n, 2)`` block per pair. Continuous
        only -- a discrete edge's fit reuses a compiled per-pair draw that
        has no batch axis.
    controls : FitControlsTorchBicop, optional
        Fit controls, shared by every pair. Defaults to
        ``FitControlsTorchBicop()``.
    cache_integrals : bool, default=True
        Precompute each returned pair's integral tables.
    device : torch.device, optional
        Placement; defaults to ``u``'s.
    dtype : torch.dtype, default=torch.float64
        Working precision.

    Returns
    -------
    list of TorchBicop
        One fitted pair copula per leading index, in order.

    Raises
    ------
    ValueError
        If ``u`` is not 3-d with two value columns.

    Notes
    -----
    A pair's fit is unaffected by *which* other pairs share its call: each
    lane's bandwidth search freezes as it converges, so the iterations a
    pair takes are the ones its own data earns.

    It is affected by *how many*, in the last bits. Torch selects
    elementwise kernels by element count, and the bandwidth search's
    ``pow`` takes a vectorized path past
    ``2 * Vectorized<double>::size()`` lanes -- 8 on AVX2, 4 on NEON -- so
    stacking ``P`` pairs agrees with fitting them one by one to floating
    point rather than bit for bit, on every device. The gap is around
    ``1e-15`` on a grid value. Where that matters, fit the pair alone.

    See Also
    --------
    TorchBicop.from_data : The single-pair entry point.
    """
    if controls is None:
      controls = FitControlsTorchBicop()
    u_t = torch.as_tensor(u, dtype=dtype, device=device)
    if u_t.ndim != 3 or u_t.shape[-1] != 2:
      raise ValueError(f"u must have shape (P, n, 2); got {tuple(u_t.shape)}")
    u_t = _trim(u_t)

    from ._fit_tll import fit_tll_constant

    grid_points, values = fit_tll_constant(
      u_t,
      grid_size=controls.grid_size,
      mult=controls.mult,
      grid_type=controls.grid_type,
      compile_fit=controls.compile_fit,
    )
    # `values` is `(P, m, m)`; each pair takes its own slice. The slice is
    # contiguous, so the grid each pair holds is its own buffer rather than
    # a view into a tensor the others share.
    return [
      cls(
        grid_points=grid_points,
        # A slice of the stack is a view; `contiguous` gives each pair a
        # buffer of its own so the stack can be freed and `state_dict` does
        # not carry P-1 unrelated grids per pair.
        values=values[i].contiguous(),
        cache_integrals=cache_integrals,
        norm_maxiter=25,
        is_linear=(controls.grid_type == "linear"),
        device=device,
        dtype=dtype,
      )
      for i in range(values.shape[0])
    ]

  @classmethod
  def from_data(
    cls,
    u,
    controls: Optional[FitControlsTorchBicop] = None,
    *,
    cache_integrals: bool = True,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
    var_types: Optional[list[str]] = None,
  ) -> "TorchBicop":
    """Fits a bicop on pseudo-observations and wraps in a ``TorchBicop``.

    Dispatches on ``controls.method`` (``"tll"``: pure-torch
    Transformed Local Likelihood, matching the ``Bicop``
    TLL fit to machine precision).

    Parameters
    ----------
    u : ndarray or Tensor, shape (n, 2) or (n, 4), dtype float
        Pseudo-observations, ``[u1, u2]`` — or ``[u1, u2, u1^-, u2^-]`` when
        ``var_types`` marks an argument discrete, of which only the values are
        fitted on.
    controls : FitControlsTorchBicop or None, default=None
        Fit-time controls. `None` defaults to TLL with
        ``grid_size=30`` on the normal-spaced grid.
    cache_integrals : bool, default=True
        See ``TorchBicop.__init__``.
    device : torch.device or None, default=None
        See ``TorchBicop.__init__``.
    dtype : torch.dtype, default=torch.float64
        See ``TorchBicop.__init__``.
    var_types : list of str or None, default=None
        The two arguments' types, ``"c"`` or ``"d"``; ``None`` means both
        continuous. A discrete argument changes only how ties are broken when
        ranking, and the fitted grid is continuous either way — the mixed-
        discrete surface is supplied by
        :class:`~pyvinecopulib.core.DiscretePair`.

    Returns
    -------
    TorchBicop
        A fitted ``TorchBicop``.

    Raises
    ------
    ValueError
        If ``u``'s column count does not match ``var_types``.
    """
    if controls is None:
      controls = FitControlsTorchBicop()
    types = ("c", "c") if var_types is None else tuple(var_types)
    discrete = "d" in types
    expected = 4 if discrete else 2

    u_t = torch.as_tensor(u, dtype=dtype, device=device)
    if u_t.ndim != 2 or u_t.shape[1] != expected:
      raise ValueError(
        f"u must have shape (n, {expected}) for var_types={list(types)}; "
        f"got {tuple(u_t.shape)}"
      )
    # ``Bicop.select`` trims before ``TllBicop::fit``, so two values above
    # ``1 - 1e-10`` are one tie group there and would be two here.
    u_t = _trim(u_t)
    values_only = u_t[:, :2]
    # An atom repeats its distribution-function value, so the ranks have ties.
    # ``TllBicop::fit`` breaks them at random from a fixed seed; reuse that draw
    # rather than reimplementing it, as structure selection reuses ``wdm``. On a
    # discrete edge those ranks only select the bandwidth: the fit itself runs on
    # the latent sample drawn with it, which ``fit_tll_constant`` handles.
    pseudo_obs = None
    if discrete:
      from ..utils import to_pseudo_obs

      pseudo_obs = torch.as_tensor(
        to_pseudo_obs(
          values_only.detach().cpu().numpy(),
          ties_method="random",
          seeds=[5],
        ),
        dtype=u_t.dtype,
        device=u_t.device,
      )

    # ``method`` is validated to be "tll" by FitControlsTorchBicop; it is
    # kept as the dispatch seam for future torch fitters.
    from ._fit_tll import fit_tll_constant

    grid_points, values = fit_tll_constant(
      values_only,
      grid_size=controls.grid_size,
      mult=controls.mult,
      grid_type=controls.grid_type,
      compile_fit=controls.compile_fit,
      pseudo_obs=pseudo_obs,
      discrete_data=u_t if discrete else None,
    )
    return cls(
      grid_points=grid_points,
      values=values,
      cache_integrals=cache_integrals,
      norm_maxiter=25,
      is_linear=(controls.grid_type == "linear"),
      device=device,
      dtype=dtype,
    )

  def flip(self) -> "TorchBicop":
    """Return the copula with its two arguments swapped (``c'(u,v)=c(v,u)``).

    Transposes the (already margin-normalized) interpolation grid — the same
    reorientation :class:`~pyvinecopulib.core.Bicop` applies to a ``tll`` pair
    — so the result is a fresh :class:`~pyvinecopulib.torch.TorchBicop`, an
    ``nn.Module`` that keeps the batched fast path. Used by
    :meth:`~pyvinecopulib.core.VinecopBase.select` to reorient a selection-time
    pair onto its slot in the finalized structure.

    Returns
    -------
    TorchBicop
        The argument-swapped copula; the object itself is left unchanged.
    """
    device = self.interp_grid.values.device
    dtype = self.interp_grid.values.dtype
    if self.is_indep:
      return TorchBicop(
        cache_integrals=self._cache_integrals, device=device, dtype=dtype
      )
    return TorchBicop(
      grid_points=self.interp_grid.grid_points,
      values=self.interp_grid.values.transpose(0, 1).contiguous(),
      cache_integrals=self._cache_integrals,
      norm_maxiter=0,
      is_linear=self.interp_grid._is_linear,
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
    return _trim(u)

  def pdf(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
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

  def cdf(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
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
      return _trim(u[:, 0] * u[:, 1])
    if self._sy is not None:
      sy, sx, pref = self._tables()
      return self.interp_grid.cdf_cached(u, sy, sx, pref)
    return self.interp_grid.integrate_2d(u)

  def rect_mass(self, a1: Tensor, b1: Tensor, a2: Tensor, b2: Tensor) -> Tensor:
    """The probability of ``(a1, b1] x (a2, b2]``, without the cancellation.

    An optional capability a discrete edge discovers with ``getattr``: it lets
    the mixed-discrete density read the atom's probability directly instead of
    differencing four :meth:`cdf` values, which amplifies any absolute error by
    ``~4 / (w1 w2)`` in the atom widths where this route amplifies by
    ``1 / w2`` alone. Available in both cache modes -- it reads the grid, not
    the prefix tables.

    Parameters
    ----------
    a1, b1, a2, b2 : Tensor, shape (n,), dtype float
        Rectangle bounds per query. An empty or inverted interval gives zero.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Rectangle probabilities. Independence returns
        ``(b1 - a1) * (b2 - a2)``.
    """
    # Coerced but *not* trimmed: `[1e-10, 1-1e-10]` is the guard `cdf` applies
    # to a query point, and a rectangle's lower bound is legitimately zero --
    # trimming it drops a sliver of width `1e-10` per axis, which is exactly
    # what a corner rectangle would then disagree with `cdf` by. The weights
    # clamp to `[0, 1]` themselves.
    ref = self.interp_grid.values
    a1, b1, a2, b2 = (
      torch.as_tensor(t, dtype=ref.dtype, device=ref.device)
      for t in (a1, b1, a2, b2)
    )
    if self.is_indep:
      return (b1 - a1).clamp_min(0.0) * (b2 - a2).clamp_min(0.0)
    return self.interp_grid.rect_mass(a1, b1, a2, b2)

  def _tables(self) -> tuple[Tensor, Tensor, Tensor]:
    """The prefix tables, rebuilt in-graph when a gradient is being taken.

    They are cumulative sums of ``values``, so the closed forms that read them
    are exact functions of the grid -- but the stored copies are detached
    buffers, and differentiating through those would silently drop most of the
    gradient rather than all of it. They are therefore recomputed whenever grad
    is live: three cumulative sums over ``(m, m)``, and never on the
    no-gradient path.
    """
    assert self._sy is not None
    assert self._sx is not None and self._prefix is not None
    if torch.is_grad_enabled() and self.interp_grid.values.requires_grad:
      return self.interp_grid.build_caches()
    return self._sy, self._sx, self._prefix

  # --------------------------------------------------------------------- #
  # h-functions                                                            #
  # --------------------------------------------------------------------- #

  def _hfunc_raw(self, u: Tensor, cond_var: int) -> Tensor:
    if self._sy is not None:
      sy, sx, _ = self._tables()
      # Conditioning on the second argument reads the cumulative along the
      # first, which is `sx` laid out the way `hfunc_cached` indexes it.
      return self.interp_grid.hfunc_cached(
        u, cond_var, sy if cond_var == 1 else sx.t()
      )
    return self.interp_grid.integrate_1d(u, cond_var=cond_var)

  def hfunc1(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
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
      return _trim(u[:, 1])
    return self._hfunc_raw(u, 1).clamp(0.0, 1.0)

  def hfunc2(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
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
      return _trim(u[:, 0])
    return self._hfunc_raw(u, 2).clamp(0.0, 1.0)

  # --------------------------------------------------------------------- #
  # Inverse h-functions (closed-form conditional quantiles).               #
  # --------------------------------------------------------------------- #

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
    # Always the closed-form inversion. Tabulating the inverse at the grid
    # nodes and interpolating between them was the one cached member whose
    # error was not merely a tolerance -- up to 1e-2 -- and unlike `cdf` and
    # the h-functions it has no O(1) exact reconstruction: locating the
    # bracketing cell needs the conditional cumulative along the whole free
    # axis, which is O(m) to assemble whatever is cached. So the cache buys
    # O(1) where it can, and consistency where it cannot.
    # The prefix table is the conditional cumulative the inversion needs to
    # locate its bracketing cell, so pass it when there is one: the same
    # quantity, a gather instead of a trapezoid and a scan.
    cum = None if self._sy is None else self._tables()[cond_var - 1]
    return self.interp_grid.inverse_integrate_1d(u, cond_var, cum).clamp(
      0.0, 1.0
    )

  def hinv1(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
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
      return _trim(u[:, 1])
    return self._hinv_raw(u, 1)

  def hinv2(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
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
      return _trim(u[:, 0])
    return self._hinv_raw(u, 2)

  # --------------------------------------------------------------------- #
  # Sampling                                                               #
  # --------------------------------------------------------------------- #

  @torch.no_grad()
  def sample(
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
      # A device-local generator: ``torch.manual_seed`` reseeds the global
      # CPU generator and every CUDA one, so seeding a single pair copula
      # would perturb every other RNG consumer in the process.
      gen = (
        torch.Generator(device=device).manual_seed(seed)
        if seed is not None
        else None
      )
      u = torch.rand(n, 2, generator=gen, dtype=dtype, device=device)
    if self.is_indep:
      return u
    u2 = self.hinv1(u).unsqueeze(-1)
    return torch.cat([u[:, 0:1], u2], dim=-1)
