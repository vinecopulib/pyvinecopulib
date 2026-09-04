"""PyTorch pair copula on a bilinear interpolation grid.

``TorchBicop``, the only class here, holds a bivariate copula density on an
``m x m`` grid over the unit square and is bilinear between its nodes. It is a
``BicopBase`` over tensors and a ``torch.nn.Module``, so ``pdf`` / ``cdf`` /
the two h-functions / their inverses / ``sample`` stay in PyTorch throughout:
the pair copula moves to an accelerator with ``.to("cuda")``, joins a
``state_dict``, and carries autograd in its arguments and -- on request -- in
its density grid.

One is built from a grid directly, lifted from a fitted ``Bicop`` with
``TorchBicop.from_bicop()``, or fitted from pseudo-observations with
``TorchBicop.from_data()`` -- and ``TorchBicop.from_data_batched()`` fits
``P`` pairs at once. Fitting reproduces the ``tll`` fit of ``Bicop`` to
floating-point tolerance.

The grid arithmetic all of them share -- interpolation, the integral tables, a
rectangle's mass -- lives in the internal ``_interp.py``, and the TLL fitter
in ``_fit_tll.py``.

See Also
--------
pyvinecopulib.core.Bicop : Reference pair copula.
pyvinecopulib.core.BicopBase : The contract this completes.
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
  """Bivariate copula held as a density grid on the unit square, in PyTorch.

  The torch counterpart of the nonparametric ``tll`` (*transformed local
  likelihood*) pair copula of ``Bicop``: the density sits on an ``m x m``
  tensor-product grid and is bilinear between its nodes. There is no
  rotation -- a ``tll`` pair copula always has ``rotation=0``.

  Three ways to get one, each building on the one before: hand it a grid with
  the parameters below (omit both for the independence copula), lift a fitted
  ``Bicop`` with ``TorchBicop.from_bicop()``, or fit pseudo-observations with
  ``TorchBicop.from_data()`` -- ``TorchBicop.from_data_batched()`` for ``P``
  pairs in one call. A grid arrives at construction, which is why those
  factories are the fitting entry points: the in-place ``fit`` inherited from
  ``BicopBase`` is not implemented here, and ``select``, which defaults to it,
  would have nothing to choose in any case -- a density grid has no family.

  Beyond ``BicopBase`` it brings its own exact ``cdf``, h-functions and their
  inverses, ``sample``, ``flip``, and the extra
  ``TorchBicop.rect_mass()``; ``loglik`` and ``plot`` are inherited. The grid
  is a buffer rather than a parameter -- a fitted density, not a learned one --
  so optimizing it is opt-in: set ``requires_grad_(True)`` on the grid values
  and every method carries a gradient in them.

  Parameters
  ----------
  grid_points : Tensor, shape (m,), or None, default=None
      Strictly increasing grid points in ``[0, 1]``, shared by both axes; the
      first and last are taken as exactly ``0`` and ``1``, so no evaluation
      extrapolates. ``None`` together with ``values`` gives the independence
      copula.
  values : Tensor, shape (m, m), or None, default=None
      Nonnegative density values at the tensor-product grid. The
      nonnegativity is a precondition rather than a convention: the integral
      tables and ``TorchBicop.rect_mass()`` are free of cancellation only
      because every term they sum is nonnegative.
  cache_integrals : bool, default=True
      Precompute the three ``(m, m)`` cumulative-trapezoid prefix tables that
      ``cdf`` / ``hfunc1`` / ``hfunc2`` then read their value from in closed
      form. The reconstruction is exact rather than approximate -- the
      interpolant is bilinear, so its integral along a grid line is piecewise
      linear across cells -- and it is differentiable in ``values`` as well as
      in ``u``, so the two modes agree to summation order and carry the same
      gradients. ``hinv1`` / ``hinv2`` invert the same closed form either way,
      reading a table only to locate the cell they invert in, and are the one
      member whose two modes differ by rounding. ``pdf`` and
      ``TorchBicop.rect_mass()`` do not depend on the setting.
  norm_maxiter : int, default=25
      Cap on the passes that rescale ``values`` until both margins integrate
      to 1; they stop as soon as both do, to within ``1e-10``. ``25`` is the
      cap a ``tll`` fit of ``Bicop`` allows. Pass ``0`` for a grid that is
      already normalized, which leaves the values exactly as given.
  is_linear : bool, default=False
      Declare ``grid_points`` uniform on ``[0, 1]``, which makes locating a
      cell ``O(1)`` instead of a search. The fitting factories set it from
      ``controls.grid_type``; pass it yourself only for a grid that really is
      uniform, since it is taken at its word.
  device : torch.device or None, default=None
      Where the tensors live; ``None`` keeps them on the inputs' device.
  dtype : torch.dtype, default=torch.float64
      Precision the grid is held and evaluated in. Agreement with ``Bicop``
      is a ``float64`` statement.

  Raises
  ------
  ValueError
      If one of ``grid_points`` / ``values`` is given without the other, if
      ``values`` is not square or has a negative entry, or if ``grid_points``
      is not strictly increasing with at least two points.

  See Also
  --------
  pyvinecopulib.core.Bicop : Reference pair copula.
  pyvinecopulib.core.BicopBase : The contract this class completes.
  FitControlsTorchBicop : Fit-time controls.
  """

  is_indep: bool
  # Class-level hints so the prefix tables registered in __init__ are
  # statically typed as Tensors instead of nn.Module. They are set together:
  # all three are present or all three are None.
  _sy: Tensor | None
  _sx: Tensor | None
  _prefix: Tensor | None
  #: Declares the grid and cache internals ``TorchVinecop``'s batched
  #: cascades read, so a vine of these pairs can take the stacked path.
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
    """Lift a fitted ``Bicop`` into a ``TorchBicop``.

    The density values are taken as fitted, so the lifted copula evaluates
    the same grid as ``cop`` -- in PyTorch, on a device, with autograd.

    Parameters
    ----------
    cop : Bicop
        A fitted pair copula of the ``tll`` family at ``rotation=0``, the one
        shape a ``TorchBicop`` represents.
    cache_integrals : bool, default=True
        As on ``TorchBicop``.
    device : torch.device or None, default=None
        As on ``TorchBicop``.
    dtype : torch.dtype, default=torch.float64
        As on ``TorchBicop``.

    Returns
    -------
    TorchBicop
        A pair copula mirroring ``cop``.

    Raises
    ------
    ValueError
        If ``cop`` is not of the ``tll`` family, is rotated, or does not carry
        a square grid of density values.
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

    What comes back is ``P`` ordinary ``TorchBicop`` objects; only the fit is
    batched. Nothing here reads the pairs as related -- the leading axis is
    ``P`` independent pairs on shared rows -- so one vine's tree level and
    several vines' levels stack alike.

    Parameters
    ----------
    u : Tensor, shape (P, n, 2), dtype float
        Copula-scale observations, one ``(n, 2)`` block per pair. Continuous
        arguments only; an edge with atoms is fitted on its own through
        ``TorchBicop.from_data()``.
    controls : FitControlsTorchBicop or None, default=None
        Fit controls, shared by every pair. ``None`` is
        ``FitControlsTorchBicop()``.
    cache_integrals : bool, default=True
        As on ``TorchBicop``, for each pair returned.
    device : torch.device or None, default=None
        As on ``TorchBicop``.
    dtype : torch.dtype, default=torch.float64
        As on ``TorchBicop``.

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
    lane's bandwidth search freezes as it converges, so the iterations a pair
    takes are the ones its own data earns.

    It is affected by *how many*, in the last bits. Torch selects elementwise
    kernels by element count, so stacking ``P`` pairs agrees with fitting them
    one by one to floating point rather than bit for bit, on every device --
    around ``1e-15`` on a grid value. Fit the pair alone where that matters.

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
    /,
    controls: Optional[FitControlsTorchBicop] = None,
    var_types: Optional[list[str]] = None,
    *,
    cache_integrals: bool = True,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchBicop":
    """Fit a pair copula on pseudo-observations, in pure PyTorch.

    Reproduces the ``tll`` fit of ``Bicop`` to floating-point tolerance, with
    ``controls.method`` naming the fitter (``"tll"`` is the one shipped).

    Parameters
    ----------
    u : ndarray or Tensor, shape (n, 2) or (n, 4), dtype float
        Pseudo-observations ``[u1, u2]``, or ``[u1, u2, u1^-, u2^-]`` when
        ``var_types`` marks an argument discrete.
    controls : FitControlsTorchBicop or None, default=None
        Fit controls. ``None`` is ``FitControlsTorchBicop()``: a TLL fit on a
        30-point normal-spaced grid.
    var_types : list of str or None, default=None
        The two arguments' types, ``"c"`` or ``"d"``; ``None`` means both
        continuous, and a ``"d"`` is what asks for the four-column layout.
        Either way the fitted grid is a continuous density -- the
        mixed-discrete surface an atom needs comes from
        :class:`~pyvinecopulib.core.DiscretePair`.
    cache_integrals : bool, default=True
        As on ``TorchBicop``.
    device : torch.device or None, default=None
        As on ``TorchBicop``.
    dtype : torch.dtype, default=torch.float64
        As on ``TorchBicop``.

    Returns
    -------
    TorchBicop
        The fitted pair copula.

    Raises
    ------
    ValueError
        If ``u``'s column count does not match ``var_types``, or if
        ``controls`` asks for a grid of fewer than two points or a
        non-positive bandwidth multiplier.

    See Also
    --------
    TorchBicop.from_data_batched : Fit ``P`` pairs in one call.
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

  def fit(
    self,
    u,
    /,
    controls: Optional[FitControlsTorchBicop] = None,
    var_types: Optional[list[str]] = None,
  ) -> "TorchBicop":
    """Refit this pair copula's density grid, in place.

    The in-place counterpart of :meth:`from_data`, keeping this module's
    identity: its device, its dtype and its cache mode are what the refitted
    grid gets, and any reference held to it -- a vine's stored pair, an
    optimizer's parameter list -- sees the new density.

    Parameters
    ----------
    u : Tensor, shape (n, 2) or (n, 4)
        Pseudo-observations, with the two left-limit columns when
        ``var_types`` marks an argument discrete.
    controls : FitControlsTorchBicop, or None, optional
        Fit configuration; defaults are used when ``None``.
    var_types : list of str, or None, optional
        The two variable types of the edge this pair sits on. ``None`` means
        both are continuous.

    Returns
    -------
    TorchBicop
        ``self``, so the call chains.

    See Also
    --------
    TorchBicop.from_data : Construct and fit in one call.
    """
    values = self.interp_grid.values
    fitted = type(self).from_data(
      u,
      controls,
      var_types,
      cache_integrals=self._cache_integrals,
      device=values.device,
      dtype=values.dtype,
    )
    self._adopt(fitted)
    return self

  def _adopt(self, other: "TorchBicop") -> None:
    """Take over another pair copula's grid and caches.

    Parameters
    ----------
    other : TorchBicop
        The pair copula whose density to adopt.

    Returns
    -------
    None
    """
    self.is_indep = other.is_indep
    self.interp_grid = other.interp_grid
    self._cache_integrals = other._cache_integrals
    for name in ("_sy", "_sx", "_prefix"):
      table = getattr(other, name)
      # A buffer cannot be reassigned by attribute once registered, so the
      # registry entry itself is replaced.
      self._buffers.pop(name, None)
      if table is None:
        setattr(self, name, None)
      else:
        self.__dict__.pop(name, None)
        self.register_buffer(name, table)

  def flip(self) -> "TorchBicop":
    """Return the copula with its two arguments swapped (``c'(u,v)=c(v,u)``).

    Returns
    -------
    TorchBicop
        A new argument-swapped copula; this one is left unchanged.
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
    """Evaluate the copula density ``c(u1, u2)``.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Pseudo-observations, clamped strictly inside the unit square before
        evaluation (``1e-10`` from each end in ``float64``); a ``NaN`` comes
        back as ``NaN``.
    x : Tensor, shape (n, k), or None, optional
        Unused: a ``TorchBicop`` is unconditional. Accepted so the class
        satisfies ``BicopLike``.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Density values, floored at ``1e-20``.
    """
    u = self._prep(u)
    if self.is_indep:
      return torch.ones(u.shape[0], dtype=u.dtype, device=u.device)
    return self.interp_grid.interpolate(u).clamp_min(1e-20)

  def cdf(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Evaluate the copula distribution function.

    .. math::

       C(u_1, u_2) = \\int_0^{u_1} \\int_0^{u_2} c(s, t)\\, ds\\, dt.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Pseudo-observations, clamped as in ``TorchBicop.pdf()``.
    x : Tensor, shape (n, k), or None, optional
        Unused: a ``TorchBicop`` is unconditional. Accepted so the class
        satisfies ``BicopLike``.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Distribution values, clamped strictly inside ``(0, 1)``.

    Notes
    -----
    Each grid line is rescaled by its own total, so ``C(1, u2) = u2`` holds
    exactly.
    """
    u = self._prep(u)
    if self.is_indep:
      return _trim(u[:, 0] * u[:, 1])
    if self._sy is not None:
      sy, sx, pref = self._tables()
      return self.interp_grid.cdf_cached(u, sy, sx, pref)
    return self.interp_grid.integrate_2d(u)

  def rect_mass(self, a1: Tensor, b1: Tensor, a2: Tensor, b2: Tensor) -> Tensor:
    """Probability of the rectangle ``(a1, b1] x (a2, b2]``.

    The value a four-corner difference of ``TorchBicop.cdf()`` defines,
    arranged so that almost none of it cancels: differencing amplifies an
    absolute error by ``~4 / (w1 w2)`` in the rectangle's widths, where this
    route amplifies by ``1 / w2`` alone. Available in both cache modes -- it
    reads the density grid, not the prefix tables -- and cancellation-free
    only because ``values`` is nonnegative.

    Parameters
    ----------
    a1, b1, a2, b2 : Tensor, shape (n,), dtype float
        Rectangle bounds per query. An empty or inverted interval gives zero.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Rectangle probabilities. The independence copula gives
        ``(b1 - a1) * (b2 - a2)``.

    Notes
    -----
    A probability, not the density grid's own mass over the rectangle: the two
    differ by the rescaling ``TorchBicop.cdf()`` applies.
    :class:`~pyvinecopulib.core.DiscretePair` deliberately leaves this
    accuracy on the table and differences ``cdf`` instead, which is what keeps
    a discrete torch vine in step with ``Vinecop``.
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

    The stored copies are detached buffers, and differentiating through those
    would silently drop most of the gradient in ``values`` rather than all of
    it. They are therefore recomputed whenever grad is live: three cumulative
    sums over ``(m, m)``, and never on the no-gradient path.
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
    """Evaluate the first h-function.

    .. math::

       h_1(u_1, u_2) = \\mathbb{P}(U_2 \\le u_2 \\mid U_1 = u_1).

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Pseudo-observations; column 0 is the conditioning value.
    x : Tensor, shape (n, k), or None, optional
        Unused: a ``TorchBicop`` is unconditional. Accepted so the class
        satisfies ``BicopLike``.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Conditional distribution values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return _trim(u[:, 1])
    return self._hfunc_raw(u, 1).clamp(0.0, 1.0)

  def hfunc2(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Evaluate the second h-function.

    .. math::

       h_2(u_1, u_2) = \\mathbb{P}(U_1 \\le u_1 \\mid U_2 = u_2).

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Pseudo-observations; column 1 is the conditioning value.
    x : Tensor, shape (n, k), or None, optional
        Unused: a ``TorchBicop`` is unconditional. Accepted so the class
        satisfies ``BicopLike``.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Conditional distribution values in ``[0, 1]``.
    """
    u = self._prep(u)
    if self.is_indep:
      return _trim(u[:, 0])
    return self._hfunc_raw(u, 2).clamp(0.0, 1.0)

  # --------------------------------------------------------------------- #
  # Inverse h-functions (closed-form conditional quantiles).               #
  # --------------------------------------------------------------------- #

  def _hinv_raw(self, u: Tensor, cond_var: int) -> Tensor:
    """Shared inverse-h-function body for ``hinv1`` / ``hinv2``.

    ``cond_var=1`` solves ``H1(u1, x) = p`` for ``x`` (free column 1,
    ``u = [u1, p]``); ``cond_var=2`` solves ``H2(x, u2) = p`` (free column 0,
    ``u = [p, u2]``). Both cache modes run the same closed-form inversion of
    the piecewise-quadratic conditional distribution function, as
    ``inverse_integrate_1d`` implements it (vinecopulib#691).
    """
    # Always the closed-form inversion: unlike `cdf` and the h-functions there
    # is no O(1) exact reconstruction to cache here, because locating the
    # bracketing cell needs the conditional cumulative along the whole free
    # axis, which is O(m) to assemble whatever is cached. That cumulative is
    # exactly what a prefix table holds, so pass it when there is one: the
    # same quantity, a gather instead of a trapezoid and a scan.
    cum = None if self._sy is None else self._tables()[cond_var - 1]
    return self.interp_grid.inverse_integrate_1d(u, cond_var, cum).clamp(
      0.0, 1.0
    )

  def hinv1(self, u: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Invert ``TorchBicop.hfunc1()`` in its second argument.

    Given ``u = [u1, p]``, returns the ``u2`` at which
    ``h_1(u_1, u_2) = p``.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Column 0 is ``u1``; column 1 is the target probability ``p``.
    x : Tensor, shape (n, k), or None, optional
        Unused: a ``TorchBicop`` is unconditional. Accepted so the class
        satisfies ``BicopLike``.

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
    """Invert ``TorchBicop.hfunc2()`` in its first argument.

    Given ``u = [p, u2]``, returns the ``u1`` at which
    ``h_2(u_1, u_2) = p``.

    Parameters
    ----------
    u : Tensor, shape (n, 2), dtype float
        Column 0 is the target probability ``p``; column 1 is ``u2``.
    x : Tensor, shape (n, k), or None, optional
        Unused: a ``TorchBicop`` is unconditional. Accepted so the class
        satisfies ``BicopLike``.

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
    """Draw ``n`` samples from the pair copula.

    Parameters
    ----------
    n : int, default=100
        Number of samples to draw; must be positive.
    x : Tensor, shape (n, k), or None, optional
        Unused: a ``TorchBicop`` is unconditional. Accepted so the class
        satisfies ``BicopLike``.
    qrng : bool, default=False
        Draw the base uniforms from a scrambled Sobol sequence instead of
        pseudo-random ones.
    seeds : list of int or None, optional
        The first entry seeds the draw: the Sobol scramble when ``qrng=True``,
        and otherwise a generator private to this call, so seeding one pair
        copula leaves every other consumer of PyTorch's RNG alone. ``None``
        draws from the default generator without reseeding it.

    Returns
    -------
    Tensor, shape (n, 2), dtype float
        Samples in the unit square.

    Raises
    ------
    ValueError
        If ``n`` is not positive.

    Notes
    -----
    Drawn under ``torch.no_grad``, so the samples carry no gradient.
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
