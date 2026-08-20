"""A one-dimensional kernel density margin evaluated in pure PyTorch.

Fitting delegates to the compiled :class:`pyvinecopulib.core.Kde1d`; every
evaluation runs on tensors, on device, under autograd. The split is not a
compromise -- ``grid_points``, ``values``, ``type`` and ``prob0`` are the *only*
state the compiled ``pdf`` / ``cdf`` / ``icdf`` read, so a lifted grid is a
complete model rather than an approximation of one. What a torch fitter would
add is bandwidth selection and the local-likelihood fit -- a separate piece of
work, which would attach as a ``method`` on a fit-controls dataclass beside the
existing ones.
"""

from __future__ import annotations

from typing import Any, Optional

import torch
from torch import Tensor

from ..core import Kde1d, MarginBase
from ..core.margin_base import _reject_covariates
from . import _kde1d_interp as interp

#: The compiled spellings of the variable type, and the contract's.
_VAR_TYPE_OF = {"continuous": "c", "discrete": "d", "zero-inflated": "zi"}


class TorchKde1d(MarginBase[Tensor], torch.nn.Module):
  """A kernel density margin on tensors, with the compiled fitter behind it.

  Satisfies the ``MarginLike`` contract, so it drops into
  :class:`pyvinecopulib.core.Vinedist` or
  :class:`pyvinecopulib.torch.TorchVinedist` exactly like any other margin, and
  is an ``nn.Module``, so ``.to(device)``, ``state_dict`` and autograd all work.
  Unlike :class:`pyvinecopulib.torch.TorchMargin` it handles **discrete** and
  **zero-inflated** variables, which is what lets a torch vine distribution be
  fitted to data with atoms at all.

  ``grid_points``, ``values`` and ``prob0`` are registered as buffers rather
  than parameters: the density is fitted, not learned. A caller who wants to
  optimize it calls ``values.requires_grad_(True)`` -- the same opt-in
  ``TorchBicop`` uses for its grid.

  Two of ``Kde1d``'s attribute names cannot be reused here, because the base
  classes already own them: ``type`` is ``nn.Module``'s legacy dtype cast, and
  ``loglik`` is the contract's *method*. The constructor argument stays
  ``type=``; read it back as :attr:`kde_type`, and read the fitted
  log-likelihood as ``loglik()``.

  Parameters
  ----------
  xmin : float or None, optional
      Lower bound of the support, or ``None`` for unbounded.
  xmax : float or None, optional
      Upper bound of the support, or ``None`` for unbounded.
  type : {"continuous", "discrete", "zero-inflated"}, optional
      The variable type, in the compiled spelling. Note the hyphen.
  multiplier : float, optional
      Bandwidth multiplier.
  bandwidth : float or None, optional
      Fixed bandwidth, or ``None`` to select one.
  degree : int, optional
      Local-polynomial degree: ``0``, ``1`` or ``2``.
  grid_size : int, optional
      Number of interpolation grid points.
  device : torch.device or None, optional
      Where the buffers live.
  dtype : torch.dtype, optional
      Buffer precision; ``float64`` by default, since the copula scale is a
      distribution function and ``float32`` costs three digits of it.

  See Also
  --------
  pyvinecopulib.core.Kde1d : The compiled estimator, and the fitter used here.
  pyvinecopulib.torch.TorchMargin : Parametric torch families, continuous only.

  Notes
  -----
  Evaluation reproduces the compiled implementation rather than improving on
  it, including where the C++ is quirky: the interpolant drops by
  ``exp(-0.5)`` exactly at the grid's right end, the unnormalized integral
  carries no Gaussian-tail mass, and the discrete normalization divides by the
  *raw* interpolation rather than the clamped density. Each is pinned by a
  parity test; a divergence there is a defect in this class, not a fix.
  """

  supports_weights: bool = True
  supported_var_types: tuple[str, ...] = ("c", "d", "zi")

  def __init__(
    self,
    *,
    xmin: Optional[float] = None,
    xmax: Optional[float] = None,
    type: str = "continuous",
    multiplier: float = 1.0,
    bandwidth: Optional[float] = None,
    degree: int = 2,
    grid_size: int = 400,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> None:
    # Initialize nn.Module explicitly: this also subclasses MarginBase (a
    # Protocol-derived ABC), whose __init__ chain would otherwise shadow
    # nn.Module's under super().
    torch.nn.Module.__init__(self)
    if type not in _VAR_TYPE_OF:
      raise ValueError(
        f"unknown type={type!r}; expected one of {list(_VAR_TYPE_OF)}"
      )
    self.xmin = xmin
    self.xmax = xmax
    self._type = type
    self.multiplier = multiplier
    self.bandwidth = bandwidth
    self.degree = degree
    self.grid_size = grid_size
    self._loglik: Optional[float] = None
    self.edf: Optional[float] = None
    self._dtype = dtype
    self._device = device
    self.register_buffer(
      "grid_points", torch.empty(0, dtype=dtype, device=device)
    )
    self.register_buffer("values", torch.empty(0, dtype=dtype, device=device))
    self.register_buffer("prob0", torch.zeros((), dtype=dtype, device=device))

  # --- construction --------------------------------------------------------- #

  def fit(
    self,
    y: Tensor,
    /,
    *,
    x: Optional[Tensor] = None,
    weights: Optional[Tensor] = None,
  ) -> "TorchKde1d":
    """Fit the density with the compiled estimator, then lift its grid.

    Parameters
    ----------
    y : Tensor, shape (n,)
        Observations on the original scale.
    x : Tensor or None, optional
        Not supported; a kernel density reads no covariates, so passing them
        raises rather than fitting an unconditional margin silently.
    weights : Tensor, shape (n,), or None, optional
        Observation weights.

    Returns
    -------
    TorchKde1d
        ``self``, so the call chains.
    """
    _reject_covariates(self, x)
    kde = Kde1d(
      xmin=self.xmin,
      xmax=self.xmax,
      type=self._type,
      multiplier=self.multiplier,
      bandwidth=self.bandwidth,
      degree=self.degree,
      grid_size=self.grid_size,
    )
    data = torch.as_tensor(y).detach().reshape(-1).cpu().numpy()
    if weights is None:
      kde.fit(data)
    else:
      kde.fit(data, torch.as_tensor(weights).detach().reshape(-1).cpu().numpy())
    return self._adopt(kde)

  @classmethod
  def from_kde1d(
    cls,
    kde: Any,
    *,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchKde1d":
    """Lift a fitted compiled estimator onto tensors.

    Parameters
    ----------
    kde : Kde1d
        A fitted estimator.
    device : torch.device or None, optional
        Where the buffers live.
    dtype : torch.dtype, optional
        Buffer precision.

    Returns
    -------
    TorchKde1d
        The same density, evaluated in torch.

    Raises
    ------
    ValueError
        If ``kde`` is not fitted.
    """
    if not kde.is_fitted:
      raise ValueError(
        "Kde1d is not fitted; call fit(y) before lifting it onto tensors"
      )
    out = cls(
      xmin=kde.xmin,
      xmax=kde.xmax,
      type=kde.type,
      multiplier=kde.multiplier,
      bandwidth=kde.bandwidth,
      degree=kde.degree,
      grid_size=kde.grid_size,
      device=device,
      dtype=dtype,
    )
    return out._adopt(kde)

  @classmethod
  def from_data(
    cls,
    y: Tensor,
    *,
    weights: Optional[Tensor] = None,
    **kwargs: Any,
  ) -> "TorchKde1d":
    """Construct and fit in one call.

    Parameters
    ----------
    y : Tensor, shape (n,)
        Observations on the original scale.
    weights : Tensor, shape (n,), or None, optional
        Observation weights.
    **kwargs
        Forwarded to the constructor.

    Returns
    -------
    TorchKde1d
        The fitted margin.
    """
    return cls(**kwargs).fit(y, weights=weights)

  @classmethod
  def from_grid(
    cls,
    grid_points: Tensor,
    values: Tensor,
    *,
    prob0: float = 0.0,
    **kwargs: Any,
  ) -> "TorchKde1d":
    """Build a margin directly from a grid, with no fit.

    The fit-free injection point, mirroring ``Kde1d.from_grid``: a density
    obtained some other way -- optimized, transferred, hand-built -- becomes a
    margin without going through the compiled estimator.

    Parameters
    ----------
    grid_points : Tensor, shape (m,)
        Ascending grid.
    values : Tensor, shape (m,)
        Density values on the grid.
    prob0 : float, optional
        Point mass at zero, for a zero-inflated margin.
    **kwargs
        Forwarded to the constructor.

    Returns
    -------
    TorchKde1d
        A margin ready to evaluate.

    Raises
    ------
    ValueError
        If the two tensors disagree in length, or the grid is not ascending.
    """
    out = cls(**kwargs)
    g = torch.as_tensor(grid_points, dtype=out._dtype, device=out._device)
    v = torch.as_tensor(values, dtype=out._dtype, device=out._device)
    if g.ndim != 1 or v.shape != g.shape:
      raise ValueError(
        "grid_points and values must be one-dimensional and the same length; "
        f"got {tuple(g.shape)} and {tuple(v.shape)}"
      )
    if g.numel() < 2 or bool(torch.any(g[1:] <= g[:-1])):
      raise ValueError("grid_points must be strictly ascending")
    out.grid_points = g
    out.values = v
    out.prob0 = torch.as_tensor(
      float(prob0), dtype=out._dtype, device=out._device
    )
    return out

  def _adopt(self, kde: Any) -> "TorchKde1d":
    """Copy a fitted compiled estimator's state onto this module's buffers."""
    ref = self.grid_points
    self.grid_points = torch.as_tensor(
      kde.grid_points, dtype=ref.dtype, device=ref.device
    )
    self.values = torch.as_tensor(
      kde.values, dtype=ref.dtype, device=ref.device
    )
    self.prob0 = torch.as_tensor(
      float(kde.prob0), dtype=ref.dtype, device=ref.device
    )
    self._type = kde.type
    self.xmin = kde.xmin
    self.xmax = kde.xmax
    self.bandwidth = kde.bandwidth
    self.degree = kde.degree
    self.multiplier = kde.multiplier
    self._loglik = float(kde.loglik())
    self.edf = float(kde.edf)
    return self

  # --- declared capabilities ------------------------------------------------ #

  @property
  def kde_type(self) -> str:
    """The variable type in ``Kde1d``'s spelling.

    Not exposed as ``type``: ``nn.Module.type`` is the legacy dtype cast, and
    shadowing it would break ``module.type(torch.float32)`` on this class alone.
    The constructor argument keeps the familiar name.

    Returns
    -------
    str
        ``"continuous"``, ``"discrete"`` or ``"zero-inflated"`` -- hyphenated,
        as the compiled class spells it.
    """
    return self._type

  @property
  def var_type(self) -> str:
    """Variable type in the contract's spelling.

    Returns
    -------
    str
        ``"c"``, ``"d"`` or ``"zi"``.
    """
    return _VAR_TYPE_OF[self._type]

  @property
  def support(self) -> tuple[float, float]:
    """Closed bounds of the support.

    Returns
    -------
    tuple of float
        ``(xmin, xmax)``, with an unset bound reported as infinite.
    """
    lo = float("-inf") if self.xmin is None else float(self.xmin)
    hi = float("inf") if self.xmax is None else float(self.xmax)
    return (lo, hi)

  @property
  def is_fitted(self) -> bool:
    """Whether a density is available to evaluate.

    Returns
    -------
    bool
        ``True`` once a grid has been fitted or supplied.
    """
    return int(self.grid_points.numel()) > 0

  @property
  def n_parameters(self) -> float:
    """Effective degrees of freedom of the fit.

    Returns
    -------
    float
        ``edf``, the trace-of-smoother convention the library uses on both the
        marginal and the copula side. ``nan`` for a grid supplied directly,
        which was not fitted here.
    """
    return float("nan") if self.edf is None else self.edf

  @property
  def family_name(self) -> str:
    """Name of the family, for a report.

    Returns
    -------
    str
        Always ``"kde1d"``.
    """
    return "kde1d"

  def _fitted_loglik(self) -> float:
    """The log-likelihood attained at the fit."""
    if self._loglik is None:
      raise RuntimeError(
        "TorchKde1d has no fitted log-likelihood: its grid was supplied "
        "directly rather than fitted. Pass y to loglik() to evaluate one."
      )
    return self._loglik

  # --- evaluation ----------------------------------------------------------- #

  def _check_fitted(self) -> None:
    if not self.is_fitted:
      raise RuntimeError("TorchKde1d is not fitted; call fit(y)")

  def _as_tensor(self, values: Any) -> Tensor:
    """Coerce a query onto the buffers' dtype and device, flattened."""
    return torch.as_tensor(
      values, dtype=self.grid_points.dtype, device=self.grid_points.device
    ).reshape(-1)

  def _levels(self) -> Tensor:
    """The integer lattice the discrete branches live on.

    Derived from the grid rather than stored, as the C++ does, so a grid that
    moves takes its levels with it.
    """
    g = self.grid_points
    lo = torch.floor(g[0])
    hi = torch.ceil(g[-1])
    n = int(round(float(hi - lo))) + 1
    return lo + torch.arange(n, dtype=g.dtype, device=g.device)

  def _pdf_continuous(self, y: Tensor) -> Tensor:
    raw = interp.interpolate(self.grid_points, self.values, y)
    return torch.clamp(raw, min=0.0)

  def _pdf_discrete(self, y: Tensor) -> Tensor:
    lvs = self._levels()
    keep = (y >= lvs[0]) & (y <= lvs[-1]) & (y == torch.round(y))
    # The denominator is the *raw* interpolation over the levels, tails and the
    # endpoint drop included -- not the clamped density. Using the clamped one
    # here would change every discrete mass in the fourth digit.
    norm = interp.interpolate(self.grid_points, self.values, lvs).sum()
    return self._pdf_continuous(y) * keep / norm

  def _level_cdf(self) -> tuple[Tensor, Tensor]:
    lvs = self._levels()
    return lvs, torch.cumsum(self._pdf_discrete(lvs), dim=0)

  def _cdf_continuous(self, y: Tensor) -> Tensor:
    return interp.integrate(self.grid_points, self.values, y, normalize=True)

  def _cdf_discrete(self, y: Tensor) -> Tensor:
    lvs, f_cum = self._level_cdf()
    lo, hi = lvs[0], lvs[-1]
    idx = torch.clamp((y - lo).to(torch.int64), min=0, max=lvs.numel() - 1)
    out = torch.where(
      y < lo,
      torch.zeros_like(y),
      torch.where(y >= hi, torch.ones_like(y), f_cum[idx]),
    )
    return torch.where(torch.isnan(y), y, out)

  def pdf(self, y: Tensor, /, *, x: Optional[Tensor] = None) -> Tensor:
    """Density with respect to this margin's own reference measure.

    Parameters
    ----------
    y : Tensor, shape (n,)
        Evaluation points.
    x : Tensor or None, optional
        Ignored; a kernel density reads no covariates.

    Returns
    -------
    Tensor, shape (n,)
        A Lebesgue density for a continuous margin, a probability mass on the
        lattice for a discrete one, and whichever applies pointwise for a
        zero-inflated one.
    """
    self._check_fitted()
    ya = self._as_tensor(y)
    if self._type == "discrete":
      return self._pdf_discrete(ya)
    if self._type == "zero-inflated":
      return torch.where(
        ya == 0.0,
        self.prob0.expand_as(ya),
        (1.0 - self.prob0) * self._pdf_continuous(ya),
      )
    return self._pdf_continuous(ya)

  def cdf(self, y: Tensor, /, *, x: Optional[Tensor] = None) -> Tensor:
    """Distribution function.

    Parameters
    ----------
    y : Tensor, shape (n,)
        Evaluation points.
    x : Tensor or None, optional
        Ignored; a kernel density reads no covariates.

    Returns
    -------
    Tensor, shape (n,)
        ``F(y)``.
    """
    self._check_fitted()
    ya = self._as_tensor(y)
    if self._type == "discrete":
      return self._cdf_discrete(ya)
    if self._type == "zero-inflated":
      atom = (ya >= 0.0).to(ya.dtype)
      tail = (
        torch.zeros_like(ya)
        if float(self.prob0) >= 1.0
        else self._cdf_continuous(ya)
      )
      return self.prob0 * atom + (1.0 - self.prob0) * tail
    return self._cdf_continuous(ya)

  def _icdf_continuous(self, p: Tensor) -> Tensor:
    """Bisect as the C++ does, then reattach an exact gradient.

    The forward value is the bisection's, bit for bit. The gradient comes from
    the implicit function theorem -- ``dq/dtheta = -(dF/dtheta) / f(q)`` -- which
    a first Newton step expresses exactly, so differentiating the correction
    while returning the bisection gives both.
    """
    q = interp.invert_integral(self.grid_points, self.values, p)
    if not self.values.requires_grad:
      return q
    residual = interp.integrate(self.grid_points, self.values, q) - p
    density = interp.interpolate(self.grid_points, self.values, q)
    tiny = torch.finfo(self.values.dtype).tiny
    corrected = q - residual / torch.clamp(density, min=tiny)
    return corrected + (q - corrected).detach()

  def icdf(self, p: Tensor, /, *, x: Optional[Tensor] = None) -> Tensor:
    """Quantile function.

    Parameters
    ----------
    p : Tensor, shape (n,)
        Probabilities in ``[0, 1]``.
    x : Tensor or None, optional
        Ignored; a kernel density reads no covariates.

    Returns
    -------
    Tensor, shape (n,)
        ``F^{-1}(p)``, on the lattice for a discrete margin and exactly ``0``
        on the atom of a zero-inflated one.

    Raises
    ------
    ValueError
        If any probability lies outside ``[0, 1]``.
    """
    self._check_fitted()
    pa = self._as_tensor(p)
    finite = pa[~torch.isnan(pa)]
    if finite.numel() and (
      bool(torch.any(finite < 0.0)) or bool(torch.any(finite > 1.0))
    ):
      raise ValueError("probabilities must lie in [0, 1]")
    if self._type == "discrete":
      lvs, f_cum = self._level_cdf()
      idx = torch.clamp(
        torch.searchsorted(f_cum.contiguous(), pa.contiguous(), right=True),
        max=lvs.numel() - 1,
      )
      return torch.where(torch.isnan(pa), pa, lvs[idx])
    if self._type == "zero-inflated":
      zero = torch.zeros(1, dtype=pa.dtype, device=pa.device)
      p0 = self.cdf(zero)[0]
      below = pa <= p0 - self.prob0
      rescaled = torch.where(
        below,
        pa / (1.0 - self.prob0),
        torch.clamp(pa - self.prob0, min=0.0) / (1.0 - self.prob0),
      )
      q = self._icdf_continuous(rescaled)
      on_atom = (pa > p0 - self.prob0) & (pa <= p0)
      return torch.where(on_atom, torch.zeros_like(q), q)
    return self._icdf_continuous(pa)

  # --- sampling ------------------------------------------------------------- #

  def _sample_uniform(self, n: int, seeds: list[int]) -> Tensor:
    """Draw ``n`` uniforms on the buffers' dtype and device."""
    ref = self.grid_points
    generator: Optional[torch.Generator] = None
    if seeds:
      generator = torch.Generator(device=ref.device).manual_seed(int(seeds[0]))
    return torch.rand(
      n, generator=generator, dtype=ref.dtype, device=ref.device
    )

  def __repr__(self) -> str:
    """Return a structural representation of the margin.

    Returns
    -------
    str
        The type, the grid size and the support.
    """
    if not self.is_fitted:
      return f"TorchKde1d(type={self._type!r}, unfitted)"
    lo, hi = self.support
    return (
      f"TorchKde1d(type={self._type!r}, grid_size={self.grid_points.numel()}, "
      f"support=({lo}, {hi}))"
    )
