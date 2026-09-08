"""A one-dimensional kernel density margin evaluated in pure PyTorch.

Home of ``TorchKde1d``, the torch marginal estimator, and the only torch
margin that serves discrete and zero-inflated variables.

Fitting delegates to ``Kde1d``; every evaluation runs on tensors, on device,
under autograd. The split is not a compromise -- ``grid_points``, ``values``,
the variable type, ``prob0`` and the declared bounds are all that
``Kde1d``'s ``pdf`` / ``cdf`` / ``icdf`` read, the bounds among them because
for a discrete variable they *are* the integer support, so a lifted grid is a
complete model rather than an approximation of one. Bandwidth selection and
the local-likelihood fit stay with ``Kde1d``, which is why a fit here takes no
controls.

The interpolation every evaluation runs on lives in ``_kde1d_interp``, a port
of kde1d's own interpolation grid whose contract is fidelity to it rather than
improvement on it.
"""

from __future__ import annotations

import math
from typing import Any, Optional, cast

import torch
from torch import Tensor

from ..core import ControlsLike, Kde1d, MarginBase
from ..core._validation import (
  reject_array_controls,
  reject_covariates,
  validate_univariate,
  validate_weights,
)
from . import _kde1d_interp as interp


def _bound(value: Optional[float], unbounded: float) -> float:
  """One end of a support, with ``Kde1d``'s ``nan`` normalized.

  Parameters
  ----------
  value : float, or None, optional
      A bound as ``Kde1d`` reports it: a number, ``None``, or ``nan``.
  unbounded : float
      What an unset bound means at this end.

  Returns
  -------
  float
      A finite bound, or ``unbounded``.
  """
  if value is None:
    return unbounded
  out = float(value)
  return unbounded if out != out else out


#: ``Kde1d``'s spellings of the variable type, and the contract's.
_VAR_TYPE_OF = {"continuous": "c", "discrete": "d", "zero-inflated": "zi"}


class TorchKde1d(MarginBase[Tensor], torch.nn.Module):
  """A kernel density margin on tensors, fitted by ``Kde1d``.

  Evaluates ``pdf`` / ``cdf`` / ``icdf`` on tensors and inherits the rest of
  the ``MarginLike`` surface -- ``logpdf`` / ``cdf_left`` / ``loglik`` /
  ``sample`` -- from ``MarginBase``. It is also an ``nn.Module``, so
  ``.to(device)``, ``state_dict`` and autograd reach the whole margin. Unlike
  ``TorchDistributionMargin`` it serves **discrete** and **zero-inflated** variables as
  well as continuous ones, which is what lets a torch vine distribution be
  fitted to data with atoms at all; ``Vinedist``, evaluating on NumPy, refuses
  a torch margin, so this one belongs to ``TorchVinedist``.

  The variable type, the bounds and the bandwidth are named at construction.
  Everything else follows from that:

  - ``fit(y, weights=...)`` estimates the density from one column and returns
    ``self``. The inherited ``select`` is that same fit: a kernel density has
    no family to choose.
  - ``from_data(y, ...)`` constructs a margin with the defaults and fits it,
    for the case where the class is the whole specification.
  - ``from_kde1d(kde)`` lifts a ``Kde1d`` that is already fitted, exactly.
  - ``from_grid(grid_points, values)`` installs a density that no fit here
    produced.

  ``grid_points``, ``values`` and ``prob0`` are registered as buffers rather
  than parameters: the density is fitted, not learned. A caller who wants to
  optimize it calls ``values.requires_grad_(True)`` -- the same opt-in
  ``TorchTllBicop`` uses for its grid.

  Two of ``Kde1d``'s attribute names cannot be reused here, because the base
  classes already own them: ``type`` is ``nn.Module``'s legacy dtype cast, and
  ``loglik`` is the contract's *method*. The constructor argument stays
  ``type=``; read it back as :attr:`kde_type`, and read the fitted
  log-likelihood as ``loglik()``.

  Parameters
  ----------
  xmin : float, or None, optional
      Lower bound of the support, or ``None`` for unbounded. What a bound
      means depends on the variable type: for a discrete variable it is the
      smallest integer the variable can take. See the
      ``concepts-kde-margins`` section of the concepts page.
  xmax : float, or None, optional
      Upper bound of the support, or ``None`` for unbounded; read as ``xmin``
      is, so the largest integer for a discrete variable.
  type : {"continuous", "discrete", "zero-inflated"}, default="continuous"
      The variable type, spelled as ``Kde1d`` reports it. The hyphen is not
      optional here.
  multiplier : float, default=1.0
      Bandwidth multiplier: the bandwidth used is ``bandwidth * multiplier``.
  bandwidth : float, or None, optional
      Fixed bandwidth, or ``None`` to select one at every fit.
  degree : int, default=2
      Local-polynomial degree -- ``0``, ``1`` or ``2``, for a log-constant,
      log-linear or log-quadratic fit.
  grid_size : int, default=400
      Number of interpolation grid points.
  boundary_repair : bool, default=True
      Whether a finite bound may be fitted with a dedicated boundary
      estimator instead of the transformed bulk fit; eligibility rather than
      a guarantee, and no effect when neither bound is set. Carried through to
      the fit, and preserved when a fitted estimator is lifted.
  device : torch.device, or None, optional
      Where the buffers live.
  dtype : torch.dtype, default=torch.float64
      Buffer precision. ``float64``, since the copula scale is a distribution
      function and ``float32`` costs three digits of it.

  Raises
  ------
  ValueError
      If ``type`` is none of the three spellings above.

  See Also
  --------
  pyvinecopulib.core.Kde1d : The estimator behind the fit, and the reference.
  pyvinecopulib.torch.TorchDistributionMargin : Parametric torch families, continuous only.
  pyvinecopulib.torch.TorchVinedist : The distribution these margins compose.

  Notes
  -----
  Evaluation reproduces ``Kde1d`` rather than improving on it, including where
  ``Kde1d`` is quirky: the unnormalized integral carries no Gaussian-tail mass
  even though the density beyond the grid does. It is pinned by a parity test;
  a divergence there is a defect in this class, not a fix.

  That parity is an equality everywhere but the quantile of a continuous or
  zero-inflated margin, which is an iteration whose last bits follow the
  instruction set ``Kde1d`` was built for -- rebuilding kde1d with
  ``-march=native`` and nothing else moves it 19 ULPs, so no port can equal
  every build of it. The tolerance is ``_QUANTILE_RTOL`` in
  ``tests/test_torch_kde1d.py``.
  """

  supports_weights: bool = True
  #: The bandwidth, bounds and variable type are named at construction, so
  #: `fit` reads no controls. Declaring it is what makes a `family_set` a
  #: refusal rather than a kernel density fitted in silence.
  supports_controls: bool = False
  #: The variable type in ``Kde1d``'s spelling; read back as
  #: :attr:`kde_type`, since ``type`` is ``nn.Module``'s dtype cast.
  _type: str

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
    boundary_repair: bool = True,
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
    self._bandwidth_spec = bandwidth
    self.degree = degree
    self.grid_size = grid_size
    self.boundary_repair = boundary_repair
    self._loglik: Optional[float] = None
    self.edf: Optional[float] = None
    self._selected_bandwidth: Optional[float] = None
    self._dtype = dtype
    self._device = device
    self.register_buffer(
      "grid_points", torch.empty(0, dtype=dtype, device=device)
    )
    self.register_buffer("values", torch.empty(0, dtype=dtype, device=device))
    self.register_buffer("prob0", torch.zeros((), dtype=dtype, device=device))

  def _load_from_state_dict(
    self,
    state_dict: dict[str, Any],
    prefix: str,
    local_metadata: dict[str, Any],
    strict: bool,
    missing_keys: list[str],
    unexpected_keys: list[str],
    error_msgs: list[str],
  ) -> None:
    """Resize the buffers before loading, since a fresh module has none.

    ``grid_points`` and ``values`` start empty -- their length is a property of
    the fit, not of the hyperparameters -- and ``load_state_dict`` checks shapes
    against what is registered. Adopting the incoming shapes first is what makes
    a fitted state loadable into a freshly constructed margin.

    Parameters
    ----------
    state_dict : dict
        The state being loaded.
    prefix : str
        Key prefix for this module.
    local_metadata : dict
        Version metadata, unused.
    strict : bool
        Whether to require an exact key match.
    missing_keys, unexpected_keys, error_msgs : list of str
        Accumulators ``nn.Module`` passes down.

    Returns
    -------
    None
    """
    for name in ("grid_points", "values", "prob0"):
      incoming = state_dict.get(prefix + name)
      current = getattr(self, name, None)
      if incoming is not None and current is not None:
        if incoming.shape != current.shape:
          setattr(
            self,
            name,
            torch.empty(
              incoming.shape, dtype=current.dtype, device=current.device
            ),
          )
    super()._load_from_state_dict(
      state_dict,
      prefix,
      local_metadata,
      strict,
      missing_keys,
      unexpected_keys,
      error_msgs,
    )

  # --- construction --------------------------------------------------------- #

  def fit(
    self,
    y: Tensor,
    /,
    controls: Optional[ControlsLike] = None,
    *,
    x: Optional[Tensor] = None,
    weights: Optional[Tensor] = None,
  ) -> "TorchKde1d":
    """Estimate the density from one column of data, in place.

    The fit is ``Kde1d``'s, and the grid it settles on is what this margin
    evaluates from then on -- carried onto the device and dtype named at
    construction.

    Parameters
    ----------
    y : Tensor, shape (n,)
        Observations on the original scale.
    controls : ControlsLike, or None, optional
        Unused; the bandwidth, bounds and variable type are named at
        construction, so a margin fitted differently is constructed
        differently. Accepted because every margin's fit takes one.
    x : Tensor, or None, optional
        Not supported; a kernel density reads no covariates, so passing them
        raises rather than fitting an unconditional margin silently.
    weights : Tensor, shape (n,), or None, optional
        Observation weights, one per observation.

    Returns
    -------
    TorchKde1d
        ``self``, so the call chains.

    Raises
    ------
    ValueError
        If ``y`` or ``weights`` is not one-dimensional, if the two have
        different lengths, or if ``x`` was supplied.

    See Also
    --------
    from_data : Construct and fit in one call.
    from_kde1d : Lift a ``Kde1d`` that is already fitted.
    """
    reject_covariates(self, x)
    # `kde.fit(x, w)` is the compiled `Kde1d`'s spelling, and here it would
    # bind the weights to `controls` and fit unweighted.
    reject_array_controls(self, controls)
    kde = Kde1d(
      xmin=self.xmin,
      xmax=self.xmax,
      type=self._type,
      multiplier=self.multiplier,
      # The construction spec, not `self.bandwidth`: a previous fit overwrote
      # that with the bandwidth it selected, and passing a value back in pins
      # it, so a refit on different data reused the first fit's bandwidth.
      bandwidth=self._bandwidth_spec,
      degree=self.degree,
      grid_size=self.grid_size,
      boundary_repair=self.boundary_repair,
    )
    y_tensor = validate_univariate(torch.as_tensor(y))
    # The shared validators, not a local pair of shape checks: their dtype,
    # finiteness, nonnegativity and positive-sum rules are what keep a weight
    # array out of `Kde1d`'s bandwidth selection, which divides by the sum.
    weight_tensor = validate_weights(weights, y_tensor)
    data = y_tensor.detach().cpu().numpy()
    if weight_tensor is None:
      kde.fit(data)
    else:
      kde.fit(data, weight_tensor.detach().cpu().numpy())
    fitted = self._adopt(kde)
    # The retained rows, not the input length: `Kde1d` drops a NaN observation
    # and a NaN or zero weight, and the log-likelihood `_adopt` just read is
    # over what is left, so that is what the criteria penalize against.
    kept = ~torch.isnan(y_tensor)
    if weight_tensor is not None:
      kept = kept & ~torch.isnan(weight_tensor) & (weight_tensor > 0)
    fitted._nobs = int(kept.sum())
    return fitted

  @classmethod
  def from_kde1d(
    cls,
    kde: Kde1d,
    *,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchKde1d":
    """Lift a fitted ``Kde1d`` onto tensors.

    An exact transfer rather than a second fit: the same grid, bounds,
    variable type and diagnostics, evaluated in torch. A later refit on other
    data selects a bandwidth again, exactly as the estimator handed over here
    would have.

    Parameters
    ----------
    kde : Kde1d
        A fitted estimator.
    device : torch.device, or None, optional
        Where the buffers live.
    dtype : torch.dtype, default=torch.float64
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
      # `bandwidth_spec` is what the compiled object was *asked* for; its
      # `bandwidth` is what it selected. Adopting the latter as the spec would
      # stop a lifted margin from ever re-selecting on a refit.
      bandwidth=kde.bandwidth_spec,
      degree=kde.degree,
      grid_size=kde.grid_size,
      boundary_repair=kde.boundary_repair,
      device=device,
      dtype=dtype,
    )
    return out._adopt(kde)

  @classmethod
  def from_grid(
    cls,
    grid_points: Tensor,
    values: Tensor,
    *,
    prob0: float = 0.0,
    **kwargs: Any,  # noqa: ANN401 - forwarded to `__init__`
  ) -> "TorchKde1d":
    """Build a margin directly from a grid, with no fit.

    The fit-free injection point, mirroring ``Kde1d.from_grid``: a density
    obtained some other way -- optimized, transferred, hand-built -- becomes a
    margin. Nothing was estimated, so ``loglik()`` and ``n_parameters`` have no
    fitted value to report.

    Parameters
    ----------
    grid_points : Tensor, shape (m,)
        Ascending grid.
    values : Tensor, shape (m,)
        Density values on the grid.
    prob0 : float, default=0.0
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
        If the two tensors disagree in length, or the grid is not strictly
        ascending.
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

  def _adopt(self, kde: Kde1d) -> "TorchKde1d":
    """Copy a fitted ``Kde1d``'s state onto this module's buffers."""
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
    self.bandwidth = float(kde.bandwidth)
    self._selected_bandwidth = float(kde.bandwidth)
    self.degree = kde.degree
    self.multiplier = kde.multiplier
    self.boundary_repair = kde.boundary_repair
    self._loglik = float(kde.loglik())
    self.edf = float(kde.edf)
    return self

  def to_json(self) -> dict[str, Any]:
    """Return this margin's JSON payload.

    Returns
    -------
    dict
        A JSON-serializable mapping that
        :func:`~pyvinecopulib.core.margin_from_json` reads back.

    Raises
    ------
    ValueError
        If the density has not been fitted, so there is no grid to store.
    """
    if not self.is_fitted:
      raise ValueError(
        "an unfitted TorchKde1d cannot be serialized; call fit(y) first"
      )
    return {
      "kind": "TorchKde1d",
      "state": self.get_extra_state(),
      "grid_points": [float(v) for v in self.grid_points.tolist()],
      "values": [float(v) for v in self.values.tolist()],
      "prob0": float(self.prob0),
    }

  @classmethod
  def from_json_payload(cls, payload: dict[str, Any]) -> "TorchKde1d":
    """Rebuild a margin from the payload :meth:`to_json` produced.

    Parameters
    ----------
    payload : dict
        The mapping :meth:`to_json` returned.

    Returns
    -------
    TorchKde1d
        The reconstructed margin, on the default device and dtype.
    """
    state = dict(payload["state"])
    out = cls.from_grid(
      torch.as_tensor(payload["grid_points"], dtype=torch.float64),
      torch.as_tensor(payload["values"], dtype=torch.float64),
      prob0=float(payload.get("prob0", 0.0)),
      xmin=state.get("xmin"),
      xmax=state.get("xmax"),
      type=str(state.get("type", "continuous")),
      multiplier=float(state.get("multiplier", 1.0)),
      degree=int(state.get("degree", 2)),
      boundary_repair=bool(state.get("boundary_repair", True)),
    )
    # The diagnostics `get_extra_state` carries and `from_grid` cannot take:
    # a grid supplied directly has no fit behind it, so these are restored
    # rather than recomputed.
    out.set_extra_state(state)
    return out

  def get_extra_state(self) -> dict[str, Any]:
    """Return non-tensor fitted state for ``state_dict`` round-trips.

    Returns
    -------
    dict
        Fitted configuration and diagnostics not stored as tensors.
    """
    return {
      "version": 1,
      "xmin": self.xmin,
      "xmax": self.xmax,
      "type": self._type,
      "multiplier": self.multiplier,
      "bandwidth": self.bandwidth,
      "bandwidth_spec": self._bandwidth_spec,
      "selected_bandwidth": self._selected_bandwidth,
      "degree": self.degree,
      "grid_size": self.grid_size,
      "boundary_repair": self.boundary_repair,
      "loglik": self._loglik,
      "edf": self.edf,
      "nobs": self._nobs,
    }

  def set_extra_state(self, state: object) -> None:
    """Restore non-tensor fitted state saved by :meth:`get_extra_state`.

    Parameters
    ----------
    state : object
        State returned by :meth:`get_extra_state`; anything else is refused.

    Raises
    ------
    RuntimeError
        If the state was not written by this version of the class.
    """
    if not isinstance(state, dict) or state.get("version") != 1:
      raise RuntimeError("unsupported TorchKde1d state-dict version")
    self.xmin = state["xmin"]
    self.xmax = state["xmax"]
    self._type = state["type"]
    self.multiplier = state["multiplier"]
    self.bandwidth = state["bandwidth"]
    self._bandwidth_spec = state["bandwidth_spec"]
    self._selected_bandwidth = state["selected_bandwidth"]
    self.degree = state["degree"]
    self.grid_size = state["grid_size"]
    self.boundary_repair = state["boundary_repair"]
    self._loglik = state["loglik"]
    self.edf = state["edf"]
    # The payload is an opaque ``object``, and the retained sample size is
    # declared on ``MarginBase`` as what the `nobs` property answers.
    self._nobs = cast("Optional[int]", state["nobs"])

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
        ``"continuous"``, ``"discrete"`` or ``"zero-inflated"``, hyphenated as
        ``Kde1d`` spells it.
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
        ``(xmin, xmax)``, with an unset bound reported as infinite. ``Kde1d``
        spells an unset bound ``nan``, which is neither ordered nor equal to
        itself, so it is normalized here: every evaluation compares against
        this pair, and for a discrete margin it is the integer support the
        masses sit on -- the fitted grid runs half a unit wider at each end.
    """
    return (_bound(self.xmin, float("-inf")), _bound(self.xmax, float("inf")))

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
  def nobs(self) -> Optional[int]:
    """Number of observations the fit **retained**.

    ``Kde1d`` drops a NaN observation and a NaN or zero weight, so this falls
    below the input length whenever one was dropped.

    Returns
    -------
    int or None
        The retained sample size, or ``None`` before a fit.
    """
    return self._nobs

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

  @property
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

  def _as_tensor(self, values: Any) -> Tensor:  # noqa: ANN401 - as_tensor input
    """Coerce a query onto the buffers' dtype and device, flattened."""
    return torch.as_tensor(
      values, dtype=self.grid_points.dtype, device=self.grid_points.device
    ).reshape(-1)

  def _levels(self) -> Tensor:
    """The integer support the discrete branches live on.

    A declared bound *is* the support endpoint; the grid runs half a unit
    wider, since the jittered observations fill the boundary cells. Where no
    bound was declared the grid is all there is to go on, so the support is
    read off it. Derived rather than stored, as ``Kde1d`` derives it, so a grid
    that moves takes its support with it -- and matched by a parity test
    against the levels ``Kde1d`` itself gives mass to, since the two are
    separate copies of one rule.
    """
    g = self.grid_points
    lo_b, hi_b = self.support
    lo = (
      torch.floor(g[0])
      if lo_b == -math.inf
      else torch.as_tensor(lo_b, dtype=g.dtype, device=g.device)
    )
    hi = (
      torch.maximum(lo, torch.ceil(g[-1]))
      if hi_b == math.inf
      else torch.as_tensor(hi_b, dtype=g.dtype, device=g.device)
    )
    n = int(round(float(hi - lo))) + 1
    return lo + torch.arange(n, dtype=g.dtype, device=g.device)

  def _pdf_continuous(self, y: Tensor) -> Tensor:
    raw = interp.interpolate(self.grid_points, self.values, y)
    out = torch.clamp(raw, min=0.0)
    # A declared bound is a hard edge of the support, not merely where the grid
    # happens to stop -- on a discrete margin the grid overhangs it by half a
    # cell, so the two are no longer the same question.
    lo, hi = self.support
    if lo != -math.inf:
      out = torch.where(y < lo, torch.zeros_like(out), out)
    if hi != math.inf:
      out = torch.where(y > hi, torch.zeros_like(out), out)
    return out

  def _pdf_discrete(self, y: Tensor) -> Tensor:
    lvs = self._levels()
    keep = (y >= lvs[0]) & (y <= lvs[-1]) & (y == torch.round(y))
    # The jitter density between the levels, and outside the support, is not
    # probability mass of the discrete model: only the ordinates at integer
    # centers are, so the normalizer is the density at the levels.
    norm = self._pdf_continuous(lvs).sum()
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
      torch.where(y >= hi, torch.ones_like(y), f_cum[idx].clamp(0.0, 1.0)),
    )
    return torch.where(torch.isnan(y), y, out)

  def pdf(self, y: Tensor, /, *, x: Optional[Tensor] = None) -> Tensor:
    """Density with respect to this margin's own reference measure.

    Parameters
    ----------
    y : Tensor, shape (n,)
        Evaluation points.
    x : Tensor, or None, optional
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
    x : Tensor, or None, optional
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
    """Invert as ``Kde1d`` does -- cell by cell -- then reattach a gradient.

    The forward value is the iteration's, bit for bit. The gradient comes from
    the implicit function theorem -- ``dq/dtheta = -(dF/dtheta) / f(q)``, and
    ``dq/dp = 1 / f(q)`` -- which a first Newton step expresses exactly, so
    differentiating the correction while returning the iterate gives both.

    The residual is written in units of mass rather than probability so that
    the total mass, which is itself a function of ``values``, carries its share
    of ``dq/dtheta``.

    The correction is skipped only when no gradient is wanted at all. Gating it
    on ``values.requires_grad`` alone would kill ``dq/dp`` for a fitted, fixed
    grid, which is the common case: the density is fitted, not learned, and the
    quantile is still a differentiable function of its probability.
    """
    q = interp.invert_integral(self.grid_points, self.values, p)
    if not torch.is_grad_enabled() or not (
      self.values.requires_grad or p.requires_grad
    ):
      return q
    total = interp.total_mass(self.grid_points, self.values)
    residual = interp.integrate(self.grid_points, self.values, q) - p * total
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
    x : Tensor, or None, optional
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

    Notes
    -----
    On a continuous or zero-inflated margin the quantile is differentiable in
    ``p``, and in the grid where a caller opted into that with
    ``values.requires_grad_(True)``; neither costs accuracy, since the value
    returned is the one a gradient-free call gives. A discrete quantile is a
    step function and carries no gradient.
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
      # The lowest level whose cdf has reached `p`, not the lowest that has
      # passed it -- `lower_bound`, as the compiled quantile does. The two
      # differ exactly where `p` lands on a level's cumulative probability,
      # which is every point of an `icdf(cdf(k))` round trip.
      idx = torch.clamp(
        torch.searchsorted(f_cum.contiguous(), pa.contiguous()),
        max=lvs.numel() - 1,
      )
      return torch.where(torch.isnan(pa), pa, lvs[idx])
    if self._type == "zero-inflated":
      zero = torch.zeros(1, dtype=pa.dtype, device=pa.device)
      if float(self.prob0) >= 1.0:
        # An all-zero column: the mass is the whole distribution and there is
        # no continuous part to invert. `cdf` already answers this way.
        return torch.where(torch.isnan(pa), pa, torch.zeros_like(pa))
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
