"""PyTorch vine distribution: a torch vine copula with torch margins.

:class:`~pyvinecopulib.core.VinedistBase` is already array-agnostic, so the joint
density, the Rosenblatt transform and simulation need no porting. What this
subclass adds is the :class:`torch.nn.Module` half: the copula and every margin
are registered children, so one ``.to(device)`` moves the whole distribution,
``state_dict()`` captures all of it, and an optimizer sees the margins'
parameters.

The parts have to be torch parts for that to mean anything, which is why both
the copula and the margins are checked at construction rather than left to fail
numerically later.

See Also
--------
pyvinecopulib.core.VinedistBase : The array-agnostic base.
TorchMargin : The margins this holds.
TorchVinecop : The copula this holds.
"""

from __future__ import annotations

import math

import dataclasses
from itertools import chain
from typing import Any, ClassVar, Optional, Sequence, cast

import torch
from torch import Tensor

from ..core import MarginLike, VinedistBase
from .controls import FitControlsTorchVinecop
from .kde1d import TorchKde1d
from .vinecop import TorchVinecop

__all__ = ["TorchVinedist"]


def _check_margin(margin: Any, name: str) -> None:
  """Reject a margin the torch lane cannot hold.

  Parameters
  ----------
  margin : object
      The candidate margin, before coercion.
  name : str
      How to refer to it in an error message.

  Returns
  -------
  None

  Raises
  ------
  TypeError
      If the margin is not a :class:`torch.nn.Module`.
  NotImplementedError
      If the margin declares atoms but supplies no left limit.
  """
  if not isinstance(margin, torch.nn.Module):
    raise TypeError(
      f"TorchVinedist requires every margin to be a torch.nn.Module, so that "
      f"`.to(device)`, `state_dict()` and autograd reach its parameters; "
      f"{name} is a {type(margin).__name__}. Wrap a torch.distributions "
      "family with TorchMargin.from_distribution(...), or use TorchKde1d, "
      "which fits any of the three variable types. A NumPy margin such as "
      "Kde1d belongs in pyvinecopulib.core.Vinedist instead."
    )
  # Atoms are fine now that the copula half has a discrete cascade -- what a
  # margin with atoms must supply is the left limit the cascade differences.
  # `TorchKde1d` inherits `cdf_left` from `MarginBase`; `TorchMargin` cannot
  # have one, because `torch.distributions`' discrete families implement
  # neither `cdf` nor `icdf`.
  if getattr(margin, "var_type", "c") != "c" and (
    getattr(margin, "cdf_left", None) is None
  ):
    raise NotImplementedError(
      f"{name} declares var_type={getattr(margin, 'var_type')!r} but has no "
      "`cdf_left`, so the copula cannot difference its atoms. Use TorchKde1d, "
      "which supplies one, or subclass MarginBase, which derives it."
    )


def _check_copula(copula: Any) -> None:
  """Reject a copula the torch lane cannot evaluate.

  Parameters
  ----------
  copula : object
      The candidate copula.

  Returns
  -------
  None

  Raises
  ------
  TypeError
      If the copula is the compiled ``Vinecop``, which evaluates on NumPy.
  """
  from ..pyvinecopulib_ext import Vinecop

  if isinstance(copula, Vinecop):
    raise TypeError(
      "TorchVinedist cannot hold a compiled Vinecop: it evaluates on NumPy "
      "arrays, so it would detach every gradient and ignore `.to(device)`. "
      "Lift it first with TorchVinecop.from_vinecop(copula)."
    )


def _declared_kwargs(controls: Optional[Any]) -> dict[str, Any]:
  """Translate one margin's declared type and support into constructor kwargs.

  A kernel density takes both at construction, so a declaration has to reach
  it before the fit: a grid fitted unbounded is already padded past the data.

  Parameters
  ----------
  controls : FitControlsMargin, or None
      This variable's marginal configuration.

  Returns
  -------
  dict
      Keyword arguments for :class:`~pyvinecopulib.torch.TorchKde1d`.
  """
  kwargs: dict[str, Any] = {}
  var_type = getattr(controls, "var_type", None)
  if var_type == "d":
    kwargs["type"] = "discrete"
  elif var_type == "zi":
    kwargs["type"] = "zero_inflated"
  support = getattr(controls, "support", None)
  if support is not None:
    lo, hi = support
    if lo is not None and math.isfinite(lo):
      kwargs["xmin"] = float(lo)
    if hi is not None and math.isfinite(hi):
      kwargs["xmax"] = float(hi)
  return kwargs


class TorchVinedist(VinedistBase[Tensor], torch.nn.Module):
  """A vine distribution whose copula and margins are all PyTorch modules.

  Same surface as :class:`~pyvinecopulib.core.VinedistBase` — ``logpdf`` / ``pdf`` /
  ``cdf`` / ``loglik`` / ``sample`` / ``rosenblatt`` /
  ``inverse_rosenblatt`` / ``marginal_cdf`` / ``marginal_icdf`` — evaluated
  entirely in PyTorch, and additionally a :class:`torch.nn.Module`, so the
  copula and the margins are registered children.

  Assembled or fitted: pair a copula from :meth:`TorchVinecop.from_data` or
  :meth:`TorchVinecop.from_vinecop` with one margin per variable, or call
  :meth:`from_data` to fit both halves end to end in torch.

  Discrete and mixed data are supported, through :class:`TorchKde1d` margins:
  a margin that declares atoms must supply the left limit ``cdf_left`` the
  copula's discrete cascade differences, which :class:`TorchMargin` cannot
  (``torch.distributions``' discrete families implement neither ``cdf`` nor
  ``icdf``).

  Parameters
  ----------
  vinecop : VinecopLike
      A fitted vine copula on ``[0, 1]^d`` that evaluates in PyTorch —
      normally a :class:`~pyvinecopulib.torch.TorchVinecop`. A
      :class:`~pyvinecopulib.core.Vinecop` is **refused**, not merely
      discouraged: it evaluates on NumPy, so it would detach every gradient.
  margins : sequence of TorchMargin, or TorchMargin
      One margin per variable, each a :class:`torch.nn.Module`. A single margin
      is accepted when it carries array-valued parameters and so already
      represents all ``d``.

  Raises
  ------
  TypeError
      If a margin is not a :class:`torch.nn.Module`, or the copula is a
      :class:`~pyvinecopulib.core.Vinecop`.
  NotImplementedError
      If a margin declares atoms but does not provide ``cdf_left``.

  See Also
  --------
  pyvinecopulib.core.VinedistBase : The array-agnostic base both routes share.
  pyvinecopulib.core.Vinedist : The NumPy route.
  TorchMargin : The margins this holds.
  TorchVinecop : The copula this holds.

  Examples
  --------
  A three-variable distribution with learnable normal margins::

      import torch
      import pyvinecopulib as pv
      from pyvinecopulib.torch import TorchMargin, TorchVinecop, TorchVinedist

      u = pv.utils.to_pseudo_obs(x)
      dist = TorchVinedist(
        TorchVinecop.from_data(torch.as_tensor(u)),
        [
          TorchMargin(torch.distributions.Normal, {"loc": 0.0, "scale": 1.0})
          for _ in range(3)
        ],
      )
      loss = -dist.logpdf(torch.as_tensor(x)).mean()
      loss.backward()
  """

  def __init__(
    self,
    vinecop: Any,
    margins: Sequence[Any] | Any,
  ) -> None:
    # Initialize nn.Module explicitly, then hand over to the Vinedist seam:
    # TorchVinedist also subclasses Vinedist, whose __init__ chain would
    # otherwise shadow nn.Module's under super(). Same shape as TorchVinecop.
    torch.nn.Module.__init__(self)
    self._bind_dist(vinecop, margins)

  def _bind_dist(
    self,
    vinecop: Any,
    margins: Sequence[Any] | Any,
  ) -> None:
    """Validate the parts, then install them as registered children.

    Parameters
    ----------
    vinecop : VinecopLike
        The vine copula.
    margins : sequence of TorchMargin, or TorchMargin
        The margins.

    Returns
    -------
    None
    """
    _check_copula(vinecop)
    if isinstance(margins, (list, tuple)):
      for j, margin in enumerate(margins):
        _check_margin(margin, f"margins[{j}]")
    else:
      _check_margin(margins, "margins")

    super()._bind_dist(vinecop, margins)
    # `nn.Module` tracks a `ModuleList`, not the plain tuple the base stores, so
    # rebind through one: without it the margins' parameters are invisible to
    # `state_dict`, `.to()` and every optimizer. The cascades only iterate and
    # index the attribute, both of which a `ModuleList` supports.
    registered = cast("list[torch.nn.Module]", list(self._margins))
    self._margins = cast(Any, torch.nn.ModuleList(registered))

  # The vine copula this route fits; the margins need a placement, so
  # `_default_margins` is overridden rather than declared.
  vinecop_class: ClassVar[Optional[type]] = TorchVinecop

  #: The torch TLL fitter and the tree criterion are both unweighted, so a
  #: weighted request is refused rather than applied to the margins alone.
  supports_weighted_copula: bool = False

  #: No torch margin reads covariates, so a conditional fit is refused outright
  #: rather than answered with an unconditional one.
  supports_fit_covariates: bool = False

  @classmethod
  def _coerce_fit_data(
    cls,
    y: Any,
    weights: Optional[Any],
    controls: Optional[Any],
  ) -> tuple[Tensor, Optional[Tensor]]:
    """Put the fit inputs on one device, in one dtype.

    The controls decide the placement and the data follow, as documented.
    Reading it off ``y`` instead left the margins wherever the caller's data
    happened to be while the copula went to ``controls.device``, so
    ``state_dict`` spanned two devices and ``logpdf`` raised; an integer ``y``
    would likewise have given the margins an integer grid.

    Parameters
    ----------
    y : object
        Observations on the original scale.
    weights : object, or None
        Observation weights.
    controls : FitControlsTorchVinecop, or None
        Carries the device and dtype when the caller set them.

    Returns
    -------
    tuple
        The observations and weights, both on the resolved placement.
    """
    ya = torch.as_tensor(y)
    resolved = controls or FitControlsTorchVinecop()
    device = resolved.device if resolved.device is not None else ya.device
    if resolved.dtype is not None:
      dtype = resolved.dtype
    elif ya.dtype.is_floating_point:
      dtype = ya.dtype
    else:
      dtype = torch.get_default_dtype()
    ya = ya.to(device=device, dtype=dtype)
    if weights is not None:
      weights = torch.as_tensor(weights).to(device=device, dtype=dtype)
    return ya, weights

  @classmethod
  def _default_margins(
    cls,
    d: int,
    controls: Optional[Any] = None,
    margin_controls: Optional[Sequence[Any]] = None,
  ) -> Sequence[Any]:
    """One :class:`TorchKde1d` per variable, on the resolved placement.

    Parameters
    ----------
    d : int
        Number of variables.
    controls : FitControlsTorchVinecop, or None, optional
        Carries the device and dtype the margins share with the copula.
    margin_controls : sequence, or None, optional
        One marginal configuration per variable, whose declared variable type
        and support each margin is built with.

    Returns
    -------
    sequence of TorchKde1d
        One unfitted margin per variable.
    """
    resolved = controls or FitControlsTorchVinecop()
    # `TorchKde1d` fixes its own default dtype, so name one only when the
    # controls actually carry it.
    placement: dict[str, Any] = {"device": resolved.device}
    if resolved.dtype is not None:
      placement["dtype"] = resolved.dtype
    per_variable = margin_controls or [None] * d
    return [
      TorchKde1d(**placement, **_declared_kwargs(mc)) for mc in per_variable
    ]

  @classmethod
  def _fit_copula(
    cls,
    u: Any,
    *,
    var_types: list[str],
    controls: Optional[Any] = None,
    structure: Optional[Any] = None,
    weights: Optional[Any] = None,
  ) -> Any:
    """Fit a :class:`TorchVinecop` on the pseudo-observations.

    Parameters
    ----------
    u : Tensor, shape (n, d + k)
        The copula-scale layout the margins produced.
    var_types : list of str
        One ``"c"`` or ``"d"`` per variable.
    controls : FitControlsTorchVinecop, or None, optional
        Copula fit controls.
    structure : RVineStructure, or None, optional
        A fixed structure, or ``None`` to select one.
    weights : None, optional
        Always ``None`` here: ``supports_weighted_copula`` is ``False``, so
        ``from_data`` refuses a weighted request before reaching this hook.

    Returns
    -------
    TorchVinecop
        The fitted copula.
    """
    del weights
    resolved = controls or FitControlsTorchVinecop()
    return TorchVinecop.from_data(
      u,
      structure=structure,
      # The same resolved placement the margins got, so the copula cannot
      # default to a different dtype than they did.
      controls=dataclasses.replace(resolved, device=u.device, dtype=u.dtype),
      var_types=var_types,
    )

  @property
  def margins(self) -> tuple[MarginLike, ...]:
    """The margins, in variable order.

    Returns
    -------
    tuple of MarginLike
        One margin per variable, read off the registered ``ModuleList``.
    """
    return cast("tuple[MarginLike, ...]", tuple(self._margins))

  def _ref_tensor(self) -> Tensor:
    """A registered tensor to crib dtype and device from.

    Returns
    -------
    Tensor
        The first registered parameter or buffer, or an empty CPU tensor when
        the parts register neither.
    """
    for tensor in chain(self.parameters(), self.buffers()):
      return tensor
    return torch.empty(0, dtype=torch.float64)

  def _prep(self, a: Any) -> Tensor:
    """Bring one input array onto this distribution's dtype and device.

    ``as_tensor`` rather than ``tensor`` or ``detach``, so a tensor that already
    matches is returned untouched and a gradient-carrying one stays in the
    graph. It is what lets an estimator hand this object plain NumPy while the
    cascade and every margin stay on device.

    Parameters
    ----------
    a : array
        An input array, tensor or NumPy.

    Returns
    -------
    Tensor
        The same values, as a tensor placed like the registered ones.
    """
    ref = self._ref_tensor()
    return torch.as_tensor(a, dtype=ref.dtype, device=ref.device)

  def log_prob(self, y: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Alias of ``logpdf``, the ``torch.distributions`` spelling.

    Parameters
    ----------
    y : Tensor, shape (n, d), dtype float
        Observations on the original scale.
    x : Tensor, shape (n, k), or None, optional
        Exogenous covariates, forwarded as ``logpdf`` forwards them.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Joint log-density values.
    """
    return self.logpdf(y, x=x)
