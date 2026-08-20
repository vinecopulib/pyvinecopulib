"""PyTorch vine distribution: a torch vine copula with torch margins.

:class:`~pyvinecopulib.core.Vinedist` is already array-agnostic, so the joint
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
pyvinecopulib.core.Vinedist : The array-agnostic base.
TorchMargin : The margins this holds.
TorchVinecop : The copula this holds.
"""

from __future__ import annotations

from itertools import chain
from typing import Any, Optional, Sequence, cast

import torch
from torch import Tensor

from ..core import MarginLike, Vinedist
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
      If the margin declares atoms.
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
  # `TorchKde1d` models atoms, but `TorchVinecop` refuses a non-continuous
  # `var_types`, so the pair would never see the left-limit columns this margin
  # supplies. Refuse here, where the message can name the whole distribution,
  # rather than one tree level down.
  if getattr(margin, "var_type", "c") != "c":
    raise NotImplementedError(
      f"TorchVinedist is continuous-only; {name} declares "
      f"var_type={getattr(margin, 'var_type')!r}. TorchKde1d models atoms, but "
      "TorchVinecop has no discrete cascade yet; use "
      "pyvinecopulib.core.Vinedist for discrete or mixed data."
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


class TorchVinedist(Vinedist[Tensor], torch.nn.Module):
  """A vine distribution whose copula and margins are all PyTorch modules.

  Same surface as :class:`~pyvinecopulib.core.Vinedist` — ``logpdf`` / ``pdf`` /
  ``cdf`` / ``loglik`` / ``sample`` / ``rosenblatt`` /
  ``inverse_rosenblatt`` / ``marginal_cdf`` / ``marginal_icdf`` — evaluated
  entirely in PyTorch, and additionally a :class:`torch.nn.Module`, so the
  copula and the margins are registered children.

  Assembled rather than fitted: pair a copula from
  :meth:`TorchVinecop.from_data` or :meth:`TorchVinecop.from_vinecop` with one
  :class:`TorchMargin` per variable. There is no torch marginal estimator yet,
  so :meth:`from_data` refuses and points at
  :meth:`pyvinecopulib.core.Vinedist.from_data` for the NumPy two-step fit.

  Continuous variables only. A margin declaring atoms is refused at
  construction; the discrete cascade is on the NumPy side.

  Parameters
  ----------
  copula : VinecopLike
      A fitted vine copula on ``[0, 1]^d``, normally a
      :class:`~pyvinecopulib.torch.TorchVinecop`.
  margins : sequence of TorchMargin, or TorchMargin
      One margin per variable, each a :class:`torch.nn.Module`. A single margin
      is accepted when it carries array-valued parameters and so already
      represents all ``d``.

  Raises
  ------
  TypeError
      If a margin is not a :class:`torch.nn.Module`, or the copula is the
      compiled ``Vinecop``.
  NotImplementedError
      If a margin declares atoms.

  See Also
  --------
  pyvinecopulib.core.Vinedist : The array-agnostic base, and the NumPy route.
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
    copula: Any,
    margins: Sequence[Any] | Any,
  ) -> None:
    # Initialize nn.Module explicitly, then hand over to the Vinedist seam:
    # TorchVinedist also subclasses Vinedist, whose __init__ chain would
    # otherwise shadow nn.Module's under super(). Same shape as TorchVinecop.
    torch.nn.Module.__init__(self)
    self._bind_dist(copula, margins)

  def _bind_dist(
    self,
    copula: Any,
    margins: Sequence[Any] | Any,
  ) -> None:
    """Validate the parts, then install them as registered children.

    Parameters
    ----------
    copula : VinecopLike
        The vine copula.
    margins : sequence of TorchMargin, or TorchMargin
        The margins.

    Returns
    -------
    None
    """
    _check_copula(copula)
    if isinstance(margins, (list, tuple)):
      for j, margin in enumerate(margins):
        _check_margin(margin, f"margins[{j}]")
    else:
      _check_margin(margins, "margins")

    super()._bind_dist(copula, margins)
    # `nn.Module` tracks a `ModuleList`, not the plain tuple the base stores, so
    # rebind through one: without it the margins' parameters are invisible to
    # `state_dict`, `.to()` and every optimizer. The cascades only iterate and
    # index the attribute, both of which a `ModuleList` supports.
    registered = cast("list[torch.nn.Module]", list(self._margins))
    self._margins = cast(Any, torch.nn.ModuleList(registered))

  @classmethod
  def from_data(
    cls,
    y: Any,
    *,
    x: Optional[Tensor] = None,
    margins: Any = None,
    controls: Optional[FitControlsTorchVinecop] = None,
    structure: Optional[Any] = None,
    weights: Optional[Tensor] = None,
    names: Optional[Any] = None,
  ) -> "TorchVinedist":
    """Fit margins and a torch vine copula to data, in that order.

    End to end in torch: :class:`TorchKde1d` per column, then
    :meth:`TorchVinecop.from_data` on the copula data the margins produce. The
    result is on one device, in one dtype, and differentiable throughout --
    which the inherited NumPy route could not give, since it fits ``Kde1d``
    margins and a compiled ``Vinecop``.

    Parameters
    ----------
    y : Tensor, shape (n, d)
        Observations on the original scale.
    x : Tensor or None, optional
        Not supported; no torch margin reads covariates, and a silently
        unconditional fit is worse than a refusal.
    margins : object, optional
        Specification per :func:`pyvinecopulib.margins.resolve_margins`.
        ``None`` means one :class:`TorchKde1d` per column. Every resolved
        margin must be an ``nn.Module``.
    controls : FitControlsTorchVinecop or None, optional
        Copula fit controls.
    structure : RVineStructure or None, optional
        A fixed structure, or ``None`` to select one.
    weights : Tensor, shape (n,), or None, optional
        Observation weights, forwarded to every margin that accepts them.
    names : sequence of str or None, optional
        Variable names, used only to resolve a mapping specification.

    Returns
    -------
    TorchVinedist
        The fitted distribution.

    Raises
    ------
    NotImplementedError
        If ``x`` is given, or if any resolved margin declares atoms -- the
        torch copula cannot difference them accurately yet.
    """
    from ..margins import resolve_margins
    from ..margins._resolve import fit_margin

    if x is not None:
      raise NotImplementedError(
        "TorchVinedist.from_data takes no covariates: no torch margin reads "
        "them, so the margins would be fitted unconditionally while the call "
        "suggested otherwise. Fit the parts yourself if the copula alone is "
        "conditional."
      )
    ya = torch.as_tensor(y)
    if ya.ndim != 2:
      raise ValueError(f"y must be two-dimensional; got {tuple(ya.shape)}")
    d = int(ya.shape[1])
    default = [TorchKde1d(device=ya.device, dtype=ya.dtype) for _ in range(d)]
    specs = resolve_margins(margins, d, names=names, default=default)
    fitted = [fit_margin(specs[j], ya[:, j], weights=weights) for j in range(d)]
    u = cls.copula_data(fitted, ya)
    var_types = cls.copula_var_types(fitted)
    copula = TorchVinecop.from_data(
      u,
      structure=structure,
      controls=controls,
      var_types=var_types,
    )
    return cls(copula, fitted)

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
