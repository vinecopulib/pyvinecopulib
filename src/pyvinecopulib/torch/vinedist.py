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

from typing import Any, NoReturn, Sequence, cast

import torch
from torch import Tensor

from ..core import MarginLike, Vinedist

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
      "family with TorchMargin.from_distribution(...). A NumPy margin such as "
      "Kde1dMargin belongs in pyvinecopulib.core.Vinedist instead — which is "
      "also what Vinedist.from_data builds, since there is no torch marginal "
      "estimator yet."
    )
  if getattr(margin, "var_type", "c") != "c":
    raise NotImplementedError(
      f"TorchVinedist is continuous-only; {name} declares "
      f"var_type={getattr(margin, 'var_type')!r}. The discrete cascade lives "
      "on pyvinecopulib.core.Vinedist; use it with a VinecopBackend-style "
      "NumPy copula for discrete or mixed data."
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
  ``cdf`` / ``loglik`` / ``simulate`` / ``rosenblatt`` /
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
  def from_data(cls, *args: Any, **kwargs: Any) -> NoReturn:
    """Refuse the two-step fit: there is no torch marginal estimator.

    The inherited route fits ``Kde1dMargin`` margins and a compiled ``Vinecop``,
    neither of which this class can hold. Fit on the NumPy lane with
    :meth:`pyvinecopulib.core.Vinedist.from_data`, or assemble the torch parts
    yourself from :meth:`TorchVinecop.from_data` and :class:`TorchMargin`.

    Parameters
    ----------
    *args, **kwargs
        Accepted only so that a caller sees this explanation rather than a
        signature mismatch.

    Raises
    ------
    NotImplementedError
        Always.
    """
    del args, kwargs
    raise NotImplementedError(
      "TorchVinedist is assembled, not fitted: there is no torch marginal "
      "estimator, so the two-step fit would produce Kde1dMargin margins and a "
      "compiled Vinecop. Use pyvinecopulib.core.Vinedist.from_data for the "
      "NumPy fit, or build TorchVinedist(TorchVinecop.from_data(u), "
      "[TorchMargin(...), ...]) yourself."
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

  def log_prob(self, x: Tensor) -> Tensor:
    """Alias of ``logpdf``, the ``torch.distributions`` spelling.

    Parameters
    ----------
    x : Tensor, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Joint log-density values.
    """
    return self.logpdf(x)
