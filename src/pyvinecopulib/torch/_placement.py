"""Placement for the torch lane: the tensor a module evaluates on.

``pyvinecopulib.core._placement`` answers the same question for any array
namespace, by asking ``array_api_compat`` which namespace each candidate array
belongs to. Here the namespace is known and the question is on a
per-evaluation path -- every ``pdf`` call places its argument -- so the search
reads the tensors directly: 12.1 us per call against the array-API route's
64.9 us.

The rule about *which* tensor is the same one, and it is not "the first": a
**floating-point** tensor wins wherever the module holds one. A module may
register an integer tensor as readily as a float one -- a family parameterized
by a count, a variable-type code, an index table -- and adopting its dtype
places every copula argument and every uniform draw at zero.
"""

from __future__ import annotations

from itertools import chain
from typing import TYPE_CHECKING, Any, Optional

import torch
from torch import Tensor

__all__ = ["TensorPlacementMixin", "reference_tensor"]

# The mixin reads a module's registered tensors *when there are any*, so the
# `nn.Module` declaration is for the type checker -- it lets `reference_tensor`
# name the type it takes -- and costs nothing at run time, `class X(object)`
# being `class X`. A real base here would put `nn.Module` ahead of the
# canonical base in every subclass's MRO. A host that is not a module is
# therefore fine, and resolves its placement from a declaration instead.
if TYPE_CHECKING:
  _ModuleBase = torch.nn.Module
else:
  _ModuleBase = object


def reference_tensor(module: torch.nn.Module) -> Optional[Tensor]:
  """A floating-point tensor ``module`` holds, naming where its numerics run.

  ``parameters()`` and ``buffers()`` recurse, so a grid held by a submodule
  counts, while a plain attribute does not -- registering the tensors it
  evaluates with is what makes a module movable by ``.to(device)`` at all.

  An integer tensor is **not** a fallback here, unlike in
  ``pyvinecopulib.core._placement``, where one still names a namespace and a
  device: every caller of this has a floating default of its own.

  Parameters
  ----------
  module : torch.nn.Module
      The module to read a placement from.

  Returns
  -------
  Tensor, or None
      A registered floating-point parameter or buffer, or ``None`` when the
      module registers none.
  """
  for tensor in chain(module.parameters(), module.buffers()):
    if tensor.is_floating_point():
      return tensor
  return None


class TensorPlacementMixin(_ModuleBase):
  """The ``_prep`` hook for a class whose numerics run on tensors.

  The torch counterpart of
  ``pyvinecopulib.core._placement.PlacementMixin``, which it shadows: a
  subclass mixes this in **ahead** of its canonical base, so ``_prep``
  resolves here rather than to the array-API inference the other classes use.
  Getting that order wrong is the failure to watch for -- the base's
  ``PlacementMixin`` linearizes first and ``_prep`` silently becomes the
  array-API inference again.

  A host does **not** have to be an ``nn.Module``. Placement resolves in three
  steps: a registered floating-point tensor if the host has any, then a
  ``device`` and ``dtype`` the host *declares*, then an empty CPU ``float64``
  tensor. So a class that keeps its device as a handle rather than as a tensor
  -- no parameters, no buffers -- says so once by exposing those two
  attributes, instead of holding a dummy tensor for the inference to find or
  writing the hook itself. Anything else overrides ``_ref_tensor``.

  Only ordinary private members belong on a mixin at that position. Where the
  host *is* an ``nn.Module`` this lands ahead of it in the MRO, so anything
  defined here that ``nn.Module`` also defines -- a dunder above all -- would
  silently shadow it.
  """

  def _ref_tensor(self) -> Tensor:
    """A tensor carrying the dtype and device this class evaluates on.

    Resolved in the three steps the class docstring lists. Override this where
    none of them fits; returning ``torch.empty(0, dtype=..., device=...)`` is
    all such an override needs.
    """
    # Guarded rather than called outright: `reference_tensor` reads
    # `parameters()` / `buffers()`, which a host that is not an `nn.Module`
    # does not have, and reaching for them there raised `AttributeError`
    # instead of resolving a placement.
    if hasattr(self, "parameters"):
      ref = reference_tensor(self)
      if ref is not None:
        return ref
    declared = torch.empty(
      0,
      dtype=getattr(self, "dtype", None) or torch.float64,
      device=getattr(self, "device", None) or "cpu",
    )
    return declared

  def _prep(self, a: Any) -> Tensor:  # noqa: ANN401 - any array type, placed
    """Bring one input array onto this module's dtype and device.

    ``as_tensor`` rather than ``tensor`` or ``detach``, so a tensor that
    already matches is returned untouched and a gradient-carrying one stays in
    the graph. That is a *guarantee* here and not one the array-API route
    makes: ``place`` converts through the namespace's own ``asarray``, whose
    answer for a tracked tensor is that library's and has changed between
    releases of it.
    """
    ref = self._ref_tensor()
    return torch.as_tensor(a, dtype=ref.dtype, device=ref.device)
