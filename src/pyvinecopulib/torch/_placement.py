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

# Everything the mixin does rests on being mixed into an `nn.Module`: it reads
# the module's registered tensors. Declaring that to the type checker is what
# lets `reference_tensor` name the type it actually takes, and it costs nothing
# at runtime, `class X(object)` being `class X`. A real base here would put
# `nn.Module` ahead of the canonical base in every subclass's MRO.
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
  """The ``_prep`` hook for a module placed on its own registered tensors.

  The torch counterpart of
  ``pyvinecopulib.core._placement.PlacementMixin``, which it shadows: a
  subclass mixes this in **ahead** of its canonical base, so ``_prep``
  resolves here rather than to the array-API inference the other classes use.

  Only ordinary private members belong on a mixin at that position. It lands
  ahead of ``torch.nn.Module`` in the resulting MRO, so anything defined here
  that ``nn.Module`` also defines -- a dunder above all -- would silently
  shadow it.
  """

  def _ref_tensor(self) -> Tensor:
    """A tensor carrying the dtype and device this module evaluates on.

    Falls back to an empty CPU ``float64`` tensor for a module that registers
    none, so a factory holding no parameters still places its input.
    """
    ref = reference_tensor(self)
    if ref is None:
      return torch.empty(0, dtype=torch.float64)
    return ref

  def _prep(self, a: Any) -> Tensor:  # noqa: ANN401 - any array type, placed
    """Bring one input array onto this module's dtype and device.

    ``as_tensor`` rather than ``tensor`` or ``detach``, so a tensor that
    already matches is returned untouched and a gradient-carrying one stays in
    the graph.
    """
    ref = self._ref_tensor()
    return torch.as_tensor(a, dtype=ref.dtype, device=ref.device)
