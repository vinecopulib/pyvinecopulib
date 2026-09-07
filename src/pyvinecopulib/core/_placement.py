"""Coerce an incoming array onto the array namespace an object evaluates on.

The four canonical bases each accept arrays from a caller and each hold arrays
of their own, and the two need not be the same type: a caller may reasonably
hand a NumPy matrix to a PyTorch vine, or a list of floats to either. Getting
one onto the other is *placement*, and it is one of three separable steps every
layer performs on its input:

- **placement** -- this module: which namespace, dtype and device the values
  must be in for the object's own arrays to combine with them.
- **layout** -- ``_validation``: which shapes are admissible, which the base
  knows from its own dimension and variable types.
- **domain** -- ``_trim``: clamping copula arguments into the open unit square
  at the working precision.

They are separated because the three do not always apply together. Exogenous
covariates are *placed* but never *trimmed* -- they are arbitrary reals, not
copula arguments -- and a manufactured evaluation grid needs placement without
any layout check at all.

Placement is *inferred* rather than declared, so hosting a custom pair copula,
margin or vine on PyTorch requires writing none of it: the object already holds
the tensors that answer the question, and :func:`reference_array` finds them.
A subclass whose arrays live somewhere this misses overrides ``_prep``.
"""

from __future__ import annotations

from typing import Any, Optional

from array_api_compat import array_namespace, device as _device_of

__all__ = ["place", "reference_array"]


def reference_array(obj: Any) -> Optional[Any]:
  """An array ``obj`` holds, naming where its numerics run.

  Looks in the two places an object keeps arrays, in the order that finds the
  most specific answer. ``parameters()`` and ``buffers()`` come first and are
  duck-typed rather than imported: they are how a ``torch.nn.Module`` exposes
  its tensors, and they recurse, so a grid held by a submodule counts -- while
  ``vars()`` on such a module yields the registries rather than the tensors.
  Then the instance's own attributes, for a subclass that stores a plain array.

  Duck-typing here mirrors what the plotting helper already does on the way
  *out*, where a returned density is brought to the host through ``detach()``
  and ``cpu()`` without importing PyTorch either.

  Parameters
  ----------
  obj : object
      The margin, pair copula, vine or vine distribution.

  Returns
  -------
  array, or None
      One of the object's arrays, or ``None`` when it holds none -- which is
      the right answer for a functional part that computes in whatever
      namespace it is handed.
  """
  for name in ("parameters", "buffers"):
    method = getattr(obj, name, None)
    if callable(method):
      try:
        for tensor in method():
          return tensor
      except TypeError:
        # Not the nn.Module member of that name; fall through to the next.
        continue
  for value in vars(obj).values():
    if _is_array(value):
      return value
  return None


def _is_array(value: Any) -> bool:
  """Whether ``value`` is an array the array API recognizes.

  Asked of the namespace rather than by duck-typing ``dtype`` and ``shape``,
  which a *module* also carries: ``array_api_compat.numpy`` re-exports both
  ``numpy.dtype`` and ``numpy.shape``, so an attribute sweep would mistake a
  memoized namespace for one of the object's arrays.

  Parameters
  ----------
  value : object
      The candidate.

  Returns
  -------
  bool
      ``True`` if an array namespace claims it.
  """
  if getattr(value, "dtype", None) is None:
    return False
  try:
    array_namespace(value)
  except TypeError:
    return False
  return True


def place(obj: Any, a: Any) -> Any:
  """Coerce ``a`` onto the namespace, dtype and device ``obj`` evaluates on.

  Placement only: no shape is checked and no value is clamped, so this is
  equally correct for a covariate matrix, a manufactured grid, and a block of
  pseudo-observations.

  Parameters
  ----------
  obj : object
      The object whose placement to match, read through
      :func:`reference_array`.
  a : array
      The values to place.

  Returns
  -------
  array
      ``a`` on ``obj``'s namespace, or unchanged when ``obj`` holds no array
      of its own to read a placement from.
  """
  reference = reference_array(obj)
  if reference is None:
    return a
  xp = array_namespace(reference)
  # Already there: skip the conversion, which would otherwise warn about the
  # `requires_grad` flag it inherits from a tensor that tracks one.
  if getattr(a, "dtype", None) is reference.dtype:
    try:
      if array_namespace(a) is xp:
        return a
    except TypeError:
      pass
  return xp.asarray(a, dtype=reference.dtype, device=_device_of(reference))
