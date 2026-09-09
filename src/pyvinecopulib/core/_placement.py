"""Coerce an incoming array onto the array namespace an object evaluates on.

The four canonical bases each accept arrays from a caller and each hold arrays
of their own, and the two need not be the same type: a caller may reasonably
hand a NumPy matrix to a PyTorch vine, or a list of floats to either. Getting
one onto the other is *placement*, and it is one of three separable steps every
layer performs on its input:

- **placement** -- this module: which namespace, dtype and device the values
  must be in for the object's own arrays to combine with them.
- **layout** -- the ``_layout`` hook, one per level: which shapes are
  admissible, which the base knows from its own dimension and variable types.
- **domain** -- ``_trim``: clamping copula arguments into the open unit square
  at the working precision.

They are separated because the three do not always apply together. Exogenous
covariates are *placed* but never *trimmed* -- they are arbitrary reals, not
copula arguments, and ``_covariates.prepare`` is the composite that applies
exactly those two steps to them. A manufactured evaluation grid needs placement
without any layout check at all, and a margin's argument is on the data scale,
so it is placed and checked but never clamped.

Placement is *inferred* rather than declared, so hosting a custom pair copula,
margin or vine on PyTorch requires writing none of it: the object already holds
the tensors that answer the question, and :func:`reference_array` finds them.
A subclass whose arrays live somewhere this misses overrides ``_prep``.

An array reaches these levels one other way -- drawn from the array
library's own RNG,
which is the one thing no inference can supply -- so the raising hook for that
lives here too, beside the hook that places what a caller supplied. The two
copula bases share it verbatim. The marginal one takes no ``qrng`` flag,
since ``MarginBase.sample`` offers none, and stays on ``MarginBase``.
"""

from __future__ import annotations

from typing import Any, Generic, Optional, cast

import numpy as np
from array_api_compat import array_namespace, device as _device_of

from .protocols import ArrayT

__all__ = [
  "PlacementMixin",
  "QrngUniformMixin",
  "place",
  "reference_array",
]


def to_numpy(a: ArrayT) -> np.ndarray:
  """Host NumPy view of ``a`` -- placement's return trip.

  The one direction inference cannot help with: a value has to come back to
  NumPy wherever something outside the array namespace reads it -- the
  criterion binding, a plot, a printed report. ``np.asarray`` alone raises on
  a tensor that requires grad and again on one that lives on an accelerator,
  so this detaches and transfers first, both through ``getattr`` since a
  NumPy array has neither method.

  Parameters
  ----------
  a : array
      Values in any array namespace.

  Returns
  -------
  numpy.ndarray
      The same values, on the host, detached.
  """
  v: Any = a
  detach = getattr(v, "detach", None)
  if detach is not None:
    v = detach()
  cpu = getattr(v, "cpu", None)
  if cpu is not None:
    v = cpu()
  return np.asarray(v)


def reference_array(obj: object) -> Optional[Any]:  # noqa: ANN401
  """An array ``obj`` holds, naming where its numerics run.

  Looks in the two places an object keeps arrays, in the order that finds the
  most specific answer. ``parameters()`` and ``buffers()`` come first and are
  duck-typed rather than imported: they are how a ``torch.nn.Module`` exposes
  its tensors, and they recurse, so a grid held by a submodule counts -- while
  ``vars()`` on such a module yields the registries rather than the tensors.
  Then the instance's own attributes, for a subclass that stores a plain array.

  A **floating-point** array wins wherever one is available, because what the
  answer is used for is placing copula arguments and covariates -- both real
  valued. An object may hold an integer array too (an index table, a variable
  -type code, a count buffer), and it is only the fallback: adopting its dtype
  would truncate every argument to zero.

  Duck-typing here mirrors what the plotting helper already does on the way
  *out*, where a returned density is brought to the host through ``detach()``
  and ``cpu()`` without importing PyTorch either.

  An **array** is its own reference, which is what lets a caller holding no
  object -- the fit engines, which are static -- place covariates onto the
  observations being fitted.

  Parameters
  ----------
  obj : object
      The margin, pair copula, vine or vine distribution, or an array to place
      onto directly.

  Returns
  -------
  array, or None
      One of the object's arrays, or ``None`` when it holds none -- which is
      the right answer for a functional part that computes in whatever
      namespace it is handed.
  """
  if _is_array(obj):
    return obj
  fallback: Optional[Any] = None
  for name in ("parameters", "buffers"):
    method = getattr(obj, name, None)
    if callable(method):
      try:
        for tensor in method():
          if _is_float(tensor):
            return tensor
          if fallback is None:
            fallback = tensor
      except TypeError:
        # Not the nn.Module member of that name; fall through to the next.
        continue
  for value in vars(obj).values():
    if not _is_array(value):
      continue
    if _is_float(value):
      return value
    if fallback is None:
      fallback = value
  return fallback


def _is_float(value: object) -> bool:
  """Whether ``value`` is an array with a floating-point dtype.

  Parameters
  ----------
  value : object
      The candidate.

  Returns
  -------
  bool
      ``True`` for a real floating-point array; ``False`` for an integer or
      boolean one, and for anything no array namespace claims.
  """
  if not _is_array(value):
    return False
  xp = array_namespace(value)
  return bool(xp.isdtype(cast("Any", value).dtype, "real floating"))


def _is_array(value: object) -> bool:
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


def place(obj: object, a: Any) -> Any:  # noqa: ANN401
  """Coerce ``a`` onto the namespace, dtype and device ``obj`` evaluates on.

  Placement only: no shape is checked and no value is clamped, so this is
  equally correct for a covariate matrix, a manufactured grid, and a block of
  pseudo-observations.

  Parameters
  ----------
  obj : object
      The object whose placement to match, read through
      :func:`reference_array` -- or an array to match directly.
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
  # Only a floating reference names a dtype worth adopting. An integer one
  # names a namespace and a device and nothing else: casting `0.25` to it
  # would place a zero, so the values keep their own precision instead.
  dtype = reference.dtype if _is_float(reference) else None
  # Already there: skip the conversion, which would otherwise warn about the
  # `requires_grad` flag it inherits from a tensor that tracks one. The device
  # has to agree as well -- an array whose dtype already matches can still sit
  # on the host while the reference is on an accelerator, and placing it is
  # what this function is for.
  if dtype is None or getattr(a, "dtype", None) is dtype:
    try:
      if array_namespace(a) is xp and _device_of(a) == _device_of(reference):
        return a
    except TypeError:
      pass
  return xp.asarray(a, dtype=dtype, device=_device_of(reference))


class PlacementMixin:
  """The default ``_prep`` hook, shared by all four canonical bases.

  All four need the same hook -- one array, brought onto the namespace the
  object evaluates on -- and inferring that from the arrays the object already
  holds is what makes hosting a subclass on PyTorch require writing none of
  it. The four wrote it out identically; it lives here instead, beside the
  ``place`` it delegates to.
  """

  def _prep(self, a: Any) -> Any:  # noqa: ANN401
    """Bring one input array onto the namespace this object evaluates on.

    Placement only -- no shape or layout check and no clamping -- so it is
    equally correct for exogenous covariates, which must be placed but never
    trimmed, and for an array the class manufactures itself. Each base's
    ``_prep_args`` is the composite that adds the layout and domain steps a
    copula argument needs.

    The default infers the placement from the arrays this object already
    holds. Override it where those arrays live somewhere the inference
    misses.

    Parameters
    ----------
    a : array
        An input array on any namespace.

    Returns
    -------
    array
        The same values, on this object's namespace, dtype and device.
    """
    return place(self, a)


class QrngUniformMixin(Generic[ArrayT]):
  """The uniform-draw hook ``BicopBase`` and ``VinecopBase`` share.

  Kept apart from ``PlacementMixin``, which all four bases inherit, because
  the marginal one's hook has a different signature -- see the module
  docstring.
  """

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> ArrayT:
    """Draw the base uniforms ``sample`` transforms into a copula draw.

    Raising default; override it to enable ``sample``, which is all ``sample``
    needs. NumPy and PyTorch differ on RNG, so this is the one hook with no
    array-agnostic default. Named after
    :func:`pyvinecopulib.utils.sample_uniform`.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    qrng : bool
        Whether to draw a quasi-random (low-discrepancy) sequence.
    seeds : list of int
        RNG seeds.

    Returns
    -------
    array, shape (n, 2) or (n, d), dtype float
        Base uniforms in ``[0, 1)`` -- two columns for a pair copula, one per
        variable for a vine.

    Raises
    ------
    NotImplementedError
        Unless a subclass overrides this hook.
    """
    raise NotImplementedError(
      f"{type(self).__name__} does not implement _sample_uniform; override it "
      "to enable sample()."
    )
