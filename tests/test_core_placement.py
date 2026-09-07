"""Tests for the placement seam behind every base's ``_prep``.

``pyvinecopulib.core._placement`` is reached directly here. It is the one step
of the input pipeline whose whole contract is *inference* -- a subclass writes
no conversion code, so what it infers from is the only thing that can be wrong
-- and none of that is observable through the public surface until an
evaluation returns the wrong number.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from pyvinecopulib.core._placement import place, reference_array


class _Holder:
  """An object whose arrays are plain instance attributes."""

  def __init__(self, **arrays: Any) -> None:
    for name, value in arrays.items():
      setattr(self, name, value)


_U = np.array([[0.25, 0.75], [0.5, 0.5]])


def test_a_float_array_wins_over_an_integer_one_stored_first() -> None:
  """The reference must be the array whose dtype the values can live in.

  An object may hold an integer array as well as a float one -- an index
  table, a variable-type code, a count buffer -- and insertion order is no
  guide to which describes its numerics.
  """
  obj = _Holder(codes=np.array([1, 2, 3]), grid=np.linspace(0.0, 1.0, 5))
  reference = reference_array(obj)
  assert reference is not None and reference.dtype == np.float64
  np.testing.assert_array_equal(place(obj, _U), _U)


def test_an_integer_reference_places_without_casting() -> None:
  """With no float array to read, the dtype is the input's, not the reference's.

  Adopting an integer dtype would truncate every copula argument to zero,
  which is a wrong answer rather than a failure.
  """
  obj = _Holder(codes=np.array([1, 2, 3]))
  reference = reference_array(obj)
  assert reference is not None and reference.dtype.kind in "iu"
  placed = place(obj, _U)
  assert placed.dtype.kind == "f"
  np.testing.assert_array_equal(placed, _U)


def test_an_object_holding_no_array_leaves_the_input_alone() -> None:
  """The right answer for a functional part is "wherever you are"."""
  obj = _Holder(threshold=0.5, name="pair")
  assert reference_array(obj) is None
  assert place(obj, _U) is _U


def test_a_memoized_namespace_is_not_mistaken_for_an_array() -> None:
  """``array_api_compat.numpy`` carries ``dtype`` and ``shape`` of its own."""
  from array_api_compat import array_namespace

  obj = _Holder(_xp=array_namespace(_U))
  assert reference_array(obj) is None


def test_placement_follows_a_float32_reference() -> None:
  """A float reference does name a precision, and it is honored."""
  obj = _Holder(grid=np.linspace(0.0, 1.0, 5, dtype=np.float32))
  placed = place(obj, _U)
  assert placed.dtype == np.float32


def test_a_torch_module_is_read_through_its_buffers() -> None:
  """Registered tensors are how an ``nn.Module`` exposes its placement."""
  torch = pytest.importorskip("torch")

  class _Module(torch.nn.Module):
    def __init__(self) -> None:
      super().__init__()
      self.register_buffer("codes", torch.tensor([1, 2, 3]))
      self.register_buffer(
        "grid", torch.linspace(0.0, 1.0, 5, dtype=torch.float32)
      )

  module = _Module()
  reference = reference_array(module)
  assert reference is not None and reference.dtype is torch.float32
  placed = place(module, _U)
  assert isinstance(placed, torch.Tensor)
  assert placed.dtype is torch.float32


# --- what the seam is for: the three steps, applied where they belong -------- #


def test_covariates_are_placed_but_never_trimmed() -> None:
  """`prepare` is the covariate half of the pipeline: place, do not clamp.

  Covariates are arbitrary reals, so the domain step that copula arguments get
  would corrupt them -- while the placement step is what lets a NumPy `x` meet
  the conditioning columns a PyTorch vine gathered.
  """
  from pyvinecopulib.core._covariates import prepare

  torch = pytest.importorskip("torch")
  u = torch.linspace(0.0, 1.0, 6, dtype=torch.float32).reshape(3, 2)
  x = np.array([[-3.0], [0.0], [7.5]])
  placed = prepare(u, x, 3)
  assert isinstance(placed, torch.Tensor)
  assert placed.dtype is torch.float32
  # Untouched values: no clamp into (0, 1), which is the point.
  assert float(placed.min()) == -3.0 and float(placed.max()) == 7.5


def test_prepare_still_refuses_a_misaligned_covariate_matrix() -> None:
  """Placement does not replace the layout check; it follows it."""
  from pyvinecopulib.core._covariates import prepare

  with pytest.raises(ValueError):
    prepare(np.zeros((3, 2)), np.zeros((4, 1)), 3)
  with pytest.raises(ValueError):
    prepare(np.zeros((3, 2)), np.zeros(3), 3)
  assert prepare(np.zeros((3, 2)), None, 3) is None
