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


def test_a_matching_dtype_does_not_excuse_the_wrong_device() -> None:
  """`place` promises namespace, dtype *and* device; the fast path skipped one.

  An input whose dtype already agrees took an early return that compared only
  the namespace, so it stayed where it was -- on the host, while the object
  evaluates somewhere else. Torch's ``meta`` device makes that reachable
  without an accelerator: it is a real device that every machine has.
  """
  torch = pytest.importorskip("torch")
  from array_api_compat import device as device_of

  from pyvinecopulib.core._placement import place

  class Holder:
    def __init__(self, grid: object) -> None:
      self.grid = grid

  elsewhere = torch.zeros(3, dtype=torch.float64, device="meta")
  host = torch.zeros(3, dtype=torch.float64)

  assert device_of(place(Holder(elsewhere), host)) == device_of(elsewhere)
  # And the fast path still returns the input untouched when it is already
  # there, which is what keeps a gradient-carrying tensor from being copied.
  assert place(Holder(host), host) is host


# --- the torch lane's own search, which has to agree with this one ----------- #


def test_the_torch_search_prefers_a_float_tensor_over_an_earlier_integer() -> (
  None
):
  """The two searches must rank candidates the same way.

  ``pyvinecopulib.torch._placement`` reads a module's tensors directly rather
  than through the array API, because placement is on a per-evaluation path.
  Ranking them differently is what made one object hold two placements at
  once: a first-hit search answered with an integer tensor while ``_prep``
  answered with the float one beside it.
  """
  torch = pytest.importorskip("torch")

  from pyvinecopulib.torch._placement import reference_tensor

  class _Module(torch.nn.Module):
    def __init__(self) -> None:
      super().__init__()
      self.register_buffer("codes", torch.tensor([1, 2, 3]))
      self.register_buffer(
        "grid", torch.linspace(0.0, 1.0, 5, dtype=torch.float64)
      )

  module = _Module()
  reference = reference_tensor(module)
  assert reference is not None and reference.dtype is torch.float64
  through_the_array_api = reference_array(module)
  assert through_the_array_api is not None
  assert through_the_array_api.dtype is reference.dtype
  # Nothing floating to name: the caller's own default answers, rather than an
  # integer dtype that neither `torch.as_tensor` nor `torch.rand` can use.
  assert reference_tensor(torch.nn.Module()) is None


def test_a_margin_with_an_integer_parameter_can_still_be_sampled() -> None:
  """A family parameterized by a count registers an integer tensor.

  ``Gamma(concentration=2, rate=1.0)`` keeps the count's own dtype, so the
  margin holds an integer tensor ahead of a float one, and ``Chi2(df=5)``
  holds nothing else at all. Placement always read past those; the draw did
  not, and asked for uniforms in ``int64``.
  """
  torch = pytest.importorskip("torch")

  from pyvinecopulib.torch import TorchDistributionMargin

  for factory, parameters in (
    (torch.distributions.Gamma, {"concentration": 2, "rate": 1.0}),
    (torch.distributions.Chi2, {"df": 5}),
  ):
    margin = TorchDistributionMargin(
      factory, parameters=parameters, trainable=False
    )
    placed = margin._prep(np.array([0.25, 0.75]))
    drawn = margin.sample(4, seeds=[7])
    assert placed.dtype is torch.float64
    assert drawn.dtype is placed.dtype
    assert drawn.shape == (4,)
    assert bool(torch.all(drawn > 0.0))
