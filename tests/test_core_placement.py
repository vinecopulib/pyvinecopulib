"""Tests for the placement hook behind every base's ``_prep``.

Placement is the one step of the input pipeline whose whole contract is
*inference* -- a subclass writes no conversion code, so what it infers from is
the only thing that can be wrong -- and none of that is observable through an
evaluation until it returns the wrong number.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from pyvinecopulib.core.extend import (
  place,
  prepare_covariates,
  reference_array,
  to_numpy,
  trim,
)


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


# --- what the hook is for: the three steps, applied where they belong -------- #


def test_covariates_are_placed_but_never_trimmed() -> None:
  """The covariate half of the pipeline: place, and do not clamp.

  Covariates are arbitrary reals, so the domain step that copula arguments get
  would corrupt them -- while the placement step is what lets a NumPy `x` meet
  the conditioning columns a PyTorch vine gathered.
  """
  torch = pytest.importorskip("torch")
  u = torch.linspace(0.0, 1.0, 6, dtype=torch.float32).reshape(3, 2)
  x = np.array([[-3.0], [0.0], [7.5]])
  placed = prepare_covariates(u, x, 3)
  assert isinstance(placed, torch.Tensor)
  assert placed.dtype is torch.float32
  # Untouched values: no clamp into (0, 1), which is the point.
  assert float(placed.min()) == -3.0 and float(placed.max()) == 7.5


def test_prepare_still_refuses_a_misaligned_covariate_matrix() -> None:
  """Placement does not replace the layout check; it follows it."""
  with pytest.raises(ValueError):
    prepare_covariates(np.zeros((3, 2)), np.zeros((4, 1)), 3)
  with pytest.raises(ValueError):
    prepare_covariates(np.zeros((3, 2)), np.zeros(3), 3)
  assert prepare_covariates(np.zeros((3, 2)), None, 3) is None


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


# --- the pipeline steps as a surface an extension can reach ------------------ #


def test_the_pipeline_steps_are_reachable_from_core() -> None:
  """A subclass composing ``_prep_args`` itself must not import a private module.

  The hooks it writes (``_prep``, ``_layout``) and the composite they feed
  (``_prep_args``) are public, so the steps composing them are too. Reaching
  them through `core._placement` / `core._trim` / `core._covariates` is what
  had a downstream extension run a covariate contract of its own, accepting a
  one-dimensional `x` on the same object whose base refuses one.
  """
  import pyvinecopulib.core as core

  # `getattr` because the generated stub declares no `__all__`, as the
  # neighboring surface and stub tests do for the same reason.
  import pyvinecopulib.core.extend as extend

  exported = set(getattr(extend, "__all__", ()))
  for name in (
    "collapse_data",
    "continuous_view",
    "covariate_row",
    "model_from_json",
    "place",
    "prepare_covariates",
    "reference_array",
    "reject_covariates",
    "to_numpy",
    "trim",
    "usable_observations",
    "validate_weights",
  ):
    assert name in exported, name
    assert callable(getattr(extend, name)), name
  # Not callables: the batching sentinel, the two fit-callback aliases and the
  # payload version.
  for name in ("FitEdge", "FitLevel", "MODEL_JSON_VERSION", "NotBatchable"):
    assert name in exported, name
    assert getattr(extend, name, None) is not None, name
  assert issubclass(extend.NotBatchable, Exception)
  # `core` keeps the type variable those signatures are written in, and none
  # of the rest: the two namespaces have different audiences.
  core_exported = set(getattr(core, "__all__", ()))
  assert "ArrayT" in core_exported
  assert exported.isdisjoint(core_exported)


def test_to_numpy_brings_back_what_asarray_refuses() -> None:
  """The return trip exists because ``np.asarray`` raises on this tensor."""
  torch = pytest.importorskip("torch")

  tracked = torch.ones(3, dtype=torch.float64, requires_grad=True)
  with pytest.raises(RuntimeError):
    np.asarray(tracked)
  back = to_numpy(tracked)
  assert isinstance(back, np.ndarray)
  np.testing.assert_array_equal(back, np.ones(3))
  # Neither `detach` nor `cpu` exists on a NumPy array, which passes through.
  plain = np.array([0.25, 0.75])
  assert to_numpy(plain) is plain


def test_trim_clamps_into_the_open_interval_at_its_own_precision() -> None:
  """The domain step, whose bounds have to be representable in the dtype.

  ``1 - 1e-10`` rounds to exactly ``1.0`` in ``float32``, so the historical
  ``float64`` pair would hand a downstream normal quantile an infinity.
  """
  from array_api_compat import array_namespace

  # Two dtype *names* rather than two classes: iterating the classes gives `a`
  # a union dtype, which numpy's own `.min` overloads do not accept -- visible
  # only once `trim` returns a real type instead of `Any`.
  for name in ("float64", "float32"):
    a = np.array([0.0, 1.0], dtype=name)
    clamped = trim(a)
    assert clamped.dtype == a.dtype
    assert float(clamped.min()) > 0.0
    assert float(clamped.max()) < 1.0
    # The namespace is the optional fast path, not part of the call.
    np.testing.assert_array_equal(clamped, trim(a, array_namespace(a)))


def test_a_part_holding_no_array_can_tell_that_placement_is_a_no_op() -> None:
  """Inference's third answer, and the documented way to detect it.

  ``_prep`` returns the values untouched when there is nothing to infer from:
  right for a part that computes in whatever namespace it is handed, and
  silently wrong for a torch part that keeps its device as a handle rather than
  as a tensor. ``reference_array(self) is None`` separates the two, and
  ``place`` taking an array as its own reference is the override.
  """
  torch = pytest.importorskip("torch")

  class _Deviced:
    """A device handle and a scalar; no array of its own."""

    def __init__(self) -> None:
      self.device = torch.device("cpu")
      self.threshold = 0.5

    def _prep(self, a: Any) -> Any:
      return place(self, a)

  part = _Deviced()
  u = np.array([[0.25, 0.75]])
  assert reference_array(part) is None
  # The no-op the check reports, rather than a placement or a failure.
  assert part._prep(u) is u
  reference = torch.empty(0, dtype=torch.float32, device=part.device)
  placed = place(reference, u)
  assert isinstance(placed, torch.Tensor)
  assert placed.dtype is torch.float32


def test_covariates_are_placed_through_the_hook_not_around_it() -> None:
  """``prepare_covariates`` honors an overridden ``_prep``.

  It used to place through the module-level ``place``, so an object whose
  placement is *declared* rather than inferable was honored on the argument
  path (``_prep_args`` calls the hook) and skipped on the covariate path -- the
  same object, the same call, two behaviors. Downstream that put a NumPy ``x``
  inside a pair copula whose backend then failed on ``.to(dtype=...)``.
  """

  class _Declared:
    """No array of its own, so only an override can place anything."""

    def __init__(self) -> None:
      self.seen: list[Any] = []

    def _prep(self, a: Any) -> Any:
      self.seen.append(a)
      return np.asarray(a, dtype=np.float32)

  onto = _Declared()
  out = prepare_covariates(onto, np.array([[1.0], [2.0]]), 2)
  assert onto.seen, "the hook was not consulted"
  assert out is not None and out.dtype == np.float32

  # An *array* as `onto` is the static fit engines' case: it has no hook, and
  # falls back to placing onto itself rather than raising.
  ref = np.zeros(2, dtype=np.float32)
  placed = prepare_covariates(ref, np.array([[1.0], [2.0]]), 2)
  assert placed is not None and placed.dtype == np.float32


def test_a_single_covariate_row_is_accepted_in_both_spellings() -> None:
  """``covariate_row`` is the narrower sibling, and narrower on purpose.

  ``prepare_covariates`` refuses a one-dimensional ``x`` because ``(n,)`` is
  ambiguous per row. A single row is not, which is why a plot takes ``(p,)``.
  """
  from pyvinecopulib.core.extend import covariate_row

  flat = covariate_row(np.array([1.0, 2.0, 3.0]))
  assert flat.shape == (1, 3)
  already = np.array([[1.0, 2.0, 3.0]])
  np.testing.assert_array_equal(covariate_row(already), already)

  for bad in (np.zeros((2, 3)), np.zeros((1, 2, 3))):
    with pytest.raises(ValueError, match="single covariate row"):
      covariate_row(bad)


def test_place_makes_no_promise_about_a_gradient() -> None:
  """What `place` inherits from the namespace, and the one case that raises.

  Conversion goes through the array namespace's own ``asarray``, so whether a
  tracked tensor stays tracked is that library's answer -- ``torch.asarray``
  defaulted ``requires_grad`` to ``False`` up to 2.11 and to the input's value
  from 2.13, so pinning either here would pin the installed torch rather than
  anything about this package. What *is* stable is the NumPy-reference case,
  and that it raises rather than silently detaching.
  """
  torch = pytest.importorskip("torch")

  tracked = torch.ones((2, 1), dtype=torch.float64, requires_grad=True)
  # A NumPy reference and a tracked tensor: reachable, because the static fit
  # engines pass an array as `onto`.
  with pytest.raises(RuntimeError, match="requires grad"):
    prepare_covariates(np.zeros((2, 3)), tracked, 2)
  # Detached, the same call is an ordinary placement onto NumPy.
  assert isinstance(
    prepare_covariates(np.zeros((2, 3)), tracked.detach(), 2), np.ndarray
  )
  # The torch route is the one that guarantees an answer, whatever the version.
  from pyvinecopulib.torch import TensorPlacementMixin

  class _Hooked(TensorPlacementMixin):
    device = torch.device("cpu")
    dtype = torch.float64

  assert _Hooked()._prep(tracked).requires_grad


def test_one_covariate_per_observation_is_checked_against_the_count() -> None:
  """`(n,)` is resolvable two ways, so `n` is what makes one of them safe.

  Choosing `covariate_column` is the caller stating that the axis is
  observations; `n` is what lets that statement be checked rather than
  trusted. Without it, an `x` of the wrong length reshapes to `(len(x), 1)`
  and pairs every observation with some other observation's covariate.
  """
  from pyvinecopulib.core.extend import covariate_column

  column = covariate_column(np.arange(3.0), 3)
  assert column.shape == (3, 1)
  np.testing.assert_array_equal(column[:, 0], np.arange(3.0))
  # `(n, 1)` says the same thing explicitly and survives the same check.
  already = np.zeros((3, 1))
  np.testing.assert_array_equal(covariate_column(already, 3), already)

  # The check that earns the required `n`: a length that cannot be n rows.
  with pytest.raises(ValueError, match="5 values but 3 observations"):
    covariate_column(np.arange(5.0), 3)
  # And a second column, having asked for one covariate.
  with pytest.raises(ValueError, match="one covariate per observation"):
    covariate_column(np.zeros((3, 2)), 3)


def test_the_two_readings_of_a_one_dimensional_x_stay_separate() -> None:
  """Where both readings are valid, the function named is the disambiguation.

  At `n == p` each produces a well-formed and *different* answer, which is the
  whole reason there are two names rather than a flag: no argument can tell
  them apart, so the call site has to.
  """
  from pyvinecopulib.core.extend import covariate_column, covariate_row

  x = np.array([0.1, 0.2, 0.3])
  assert covariate_column(x, 3).shape == (3, 1)
  assert covariate_row(x).shape == (1, 3)

  # Each refusal names the other, since picking the wrong one is the likely
  # mistake rather than a malformed array.
  with pytest.raises(ValueError, match="covariate_row"):
    covariate_column(np.zeros((1, 4)), 4)
  with pytest.raises(ValueError, match="covariate_column"):
    covariate_row(np.zeros((4, 1)))


def test_the_refusal_a_caller_hits_names_both_ways_out() -> None:
  """`prepare_covariates` is where a one-dimensional `x` is actually rejected.

  A message that only says "must have shape (n, p)" invites a reshape, which
  is the caller picking a reading by accident.
  """
  with pytest.raises(ValueError, match="covariate_column.*covariate_row"):
    prepare_covariates(np.zeros(3), np.arange(3.0), 3)
  # A two-dimensional mismatch is a different error and stays terse.
  with pytest.raises(ValueError, match="one row per observation"):
    prepare_covariates(np.zeros(3), np.zeros((2, 1)), 3)
