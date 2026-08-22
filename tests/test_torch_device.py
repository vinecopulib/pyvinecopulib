"""The torch cascade on CUDA: same numbers, same placement, no host detours.

The C++-parity suite in ``test_torch_bicop`` / ``test_torch_vinecop`` pins
torch against the compiled library on the cpu, in float64. This file pins
cuda against *that*, so the compiled-library agreement carries over to the
device without running the whole parity suite twice.

``float32`` is compared only to our own ``float64`` answer, never to the
compiled library: it is a lower precision, so the float64 tolerances do not
apply to it and loosening them would weaken the gate that does.
"""

import numpy as np
import pytest

torch = pytest.importorskip("torch")

import pyvinecopulib as pv  # noqa: E402
from pyvinecopulib.torch import (  # noqa: E402
  FitControlsTorchVinecop,
  TorchBicop,
  TorchVinecop,
)
from tests.helpers import assert_on_device, count_transfers  # noqa: E402

# cuda vs cpu at the same dtype. Tight enough that the cpu-vs-C++ 1e-10
# gates carry over to cuda by the triangle inequality.
DEVICE_TOL = {torch.float64: 1e-12, torch.float32: 2e-5}
# float32 against our own float64 result.
F32_VS_F64_TOL = 5e-5

_EVAL_OPS = ("pdf", "rosenblatt", "inverse_rosenblatt")


def _u(d: int, n: int, seed: int) -> np.ndarray:
  rng = np.random.default_rng(seed)
  a = rng.standard_normal((d, d))
  s = a @ a.T + d * np.eye(d)
  sd = np.sqrt(np.diag(s))
  return pv.to_pseudo_obs(
    rng.multivariate_normal(np.zeros(d), s / np.outer(sd, sd), size=n)
  )


@pytest.fixture
def cpp_vine() -> pv.Vinecop:
  """A fitted reference vine.

  Function-scoped like the rest of the suite: a module-scoped compiled
  object outlives pytest's own teardown and nanobind reports it as a leak
  at interpreter shutdown.
  """
  return pv.Vinecop.from_data(
    _u(5, 600, 7),
    controls=pv.FitControlsVinecop(family_set=[pv.BicopFamily.tll]),
  )


@pytest.fixture
def u_eval() -> np.ndarray:
  return _u(5, 400, 8)


def _np(t: "torch.Tensor") -> np.ndarray:
  return t.detach().cpu().numpy()


@pytest.mark.parametrize("op", _EVAL_OPS)
@pytest.mark.parametrize("cache", [False, True])
def test_eval_matches_cpu(
  device: str, op: str, cache: bool, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """Every evaluation op returns the cpu answer, on the device."""
  ref = TorchVinecop.from_vinecop(
    cpp_vine, cache_integrals=cache, device=torch.device("cpu")
  )
  got = TorchVinecop.from_vinecop(
    cpp_vine, cache_integrals=cache, device=torch.device(device)
  )
  a = getattr(ref, op)(torch.as_tensor(u_eval))
  b = getattr(got, op)(torch.as_tensor(u_eval, device=device))
  assert b.device.type == torch.device(device).type
  np.testing.assert_allclose(
    _np(b),
    _np(a),
    rtol=DEVICE_TOL[torch.float64],
    atol=DEVICE_TOL[torch.float64],
  )


@pytest.mark.parametrize("op", ("pdf", "rosenblatt"))
def test_batched_matches_unbatched_on_device(
  device: str, op: str, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """The batched fast path agrees with the per-edge cascade on the device."""
  vine = TorchVinecop.from_vinecop(cpp_vine, device=torch.device(device))
  ut = torch.as_tensor(u_eval, device=device)
  np.testing.assert_allclose(
    _np(getattr(vine, op)(ut, batched=True)),
    _np(getattr(vine, op)(ut, batched=False)),
    rtol=1e-12,
    atol=1e-13,
  )


def test_float32_tracks_float64(
  device: str, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """float32 reproduces the float64 answer to float32 precision.

  Pinned against our own float64 result rather than the compiled library:
  the C++ tolerances describe float64 and do not transfer.
  """
  f64 = TorchVinecop.from_vinecop(cpp_vine, device=torch.device(device))
  f32 = TorchVinecop.from_vinecop(
    cpp_vine, device=torch.device(device), dtype=torch.float32
  )
  a = _np(f64.pdf(torch.as_tensor(u_eval, device=device)))
  b = _np(f32.pdf(torch.as_tensor(u_eval, device=device, dtype=torch.float32)))
  assert np.isfinite(b).all()
  np.testing.assert_allclose(b, a, rtol=F32_VS_F64_TOL, atol=F32_VS_F64_TOL)


def test_float32_keeps_arguments_inside_the_unit_square(device: str) -> None:
  """The trim survives float32.

  ``1 - 1e-10`` rounds to exactly ``1.0`` in float32, so a bound written as
  that literal would admit the value it exists to exclude -- and ``ndtri(1)``
  is an infinity for whatever refits on the result.
  """
  u = pv.to_pseudo_obs(
    np.random.default_rng(3).multivariate_normal(
      [0, 0], [[1, 0.9], [0.9, 1]], size=1500
    )
  )
  cpp = pv.Bicop.from_data(
    u, controls=pv.FitControlsBicop(family_set=[pv.BicopFamily.tll])
  )
  bc = TorchBicop.from_bicop(
    cpp, device=torch.device(device), dtype=torch.float32
  )
  edge = torch.tensor(
    [[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0], [1.0, 1.0], [0.0, 0.0]],
    dtype=torch.float32,
    device=device,
  )
  for name in ("hfunc1", "hfunc2", "cdf"):
    out = getattr(bc, name)(edge)
    assert (out > 0.0).all(), f"{name} returned 0"
    assert (out < 1.0).all(), f"{name} returned 1"
    assert torch.isfinite(torch.special.ndtri(out)).all()


def test_every_buffer_follows_to_device(
  device: str, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """``.to(device)`` moves the whole object, batched cache included."""
  vine = TorchVinecop.from_vinecop(cpp_vine)
  vine.to(device)
  ut = torch.as_tensor(u_eval, device=device)
  out = vine.pdf(ut, batched=True)
  assert_on_device(vine, device, out, extra=(vine._batched,))


@pytest.mark.parametrize("op", _EVAL_OPS + ("sample",))
def test_evaluation_does_not_round_trip_through_the_host(
  device: str, op: str, cpp_vine: pv.Vinecop, u_eval: np.ndarray
) -> None:
  """An evaluation call must not move data back to the host.

  ``fit`` / ``select`` legitimately do -- Kendall's tau goes through the
  compiled ``wdm`` -- which is why this is scoped to evaluation. The first
  call is untimed: it bakes the batched cache, which reads the structure
  from the compiled extension.
  """
  vine = TorchVinecop.from_vinecop(cpp_vine, device=torch.device(device))
  ut = torch.as_tensor(u_eval, device=device)
  call = (
    (lambda: vine.sample(64, seeds=[1]))
    if op == "sample"
    else (lambda: getattr(vine, op)(ut))
  )
  call()
  with count_transfers(device) as c:
    call()
  c.assert_no_d2h(f"TorchVinecop.{op}")


def test_fit_and_select_run_on_device(device: str) -> None:
  """Fitting and structure selection work with device-resident data."""
  u = _u(4, 600, 11)
  ut = torch.as_tensor(u, device=device)
  ctl = FitControlsTorchVinecop(device=torch.device(device))
  fixed = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.BicopFamily.tll])
  )
  vine = TorchVinecop.from_data(ut, fixed.structure, controls=ctl)
  assert vine.pdf(ut).device.type == torch.device(device).type
  selected = TorchVinecop.from_data(ut, controls=ctl)
  assert selected.pdf(ut).device.type == torch.device(device).type
