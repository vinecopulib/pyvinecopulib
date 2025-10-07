# tests/test_gil_release.py

import threading
import time
from typing import Callable

import numpy as np
import pytest

import pyvinecopulib as pv

# Tunable parameters for "heavy enough" tests
N_BICOP = 8_000
N_VINECOP = 6_000
N_VINECOP_SIM = 100
N_RVINE_SIM = 200
D_VINECOP = 3
N_KDE = 100_000
T_PROBE = 1.0


@pytest.fixture
def bicop_factory() -> Callable[[], pv.Bicop]:
  def _make() -> pv.Bicop:
    data = np.random.rand(N_BICOP, 2)
    controls = pv.FitControlsBicop(family_set=[pv.tll])
    return pv.Bicop.from_data(data, controls)

  return _make


@pytest.fixture
def vinecop_factory() -> Callable[[], pv.Vinecop]:
  def _make() -> pv.Vinecop:
    data = np.random.rand(N_VINECOP, D_VINECOP)
    controls = pv.FitControlsVinecop(family_set=[pv.tll])
    return pv.Vinecop.from_data(data, controls=controls)

  return _make


@pytest.fixture
def rvine_factory() -> Callable[[], pv.RVineStructure]:
  def _make() -> pv.RVineStructure:
    return pv.RVineStructure.simulate(D_VINECOP)

  return _make


@pytest.fixture
def kde_factory() -> Callable[[], pv.Kde1d]:
  def _make() -> pv.Kde1d:
    kde = pv.Kde1d()
    x = np.random.normal(0, 1, N_KDE)
    kde.fit(x)
    return kde

  return _make


@pytest.mark.flaky(reruns=2)
@pytest.mark.parametrize(
  "factory, method_name, args_factory",
  [
    ("bicop_factory", "hinv1", lambda: (np.random.rand(N_BICOP, 2),)),
    ("bicop_factory", "fit", lambda: (np.random.rand(N_BICOP, 2),)),
    (
      "vinecop_factory",
      "pdf",
      lambda: (np.random.rand(N_VINECOP, D_VINECOP),),
    ),
    ("vinecop_factory", "simulate", lambda: (N_VINECOP_SIM,)),
    ("rvine_factory", "simulate", lambda: (N_RVINE_SIM,)),
    ("kde_factory", "fit", lambda: (np.random.normal(0, 1, N_KDE),)),
    ("kde_factory", "quantile", lambda: (np.linspace(1e-3, 1 - 1e-3, N_KDE),)),
  ],
)
@pytest.mark.flaky(reruns=3)
def test_python_progress_during_cpp(
  request: pytest.FixtureRequest,
  factory: str,
  method_name: str,
  args_factory: Callable[[], tuple[object, ...]],
) -> None:
  """
  Ensure that while a C++ function executes, Python bytecode continues running.
  It simply verifies that a Python thread can make measurable progress while
  the C++ method runs, i.e. the interpreter is not blocked.
  """
  obj = request.getfixturevalue(factory)()
  method = getattr(obj, method_name)

  progress_counter = 0
  ready = threading.Event()
  done = threading.Event()

  def cpp_runner() -> None:
    ready.set()
    _ = method(*args_factory())
    done.set()

  def python_probe() -> None:
    nonlocal progress_counter
    ready.wait()
    start = time.perf_counter()
    while not done.is_set() and (time.perf_counter() - start) < 2.0:
      progress_counter += 1
      time.sleep(0.005)  # yield frequently

  t_cpp = threading.Thread(target=cpp_runner)
  t_py = threading.Thread(target=python_probe)

  t_cpp.start()
  t_py.start()
  t_cpp.join()
  t_py.join()

  # The probe must have looped at least once -> Python wasn't blocked
  assert progress_counter > 0, f"Python thread was blocked during {method_name}"
