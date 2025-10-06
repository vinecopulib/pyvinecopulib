# tests/test_gil_release.py

import platform
import threading
import time
from typing import Callable

import numpy as np
import pytest

import pyvinecopulib as pv

if "musl" in platform.libc_ver()[0]:
  pytest.skip(
    "Skip GIL release timing tests on musllinux (unstable timing)",
    allow_module_level=True,
  )

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


def _python_probe(duration: float = T_PROBE) -> None:
  """Sleep-based probe to detect if Python runs concurrently with C++."""
  time.sleep(duration)


@pytest.mark.flaky(reruns=3)
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
def test_gil_release_parallel(
  request: pytest.FixtureRequest,
  factory: str,
  method_name: str,
  args_factory: Callable[[], tuple[object, ...]],
) -> None:
  """Check that two concurrent calls run faster than sequential ones."""
  obj = request.getfixturevalue(factory)()
  method = getattr(obj, method_name)

  # Warm up
  _ = method(*args_factory())

  # Sequential
  t0 = time.perf_counter()
  method(*args_factory())
  method(*args_factory())
  elapsed_seq = time.perf_counter() - t0

  # Concurrent
  def runner() -> None:
    method(*args_factory())

  t1 = threading.Thread(target=runner)
  t2 = threading.Thread(target=runner)
  t1.start()
  t2.start()
  t0 = time.perf_counter()
  t1.join()
  t2.join()
  elapsed_conc = time.perf_counter() - t0

  # With GIL release, expect at least some improvement
  assert elapsed_conc < 0.9 * elapsed_seq, (
    f"GIL not released for {method_name}: "
    f"concurrent {elapsed_conc:.3f}s vs sequential {elapsed_seq:.3f}s"
  )


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
  """Ensure Python bytecode executes while C++ method runs."""
  obj = request.getfixturevalue(factory)()
  method = getattr(obj, method_name)

  def runner() -> None:
    _ = method(*args_factory())

  # Time cpp alone
  t0 = time.perf_counter()
  _ = method(*args_factory())
  cpp_time = time.perf_counter() - t0

  # Choose probe duration dynamically
  probe_duration = min(max(0.5, cpp_time), 2.0)

  # Time probe alone
  t0 = time.perf_counter()
  _python_probe(probe_duration)
  probe_time = time.perf_counter() - t0

  # Run probe + C++ together
  t_py = threading.Thread(target=_python_probe, args=(probe_duration,))
  t_cpp = threading.Thread(target=runner)
  t_py.start()
  t_cpp.start()

  t0 = time.perf_counter()
  t_cpp.join()
  t_py.join()
  elapsed = time.perf_counter() - t0

  # Expect overlap: elapsed closer to max(probe, cpp) than sum
  assert elapsed < 1.2 * max(probe_time, cpp_time), (
    f"Python blocked by {method_name}: "
    f"elapsed {elapsed:.3f}s vs probe {probe_time:.3f}s + cpp {cpp_time:.3f}s"
  )
