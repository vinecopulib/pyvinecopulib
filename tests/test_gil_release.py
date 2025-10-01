# tests/test_gil_release.py

import threading
import time
from typing import Callable

import numpy as np
import pytest

import pyvinecopulib as pv

# Tunable parameters for "heavy enough" tests
N_BICOP = 20_000
N_VINECOP = 10_000
D_VINECOP = 6
N_KDE = 80_000
N_BUSYWORK = 5_000_000


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


def _python_busywork(n: int = N_BUSYWORK) -> int:
  """Pure Python busy loop, should make progress while C++ holds no GIL."""
  s = 0
  for i in range(n):
    s += i
  return s


@pytest.mark.flaky(reruns=2)
@pytest.mark.parametrize(
  "factory, method_name, args_factory",
  [
    ("bicop_factory", "pdf", lambda: (np.random.rand(N_BICOP, 2),)),
    ("bicop_factory", "fit", lambda: (np.random.rand(N_BICOP, 2),)),
    (
      "vinecop_factory",
      "pdf",
      lambda: (np.random.rand(N_VINECOP, D_VINECOP),),
    ),
    ("vinecop_factory", "simulate", lambda: (500,)),
    ("rvine_factory", "simulate", lambda: (100,)),
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
  t0 = time.perf_counter()
  t1.start()
  t2.start()
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
    ("bicop_factory", "pdf", lambda: (np.random.rand(N_BICOP, 2),)),
    ("bicop_factory", "fit", lambda: (np.random.rand(N_BICOP, 2),)),
    (
      "vinecop_factory",
      "pdf",
      lambda: (np.random.rand(N_VINECOP, D_VINECOP),),
    ),
    ("vinecop_factory", "simulate", lambda: (500,)),
    ("rvine_factory", "simulate", lambda: (100,)),
    ("kde_factory", "fit", lambda: (np.random.normal(0, 1, N_KDE),)),
    ("kde_factory", "quantile", lambda: (np.linspace(1e-3, 1 - 1e-3, N_KDE),)),
  ],
)
@pytest.mark.flaky(reruns=2)
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

  # Time busywork alone
  t0 = time.perf_counter()
  _python_busywork()
  busy_alone = time.perf_counter() - t0

  # Run busywork + C++ together
  t_py = threading.Thread(target=_python_busywork)
  t_cpp = threading.Thread(target=runner)

  t0 = time.perf_counter()
  t_py.start()
  t_cpp.start()
  t_cpp.join()
  t_py.join()
  elapsed = time.perf_counter() - t0

  # Expect overlap: elapsed closer to max(busy, cpp) than sum
  assert elapsed < 1.5 * busy_alone, (
    f"Busywork blocked by {method_name}: "
    f"elapsed {elapsed:.3f}s vs busy_alone {busy_alone:.3f}s"
  )
