"""Benchmark: ``Vinecop`` structure selection vs the array-agnostic
``VinecopBase.select`` port.

Both use the *same* pair-copula fits (a ``Bicop`` with the Gaussian family —
fast, and it isolates the Python <-> C++ orchestration overhead of the
selector from the pair-fit cost) and produce the same vine (identical matrix
encoding and reused pairs). We report:

- ``cpp``    : ``pv.Vinecop.from_data(u, gaussian)`` (select + fit).
- ``select`` : ``VinecopBase.select(u, fit_edge)`` — the full selection path
               (select edges, fit each pair once, reorient the fitted pairs
               onto the finalized structure), the analog of a torch
               ``from_data(structure=None)``.

Run on a quiet machine: wall-clock ratios are noisy under load.

    uv run --no-sync python scripts/bench_structure_selection.py
"""

from __future__ import annotations

import time
from typing import Any, Callable, Optional

import numpy as np

import pyvinecopulib as pv
from pyvinecopulib.core import VinecopBase

_GAUSSIAN = pv.FitControlsBicop(family_set=[pv.families.gaussian])
_GAUSSIAN_VINE = pv.FitControlsVinecop(
  family_set=[pv.families.gaussian], trunc_lvl=20, num_threads=1
)


class _CppBicopLike:
  """Adapt a ``Bicop`` to the ``(u, x)`` ``BicopLike`` signature."""

  def __init__(self, bicop: pv.Bicop) -> None:
    self._b = bicop

  def pdf(self, u: Any, x: Any = None) -> Any:
    return self._b.pdf(np.asarray(u))

  def cdf(self, u: Any, x: Any = None) -> Any:
    return self._b.cdf(np.asarray(u))

  def hfunc1(self, u: Any, x: Any = None) -> Any:
    return self._b.hfunc1(np.asarray(u))

  def hfunc2(self, u: Any, x: Any = None) -> Any:
    return self._b.hfunc2(np.asarray(u))

  def hinv1(self, u: Any, x: Any = None) -> Any:
    return self._b.hinv1(np.asarray(u))

  def hinv2(self, u: Any, x: Any = None) -> Any:
    return self._b.hinv2(np.asarray(u))

  def simulate(
    self,
    n: int,
    *,
    x: Any = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> Any:
    return self._b.simulate(n, qrng=qrng, seeds=seeds or [])

  def flip(self) -> "_CppBicopLike":
    return _CppBicopLike(self._b.flip())


def _fit_edge(tree: int, edge: int, u_e: Any, x_e: Any) -> _CppBicopLike:
  return _CppBicopLike(pv.Bicop.from_data(np.asarray(u_e), controls=_GAUSSIAN))


def _best(fn: Callable[[], Any], reps: int) -> float:
  best = float("inf")
  for _ in range(reps):
    t0 = time.perf_counter()
    fn()
    best = min(best, time.perf_counter() - t0)
  return best


def _data(seed: int, d: int, n: int) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return pv.to_pseudo_obs(
    rng.standard_normal((n, d)) @ rng.standard_normal((d, d))
  )


def _run(n: int, dims: tuple[int, ...], reps: int) -> None:
  print(f"\nn={n}, reps={reps} (best-of); Gaussian pair fits")
  header = f"{'d':>3} {'cpp (ms)':>10} {'select (ms)':>12} {'select/cpp':>11}"
  print(header)
  print("-" * len(header))
  for d in dims:
    u = _data(0, d, n)

    def do_cpp(u: np.ndarray = u) -> Any:
      return pv.Vinecop.from_data(u, controls=_GAUSSIAN_VINE)

    def do_select(u: np.ndarray = u) -> Any:
      return VinecopBase.select(u, _fit_edge, tree_criterion="tau")

    t_cpp = _best(do_cpp, reps) * 1e3
    t_sel = _best(do_select, reps) * 1e3
    print(f"{d:>3} {t_cpp:>10.2f} {t_sel:>12.2f} {t_sel / t_cpp:>10.1f}x")


def main() -> None:
  reps = 3
  dims = (5, 8, 12, 16)
  for n in (100, 500, 1000, 5000):
    _run(n, dims, reps)


if __name__ == "__main__":
  main()
