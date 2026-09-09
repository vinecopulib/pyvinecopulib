"""Benchmark: the array-agnostic ``VinecopBase`` against the compiled ``Vinecop``.

Both sides use the *same* pair copulas -- a ``Bicop`` with the Gaussian family,
fast and orientation-symmetric -- so what the ratios measure is the Python
orchestration of the tree-by-tree walk, not the pair-copula cost. Each row is one
operation:

- ``pdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` : the evaluation cascades.
- ``fit``    : pair copulas fitted along a fixed structure.
- ``select`` : structure selection plus one fit per edge (the analog of a torch
  ``from_data(structure=None)``).
- ``sample_conditional`` : the two cascade passes conditional sampling makes.

A second table splits the relabeling that ``conditioning_set`` rests on, since it
is the one place this layer hands work back to C++: building a throwaway
independence vine, running the compiled peel through ``Vinecop.reorient``, and
matching the slots back up in Python. It is reported against the cascade passes
it precedes, which is what says whether the round trip is worth replacing with a
scoped ``RVineStructure`` binding.

Each is reported for a continuous vine and for one whose first half of the
variables are discrete, where every discrete edge carries a parallel left-limit
cascade and costs up to two extra pair-copula calls.

Run on a quiet machine: wall-clock ratios are noisy under load.

    uv run --no-sync python scripts/bench_vinecop_base.py
"""

from __future__ import annotations

import copy
import itertools
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

  def sample(
    self,
    n: int,
    *,
    x: Any = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> Any:
    return self._b.sample(n, qrng=qrng, seeds=seeds or [])

  def flip(self) -> "_CppBicopLike":
    return _CppBicopLike(self._b.flip())

  def as_continuous(self) -> "_CppBicopLike":
    # The inverse cascade asks for this; forward what the wrapped copula has.
    return _CppBicopLike(self._b.as_continuous())


def _fit_edge(
  tree: int,
  edge: int,
  u_e: Any,
  x_e: Any,
  var_types: Any = ("c", "c"),
) -> _CppBicopLike:
  return _CppBicopLike(
    pv.Bicop.from_data(
      np.asarray(u_e), controls=_GAUSSIAN, var_types=list(var_types)
    )
  )


def _best(fn: Callable[[], Any], reps: int) -> float:
  best = float("inf")
  for _ in range(reps):
    t0 = time.perf_counter()
    fn()
    best = min(best, time.perf_counter() - t0)
  return best


def _data(
  seed: int, d: int, n: int, n_disc: int = 0
) -> tuple[np.ndarray, list[str]]:
  """Pseudo-observations in the compact layout, plus the variable types.

  The first ``n_disc`` variables are collapsed onto five ordered categories, so
  each contributes an ``F(x)`` column and an ``F(x^-)`` one.
  """
  rng = np.random.default_rng(seed)
  u = pv.to_pseudo_obs(
    rng.standard_normal((n, d)) @ rng.standard_normal((d, d))
  )
  var_types = ["d"] * n_disc + ["c"] * (d - n_disc)
  if n_disc == 0:
    return u, var_types
  bounds = np.linspace(0.04, 0.96, 6)
  values, limits = u.copy(), []
  for j in range(n_disc):
    level = np.searchsorted(bounds[1:-1], u[:, j])
    values[:, j] = bounds[level + 1]
    limits.append(bounds[level])
  return np.column_stack([values, *limits]), var_types


class _ListVinecop(VinecopBase[Any]):
  """A vine over a nested list of pair copulas, declaring its variable types."""

  def __init__(
    self,
    pairs: list[list[Any]],
    structure: pv.RVineStructure,
    var_types: list[str],
  ) -> None:
    self._pairs = pairs
    self._bind_vine(structure, var_types=var_types)

  def get_pair_copula(self, tree: int, edge: int) -> Any:
    return self._pairs[tree][edge]

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> Any:
    return pv.utils.sample_uniform(n, self.d, qrng, seeds)


def _hosts(
  u: np.ndarray, var_types: list[str]
) -> tuple[_ListVinecop, pv.Vinecop]:
  """One selected structure with one set of fits, hosted in both evaluators."""
  cpp = pv.Vinecop.from_data(u, var_types=var_types, controls=_GAUSSIAN_VINE)
  pairs = [
    [_CppBicopLike(cpp.get_pair_copula(t, e)) for e in range(cpp.dim - 1 - t)]
    for t in range(cpp.dim - 1)
  ]
  return _ListVinecop(pairs, cpp.structure, var_types), cpp


_OPS = (
  "pdf",
  "rosenblatt",
  "inverse_rosenblatt",
  "fit",
  "select",
  "sample_conditional",
)


def _timings(
  u: np.ndarray, var_types: list[str], reps: int
) -> dict[str, tuple[float, float]]:
  mine, cpp = _hosts(u, var_types)
  d = cpp.dim
  # `inverse_rosenblatt` consumes independent uniforms, not observations.
  w = u[:, :d]
  # One conditioning point per output row, on the variable at the tail of the
  # order -- the set both samplers condition on by default. Taken from `u`'s own
  # columns so the layout is valid: a discrete conditioner needs its left limit,
  # and `_data` puts variable i's at column d + i.
  tail = int(cpp.structure.order[-1]) - 1
  cond_cols = [tail] + ([d + tail] if var_types[tail] == "d" else [])
  u_cond = u[:, cond_cols]
  cases: dict[str, tuple[Callable[[], Any], Callable[[], Any]]] = {
    "pdf": (lambda: mine.pdf(u), lambda: cpp.pdf(u)),
    "rosenblatt": (
      lambda: mine.rosenblatt(u, randomize_discrete=False),
      lambda: cpp.rosenblatt(u, randomize_discrete=False),
    ),
    "inverse_rosenblatt": (
      lambda: mine.inverse_rosenblatt(w),
      lambda: cpp.inverse_rosenblatt(w),
    ),
    "fit": (
      lambda: _ListVinecop.from_data(
        u, structure=cpp.structure, var_types=var_types, fit_edge=_fit_edge
      ),
      lambda: pv.Vinecop.from_data(
        u,
        structure=cpp.structure,
        var_types=var_types,
        controls=_GAUSSIAN_VINE,
      ),
    ),
    "select": (
      lambda: _ListVinecop.from_data(
        u, var_types=var_types, fit_edge=_fit_edge
      ),
      lambda: pv.Vinecop.from_data(
        u, var_types=var_types, controls=_GAUSSIAN_VINE
      ),
    ),
    "sample_conditional": (
      lambda: mine.sample_conditional(u_cond, seeds=[1]),
      lambda: cpp.sample_conditional(u_cond, seeds=[1]),
    ),
  }
  return {
    op: (_best(cases[op][1], reps) * 1e3, _best(cases[op][0], reps) * 1e3)
    for op in _OPS
  }


def _run(n: int, dims: tuple[int, ...], reps: int) -> None:
  print(f"\nn={n}, reps={reps} (best-of); Gaussian pair copulas")
  head = f"{'d':>3} {'types':>10} " + " ".join(f"{op:>21}" for op in _OPS)
  print(head)
  print(
    f"{'':>3} {'':>10} " + " ".join(f"{'cpp/py (ms)  ratio':>21}" for _ in _OPS)
  )
  print("-" * len(head))
  for d in dims:
    for label, n_disc in (("continuous", 0), ("half disc.", d // 2)):
      u, var_types = _data(0, d, n, n_disc)
      timings = _timings(u, var_types, reps)
      cells = []
      for op in _OPS:
        t_cpp, t_mine = timings[op]
        cells.append(f"{t_cpp:>7.1f}/{t_mine:<7.1f}{t_mine / t_cpp:>5.1f}x")
      print(f"{d:>3} {label:>10} " + " ".join(cells))


def _run_reorient(dims: tuple[int, ...], reps: int) -> None:
  """Split the compiled round trip the ``conditioning_set`` views rest on."""
  from pyvinecopulib.core._reorient import reorientation

  print(f"\nrelabeling breakdown, reps={reps} (best-of)")
  head = (
    f"{'d':>3} {'placeholder (ms)':>18} {'compiled peel (ms)':>20} "
    f"{'slot map (ms)':>15} {'one pdf pass (ms)':>19}"
  )
  print(head)
  print("-" * len(head))
  for d in dims:
    u, var_types = _data(0, d, 2000)
    mine, cpp = _hosts(u, var_types)
    structure = cpp.structure
    order = [int(v) for v in structure.order]
    tail = set(order[-2:])
    cond = None
    for candidate in itertools.combinations(order, 2):
      if set(candidate) == tail:
        continue
      try:
        probe = reorientation(structure, list(candidate))
      except RuntimeError:
        continue
      if not probe.identity:
        cond = list(candidate)
        break
    if cond is None:
      raise RuntimeError(
        f"could not find a non-identity admissible conditioning set for d={d}"
      )
    assert not reorientation(structure, cond).identity
    indep = pv.Bicop()
    store = [[indep] * (d - 1 - t) for t in range(d - 1)]

    def build(structure: Any = structure, store: Any = store) -> Any:
      return pv.Vinecop.from_structure(structure=structure, pair_copulas=store)

    placeholder = build()

    def peel(placeholder: Any = placeholder, cond: Any = cond) -> Any:
      # `reorient` mutates, so work on a fresh copy each rep.
      vc = copy.deepcopy(placeholder)
      vc.reorient(cond)
      return vc.structure

    def whole(structure: Any = structure, cond: Any = cond) -> Any:
      return reorientation(structure, cond)

    t_build = _best(build, reps) * 1e3
    t_peel = _best(peel, reps) * 1e3
    t_whole = _best(whole, reps) * 1e3
    t_pdf = _best(lambda mine=mine, u=u: mine.pdf(u), reps) * 1e3
    print(
      f"{d:>3} {t_build:>18.2f} {t_peel:>20.2f} "
      f"{max(t_whole - t_build - t_peel, 0.0):>15.2f} {t_pdf:>19.2f}"
    )


def main() -> None:
  reps = 3
  dims = (5, 10, 20)
  for n in (1000, 10000):
    _run(n, dims, reps)
  _run_reorient(dims, reps)


if __name__ == "__main__":
  main()
