#!/usr/bin/env python
"""Single-pass perf benchmark for pyvinecopulib (bicop + vine fits).

Mirrors examples/05_benchmark.ipynb. Run twice, once per CMake `-march=`
variant, and the companion driver script summarizes the gap.
"""

import argparse
import json
import timeit
from collections import defaultdict

import numpy as np

import pyvinecopulib as pv


def gen_bicop(n: int, seed: int) -> np.ndarray:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  return cop.sample(n, seeds=[seed])


def gen_vine(n: int, d: int, seed: int) -> np.ndarray:
  np.random.seed(seed)
  x = np.random.normal(size=n)[:, None] * np.ones(
    (n, d)
  ) + 0.5 * np.random.normal(size=(n, d))
  return pv.to_pseudo_obs(x)


def bench_once(fn, u) -> float:
  t = timeit.default_timer()
  fn(u)
  return timeit.default_timer() - t


def main() -> None:
  ap = argparse.ArgumentParser()
  ap.add_argument("--label", required=True)
  ap.add_argument("--output", required=True)
  ap.add_argument("--repeats", type=int, default=5)
  ap.add_argument("--n", type=int, default=1000)
  ap.add_argument("--d", type=int, default=5)
  args = ap.parse_args()

  itau_set = pv.families.itau
  tll_set = [pv.families.tll]

  def b_itau(u):
    pv.Bicop.from_data(
      u, controls=pv.FitControlsBicop(family_set=itau_set, num_threads=1)
    )

  def b_itau_par(u):
    pv.Bicop.from_data(
      u,
      controls=pv.FitControlsBicop(
        family_set=itau_set, parametric_method="itau", num_threads=1
      ),
    )

  def b_tll(u):
    pv.Bicop.from_data(
      u, controls=pv.FitControlsBicop(family_set=tll_set, num_threads=1)
    )

  ctrl_v_itau = pv.FitControlsVinecop(family_set=itau_set, num_threads=1)
  ctrl_v_itau_par = pv.FitControlsVinecop(
    family_set=itau_set, parametric_method="itau", num_threads=1
  )
  ctrl_v_tll = pv.FitControlsVinecop(family_set=tll_set, num_threads=1)

  def v_itau(u):
    pv.Vinecop.from_data(u, controls=ctrl_v_itau)

  def v_itau_par(u):
    pv.Vinecop.from_data(u, controls=ctrl_v_itau_par)

  def v_tll(u):
    pv.Vinecop.from_data(u, controls=ctrl_v_tll)

  bicop_fns = {"itau": b_itau, "itau_par": b_itau_par, "tll": b_tll}
  vine_fns = {"itau": v_itau, "itau_par": v_itau_par, "tll": v_tll}

  results: dict[str, dict[str, list[float]]] = {
    "bicop": defaultdict(list),
    "vine": defaultdict(list),
  }
  for seed in range(args.repeats):
    u_b = gen_bicop(args.n, seed)
    u_v = gen_vine(args.n, args.d, seed)
    for name, fn in bicop_fns.items():
      results["bicop"][name].append(bench_once(fn, u_b))
    for name, fn in vine_fns.items():
      results["vine"][name].append(bench_once(fn, u_v))

  out = {
    "label": args.label,
    "n": args.n,
    "d": args.d,
    "repeats": args.repeats,
    "results": {k: dict(v) for k, v in results.items()},
  }
  with open(args.output, "w") as f:
    json.dump(out, f, indent=2)

  print(f"\n=== {args.label} ===")
  for kind in ["bicop", "vine"]:
    print(f"  {kind}:")
    for name, times in results[kind].items():
      ms = np.array(times) * 1000
      print(
        f"    {name:10s} min={ms.min():8.2f} mean={ms.mean():8.2f} "
        f"max={ms.max():8.2f}  (ms, n={len(times)})"
      )


if __name__ == "__main__":
  main()
