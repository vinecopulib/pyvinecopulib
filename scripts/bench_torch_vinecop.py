#!/usr/bin/env python
"""Bench TorchVinecop vs pv.Vinecop on pdf / rosenblatt / inverse_rosenblatt.

For each (n, d) cell in the grid:
  - generate n correlated pseudo-obs of dim d (Gaussian latent),
  - fit a pv.Vinecop with the TLL family (single thread),
  - wrap with TorchVinecop,
  - time pv.Vinecop.<method>(u) vs bc.<method>(u, impl="legacy")
    vs bc.<method>(u, impl="lazy") across `--repeats` runs,
  - report the median wall-clock per call (in milliseconds).

Both torch and pv are pinned to a single thread so the comparison is
apples-to-apples.

Outputs a CSV table to --output (default: stdout) with columns
  method, n, d, t_cpp_ms, t_legacy_ms, t_lazy_ms.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from statistics import median

import numpy as np
import torch

import pyvinecopulib as pv
from pyvinecopulib.torch import TorchVinecop


def _parse_int_list(s: str) -> list[int]:
  return [int(x) for x in s.split(",") if x.strip()]


def _simulate(n: int, d: int, seed: int) -> np.ndarray:
  rng = np.random.default_rng(seed)
  base = rng.standard_normal(size=(n, 1))
  noise = rng.standard_normal(size=(n, d))
  return pv.to_pseudo_obs(0.6 * base + 0.4 * noise)


def _time_repeats(fn, repeats: int) -> float:
  """Median wall-clock per call, in milliseconds."""
  fn()  # warm-up
  times = []
  for _ in range(repeats):
    t0 = time.perf_counter()
    fn()
    times.append(time.perf_counter() - t0)
  return 1000.0 * median(times)


def _bench_cell(
  n: int, d: int, repeats: int, seed: int
) -> list[tuple[str, int, int, float, float, float]]:
  u_fit = _simulate(n=n, d=d, seed=seed)
  ctl = pv.FitControlsVinecop(family_set=[pv.families.tll], num_threads=1)
  cop = pv.Vinecop.from_data(u_fit, controls=ctl)
  bc = TorchVinecop.from_vinecop(cop)

  rng = np.random.default_rng(seed + 1)
  u_eval = rng.uniform(0.05, 0.95, size=(n, d))
  u_t = torch.from_numpy(u_eval)

  rows: list[tuple[str, int, int, float, float, float]] = []
  for name, fn_cpp, fn_torch in [
    ("pdf", cop.pdf, bc.pdf),
    ("rosenblatt", cop.rosenblatt, bc.rosenblatt),
    ("inverse_rosenblatt", cop.inverse_rosenblatt, bc.inverse_rosenblatt),
  ]:
    t_cpp = _time_repeats(lambda: fn_cpp(u_eval), repeats)
    t_legacy = _time_repeats(lambda: fn_torch(u_t, impl="legacy"), repeats)
    t_lazy = _time_repeats(lambda: fn_torch(u_t, impl="lazy"), repeats)
    rows.append((name, n, d, t_cpp, t_legacy, t_lazy))
  return rows


def main() -> None:
  ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  ap.add_argument("--n", default="500,2000,10000", type=_parse_int_list)
  ap.add_argument("--d", default="5,10,20,40", type=_parse_int_list)
  ap.add_argument("--repeats", default=5, type=int)
  ap.add_argument("--seed", default=42, type=int)
  ap.add_argument(
    "--output",
    default="-",
    help="CSV output path; '-' for stdout (default).",
  )
  args = ap.parse_args()

  torch.set_num_threads(1)

  out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
  writer = csv.writer(out)
  writer.writerow(["method", "n", "d", "t_cpp_ms", "t_legacy_ms", "t_lazy_ms"])
  out.flush()
  for n in args.n:
    for d in args.d:
      print(f"# cell n={n} d={d}", file=sys.stderr, flush=True)
      rows = _bench_cell(n=n, d=d, repeats=args.repeats, seed=args.seed)
      for row in rows:
        writer.writerow(
          [
            row[0],
            row[1],
            row[2],
            f"{row[3]:.3f}",
            f"{row[4]:.3f}",
            f"{row[5]:.3f}",
          ]
        )
        out.flush()
  if args.output != "-":
    out.close()


if __name__ == "__main__":
  main()
