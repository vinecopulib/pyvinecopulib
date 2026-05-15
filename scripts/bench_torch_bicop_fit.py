#!/usr/bin/env python
"""Bench TorchBicop.from_data vs pv.Bicop.from_data on the TLL constant fit.

For each ``n`` in the grid:
  - simulate n correlated pseudo-obs via a Gaussian copula,
  - time pv.Bicop.from_data with TLL family for every t in --threads,
  - time TorchBicop.from_data on every device in --devices,
  - report the median wall-clock per fit (in milliseconds).

Outputs a long-format CSV (one row per timed configuration) to --output
(default: stdout) with columns:
    n, backend, threads, device, time_ms

For C++ rows, device is empty. For torch rows, threads is empty.
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
from pyvinecopulib.torch import TorchBicop


def _parse_int_list(s: str) -> list[int]:
  return [int(x) for x in s.split(",") if x.strip()]


def _parse_str_list(s: str) -> list[str]:
  return [x.strip() for x in s.split(",") if x.strip()]


def _simulate(n: int, seed: int) -> np.ndarray:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  return cop.simulate(n, seeds=[seed, seed + 1, seed + 2])


def _time_repeats(fn, repeats: int, sync=None) -> float:
  if sync is not None:
    sync()
  fn()  # warm-up
  if sync is not None:
    sync()
  times = []
  for _ in range(repeats):
    if sync is not None:
      sync()
    t0 = time.perf_counter()
    fn()
    if sync is not None:
      sync()
    times.append(time.perf_counter() - t0)
  return 1000.0 * median(times)


def _bench_cell(
  n: int, threads: list[int], devices: list[str], repeats: int, seed: int
) -> list[dict]:
  u_np = _simulate(n=n, seed=seed)
  rows: list[dict] = []

  for t in threads:
    ctl = pv.FitControlsBicop(family_set=[pv.families.tll], num_threads=t)
    ms = _time_repeats(lambda: pv.Bicop.from_data(u_np, controls=ctl), repeats)
    rows.append(
      {
        "n": n,
        "backend": "cpp",
        "threads": t,
        "device": "",
        "time_ms": ms,
      }
    )

  for device in devices:
    sync = torch.cuda.synchronize if device.startswith("cuda") else None
    u_t = torch.from_numpy(u_np).to(device)
    ms = _time_repeats(
      lambda u=u_t: TorchBicop.from_data(u), repeats, sync=sync
    )
    rows.append(
      {
        "n": n,
        "backend": "torch",
        "threads": "",
        "device": device,
        "time_ms": ms,
      }
    )
  return rows


def main() -> None:
  ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  ap.add_argument("--n", default="500,2000,10000", type=_parse_int_list)
  ap.add_argument(
    "--threads",
    default="1,16",
    type=_parse_int_list,
    help="C++ thread counts to sweep (default: 1,16).",
  )
  ap.add_argument(
    "--devices",
    default="cpu,cuda",
    type=_parse_str_list,
    help="Torch devices to sweep (default: cpu,cuda).",
  )
  ap.add_argument("--repeats", default=3, type=int)
  ap.add_argument("--seed", default=42, type=int)
  ap.add_argument("--output", default="-")
  args = ap.parse_args()

  torch.set_num_threads(1)

  devices = list(args.devices)
  if "cuda" in devices and not torch.cuda.is_available():
    print(
      "# WARNING: cuda requested but not available; skipping.",
      file=sys.stderr,
    )
    devices = [d for d in devices if not d.startswith("cuda")]

  out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
  writer = csv.DictWriter(
    out, fieldnames=["n", "backend", "threads", "device", "time_ms"]
  )
  writer.writeheader()
  out.flush()
  for n in args.n:
    print(f"# cell n={n}", file=sys.stderr, flush=True)
    rows = _bench_cell(
      n=n,
      threads=args.threads,
      devices=devices,
      repeats=args.repeats,
      seed=args.seed,
    )
    for row in rows:
      row["time_ms"] = f"{row['time_ms']:.3f}"
      writer.writerow(row)
      out.flush()
  if args.output != "-":
    out.close()


if __name__ == "__main__":
  main()
