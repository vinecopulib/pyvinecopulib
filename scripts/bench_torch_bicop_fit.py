#!/usr/bin/env python
"""Bench TorchBicop / pv.Bicop on the TLL constant family.

Two modes selected via ``--mode``:

* ``fit`` (default) — time the from_data fit:
    - simulate n pseudo-obs via a Gaussian copula,
    - time ``pv.Bicop.from_data`` with TLL family per thread count,
    - time ``TorchBicop.from_data`` per device and grid_type.
  Output columns:
    mode, n, backend, threads, device, grid_type, time_ms

* ``eval`` — fit once, time the eval ops (pdf, cdf, hfunc1, hfunc2,
  hinv1, hinv2) on a separate ``n_eval`` sample:
    - C++ baseline: same fit, no grid_type axis (TLL only has normal).
    - Torch sweep: device x grid_type x cache_integrals.
  Output columns:
    mode, op, n_fit, n_eval, backend, threads, device,
    cache_integrals, grid_type, time_ms

For ``cpp`` rows, ``device`` / ``cache_integrals`` / ``grid_type`` are
empty. For ``torch`` rows, ``threads`` is empty.
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


def _parse_bool_list(s: str) -> list[bool]:
  out: list[bool] = []
  for raw in s.split(","):
    tok = raw.strip().lower()
    if not tok:
      continue
    if tok in ("1", "true", "yes", "y", "t"):
      out.append(True)
    elif tok in ("0", "false", "no", "n", "f"):
      out.append(False)
    else:
      raise argparse.ArgumentTypeError(f"not a bool: {raw!r}")
  return out


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


# --------------------------------------------------------------------------- #
# Mode = fit                                                                    #
# --------------------------------------------------------------------------- #


def _bench_fit(
  n: int,
  threads: list[int],
  devices: list[str],
  grid_types: list[str],
  repeats: int,
  seed: int,
) -> list[dict]:
  u_np = _simulate(n=n, seed=seed)
  rows: list[dict] = []

  for t in threads:
    ctl = pv.FitControlsBicop(family_set=[pv.families.tll], num_threads=t)
    ms = _time_repeats(lambda: pv.Bicop.from_data(u_np, controls=ctl), repeats)
    rows.append(
      {
        "mode": "fit",
        "n": n,
        "backend": "cpp",
        "threads": t,
        "device": "",
        "grid_type": "",
        "time_ms": ms,
      }
    )

  for device in devices:
    sync = torch.cuda.synchronize if device.startswith("cuda") else None
    u_t = torch.from_numpy(u_np).to(device)
    for grid_type in grid_types:
      ms = _time_repeats(
        lambda u=u_t, g=grid_type: TorchBicop.from_data(u, grid_type=g),
        repeats,
        sync=sync,
      )
      rows.append(
        {
          "mode": "fit",
          "n": n,
          "backend": "torch",
          "threads": "",
          "device": device,
          "grid_type": grid_type,
          "time_ms": ms,
        }
      )
  return rows


# --------------------------------------------------------------------------- #
# Mode = eval                                                                  #
# --------------------------------------------------------------------------- #


_EVAL_OPS = ("pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2")


def _bench_eval(
  n_fit: int,
  n_eval: int,
  threads: list[int],
  devices: list[str],
  grid_types: list[str],
  caches: list[bool],
  repeats: int,
  seed: int,
) -> list[dict]:
  u_fit_np = _simulate(n=n_fit, seed=seed)
  rng = np.random.default_rng(seed + 100)
  u_eval_np = rng.uniform(0.05, 0.95, size=(n_eval, 2))
  rows: list[dict] = []

  # C++ baseline: fit once per thread count (the threading affects fit, not
  # eval), then time eval ops. TLL only has the Phi-spaced grid in C++.
  for t in threads:
    ctl = pv.FitControlsBicop(family_set=[pv.families.tll], num_threads=t)
    cop = pv.Bicop.from_data(u_fit_np, controls=ctl)
    for op in _EVAL_OPS:
      cpp_fn = getattr(cop, op)
      ms = _time_repeats(lambda fn=cpp_fn: fn(u_eval_np), repeats)
      rows.append(
        {
          "mode": "eval",
          "op": op,
          "n_fit": n_fit,
          "n_eval": n_eval,
          "backend": "cpp",
          "threads": t,
          "device": "",
          "cache_integrals": "",
          "grid_type": "",
          "time_ms": ms,
        }
      )

  for device in devices:
    sync = torch.cuda.synchronize if device.startswith("cuda") else None
    u_fit_t = torch.from_numpy(u_fit_np).to(device)
    u_eval_t = torch.from_numpy(u_eval_np).to(device)
    for grid_type in grid_types:
      for cache in caches:
        bc = TorchBicop.from_data(
          u_fit_t, grid_type=grid_type, cache_integrals=cache
        )
        for op in _EVAL_OPS:
          torch_fn = getattr(bc, op)
          ms = _time_repeats(
            lambda fn=torch_fn, u=u_eval_t: fn(u), repeats, sync=sync
          )
          rows.append(
            {
              "mode": "eval",
              "op": op,
              "n_fit": n_fit,
              "n_eval": n_eval,
              "backend": "torch",
              "threads": "",
              "device": device,
              "cache_integrals": str(cache).lower(),
              "grid_type": grid_type,
              "time_ms": ms,
            }
          )
        del bc
        if sync is not None:
          torch.cuda.empty_cache()
  return rows


# --------------------------------------------------------------------------- #
# Driver                                                                        #
# --------------------------------------------------------------------------- #


def main() -> None:
  ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  ap.add_argument(
    "--mode",
    default="fit",
    choices=("fit", "eval"),
    help="Which thing to time (default: fit).",
  )
  ap.add_argument("--n", default="500,2000,10000", type=_parse_int_list)
  ap.add_argument(
    "--n-eval",
    default=10000,
    type=int,
    help="Number of eval-time queries (eval mode only; default 10000).",
  )
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
  ap.add_argument(
    "--grid-types",
    default="normal,linear",
    type=_parse_str_list,
    help="TorchBicop grid types to sweep (default: normal,linear).",
  )
  ap.add_argument(
    "--cache",
    default="false,true",
    type=_parse_bool_list,
    help="cache_integrals values to sweep (eval mode only; default false,true).",
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

  if args.mode == "fit":
    fieldnames = [
      "mode",
      "n",
      "backend",
      "threads",
      "device",
      "grid_type",
      "time_ms",
    ]
  else:
    fieldnames = [
      "mode",
      "op",
      "n_fit",
      "n_eval",
      "backend",
      "threads",
      "device",
      "cache_integrals",
      "grid_type",
      "time_ms",
    ]

  out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
  writer = csv.DictWriter(out, fieldnames=fieldnames)
  writer.writeheader()
  out.flush()
  for n in args.n:
    print(f"# {args.mode} cell n={n}", file=sys.stderr, flush=True)
    if args.mode == "fit":
      rows = _bench_fit(
        n=n,
        threads=args.threads,
        devices=devices,
        grid_types=args.grid_types,
        repeats=args.repeats,
        seed=args.seed,
      )
    else:
      rows = _bench_eval(
        n_fit=n,
        n_eval=args.n_eval,
        threads=args.threads,
        devices=devices,
        grid_types=args.grid_types,
        caches=args.cache,
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
