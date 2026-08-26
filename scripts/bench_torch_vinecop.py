#!/usr/bin/env python
"""Bench TorchVinecop against pv.Vinecop on the whole evaluation surface.

For each (n, d) cell in the grid:
  - generate n correlated pseudo-obs of dim d (Gaussian latent),
  - fit a pv.Vinecop with the TLL family (single thread),
  - time pv.Vinecop.<method>(u, num_threads=t) for every t in --threads,
  - time TorchVinecop.<method>(u, ...) for every combination of
    --devices x --cache x --batched x --compile,
  - report the median wall-clock per call (in milliseconds).

The methods are pdf / rosenblatt / inverse_rosenblatt plus sample and cdf,
which are the inverse cascade wrapped in a base draw and in a Monte-Carlo
dominance count respectively -- so they are where a slow inverse shows up in
practice.

The TorchVinecop CPU runs are pinned to a single torch thread so they
stay apples-to-apples with the laptop numbers; the C++ side gets swept
across --threads explicitly.

Outputs a long-format CSV (one row per timed configuration) to --output
(default: stdout) with columns:
    method, n, d, backend, threads, device, cache_integrals,
    batched, compile, grid_type, dtype, time_ms

For C++ rows, device / cache_integrals / batched / compile / grid_type
are empty (C++ TLL only supports the Phi-spaced grid). For torch rows,
threads is empty. `--compile true` costs tens of seconds of Inductor per
cell and only applies where `batched` is true, so it is off by default.

The torch vines are built via
``TorchVinecop.from_data(u_fit, structure=...)`` so both grid types
share the same fit path; the structure is the one selected by the C++
``pv.Vinecop.from_data(u_fit)`` call so the comparison stays apples-to-
apples on the structure axis.
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
from pyvinecopulib.torch import (
  FitControlsTorchBicop,
  FitControlsTorchVinecop,
  TorchVinecop,
)


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


def _simulate(n: int, d: int, seed: int) -> np.ndarray:
  rng = np.random.default_rng(seed)
  base = rng.standard_normal(size=(n, 1))
  noise = rng.standard_normal(size=(n, d))
  return pv.to_pseudo_obs(0.6 * base + 0.4 * noise)


def _time_repeats(fn, repeats: int, sync=None) -> float:
  """Median wall-clock per call, in milliseconds. `sync()` is called
  around each timed run so CUDA kernels are fully drained."""
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


#: Timed operations. `sample` takes a draw count instead of data and `cdf` a
#: Monte-Carlo budget, so each backend builds its own call.
METHODS = ("pdf", "rosenblatt", "inverse_rosenblatt", "sample", "cdf")


def _cpp_call(cop, method: str, u, n: int, threads: int, mc: int):
  if method == "sample":
    return lambda: cop.sample(n, num_threads=threads)
  if method == "cdf":
    return lambda: cop.cdf(u, N=mc, num_threads=threads)
  return lambda: getattr(cop, method)(u, num_threads=threads)


def _torch_call(vine, method: str, u, n: int, batched: bool, mc: int):
  if method == "sample":
    return lambda: vine.sample(n, batched=batched)
  if method == "cdf":
    return lambda: vine.cdf(u, N=mc, batched=batched)
  return lambda: getattr(vine, method)(u, batched=batched)


def _bench_cell(
  n: int,
  d: int,
  threads: list[int],
  devices: list[str],
  caches: list[bool],
  batched_modes: list[bool],
  compile_modes: list[bool],
  grid_types: list[str],
  dtypes: list[str],
  repeats: int,
  seed: int,
  mc: int,
) -> list[dict]:
  u_fit = _simulate(n=n, d=d, seed=seed)
  ctl = pv.FitControlsVinecop(family_set=[pv.families.tll], num_threads=1)
  cop = pv.Vinecop.from_data(u_fit, controls=ctl)

  rng = np.random.default_rng(seed + 1)
  u_eval = rng.uniform(0.05, 0.95, size=(n, d))

  rows: list[dict] = []

  # ---- C++ side --------------------------------------------------------
  for method in METHODS:
    for t in threads:
      ms = _time_repeats(_cpp_call(cop, method, u_eval, n, t, mc), repeats)
      rows.append(
        {
          "method": method,
          "n": n,
          "d": d,
          "backend": "cpp",
          "threads": t,
          "device": "",
          "cache_integrals": "",
          "batched": "",
          "compile": "",
          "grid_type": "",
          "dtype": "f64",
          "time_ms": ms,
        }
      )

  # ---- Torch side ------------------------------------------------------
  for device in devices:
    sync = torch.cuda.synchronize if device.startswith("cuda") else None
    for grid_type in grid_types:
      for dtype_name in dtypes:
        torch_dtype = getattr(
          torch, "float64" if dtype_name == "f64" else "float32"
        )
        for cache in caches:
          # Fit a torch vine on the C++-selected structure so each variant
          # of the storage grid sits on the same skeleton and pair-copula
          # fits, isolating the eval-time effect.
          bc = TorchVinecop.from_data(
            torch.from_numpy(u_fit).to(device),
            cop.structure,
            controls=FitControlsTorchVinecop(
              bicop_controls=FitControlsTorchBicop(grid_type=grid_type),
              cache_integrals=cache,
              device=torch.device(device),
              dtype=torch_dtype,
            ),
          )
          u_t = torch.from_numpy(u_eval).to(device=device, dtype=torch_dtype)
          if sync is not None:
            sync()
          for method in METHODS:
            for batched in batched_modes:
              for compiled in compile_modes:
                # `compile` wraps the batched cascades, so it is a no-op
                # without `batched` -- timing it would only duplicate a row.
                if compiled and not batched:
                  continue
                bc.compile_cascades = compiled
                ms = _time_repeats(
                  _torch_call(bc, method, u_t, n, batched, mc),
                  repeats,
                  sync=sync,
                )
                rows.append(
                  {
                    "method": method,
                    "n": n,
                    "d": d,
                    "backend": "torch",
                    "threads": "",
                    "device": device,
                    "cache_integrals": str(cache).lower(),
                    "batched": str(batched).lower(),
                    "compile": str(compiled).lower(),
                    "grid_type": grid_type,
                    "dtype": dtype_name,
                    "time_ms": ms,
                  }
                )
              bc.compile_cascades = False
          # Free GPU memory before building the next variant
          del bc, u_t
          if sync is not None:
            torch.cuda.empty_cache()

  return rows


def main() -> None:
  ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  ap.add_argument("--n", default="500,2000", type=_parse_int_list)
  ap.add_argument("--d", default="5,10,20,40", type=_parse_int_list)
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
    "--cache",
    default="false,true",
    type=_parse_bool_list,
    help="cache_integrals values to sweep (default: false,true).",
  )
  ap.add_argument(
    "--batched",
    default="false,true",
    type=_parse_bool_list,
    help="batched=True/False values to sweep (default: false,true).",
  )
  ap.add_argument(
    "--compile",
    default="false",
    type=_parse_bool_list,
    help="compile_cascades values to sweep (default: false). Each compiled "
    "cell pays tens of seconds of Inductor; only applies where batched.",
  )
  ap.add_argument(
    "--grid-types",
    default="normal,linear",
    type=_parse_str_list,
    help="TorchBicop grid types to sweep (default: normal,linear).",
  )
  ap.add_argument(
    "--dtypes",
    default="f64",
    type=_parse_str_list,
    help="torch dtypes to sweep: f64, f32 (default: f64). The C++ side is "
    "always double, so its rows are reported as f64.",
  )
  ap.add_argument(
    "--mc-samples",
    default=10000,
    type=int,
    help="Monte-Carlo budget for cdf (default: 10000, the library default).",
  )
  ap.add_argument("--repeats", default=3, type=int)
  ap.add_argument("--seed", default=42, type=int)
  ap.add_argument(
    "--output",
    default="-",
    help="CSV output path; '-' for stdout (default).",
  )
  args = ap.parse_args()

  torch.set_num_threads(1)

  if any(args.compile):
    # Each cell compiles a fresh vine, and torch keeps only
    # `cache_size_limit` compiled variants of one code object (8 by default),
    # so a sweep of more than that many cells would silently start measuring
    # the eager path again. A script may raise the cap; the library may not.
    torch._dynamo.config.cache_size_limit = 256
    torch._dynamo.config.accumulated_recompile_limit = 4096

  devices = list(args.devices)
  if "cuda" in devices and not torch.cuda.is_available():
    print(
      "# WARNING: cuda requested but not available; skipping.",
      file=sys.stderr,
    )
    devices = [d for d in devices if not d.startswith("cuda")]

  out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
  fieldnames = [
    "method",
    "n",
    "d",
    "backend",
    "threads",
    "device",
    "cache_integrals",
    "batched",
    "compile",
    "grid_type",
    "dtype",
    "time_ms",
  ]
  writer = csv.DictWriter(out, fieldnames=fieldnames)
  writer.writeheader()
  out.flush()
  for n in args.n:
    for d in args.d:
      print(f"# cell n={n} d={d}", file=sys.stderr, flush=True)
      rows = _bench_cell(
        n=n,
        d=d,
        threads=args.threads,
        devices=devices,
        caches=args.cache,
        batched_modes=args.batched,
        compile_modes=args.compile,
        grid_types=args.grid_types,
        dtypes=args.dtypes,
        repeats=args.repeats,
        seed=args.seed,
        mc=args.mc_samples,
      )
      for row in rows:
        row["time_ms"] = f"{row['time_ms']:.3f}"
        writer.writerow(row)
        out.flush()
  if args.output != "-":
    out.close()


if __name__ == "__main__":
  main()
