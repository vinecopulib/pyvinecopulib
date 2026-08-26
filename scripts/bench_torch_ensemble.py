#!/usr/bin/env python
"""Bench BatchedVineEnsemble against a Python loop over the same vines.

For each (n, d, M) cell: fit M vines on one dataset, each on a different
sampled R-vine structure -- the random-search shape the ensemble exists
for -- then time `pdf` and `rosenblatt` two ways:

* ``torch_loop``     -- M sequential per-vine calls, stacked;
* ``torch_ensemble`` -- one stacked cascade.

Both arms are timed with and without ``torch.compile`` where --compile
asks for it. The loop arm is the reason this script raises
``torch._dynamo.config.cache_size_limit``: the honest comparison is one
ensemble against M *compiled* vines, and past the default cap of eight a
looped ensemble silently measures the eager path instead. A script may
raise the cap; the library may not.

Outputs a long-format CSV to --output (default stdout), one row per timed
configuration, with a ``# <hardware>`` banner and per-cell progress on
stderr. Columns:

    method, n, d, m_vines, backend, device, dtype, grid_size,
    chunk_size, compile, time_ms, launches, bake_ms, peak_mib

``launches`` is filled only under --count-launches, which profiles one
extra call per row: it perturbs the timings, so it is off by default and
its number is reported beside a median measured without it. ``bake_ms``
is the one-off ensemble construction, on the ensemble rows only.
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
  BatchedVineEnsemble,
  FitControlsTorchBicop,
  FitControlsTorchVinecop,
  TorchVinecop,
)

METHODS = ("pdf", "rosenblatt")

_DTYPES: dict[str, torch.dtype] = {
  "f64": torch.float64,
  "f32": torch.float32,
}


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


def _banner() -> str:
  import os

  parts = [f"torch {torch.__version__}", f"{os.cpu_count()} cpu threads"]
  if torch.cuda.is_available():
    p = torch.cuda.get_device_properties(0)
    parts.insert(
      0,
      f"{p.name} (sm_{p.major}{p.minor}, {p.total_memory / 2**30:.0f} GB)",
    )
  return " | ".join(parts)


def _count_launches(fn, sync=None) -> int:
  """CUDA kernels launched by one call, after warming up.

  The first calls allocate and run one-off initialization, so they are
  excluded; the count is taken per evaluation, which is the unit the
  launch-bound argument is made in.
  """
  from torch.profiler import ProfilerActivity, profile

  for _ in range(3):
    fn()
  if sync is not None:
    sync()
  with profile(activities=[ProfilerActivity.CUDA]) as prof:
    fn()
    if sync is not None:
      sync()
  total = 0
  for k in prof.key_averages():
    device_time = getattr(k, "self_device_time_total", None)
    if device_time is None:
      device_time = getattr(k, "self_cuda_time_total", 0)
    if device_time and device_time > 0:
      total += int(k.count)
  return total


def _fit_vines(
  u: np.ndarray,
  m_vines: int,
  d: int,
  grid_size: int,
  device: str,
  dtype: torch.dtype,
) -> list[TorchVinecop]:
  """M vines on one dataset, each on its own sampled structure."""
  u_t = torch.from_numpy(u).to(device=device, dtype=dtype)
  controls = FitControlsTorchVinecop(
    bicop_controls=FitControlsTorchBicop(grid_size=grid_size),
    device=torch.device(device),
    dtype=dtype,
  )
  return [
    TorchVinecop.from_data(
      u_t, pv.RVineStructure.sample(d, seeds=[k + 1]), controls=controls
    )
    for k in range(m_vines)
  ]


def _loop_call(vines: list[TorchVinecop], method: str, u):
  """The arm the ensemble replaces: M per-vine calls, stacked."""
  return lambda: torch.stack(
    [getattr(v, method)(u, batched=True) for v in vines], 0
  )


def _ens_call(ens: BatchedVineEnsemble, method: str, u):
  """One stacked call over the whole set."""
  return lambda: getattr(ens, method)(u)


def _bench_cell(
  n: int,
  d: int,
  m_vines: int,
  devices: list[str],
  dtypes: list[str],
  grid_sizes: list[int],
  chunk_sizes: list[int],
  compile_modes: list[bool],
  repeats: int,
  seed: int,
  count_launches: bool,
) -> list[dict]:
  u_fit = _simulate(n=n, d=d, seed=seed)
  rng = np.random.default_rng(seed + 1)
  u_eval = rng.uniform(0.05, 0.95, size=(n, d))
  rows: list[dict] = []

  for device in devices:
    sync = torch.cuda.synchronize if device.startswith("cuda") else None
    for dtype_name in dtypes:
      dtype = _DTYPES[dtype_name]
      for g in grid_sizes:
        vines = _fit_vines(u_fit, m_vines, d, g, device, dtype)
        u_t = torch.from_numpy(u_eval).to(device=device, dtype=dtype)

        def row(**kw) -> dict:
          base = {
            "n": n,
            "d": d,
            "m_vines": m_vines,
            "device": device,
            "dtype": dtype_name,
            "grid_size": g,
            "chunk_size": "",
            "compile": "",
            "launches": "",
            "bake_ms": "",
            "peak_mib": "",
          }
          base.update(kw)
          return base

        # ---- the per-vine loop, the arm being replaced ----------------
        for compiled in compile_modes:
          for v in vines:
            v.compile_cascades = compiled
          for method in METHODS:
            call = _loop_call(vines, method, u_t)
            ms = _time_repeats(call, repeats, sync=sync)
            rows.append(
              row(
                method=method,
                backend="torch_loop",
                compile=str(compiled).lower(),
                time_ms=ms,
                launches=(
                  _count_launches(call, sync)
                  if count_launches and sync is not None
                  else ""
                ),
              )
            )
          for v in vines:
            v.compile_cascades = False

        # ---- the ensemble --------------------------------------------
        for chunk in chunk_sizes:
          if sync is not None:
            torch.cuda.reset_peak_memory_stats()
          t0 = time.perf_counter()
          ens = BatchedVineEnsemble(
            vines, chunk_size=None if chunk <= 0 else chunk
          )
          if sync is not None:
            sync()
          bake_ms = 1000.0 * (time.perf_counter() - t0)
          for compiled in compile_modes:
            ens.compile_cascades = compiled
            for method in METHODS:
              call = _ens_call(ens, method, u_t)
              ms = _time_repeats(call, repeats, sync=sync)
              peak = (
                f"{torch.cuda.max_memory_allocated() / 2**20:.1f}"
                if sync is not None
                else ""
              )
              rows.append(
                row(
                  method=method,
                  backend="torch_ensemble",
                  chunk_size=ens.chunk_size,
                  compile=str(compiled).lower(),
                  time_ms=ms,
                  bake_ms=f"{bake_ms:.1f}",
                  peak_mib=peak,
                  launches=(
                    _count_launches(call, sync)
                    if count_launches and sync is not None
                    else ""
                  ),
                )
              )
            ens.compile_cascades = False
          if sync is not None:
            torch.cuda.empty_cache()

        if sync is not None:
          torch.cuda.empty_cache()

  return rows


_FIELDNAMES = [
  "method",
  "n",
  "d",
  "m_vines",
  "backend",
  "device",
  "dtype",
  "grid_size",
  "chunk_size",
  "compile",
  "time_ms",
  "launches",
  "bake_ms",
  "peak_mib",
]


def main() -> None:
  ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  ap.add_argument("--n", default="2000", type=_parse_int_list)
  ap.add_argument("--d", default="9", type=_parse_int_list)
  ap.add_argument(
    "--ensemble",
    default="1,8,64",
    type=_parse_int_list,
    help="Ensemble sizes M to sweep (default: 1,8,64).",
  )
  ap.add_argument(
    "--chunk-sizes",
    default="0",
    type=_parse_int_list,
    help="chunk_size values; 0 means the computed default (default: 0).",
  )
  ap.add_argument("--devices", default="cpu,cuda", type=_parse_str_list)
  ap.add_argument("--dtypes", default="f64", type=_parse_str_list)
  ap.add_argument("--grid-sizes", default="30", type=_parse_int_list)
  ap.add_argument(
    "--compile",
    default="false",
    type=_parse_bool_list,
    help="compile_cascades values to sweep (default: false). Each compiled "
    "cell pays tens of seconds of Inductor.",
  )
  ap.add_argument(
    "--count-launches",
    action="store_true",
    help="Also profile one call per row for its CUDA kernel count. CUDA "
    "only; perturbs nothing, since the timed median is measured apart.",
  )
  ap.add_argument("--repeats", default=5, type=int)
  ap.add_argument("--seed", default=42, type=int)
  ap.add_argument("--output", default="-")
  args = ap.parse_args()

  unknown = [t for t in args.dtypes if t not in _DTYPES]
  if unknown:
    ap.error(f"unknown --dtypes {unknown}; expected {sorted(_DTYPES)}")

  torch.set_num_threads(1)

  if any(args.compile):
    # The loop arm compiles one variant per vine, and torch keeps only
    # `cache_size_limit` (8) per code object before silently falling back
    # to eager -- which is the very thing the ensemble removes, so the
    # baseline has to be measured with the cap raised or the comparison
    # flatters the ensemble. A script may raise it; the library may not.
    torch._dynamo.config.cache_size_limit = 4096
    torch._dynamo.config.accumulated_recompile_limit = 16384

  devices = list(args.devices)
  if "cuda" in devices and not torch.cuda.is_available():
    print(
      "# WARNING: cuda requested but not available; skipping.",
      file=sys.stderr,
    )
    devices = [d for d in devices if not d.startswith("cuda")]

  print(f"# {_banner()}", file=sys.stderr, flush=True)

  out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
  writer = csv.DictWriter(out, fieldnames=_FIELDNAMES)
  writer.writeheader()
  out.flush()
  for n in args.n:
    for d in args.d:
      for m_vines in args.ensemble:
        print(f"# cell n={n} d={d} M={m_vines}", file=sys.stderr, flush=True)
        rows = _bench_cell(
          n=n,
          d=d,
          m_vines=m_vines,
          devices=devices,
          dtypes=args.dtypes,
          grid_sizes=args.grid_sizes,
          chunk_sizes=args.chunk_sizes,
          compile_modes=args.compile,
          repeats=args.repeats,
          seed=args.seed,
          count_launches=args.count_launches,
        )
        for r in rows:
          r["time_ms"] = f"{r['time_ms']:.3f}"
          writer.writerow(r)
          out.flush()
  if args.output != "-":
    out.close()


if __name__ == "__main__":
  main()
