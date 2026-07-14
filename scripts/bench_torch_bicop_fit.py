#!/usr/bin/env python
"""Bench TorchBicop / pv.Bicop bicop fitters on a shared Gaussian sample.

Four modes selected via ``--mode``:

* ``fit`` (default) — time the from_data fit:
    - simulate n pseudo-obs via a Gaussian copula,
    - time ``pv.Bicop.from_data`` with TLL family per thread count,
    - time ``TorchBicop.from_data`` per device and grid_type (TLL).
  Output columns:
    mode, n, backend, threads, device, grid_type, grid_size, time_ms

* ``eval`` — fit once, time the eval ops (pdf, cdf, hfunc1, hfunc2,
  hinv1, hinv2) on a separate ``n_eval`` sample:
    - C++ baseline: fit at each requested grid_size; no grid_type axis
      (TLL only has the Phi-spaced grid in C++).
    - Torch sweep: device x grid_type x grid_size x cache_integrals.
  Output columns:
    mode, op, n_fit, n_eval, backend, threads, device,
    cache_integrals, grid_type, grid_size, time_ms

* ``hinv`` — fit a TLL bicop once, then compare the three torch inverse
  h-function methods against the C++ reference on ``n_eval`` queries:
    - ``closed_form`` (no-cache, exact conditional-quantile inversion),
    - ``cached`` (``cache_integrals=True``, one bilinear interp),
    - ``itp`` (the pre-vinecopulib#691 vectorized ITP root-finder,
      kept here as a self-contained reference).
  Output columns:
    mode, op, n_fit, n_eval, method, device, grid_size, time_ms,
    max_abs_diff_vs_cpp

* ``cellidx`` — micro-bench the cell-index primitive
  ``_batched_cell_index`` on the fixed grid axis: ``searchsorted`` vs the
  vinecopulib#691 bucket table, over device x grid_type x grid_size x
  batch size. Asserts the two paths return identical indices.
  Output columns:
    mode, device, grid_type, grid_size, n, method, time_ms

For ``cpp`` rows, ``device`` / ``cache_integrals`` / ``grid_type`` are
empty. For ``torch`` rows, ``threads`` is empty.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from statistics import median
from typing import Callable

import numpy as np
import torch

import pyvinecopulib as pv
from pyvinecopulib.torch import FitControlsTorchBicop, TorchBicop
from pyvinecopulib.torch._batched import _batched_cell_index, _build_cell_lookup
from pyvinecopulib.torch._interp import InterpolationGrid2D


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
  grid_sizes: list[int],
  backends: list[str],
  repeats: int,
  seed: int,
) -> list[dict]:
  u_np = _simulate(n=n, seed=seed)
  rows: list[dict] = []

  if "cpp" in backends:
    for t in threads:
      for g in grid_sizes:
        ctl = pv.FitControlsBicop(
          family_set=[pv.families.tll],
          num_threads=t,
          nonparametric_grid_size=g,
        )
        ms = _time_repeats(
          lambda c=ctl: pv.Bicop.from_data(u_np, controls=c), repeats
        )
        rows.append(
          {
            "mode": "fit",
            "n": n,
            "backend": "cpp",
            "threads": t,
            "device": "",
            "grid_type": "",
            "grid_size": g,
            "time_ms": ms,
          }
        )

  if "torch" in backends:
    for device in devices:
      sync = torch.cuda.synchronize if device.startswith("cuda") else None
      u_t = torch.from_numpy(u_np).to(device)
      for grid_type in grid_types:
        for g in grid_sizes:
          ctl_torch = FitControlsTorchBicop(grid_type=grid_type, grid_size=g)
          ms = _time_repeats(
            lambda u=u_t, c=ctl_torch: TorchBicop.from_data(u, c),
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
              "grid_size": g,
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
  grid_sizes: list[int],
  caches: list[bool],
  backends: list[str],
  repeats: int,
  seed: int,
) -> list[dict]:
  u_fit_np = _simulate(n=n_fit, seed=seed)
  rng = np.random.default_rng(seed + 100)
  u_eval_np = rng.uniform(0.05, 0.95, size=(n_eval, 2))
  rows: list[dict] = []

  # C++ baseline: fit once per (thread, grid_size); threading affects fit, not
  # eval. TLL only has the Phi-spaced grid in C++.
  if "cpp" in backends:
    for t in threads:
      for g in grid_sizes:
        ctl = pv.FitControlsBicop(
          family_set=[pv.families.tll],
          num_threads=t,
          nonparametric_grid_size=g,
        )
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
              "grid_size": g,
              "time_ms": ms,
            }
          )

  if "torch" in backends:
    for device in devices:
      sync = torch.cuda.synchronize if device.startswith("cuda") else None
      u_fit_t = torch.from_numpy(u_fit_np).to(device)
      u_eval_t = torch.from_numpy(u_eval_np).to(device)
      for grid_type in grid_types:
        for g in grid_sizes:
          for cache in caches:
            bc = TorchBicop.from_data(
              u_fit_t,
              FitControlsTorchBicop(grid_type=grid_type, grid_size=g),
              cache_integrals=cache,
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
                  "grid_size": g,
                  "time_ms": ms,
                }
              )
            del bc
            if sync is not None:
              torch.cuda.empty_cache()
  return rows


# --------------------------------------------------------------------------- #
# Mode = hinv                                                                  #
# --------------------------------------------------------------------------- #


@torch.no_grad()
def _solve_itp(
  fun: Callable[[torch.Tensor], torch.Tensor],
  x_a: torch.Tensor,
  x_b: torch.Tensor,
  n_iter: int = 30,
  k1: float = 0.2,
  k2: float = 2.0,
  eps: float = 1e-12,
) -> torch.Tensor:
  """The pre-vinecopulib#691 vectorized ITP root-finder (Oliveira &
  Takahashi, 2020), kept here as the reference the closed form replaced."""
  x_a = x_a.clone()
  x_b = x_b.clone()
  fa = fun(x_a)
  fb = fun(x_b)
  swap = fa > 0
  x_a, x_b = torch.where(swap, x_b, x_a), torch.where(swap, x_a, x_b)
  fa, fb = torch.where(swap, fb, fa), torch.where(swap, fa, fb)
  for _ in range(n_iter):
    width = x_b - x_a
    mid = 0.5 * (x_a + x_b)
    denom = (fb - fa).clamp_min(1e-300)
    x_f = (fb * x_a - fa * x_b) / denom
    diff = mid - x_f
    sigma = torch.sign(diff)
    delta = k1 * width.abs().pow(k2)
    x_t = torch.where(delta <= diff.abs(), x_f + sigma * delta, mid)
    r = (0.5 * width.abs() - eps).clamp_min(0.0)
    x_itp = torch.where((x_t - mid).abs() <= r, x_t, mid - sigma * r)
    fx = fun(x_itp)
    below = fx < 0.0
    x_a = torch.where(below, x_itp, x_a)
    x_b = torch.where(below, x_b, x_itp)
    fa = torch.where(below, fx, fa)
    fb = torch.where(below, fb, fx)
  return 0.5 * (x_a + x_b)


def _hinv_itp(bc: TorchBicop, u: torch.Tensor, cond_var: int) -> torch.Tensor:
  """ITP inversion over the on-the-fly h-function (reference for ``hinv``)."""
  if cond_var == 1:
    fixed, p = u[:, 0:1], u[:, 1:2]
  else:
    fixed, p = u[:, 1:2], u[:, 0:1]

  def fun(x: torch.Tensor) -> torch.Tensor:
    u_e = (
      torch.cat([fixed, x], dim=-1)
      if cond_var == 1
      else torch.cat([x, fixed], dim=-1)
    )
    return bc.interp_grid.integrate_1d(u_e, cond_var=cond_var).unsqueeze(-1) - p

  return _solve_itp(fun, torch.zeros_like(p), torch.ones_like(p)).squeeze(-1)


def _bench_hinv(
  n_fit: int,
  n_eval: int,
  devices: list[str],
  grid_sizes: list[int],
  repeats: int,
  seed: int,
) -> list[dict]:
  u_fit_np = _simulate(n=n_fit, seed=seed)
  rng = np.random.default_rng(seed + 100)
  u_eval_np = rng.uniform(0.02, 0.98, size=(n_eval, 2))
  rows: list[dict] = []

  for g in grid_sizes:
    cpp = pv.Bicop.from_data(
      u_fit_np,
      controls=pv.FitControlsBicop(
        family_set=[pv.families.tll], nonparametric_grid_size=g
      ),
      var_types=["c", "c"],
    )
    ref = {"hinv1": cpp.hinv1(u_eval_np), "hinv2": cpp.hinv2(u_eval_np)}
    for device in devices:
      sync = torch.cuda.synchronize if device.startswith("cuda") else None
      u_t = torch.as_tensor(u_eval_np, dtype=torch.float64, device=device)
      bc_plain = TorchBicop.from_bicop(cpp, cache_integrals=False).to(device)
      bc_cached = TorchBicop.from_bicop(cpp, cache_integrals=True).to(device)
      for op in ("hinv1", "hinv2"):
        cv = 1 if op == "hinv1" else 2
        methods = {
          "closed_form": lambda o=op, bp=bc_plain: getattr(bp, o)(u_t),
          "cached": lambda o=op, bc=bc_cached: getattr(bc, o)(u_t),
          "itp": lambda c=cv, bp=bc_plain: _hinv_itp(bp, u_t, c),
        }
        for name, fn in methods.items():
          ms = _time_repeats(fn, repeats, sync=sync)
          diff = float(np.max(np.abs(fn().cpu().numpy() - ref[op])))
          rows.append(
            {
              "mode": "hinv",
              "op": op,
              "n_fit": n_fit,
              "n_eval": n_eval,
              "method": name,
              "device": device,
              "grid_size": g,
              "time_ms": ms,
              "max_abs_diff_vs_cpp": f"{diff:.3e}",
            }
          )
      del bc_plain, bc_cached
      if sync is not None:
        torch.cuda.empty_cache()
  return rows


# --------------------------------------------------------------------------- #
# Mode = cellidx                                                               #
# --------------------------------------------------------------------------- #


def _bench_cellidx(
  n: int,
  devices: list[str],
  grid_types: list[str],
  grid_sizes: list[int],
  repeats: int,
  seed: int,
) -> list[dict]:
  """Micro-bench the fixed-grid cell search: ``searchsorted`` vs the
  vinecopulib#691 bucket table, on each grid geometry. Both paths run with
  ``is_linear=False`` so this isolates searchsorted vs bucket (the linear
  grid's O(1) floor path is a separate, unchanged branch). Asserts index
  equality before timing."""
  rng = np.random.default_rng(seed)
  u_np = rng.uniform(0.0, 1.0, size=(n,))
  rows: list[dict] = []

  for device in devices:
    sync = torch.cuda.synchronize if device.startswith("cuda") else None
    u_t = torch.as_tensor(u_np, dtype=torch.float64, device=device)
    for grid_type in grid_types:
      for g in grid_sizes:
        grid = InterpolationGrid2D.make_grid_points(grid_type, g).to(device)
        cell_lookup, max_advance = _build_cell_lookup(grid)
        # Correctness: bucket path must equal searchsorted exactly.
        ref = _batched_cell_index(grid, u_t)
        got = _batched_cell_index(grid, u_t, False, cell_lookup, max_advance)
        if not torch.equal(got, ref):
          raise AssertionError(
            f"bucket != searchsorted for {grid_type} grid, m={g}"
          )
        methods = {
          "searchsorted": lambda: _batched_cell_index(grid, u_t),
          "bucket": lambda gp=grid, cl=cell_lookup, ma=max_advance: (
            _batched_cell_index(gp, u_t, False, cl, ma)
          ),
        }
        for name, fn in methods.items():
          ms = _time_repeats(fn, repeats, sync=sync)
          rows.append(
            {
              "mode": "cellidx",
              "device": device,
              "grid_type": grid_type,
              "grid_size": g,
              "n": n,
              "method": name,
              "time_ms": ms,
            }
          )
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
    choices=("fit", "eval", "hinv", "cellidx"),
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
    "--grid-sizes",
    default="30,64",
    type=_parse_int_list,
    help="Density-grid sizes to sweep for cpp/torch (default: 30,64).",
  )
  ap.add_argument(
    "--cache",
    default="false,true",
    type=_parse_bool_list,
    help="cache_integrals values to sweep (eval mode only; default false,true).",
  )
  ap.add_argument(
    "--backends",
    default="cpp,torch",
    type=_parse_str_list,
    help="Backends to bench (default: cpp,torch).",
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

  backends = list(args.backends)

  fieldnames_by_mode = {
    "fit": [
      "mode",
      "n",
      "backend",
      "threads",
      "device",
      "grid_type",
      "grid_size",
      "time_ms",
    ],
    "eval": [
      "mode",
      "op",
      "n_fit",
      "n_eval",
      "backend",
      "threads",
      "device",
      "cache_integrals",
      "grid_type",
      "grid_size",
      "time_ms",
    ],
    "hinv": [
      "mode",
      "op",
      "n_fit",
      "n_eval",
      "method",
      "device",
      "grid_size",
      "time_ms",
      "max_abs_diff_vs_cpp",
    ],
    "cellidx": [
      "mode",
      "device",
      "grid_type",
      "grid_size",
      "n",
      "method",
      "time_ms",
    ],
  }
  fieldnames = fieldnames_by_mode[args.mode]

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
        grid_sizes=args.grid_sizes,
        backends=backends,
        repeats=args.repeats,
        seed=args.seed,
      )
    elif args.mode == "eval":
      rows = _bench_eval(
        n_fit=n,
        n_eval=args.n_eval,
        threads=args.threads,
        devices=devices,
        grid_types=args.grid_types,
        grid_sizes=args.grid_sizes,
        caches=args.cache,
        backends=backends,
        repeats=args.repeats,
        seed=args.seed,
      )
    elif args.mode == "hinv":
      rows = _bench_hinv(
        n_fit=n,
        n_eval=args.n_eval,
        devices=devices,
        grid_sizes=args.grid_sizes,
        repeats=args.repeats,
        seed=args.seed,
      )
    else:  # cellidx: ``n`` is the query batch size.
      rows = _bench_cellidx(
        n=n,
        devices=devices,
        grid_types=args.grid_types,
        grid_sizes=args.grid_sizes,
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
