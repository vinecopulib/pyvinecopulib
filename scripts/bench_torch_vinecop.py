#!/usr/bin/env python
"""Bench TorchVinecop against pv.Vinecop, on evaluation and on fitting.

Two modes, selected with --mode:

* ``eval`` (the default) times the evaluation surface of an
  already-fitted vine -- pdf / rosenblatt / inverse_rosenblatt plus
  sample and cdf, the last two being the inverse cascade wrapped in a
  base draw and in a Monte-Carlo dominance count, so they are where a
  slow inverse shows up in practice. Columns:

      mode, method, n, d, backend, threads, device, cache_integrals,
      batched, compile, grid_type, grid_size, dtype, time_ms

* ``fit`` times the fit itself -- ``TorchVinecop.from_data`` against
  ``pv.Vinecop.from_data`` -- on two arms, ``select`` (no structure
  given, so the skeleton is selected as well) and ``fit`` (pair copulas
  fitted along the C++-selected skeleton). Columns:

      mode, n, d, backend, threads, device, cache_integrals, grid_type,
      grid_size, dtype, structure, time_ms

For each (n, d) cell both modes generate n correlated pseudo-obs of
dimension d from a Gaussian latent, and select one C++ skeleton that
every fixed-structure arm in the cell sits on.

For C++ rows, device / cache_integrals / grid_type are empty (C++ TLL
only supports the Phi-spaced grid) and dtype is always f64, since
``pv.Vinecop.from_data`` is double throughout. For torch rows, threads is
empty: ``torch.set_num_threads(1)`` pins the torch arm so that only the
C++ arm sweeps --threads.

Notes on the fit mode:

* Both sides are given ``trunc_lvl = d - 1`` so neither truncates.
  ``FitControlsVinecop`` does not truncate by default while
  ``FitControlsTorchVinecop`` caps at 20, so past d = 21 the two
  ``select`` arms would grow different numbers of trees and the ratio
  would stop meaning anything. On the fixed-structure torch arm the cap
  is ignored anyway -- truncation comes from ``structure.trunc_lvl``.
* ``cache_integrals`` is a fit-time knob as well as an eval-time one:
  ``TorchBicop`` builds the three prefix tables per fitted pair, so a
  d = 20 vine pays 190 of them.
* Timings are per-arm fault-tolerant: a fit that raises records
  ``time_ms = nan`` and a ``# FAILED`` line on stderr rather than losing
  a sweep that runs for hours.
* Grid memory: the local-likelihood KDE materializes an ``(m^2, n, 2)``
  tensor, so ``--grid-sizes 64`` at n = 5000 in float64 is ~330 MB per
  pair fit before temporaries. Mind an 8 GB card.
* ``--profile-ace`` runs one extra, separately-reported fit under a probe
  that counts and times the host synchronizations a ``tll`` fit pays.
  ``.item()`` occurs in the whole installed package in five places, all
  in ``pyvinecopulib.torch._fit_tll``, so patching ``torch.Tensor.item``
  is exact attribution rather than sampling. Reaching into a private
  module is deliberate: a bench script may do what the library may not,
  the same rule already written beside ``cache_size_limit`` below.

Outputs a long-format CSV (one row per timed configuration) to --output,
default stdout, flushed per row, with a ``# <hardware>`` banner and
``# <mode> cell n=<n> d=<d>`` progress on stderr.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from statistics import median
from typing import Any

import numpy as np
import torch

import pyvinecopulib as pv
from pyvinecopulib.torch import (
  FitControlsTorchBicop,
  FitControlsTorchVinecop,
  TorchVinecop,
  _fit_tll,
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


#: `--dtypes` tokens, spelled out rather than resolved with `getattr`, so a
#: typo is an error instead of a silent float32 run labeled with it.
_DTYPES: dict[str, torch.dtype] = {
  "f64": torch.float64,
  "f32": torch.float32,
}

#: `--structures` tokens: whether the timed fit also selects the skeleton.
_STRUCTURES: tuple[str, ...] = ("select", "fit")

#: Flags whose sensible default depends on `--mode`. An eval cell times
#: microsecond-to-millisecond calls; a fit cell times whole vine fits, which
#: are seconds on cuda and minutes on torch cpu, so the same grid is either
#: too small to resolve or too large to finish. Applied after parsing, so an
#: explicit flag always wins.
_MODE_DEFAULTS: dict[str, dict[str, Any]] = {
  "eval": {
    "n": [500, 2000],
    "d": [5, 10, 20, 40],
    "repeats": 3,
    "cache": [False, True],
    "grid_types": ["normal", "linear"],
  },
  "fit": {
    "n": [1000, 5000],
    "d": [5, 10],
    # One timed fit after the untimed warm-up: `_time_repeats` pays
    # `repeats + 1`, and a multi-second fit has little relative jitter.
    "repeats": 1,
    # `cache_integrals=None` resolves to True, so this is the library
    # default; `--cache false,true` prices the per-pair table build.
    "cache": [True],
    # The linear grid moves only the *storage* grid -- the KDE and the
    # bandwidth search are identical -- so pass both to confirm the fit
    # cost does not move rather than assuming it.
    "grid_types": ["normal"],
  },
}

#: Extra columns `--profile-ace` appends, and their value on a row that has
#: nothing to profile (every C++ row, and any arm whose timing raised).
_PROFILE_FIELDS = ["profiled_ms", "ace_calls", "item_calls", "item_ms"]
_PROFILE_BLANK: dict[str, Any] = dict.fromkeys(_PROFILE_FIELDS, "")


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


def _timed_or_nan(fn, repeats: int, sync=None, label: str = "") -> float:
  """`_time_repeats`, but one bad arm costs a row rather than the run.

  A `tll` fit is an iteration on the data, and a narrow dtype can walk it
  somewhere the arithmetic does not survive; a sweep of hours should not
  lose everything to one such arm. The row is still written, with
  `time_ms` as `nan`, so the gap is visible in the CSV rather than
  inferred from a combination missing from it.
  """
  try:
    return _time_repeats(fn, repeats, sync=sync)
  except Exception as exc:
    # Deliberately broad: the point is that no failure mode of a fit takes
    # the sweep down with it, and every one of them is reported.
    print(f"# FAILED {label}: {exc!r}", file=sys.stderr, flush=True)
    if sync is not None:
      torch.cuda.empty_cache()
    return float("nan")


class _AceProbe:
  """Counts the host syncs a `tll` fit pays inside `_ace`, and times them.

  `_fit_tll` reads a scalar back to the host in three places on the
  bandwidth-selection path: once in `_select_bandwidth_constant`, once in
  `_pairwise_mcor`, and once per iteration of each of `_ace`'s two
  convergence loops. `_ace` is entered once per pair-copula fit, so
  `ace_calls` is the pair count and `item_calls - 2 * ace_calls` is the
  total number of ACE iterations -- which is what separates "a few
  expensive syncs" from "thousands of cheap ones", the two readings of
  the profile note recorded under #305.

  On cuda `.item()` blocks until every kernel enqueued so far has run, so
  `item_ms` is the host stall: an *upper* bound on what removing the
  syncs could save, since part of it is device work that had to happen
  anyway. On cpu there is nothing to drain, which is what makes the two
  devices' `item_ms` worth reading side by side.

  It does not see the other host round trips a pair fit makes -- the
  `values < 0` precondition on the interpolation grid, and the margin
  renormalization, which runs on the host for grids this small -- so
  `item_ms` is the ACE loops' share specifically, not the total.

  Both patched attributes are restored on exit.
  """

  def __init__(self) -> None:
    self.item_calls = 0
    self.item_ms = 0.0
    self.ace_calls = 0

  def __enter__(self) -> "_AceProbe":
    probe = self
    self._orig_item = torch.Tensor.item
    self._orig_ace = _fit_tll._ace

    def timed_item(tensor):
      t0 = time.perf_counter()
      out = probe._orig_item(tensor)
      probe.item_ms += 1000.0 * (time.perf_counter() - t0)
      probe.item_calls += 1
      return out

    def counted_ace(*args, **kwargs):
      probe.ace_calls += 1
      return probe._orig_ace(*args, **kwargs)

    torch.Tensor.item = timed_item
    _fit_tll._ace = counted_ace
    return self

  def __exit__(self, *exc) -> None:
    torch.Tensor.item = self._orig_item
    _fit_tll._ace = self._orig_ace


def _profile_once(fn, sync=None) -> dict:
  """One extra fit under the probe, reported beside the timed median."""
  if sync is not None:
    sync()
  with _AceProbe() as probe:
    t0 = time.perf_counter()
    fn()
    if sync is not None:
      sync()
    total = 1000.0 * (time.perf_counter() - t0)
  return {
    "profiled_ms": f"{total:.3f}",
    "ace_calls": probe.ace_calls,
    "item_calls": probe.item_calls,
    "item_ms": f"{probe.item_ms:.3f}",
  }


def _banner() -> str:
  """The hardware / version line a PR table has to carry."""
  parts = [f"torch {torch.__version__}", f"{os.cpu_count()} cpu threads"]
  if torch.cuda.is_available():
    p = torch.cuda.get_device_properties(0)
    parts.insert(
      0,
      f"{p.name} (sm_{p.major}{p.minor}, {p.total_memory / 2**30:.0f} GB)",
    )
  return " | ".join(parts)


def _cpp_controls(
  threads: int, grid_size: int, trunc_lvl: int
) -> pv.FitControlsVinecop:
  return pv.FitControlsVinecop(
    family_set=[pv.families.tll],
    num_threads=threads,
    nonparametric_grid_size=grid_size,
    trunc_lvl=trunc_lvl,
  )


def _torch_controls(
  grid_type: str,
  grid_size: int,
  cache: bool,
  trunc_lvl: int,
  device: str,
  dtype: torch.dtype,
) -> FitControlsTorchVinecop:
  return FitControlsTorchVinecop(
    bicop_controls=FitControlsTorchBicop(
      grid_type=grid_type, grid_size=grid_size
    ),
    cache_integrals=cache,
    trunc_lvl=trunc_lvl,
    device=torch.device(device),
    dtype=dtype,
  )


def _reference_vine(
  u_fit: np.ndarray, grid_size: int, trunc_lvl: int
) -> pv.Vinecop:
  """The C++-selected vine every fixed-structure arm in a cell sits on.

  Selected once per cell, single-threaded, at the first `--grid-sizes`
  entry. Tree 0's edge weights come from the data, but trees 1 and up are
  weighted through the fitted grids' h-functions, so selecting per grid
  size would move the skeleton under the grid-size axis instead of
  leaving it the one thing every row in the cell shares.
  """
  return pv.Vinecop.from_data(
    u_fit, controls=_cpp_controls(1, grid_size, trunc_lvl)
  )


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


def _cpp_fit_call(u_fit, structure, ctl):
  """The C++ vine fit to time; `structure=None` means the `select` arm."""
  if structure is None:
    return lambda: pv.Vinecop.from_data(u_fit, controls=ctl)
  return lambda: pv.Vinecop.from_data(u_fit, structure=structure, controls=ctl)


def _torch_fit_call(u_t, structure, ctl):
  """The torch vine fit to time; `structure=None` means the `select` arm.

  `structure` is `from_data`'s second positional argument: `None` routes
  to the array-agnostic selector, a skeleton routes to the fit engine.
  """
  return lambda: TorchVinecop.from_data(u_t, structure, controls=ctl)


# ---- Mode = eval ------------------------------------------------------


def _bench_eval_cell(
  n: int,
  d: int,
  threads: list[int],
  devices: list[str],
  caches: list[bool],
  batched_modes: list[bool],
  compile_modes: list[bool],
  grid_types: list[str],
  grid_sizes: list[int],
  dtypes: list[str],
  backends: list[str],
  repeats: int,
  seed: int,
  mc: int,
) -> list[dict]:
  u_fit = _simulate(n=n, d=d, seed=seed)
  trunc_lvl = max(d - 1, 1)
  ref = _reference_vine(u_fit, grid_sizes[0], trunc_lvl)

  rng = np.random.default_rng(seed + 1)
  u_eval = rng.uniform(0.05, 0.95, size=(n, d))

  rows: list[dict] = []

  # ---- C++ side --------------------------------------------------------
  if "cpp" in backends:
    # The first grid size's model *is* the selection, so a single-entry
    # sweep pays exactly one C++ fit per cell.
    cpp_by_g = {grid_sizes[0]: ref}
    for g in grid_sizes[1:]:
      cpp_by_g[g] = pv.Vinecop.from_data(
        u_fit,
        structure=ref.structure,
        controls=_cpp_controls(1, g, trunc_lvl),
      )
    for method in METHODS:
      for t in threads:
        for g in grid_sizes:
          ms = _time_repeats(
            _cpp_call(cpp_by_g[g], method, u_eval, n, t, mc), repeats
          )
          rows.append(
            {
              "mode": "eval",
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
              "grid_size": g,
              "dtype": "f64",
              "time_ms": ms,
            }
          )

  # ---- Torch side ------------------------------------------------------
  if "torch" in backends:
    for device in devices:
      sync = torch.cuda.synchronize if device.startswith("cuda") else None
      for grid_type in grid_types:
        for dtype_name in dtypes:
          torch_dtype = _DTYPES[dtype_name]
          for g in grid_sizes:
            for cache in caches:
              # Fit a torch vine on the C++-selected structure so each
              # variant of the storage grid sits on the same skeleton and
              # pair-copula fits, isolating the eval-time effect.
              bc = TorchVinecop.from_data(
                torch.from_numpy(u_fit).to(device),
                ref.structure,
                controls=_torch_controls(
                  grid_type, g, cache, trunc_lvl, device, torch_dtype
                ),
              )
              u_t = torch.from_numpy(u_eval).to(
                device=device, dtype=torch_dtype
              )
              if sync is not None:
                sync()
              for method in METHODS:
                for batched in batched_modes:
                  for compiled in compile_modes:
                    # `compile` wraps the batched cascades, so it is a
                    # no-op without `batched` -- timing it would only
                    # duplicate a row.
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
                        "mode": "eval",
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
                        "grid_size": g,
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


# ---- Mode = fit -------------------------------------------------------


def _bench_fit_cell(
  n: int,
  d: int,
  threads: list[int],
  devices: list[str],
  caches: list[bool],
  grid_types: list[str],
  grid_sizes: list[int],
  dtypes: list[str],
  structures: list[str],
  backends: list[str],
  repeats: int,
  seed: int,
  profile: bool,
) -> list[dict]:
  u_fit = _simulate(n=n, d=d, seed=seed)
  trunc_lvl = max(d - 1, 1)
  ref = _reference_vine(u_fit, grid_sizes[0], trunc_lvl)

  rows: list[dict] = []

  # ---- C++ side --------------------------------------------------------
  if "cpp" in backends:
    for arm in structures:
      for t in threads:
        for g in grid_sizes:
          call = _cpp_fit_call(
            u_fit,
            ref.structure if arm == "fit" else None,
            _cpp_controls(t, g, trunc_lvl),
          )
          row = {
            "mode": "fit",
            "n": n,
            "d": d,
            "backend": "cpp",
            "threads": t,
            "device": "",
            "cache_integrals": "",
            "grid_type": "",
            "grid_size": g,
            "dtype": "f64",
            "structure": arm,
            "time_ms": _timed_or_nan(
              call, repeats, label=f"cpp {arm} t={t} g={g}"
            ),
          }
          if profile:
            row.update(_PROFILE_BLANK)
          rows.append(row)

  # ---- Torch side ------------------------------------------------------
  if "torch" in backends:
    for device in devices:
      sync = torch.cuda.synchronize if device.startswith("cuda") else None
      for dtype_name in dtypes:
        torch_dtype = _DTYPES[dtype_name]
        # Staged outside the timed region: what is measured is the fit,
        # not the host-to-device copy or the cast, and `from_data`'s own
        # `as_tensor` is then a no-op.
        u_t = torch.from_numpy(u_fit).to(device=device, dtype=torch_dtype)
        for grid_type in grid_types:
          for g in grid_sizes:
            for cache in caches:
              ctl = _torch_controls(
                grid_type, g, cache, trunc_lvl, device, torch_dtype
              )
              for arm in structures:
                label = (
                  f"torch {arm} {device} {dtype_name} {grid_type} "
                  f"g={g} cache={cache}"
                )
                call = _torch_fit_call(
                  u_t, ref.structure if arm == "fit" else None, ctl
                )
                ms = _timed_or_nan(call, repeats, sync=sync, label=label)
                row = {
                  "mode": "fit",
                  "n": n,
                  "d": d,
                  "backend": "torch",
                  "threads": "",
                  "device": device,
                  "cache_integrals": str(cache).lower(),
                  "grid_type": grid_type,
                  "grid_size": g,
                  "dtype": dtype_name,
                  "structure": arm,
                  "time_ms": ms,
                }
                if profile:
                  row.update(
                    _PROFILE_BLANK
                    if math.isnan(ms)
                    else _profile_once(call, sync=sync)
                  )
                rows.append(row)
                if sync is not None:
                  torch.cuda.empty_cache()
        del u_t
        if sync is not None:
          torch.cuda.empty_cache()

  return rows


def _build_parser() -> argparse.ArgumentParser:
  ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  ap.add_argument(
    "--mode",
    default="eval",
    choices=("eval", "fit"),
    help="Which thing to time (default: eval).",
  )
  ap.add_argument(
    "--n",
    default=None,
    type=_parse_int_list,
    help="Sample sizes (default: 500,2000 in eval mode, 1000,5000 in fit).",
  )
  ap.add_argument(
    "--d",
    default=None,
    type=_parse_int_list,
    help="Dimensions (default: 5,10,20,40 in eval mode, 5,10 in fit).",
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
    "--backends",
    default="cpp,torch",
    type=_parse_str_list,
    help="Backends to bench (default: cpp,torch). The C++ selection runs "
    "either way -- it is where the fixed skeleton comes from.",
  )
  ap.add_argument(
    "--structures",
    default="select,fit",
    type=_parse_str_list,
    help="Fit arms to sweep: select, fit (default: both). fit mode only.",
  )
  ap.add_argument(
    "--cache",
    default=None,
    type=_parse_bool_list,
    help="cache_integrals values to sweep (default: false,true in eval "
    "mode, true in fit mode).",
  )
  ap.add_argument(
    "--batched",
    default="false,true",
    type=_parse_bool_list,
    help="batched=True/False values to sweep (default: false,true). eval "
    "mode only.",
  )
  ap.add_argument(
    "--compile",
    default="false",
    type=_parse_bool_list,
    help="compile_cascades values to sweep (default: false). Each compiled "
    "cell pays tens of seconds of Inductor; only applies where batched. "
    "eval mode only.",
  )
  ap.add_argument(
    "--grid-types",
    default=None,
    type=_parse_str_list,
    help="TorchBicop grid types (default: normal,linear in eval mode, "
    "normal in fit mode).",
  )
  ap.add_argument(
    "--grid-sizes",
    default="30",
    type=_parse_int_list,
    help="Density-grid sizes to sweep for cpp/torch (default: 30).",
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
    help="Monte-Carlo budget for cdf (default: 10000, the library default). "
    "eval mode only.",
  )
  ap.add_argument(
    "--profile-ace",
    action="store_true",
    help="Add one probed fit per arm, reporting the number and cost of the "
    "host syncs inside the tll bandwidth search. fit mode only.",
  )
  ap.add_argument("--repeats", default=None, type=int)
  ap.add_argument("--seed", default=42, type=int)
  ap.add_argument(
    "--output",
    default="-",
    help="CSV output path; '-' for stdout (default).",
  )
  return ap


#: Per-mode CSV schemas. `csv.DictWriter` raises on an unexpected key, so a
#: mode's rows must match its list exactly.
_FIELDNAMES: dict[str, list[str]] = {
  "eval": [
    "mode",
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
    "grid_size",
    "dtype",
    "time_ms",
  ],
  "fit": [
    "mode",
    "n",
    "d",
    "backend",
    "threads",
    "device",
    "cache_integrals",
    "grid_type",
    "grid_size",
    "dtype",
    "structure",
    "time_ms",
  ],
}


def _validate(ap: argparse.ArgumentParser, args: argparse.Namespace) -> None:
  """Reject what would otherwise be a silently wrong or ignored sweep."""
  for name, value in _MODE_DEFAULTS[args.mode].items():
    if getattr(args, name) is None:
      setattr(args, name, value)

  unknown = [t for t in args.dtypes if t not in _DTYPES]
  if unknown:
    ap.error(f"unknown --dtypes {unknown}; expected {sorted(_DTYPES)}")
  bad_arms = [s for s in args.structures if s not in _STRUCTURES]
  if bad_arms:
    ap.error(f"unknown --structures {bad_arms}; expected {list(_STRUCTURES)}")
  if not args.backends:
    ap.error("--backends must name at least one of cpp, torch")

  # A flag that changes the CSV schema, or that cannot apply, is an error
  # rather than silently ignored.
  if args.profile_ace and args.mode != "fit":
    ap.error("--profile-ace only applies to --mode fit")
  if args.mode == "fit" and any(args.compile):
    ap.error("--compile has no effect on --mode fit; drop it")

  if args.profile_ace:
    with _AceProbe() as probe:
      torch.zeros(1).sum().item()
    if probe.item_calls != 1:
      ap.error("--profile-ace: Tensor.item is not interposable here")


def main() -> None:
  ap = _build_parser()
  args = ap.parse_args()
  _validate(ap, args)

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

  print(f"# {_banner()}", file=sys.stderr, flush=True)

  out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
  fieldnames = list(_FIELDNAMES[args.mode])
  if args.profile_ace:
    fieldnames += _PROFILE_FIELDS
  writer = csv.DictWriter(out, fieldnames=fieldnames)
  writer.writeheader()
  out.flush()
  for n in args.n:
    for d in args.d:
      print(f"# {args.mode} cell n={n} d={d}", file=sys.stderr, flush=True)
      if args.mode == "fit":
        rows = _bench_fit_cell(
          n=n,
          d=d,
          threads=args.threads,
          devices=devices,
          caches=args.cache,
          grid_types=args.grid_types,
          grid_sizes=args.grid_sizes,
          dtypes=args.dtypes,
          structures=args.structures,
          backends=args.backends,
          repeats=args.repeats,
          seed=args.seed,
          profile=args.profile_ace,
        )
      else:
        rows = _bench_eval_cell(
          n=n,
          d=d,
          threads=args.threads,
          devices=devices,
          caches=args.cache,
          batched_modes=args.batched,
          compile_modes=args.compile,
          grid_types=args.grid_types,
          grid_sizes=args.grid_sizes,
          dtypes=args.dtypes,
          backends=args.backends,
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
