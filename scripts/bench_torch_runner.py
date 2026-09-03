#!/usr/bin/env python
"""Run paired Torch benchmark arms with raw samples, confidence intervals, and memory.

Benchmark drivers call :func:`paired_repetitions` instead of reporting a
stand-alone median.  It alternates the two arms, synchronizes CUDA, and emits
one JSON record per repetition plus a summary suitable for release evidence.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import time
from collections.abc import Callable
from typing import Any

import torch


def _rss_bytes() -> int:
  """Return this process's resident memory in bytes."""
  with open("/proc/self/statm", encoding="utf-8") as statm:
    pages = int(statm.read().split()[1])
  return pages * os.sysconf("SC_PAGE_SIZE")


def paired_repetitions(
  first: Callable[[], Any], second: Callable[[], Any], repeats: int
) -> list[dict[str, float | int]]:
  """Time alternating arms and return raw paired measurements.

  Parameters
  ----------
  first, second : callable
      Comparable benchmark arms.
  repeats : int
      Number of paired repetitions.

  Returns
  -------
  list of dict
      Per-repetition wall time, RSS, and CUDA peak allocation measurements.
  """
  if repeats < 2:
    raise ValueError("repeats must be at least 2 for paired evidence")
  cuda = torch.cuda.is_available()
  rows: list[dict[str, float | int]] = []
  for index in range(repeats):
    arms = (("first", first), ("second", second))
    if index % 2:
      arms = tuple(reversed(arms))
    measured: dict[str, tuple[float, int, int]] = {}
    for name, fn in arms:
      if cuda:
        torch.cuda.synchronize()
        torch.cuda.reset_peak_memory_stats()
      rss_before = _rss_bytes()
      start = time.perf_counter()
      fn()
      if cuda:
        torch.cuda.synchronize()
      measured[name] = (
        1000.0 * (time.perf_counter() - start),
        max(0, _rss_bytes() - rss_before),
        torch.cuda.max_memory_allocated() if cuda else 0,
      )
    rows.append(
      {
        "repeat": index,
        "first_ms": measured["first"][0],
        "second_ms": measured["second"][0],
        "paired_delta_ms": measured["second"][0] - measured["first"][0],
        "first_rss_bytes": measured["first"][1],
        "second_rss_bytes": measured["second"][1],
        "first_cuda_peak_bytes": measured["first"][2],
        "second_cuda_peak_bytes": measured["second"][2],
      }
    )
  return rows


def paired_summary(rows: list[dict[str, float | int]]) -> dict[str, float]:
  """Return a normal-approximation 95% CI for paired time deltas."""
  deltas = [float(row["paired_delta_ms"]) for row in rows]
  mean = statistics.fmean(deltas)
  half_width = 1.96 * statistics.stdev(deltas) / len(deltas) ** 0.5
  return {
    "paired_delta_mean_ms": mean,
    "ci95_low_ms": mean - half_width,
    "ci95_high_ms": mean + half_width,
  }


def main() -> None:
  """Expose a reproducible no-op harness smoke test."""
  parser = argparse.ArgumentParser()
  parser.add_argument("--repeats", type=int, default=10)
  args = parser.parse_args()
  rows = paired_repetitions(lambda: None, lambda: None, args.repeats)
  print(json.dumps({"raw": rows, "summary": paired_summary(rows)}))


if __name__ == "__main__":
  main()
