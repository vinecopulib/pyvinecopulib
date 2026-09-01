"""Pre-import compatibility checks for deliberately non-baseline wheels."""

from __future__ import annotations

import platform
from pathlib import Path


def _linux_flags() -> set[str] | None:
  """Return Linux CPU flags, or ``None`` when the kernel exposes none."""
  try:
    for line in Path("/proc/cpuinfo").read_text().splitlines():
      key, _, value = line.partition(":")
      if key.strip() in {"flags", "Features"}:
        return set(value.split())
  except OSError:
    pass
  return None


def require_x86_64_v3() -> None:
  """Raise before loading the extension on an incompatible x86-64 CPU."""
  if platform.machine().lower() not in {"x86_64", "amd64"}:
    return
  flags = _linux_flags()
  if flags is None:
    return
  required = {"avx2", "fma"}
  missing = required - flags
  if missing:
    names = ", ".join(sorted(missing))
    raise ImportError(
      "This pyvinecopulib x86-64 wheel requires the x86-64-v3 instruction "
      f"set (missing: {names}). Install a compatible source build or run on "
      "a CPU/VM exposing AVX2 and FMA."
    )
