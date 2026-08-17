"""Pre-import compatibility checks for deliberately non-baseline wheels."""

from __future__ import annotations

import ctypes
import ctypes.util
import platform
from pathlib import Path

#: `IsProcessorFeaturePresent` codes (winnt.h).
_PF_AVX2 = 40
_PF_FMA = 39


def _requires_x86_64_v3() -> bool:
  """Report whether this build was compiled for the x86-64-v3 baseline."""
  try:
    from ._build_info import REQUIRES_X86_64_V3
  except ImportError:
    # No stamp: a source tree that was never configured by CMake. The
    # extension cannot be importable either, so there is nothing to guard.
    return False
  return bool(REQUIRES_X86_64_V3)


def _linux_has_v3() -> bool | None:
  """Report AVX2+FMA from ``/proc/cpuinfo``, or ``None`` if unreadable."""
  try:
    text = Path("/proc/cpuinfo").read_text()
  except OSError:
    return None
  for line in text.splitlines():
    key, _, value = line.partition(":")
    if key.strip() == "flags":
      flags = set(value.split())
      return {"avx2", "fma"} <= flags
  return None


def _macos_sysctl_flag(name: bytes) -> bool | None:
  """Read one boolean ``sysctl`` through libc, or ``None`` if unavailable."""
  try:
    libc = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)
    value = ctypes.c_int32(0)
    size = ctypes.c_size_t(ctypes.sizeof(value))
    rc = libc.sysctlbyname(
      name, ctypes.byref(value), ctypes.byref(size), None, ctypes.c_size_t(0)
    )
  except (AttributeError, OSError, TypeError, ValueError):
    # Not macOS, no libc to load, or no such symbol.
    return None
  if rc != 0:
    return None
  return value.value == 1


def _macos_has_v3() -> bool | None:
  """Report AVX2+FMA from ``sysctlbyname``, or ``None`` if unavailable."""
  avx2 = _macos_sysctl_flag(b"hw.optional.avx2_0")
  fma = _macos_sysctl_flag(b"hw.optional.fma")
  if avx2 is None or fma is None:
    return None
  return avx2 and fma


def _windows_has_v3() -> bool | None:
  """Report AVX2+FMA from ``IsProcessorFeaturePresent``, or ``None``."""
  # `ctypes.windll` exists only on Windows, so it is read dynamically rather
  # than referenced -- a static reference does not type-check elsewhere.
  windll = getattr(ctypes, "windll", None)
  if windll is None:
    return None
  try:
    is_present = windll.kernel32.IsProcessorFeaturePresent
  except (AttributeError, OSError):
    return None
  try:
    return bool(is_present(_PF_AVX2)) and bool(is_present(_PF_FMA))
  except OSError:
    return None


def _has_x86_64_v3() -> bool | None:
  """Report AVX2+FMA support, or ``None`` when it cannot be determined."""
  system = platform.system()
  if system == "Linux":
    return _linux_has_v3()
  if system == "Darwin":
    return _macos_has_v3()
  if system == "Windows":
    return _windows_has_v3()
  return None


def require_x86_64_v3() -> None:
  """Raise before loading the extension on an incompatible x86-64 CPU.

  Raises
  ------
  ImportError
      If this build targets the x86-64-v3 baseline and the CPU does not
      provide AVX2 and FMA. A build that does not target the baseline, and a
      CPU whose features cannot be determined, are both left alone.
  """
  if platform.machine().lower() not in {"x86_64", "amd64"}:
    return
  if not _requires_x86_64_v3():
    return
  if _has_x86_64_v3() is not False:
    # Supported, or undeterminable -- refusing on a failed probe would be
    # worse than letting the loader report an unsupported instruction.
    return
  raise ImportError(
    "This pyvinecopulib x86-64 wheel is compiled for the x86-64-v3 "
    "instruction set, and this CPU does not report AVX2 and FMA. Build from "
    "source instead -- `pip install --no-binary pyvinecopulib pyvinecopulib` "
    "-- which targets the x86-64 baseline."
  )
