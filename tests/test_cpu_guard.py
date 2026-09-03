"""Tests for the wheel ISA compatibility guard."""

from __future__ import annotations

import pytest

from pyvinecopulib import _cpu
from pyvinecopulib._cpu import require_x86_64_v3


def test_guard_explains_missing_features(monkeypatch) -> None:
  monkeypatch.setattr(_cpu.platform, "machine", lambda: "x86_64")
  monkeypatch.setattr(_cpu, "_requires_x86_64_v3", lambda: True)
  monkeypatch.setattr(_cpu, "_has_x86_64_v3", lambda: False)
  with pytest.raises(ImportError, match="AVX2 and FMA"):
    require_x86_64_v3()


def test_guard_is_silent_for_a_baseline_build(monkeypatch) -> None:
  # `-march=x86-64-v3` is applied only to a cibuildwheel build, so an ordinary
  # source build is correct on a pre-AVX2 CPU and must not be refused -- the
  # message would tell the user to do what they already did.
  monkeypatch.setattr(_cpu.platform, "machine", lambda: "x86_64")
  monkeypatch.setattr(_cpu, "_requires_x86_64_v3", lambda: False)
  monkeypatch.setattr(_cpu, "_has_x86_64_v3", lambda: False)
  require_x86_64_v3()


def test_guard_fails_open_when_the_cpu_cannot_be_probed(monkeypatch) -> None:
  monkeypatch.setattr(_cpu.platform, "machine", lambda: "x86_64")
  monkeypatch.setattr(_cpu, "_requires_x86_64_v3", lambda: True)
  monkeypatch.setattr(_cpu, "_has_x86_64_v3", lambda: None)
  require_x86_64_v3()


def test_guard_ignores_other_architectures(monkeypatch) -> None:
  monkeypatch.setattr(_cpu.platform, "machine", lambda: "arm64")

  def _unreachable() -> bool:
    raise AssertionError("must not probe a non-x86-64 machine")

  monkeypatch.setattr(_cpu, "_requires_x86_64_v3", _unreachable)
  require_x86_64_v3()


@pytest.mark.parametrize("system", ["Linux", "Darwin", "Windows", "SunOS"])
def test_every_platform_has_a_probe_that_returns_a_verdict(
  system, monkeypatch
) -> None:
  # The guard used to read /proc/cpuinfo only, so it was silently inert on the
  # macOS x86-64 and Windows wheels -- two of the four it is meant to protect.
  monkeypatch.setattr(_cpu.platform, "system", lambda: system)
  assert _cpu._has_x86_64_v3() in (True, False, None)


def test_linux_probe_reports_none_without_proc(monkeypatch) -> None:
  def _raise(*_args, **_kwargs):
    raise OSError("no /proc")

  monkeypatch.setattr(_cpu.Path, "read_text", _raise)
  assert _cpu._linux_has_v3() is None


def test_build_stamp_is_a_bool() -> None:
  assert isinstance(_cpu._requires_x86_64_v3(), bool)
