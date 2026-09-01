"""Tests for the wheel ISA compatibility guard."""

from __future__ import annotations

import pytest

from pyvinecopulib._cpu import require_x86_64_v3


def test_x86_64_v3_guard_explains_missing_features(monkeypatch) -> None:
  monkeypatch.setattr("pyvinecopulib._cpu.platform.machine", lambda: "x86_64")
  monkeypatch.setattr("pyvinecopulib._cpu._linux_flags", lambda: {"sse2"})
  with pytest.raises(ImportError, match="AVX2 and FMA"):
    require_x86_64_v3()
