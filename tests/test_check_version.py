"""Tests for the version-consistency check (`scripts/check_version.py`).

The script reads the files beside it, so each case runs a copy of it in a
directory holding the metadata that case describes.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

_SCRIPT = (
  Path(__file__).resolve().parent.parent / "scripts" / "check_version.py"
)


def _check(
  tmp_path: Path,
  version: str,
  heading: str,
  *args: str,
  citation: str | None = None,
) -> subprocess.CompletedProcess[str]:
  (tmp_path / "scripts").mkdir(parents=True)
  shutil.copy(_SCRIPT, tmp_path / "scripts" / "check_version.py")
  (tmp_path / "pyproject.toml").write_text(
    f'[project]\nname = "pyvinecopulib"\nversion = "{version}"\n'
  )
  (tmp_path / "CHANGELOG.md").write_text(
    f"# Changelog\n\n## {heading}\n\n## 1.0.0 (2026-09-18)\n"
  )
  if citation is not None:
    (tmp_path / "CITATION.cff").write_text(f"version: {citation}\n")
  return subprocess.run(
    [sys.executable, str(tmp_path / "scripts" / "check_version.py"), *args],
    capture_output=True,
    text=True,
    check=False,
  )


@pytest.mark.parametrize(
  ("version", "heading", "args"),
  [
    ("1.0.1.dev0", "1.0.1 (unreleased)", ()),
    ("1.0.1.dev3", "1.0.1 (unreleased)", ()),
    ("1.0.1", "1.0.1 (2026-10-01)", ()),
    ("1.0.1", "1.0.1 (2026-10-01)", ("--released", "--tag", "v1.0.1")),
  ],
)
def test_consistent_states_pass(
  tmp_path: Path, version: str, heading: str, args: tuple[str, ...]
) -> None:
  result = _check(tmp_path, version, heading, *args, citation="1.0.1")
  assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
  ("version", "heading", "args", "message"),
  [
    # An open cycle must not build as the release it leads to.
    ("1.0.1", "1.0.1 (unreleased)", (), "must say 1.0.1.devN"),
    # A dated heading is a release, which carries no suffix.
    ("1.0.1.dev0", "1.0.1 (2026-10-01)", (), "must say 1.0.1, not"),
    ("1.0.1.dev0", "1.0.1 (unreleased)", ("--released",), "unreleased"),
    ("1.0.2.dev0", "1.0.1 (unreleased)", (), "top heading is 1.0.1"),
    ("1.0.1rc1", "1.0.1 (unreleased)", (), "neither X.Y.Z nor X.Y.Z.devN"),
    ("1.0.1", "1.0.1 (2026-10-01)", ("--tag", "v1.0.0"), "tag is v1.0.0"),
  ],
)
def test_inconsistent_states_fail(
  tmp_path: Path,
  version: str,
  heading: str,
  args: tuple[str, ...],
  message: str,
) -> None:
  result = _check(tmp_path, version, heading, *args, citation="1.0.1")
  assert result.returncode == 1
  assert message in result.stderr


def test_citation_names_the_release_not_the_dev_build(tmp_path: Path) -> None:
  ok = _check(
    tmp_path / "a", "1.0.1.dev0", "1.0.1 (unreleased)", citation="1.0.1"
  )
  assert ok.returncode == 0, ok.stderr
  bad = _check(
    tmp_path / "b", "1.0.1.dev0", "1.0.1 (unreleased)", citation="1.0.1.dev0"
  )
  assert bad.returncode == 1
  assert "CITATION.cff says 1.0.1.dev0, expected 1.0.1" in bad.stderr
