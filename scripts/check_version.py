#!/usr/bin/env python3
"""Check that every place the version is written agrees.

The source of truth is ``[project].version`` in ``pyproject.toml``. Everything
else must match it: the top ``CHANGELOG.md`` heading, and ``CITATION.cff`` /
``.zenodo.json`` when they carry a version.

Reports every mismatch in one run, rather than stopping at the first.

Usage
-----
    python scripts/check_version.py
    python scripts/check_version.py --released
    python scripts/check_version.py --released --tag v1.0.0
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: ``## 1.0.0 (unreleased)`` or ``## 1.0.0 (2026-08-20)``; the date form is
#: what marks a version as released.
_HEADING = re.compile(
  r"^##\s+(?P<version>\d+\.\d+\.\d+)\s*(?:\((?P<state>[^)]*)\))?\s*$"
)


def _project_version() -> str:
  with (ROOT / "pyproject.toml").open("rb") as f:
    return str(tomllib.load(f)["project"]["version"])


def _changelog_heading() -> tuple[str, str] | None:
  """Return ``(version, state)`` of the newest ``CHANGELOG.md`` heading."""
  for line in (ROOT / "CHANGELOG.md").read_text().splitlines():
    match = _HEADING.match(line)
    if match:
      return match["version"], (match["state"] or "").strip()
  return None


def _citation_version() -> str | None:
  path = ROOT / "CITATION.cff"
  if not path.exists():
    return None
  for line in path.read_text().splitlines():
    if line.startswith("version:"):
      return line.split(":", 1)[1].strip().strip("\"'")
  return None


def _zenodo_version() -> str | None:
  path = ROOT / ".zenodo.json"
  if not path.exists():
    return None
  return json.loads(path.read_text()).get("version")


def main() -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--released",
    action="store_true",
    help="require the changelog heading to be dated, not '(unreleased)'",
  )
  parser.add_argument(
    "--tag", help="require the version to match this tag (e.g. v1.0.0)"
  )
  args = parser.parse_args()

  version = _project_version()
  problems: list[str] = []

  heading = _changelog_heading()
  if heading is None:
    problems.append("CHANGELOG.md has no '## X.Y.Z' heading")
  else:
    changelog_version, state = heading
    if changelog_version != version:
      problems.append(
        f"CHANGELOG.md top heading is {changelog_version}, "
        f"pyproject.toml says {version}"
      )
    elif not state:
      problems.append(
        f"CHANGELOG.md heading '## {changelog_version}' says neither "
        "'(unreleased)' nor a release date, so a released version cannot be "
        "told from an unreleased one"
      )
    elif args.released and state.lower() == "unreleased":
      problems.append(
        f"CHANGELOG.md still marks {changelog_version} as unreleased; "
        "date the heading before tagging"
      )
    elif not args.released and state.lower() != "unreleased":
      # Not fatal: the release pull request dates the heading before the tag
      # exists, so this state is expected for exactly one commit.
      print(
        f"note: CHANGELOG.md marks {changelog_version} as released "
        f"({state}); expected only in a release pull request",
        file=sys.stderr,
      )

  for name, found in (
    ("CITATION.cff", _citation_version()),
    (".zenodo.json", _zenodo_version()),
  ):
    if found is not None and found != version:
      problems.append(f"{name} says {found}, pyproject.toml says {version}")

  if args.tag is not None:
    expected = f"v{version}"
    if args.tag != expected:
      problems.append(f"tag is {args.tag}, pyproject.toml implies {expected}")

  for problem in problems:
    print(f"error: {problem}", file=sys.stderr)
  if problems:
    return 1

  print(f"version {version} is consistent")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
