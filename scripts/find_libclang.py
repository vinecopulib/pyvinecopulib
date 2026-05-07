"""Discover libclang.so for CMake's docstring-generation step.

Print the absolute path to a libclang shared library on stdout (empty if
none found). Resolution order:
  1. LIBCLANG_PATH env var (explicit override).
  2. The `libclang` PyPI package's bundled shared library, found in:
     - this Python's site-packages (works for --no-build-isolation), OR
     - any path passed via --search-paths (scikit-build-core's
       CMAKE_PREFIX_PATH, used for isolated builds where libclang lives
       in the build env, NOT in the target Python's site-packages).
  3. The active conda env at CONDA_PREFIX/lib/libclang.so (used for
     `pip install -e . --no-build-isolation` against a conda dev env).
"""

import argparse
import glob
import os
import sys
import sysconfig
from pathlib import Path

PATTERNS = [
  "libclang.so",
  "libclang.so.*",
  "libclang-*.so*",
  "libclang.dylib",
  "libclang.dll",
]


def candidates(extra_paths):
  override = os.environ.get("LIBCLANG_PATH")
  if override and Path(override).is_file():
    yield override

  search_roots = [
    sysconfig.get_paths()["purelib"],
    sysconfig.get_paths()["platlib"],
  ]
  # scikit-build-core's CMAKE_PREFIX_PATH includes the (isolated) build env's
  # site-packages — that's where libclang is installed during `pip install`
  # since libclang is in [build-system] requires.
  for p in extra_paths:
    if not p:
      continue
    search_roots.append(p)
    # Each prefix may itself be a venv-style root with lib/ and Library/bin/.
    search_roots.extend(
      [os.path.join(p, "lib"), os.path.join(p, "Library", "bin")]
    )
  # Active conda env (CONDA_PREFIX is set when the env is activated; sys.prefix
  # is set whenever the interpreter is the env's Python, even without the
  # activate hook having run).
  for root in (os.environ.get("CONDA_PREFIX"), sys.prefix):
    if not root:
      continue
    search_roots.extend(
      [os.path.join(root, "lib"), os.path.join(root, "Library", "bin")]
    )

  for root in search_roots:
    if not root:
      continue
    for pat in PATTERNS:
      yield from glob.glob(os.path.join(root, pat))
      yield from glob.glob(os.path.join(root, "**", pat), recursive=True)


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--search-paths",
    default="",
    help="Semicolon-separated extra search paths (CMAKE_PREFIX_PATH).",
  )
  args = parser.parse_args()
  extra = [p for p in args.search_paths.split(";") if p]

  for c in candidates(extra):
    if Path(c).is_file():
      sys.stdout.write(c)
      return
  # No path printed -> CMake's existence check will produce a clear error.


if __name__ == "__main__":
  main()
