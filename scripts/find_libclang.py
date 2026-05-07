"""Discover libclang.so for CMake's docstring-generation step.

Print the absolute path to a libclang shared library on stdout (empty if
none found). Resolution order:
  1. LIBCLANG_PATH env var (explicit override).
  2. The `libclang` PyPI package's bundled shared library (used in
     isolated builds where libclang is in [build-system] requires).
  3. The active conda env at CONDA_PREFIX/lib/libclang.so (used for
     `pip install -e . --no-build-isolation` against a conda dev env).
"""

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


def candidates():
  override = os.environ.get("LIBCLANG_PATH")
  if override and Path(override).is_file():
    yield override

  search_roots = [
    sysconfig.get_paths()["purelib"],
    sysconfig.get_paths()["platlib"],
  ]
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
  for c in candidates():
    if Path(c).is_file():
      sys.stdout.write(c)
      return
  # No path printed -> CMake's existence check will produce a clear error.


if __name__ == "__main__":
  main()
