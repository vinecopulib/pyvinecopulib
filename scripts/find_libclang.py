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

With ``--print-resource-dir`` the script instead prints a clang *resource
directory* (the dir holding the compiler builtin headers `stddef.h`,
`xmmintrin.h`, …). The pip `libclang` wheel ships **no** such headers, so on
platforms where libclang cannot auto-discover them (musllinux, macOS) the
docstring parse fails on `cmath` / `*mmintrin.h`. Resolution order:
  1. LIBCLANG_RESOURCE_DIR env var (explicit override).
  2. ``<clang|clang-NN|xcrun clang> -print-resource-dir`` (a system clang).
  3. System globs (`/usr/lib/clang/*`, `/usr/lib/llvm*/lib/clang/*`, …).
  4. Relative to the discovered libclang library.
An empty print means "use libclang's defaults" (correct on manylinux/Windows).
"""

import argparse
import glob
import os
import subprocess
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


def _is_resource_dir(d):
  # A usable resource dir has the builtin headers under include/.
  return bool(d) and Path(d, "include", "stddef.h").is_file()


def resource_dir_candidates(libclang_path):
  override = os.environ.get("LIBCLANG_RESOURCE_DIR")
  if override:
    yield override

  # A system clang knows its own resource dir. Try plain `clang`, versioned
  # `clang-NN`, and macOS' `xcrun clang` (AppleClang lives inside the SDK).
  clang_invocations = [["clang"], ["xcrun", "clang"]] + [
    ["clang-{}".format(v)] for v in range(20, 13, -1)
  ]
  for inv in clang_invocations:
    try:
      out = subprocess.run(
        inv + ["-print-resource-dir"],
        capture_output=True,
        text=True,
        timeout=30,
      )
    except (OSError, subprocess.SubprocessError):
      continue
    if out.returncode == 0 and out.stdout.strip():
      yield out.stdout.strip()

  # System install globs (newest first).
  globs = [
    "/usr/lib/clang/*",
    "/usr/lib/llvm*/lib/clang/*",
    "/usr/local/lib/clang/*",
    "/usr/lib64/clang/*",
  ]
  for pat in globs:
    for hit in sorted(glob.glob(pat), reverse=True):
      yield hit

  # Relative to the libclang shared library (some distributions co-locate the
  # resource dir next to the .so).
  if libclang_path:
    libdir = Path(libclang_path).resolve().parent
    for rel in ("clang", "../lib/clang", "../clang"):
      for hit in sorted(glob.glob(str(libdir / rel / "*")), reverse=True):
        yield hit


def find_libclang(extra):
  for c in candidates(extra):
    if Path(c).is_file():
      return c
  return ""


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--search-paths",
    default="",
    help="Semicolon-separated extra search paths (CMAKE_PREFIX_PATH).",
  )
  parser.add_argument(
    "--print-resource-dir",
    action="store_true",
    help="Print a clang resource dir (builtin headers) instead of the "
    "libclang library path. Empty output means 'use libclang defaults'.",
  )
  args = parser.parse_args()
  extra = [p for p in args.search_paths.split(";") if p]

  if args.print_resource_dir:
    libclang_path = find_libclang(extra)
    for d in resource_dir_candidates(libclang_path):
      if _is_resource_dir(d):
        sys.stdout.write(d)
        return
    # Nothing printed -> CMake omits -resource-dir and libclang uses its
    # own defaults (correct where auto-discovery already works).
    return

  path = find_libclang(extra)
  if path:
    sys.stdout.write(path)
  # No path printed -> CMake's existence check will produce a clear error.


if __name__ == "__main__":
  main()
