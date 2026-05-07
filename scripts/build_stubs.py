"""CMake POST_BUILD wrapper for stub generation.

Called from CMakeLists.txt after the C++ extension is built. The stub
generator (scripts/generate_stubs.py) needs to import the *fully assembled*
pyvinecopulib package (pure-Python __init__.py + the just-built .so + the
pair_copuladata helper). At POST_BUILD time the .so lives in the CMake
build dir and the Python sources live in src/pyvinecopulib/, so we stage
both into a tempdir and run the generator against that.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--source-pkg-dir",
    required=True,
    type=Path,
    help="Source pyvinecopulib/ directory (contains __init__.py).",
  )
  parser.add_argument(
    "--ext-so",
    required=True,
    type=Path,
    help="Path to the freshly built pyvinecopulib_ext shared library.",
  )
  parser.add_argument(
    "--output",
    required=True,
    type=Path,
    help="Output __init__.pyi path.",
  )
  parser.add_argument(
    "--py-typed",
    type=Path,
    default=None,
    help="Optional path to py.typed marker (defaults to <output dir>/py.typed).",
  )
  args = parser.parse_args()

  if not args.ext_so.exists():
    sys.exit(f"Extension not found at: {args.ext_so}")
  if not args.source_pkg_dir.is_dir():
    sys.exit(f"Source package dir not found: {args.source_pkg_dir}")

  py_typed = args.py_typed or (args.output.parent / "py.typed")

  with tempfile.TemporaryDirectory(prefix="pyvinecopulib-stubs-") as tmp:
    site = Path(tmp)
    pkg = site / "pyvinecopulib"
    # Copy pure-Python sources, skipping any pre-existing build artifacts.
    shutil.copytree(
      args.source_pkg_dir,
      pkg,
      ignore=shutil.ignore_patterns(
        "__init__.pyi", "py.typed", "*.so", "*.pyd", "*.dylib", "__pycache__"
      ),
    )
    # Drop the just-built .so alongside __init__.py.
    shutil.copy2(args.ext_so, pkg / args.ext_so.name)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(site) + os.pathsep + env.get("PYTHONPATH", "")
    gen = Path(__file__).resolve().parent / "generate_stubs.py"
    subprocess.run(
      [
        sys.executable,
        str(gen),
        "--output",
        str(args.output),
        "--py-typed",
        str(py_typed),
      ],
      env=env,
      check=True,
    )


if __name__ == "__main__":
  main()
