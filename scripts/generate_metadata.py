import argparse
import subprocess
from pathlib import Path


def generate_docstrings(env_name: str) -> None:
  print("Generating C++ docstrings...")
  print("-------------------")
  print("Conda environment:", env_name)
  clang_lib = Path.home() / f"miniforge3/envs/{env_name}/lib/libclang.so"
  print("Clang library path:", clang_lib)
  headers = subprocess.check_output(
    [
      "find",
      "lib/vinecopulib/include",
      "-regextype",
      "awk",
      "-regex",
      ".*(class|controls|family|structure|stats)\\.(hpp|ipp)",
      "-print",
    ],
    text=True,
  ).splitlines()
  print("Found headers:\n", "\n  ".join(headers))

  cmd = [
    "python3",
    "scripts/generate_docstring.py",
    "-I",
    "lib/vinecopulib/include",
    "-output",
    "src/include/docstr.hpp",
    "-library_file",
    str(clang_lib),
  ] + headers
  subprocess.run(cmd, check=True)


def generate_stubs(env_name: str) -> None:
  print("Generating stub files...")
  print("-------------------")
  env_python = Path.home() / f"miniforge3/envs/{env_name}/bin/python"

  site_packages = subprocess.check_output(
    [str(env_python), "-c", "import site; print(site.getsitepackages()[0])"],
    text=True,
  ).strip()

  cmd = [str(env_python), "scripts/generate_stubs.py", site_packages]
  subprocess.run(cmd, check=True)


def main():
  parser = argparse.ArgumentParser(
    description="Generate documentation and stub files."
  )
  parser.add_argument(
    "--env", type=str, default="pyvinecopulib", help="Conda environment name"
  )
  parser.add_argument(
    "--no-docstrings", action="store_true", help="Skip C++ docstring generation"
  )
  parser.add_argument(
    "--no-stubs", action="store_true", help="Skip stub file generation"
  )

  args = parser.parse_args()

  if not args.no_docstrings:
    generate_docstrings(args.env)

  if not args.no_stubs:
    generate_stubs(args.env)


if __name__ == "__main__":
  main()
