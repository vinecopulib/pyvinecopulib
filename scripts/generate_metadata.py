import argparse
import subprocess
from pathlib import Path
import nbformat


def inject_image_metadata(examples_dir: Path) -> None:
  """
  Adds metadata to embedded images in Jupyter notebooks so that
  nbsphinx can extract them during the documentation build.

  Parameters:
      examples_dir (Path): Path to the directory containing the .ipynb files.
  """
  print("Injecting image metadata into notebooks...")
  print("-------------------")

  notebook_paths = list(examples_dir.glob("*.ipynb"))
  if not notebook_paths:
    print(f"No notebooks found in: {examples_dir}")
    return

  for nb_path in notebook_paths:
    nb = nbformat.read(nb_path, as_version=4)
    image_counter = 0

    for cell_idx, cell in enumerate(nb.cells):
      if cell.cell_type != "code":
        continue
      for output_idx, output in enumerate(cell.get("outputs", [])):
        if "image/png" in output.get("data", {}):
          metadata = output.setdefault("metadata", {}).setdefault(
            "image/png", {}
          )
          if "name" not in metadata:
            metadata["name"] = f"{nb_path.stem}_{cell_idx}_{output_idx}.png"
            image_counter += 1

    if image_counter > 0:
      nbformat.write(nb, nb_path)
      print(f"{nb_path}: added metadata to {image_counter} image(s).")
    else:
      print(f"{nb_path}: no updates needed.")


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
  parser.add_argument(
    "--no-examples",
    action="store_true",
    help="Inject image metadata into Jupyter notebooks (for nbsphinx).",
  )

  args = parser.parse_args()

  if not args.no_docstrings:
    generate_docstrings(args.env)

  if not args.no_stubs:
    generate_stubs(args.env)

  if not args.no_examples:
    examples_dir = Path("examples")
    for file in examples_dir.glob("*.ipynb"):
      subprocess.run(
        [
          "jupyter",
          "nbconvert",
          "--to",
          "notebook",
          "--execute",
          "--inplace",
          file,
        ],
        check=True,
      )
    inject_image_metadata(examples_dir)


if __name__ == "__main__":
  main()
