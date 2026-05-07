"""Re-execute example notebooks and inject image metadata for nbsphinx.

Docstring (`src/include/docstr.hpp`) and stub (`src/pyvinecopulib/__init__.pyi`)
generation are now handled by CMake at build time; this script is responsible
only for the notebook flow used by docs/Phase 2.
"""

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


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--examples-dir",
    type=Path,
    default=Path("examples"),
    help="Directory containing example notebooks (default: examples).",
  )
  args = parser.parse_args()

  for file in sorted(args.examples_dir.glob("*.ipynb")):
    subprocess.run(
      [
        "jupyter",
        "nbconvert",
        "--to",
        "notebook",
        "--execute",
        "--inplace",
        str(file),
      ],
      check=True,
    )
  inject_image_metadata(args.examples_dir)


if __name__ == "__main__":
  main()
