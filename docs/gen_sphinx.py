"""
Generates documentation for `pyvinecopulib`.
"""

import argparse
import importlib.metadata
import os
import sys
import tempfile
from os import listdir, mkdir, symlink
from os.path import abspath, dirname, isabs, join
from shutil import rmtree


def get_description(package_name):
  """Retrieve the Description field from package metadata."""
  try:
    metadata = importlib.metadata.metadata(package_name)
    return metadata.get("Description")
  except importlib.metadata.PackageNotFoundError:
    raise RuntimeError(f"The package '{package_name}' is not installed.")


CLASSES = [
  "BicopFamily",
  "Bicop",
  "FitControlsBicop",
  "Vinecop",
  "FitControlsVinecop",
  "CVineStructure",
  "DVineStructure",
  "RVineStructure",
]

FUNCTIONS = ["to_pseudo_obs", "simulate_uniform", "ghalton", "sobol"]

EXCLUDE = []

EXTRA_FILES = [
  "../examples",
  "../README.md",
  "../CHANGELOG.md",
]


def get_submodules(name):
  prefix = name + "."
  out = []

  for s_name in sys.modules.keys():
    if s_name in EXCLUDE:
      continue

    if not s_name.startswith(prefix):
      continue
    sub = s_name[len(prefix) :]
    # Ensure its an immediate child.

    if "." in sub or sub.startswith("_"):
      continue
    # For some reason, things like `pyvinecopulib.common` has submodules like
    # `inspect`, etc, whose value in `sys.modules` are none. Ignore those.
    # TODO(eric.cousineau): Figure out where these come from, and remove
    # them.

    if sys.modules[s_name] is None:
      continue
    out.append(s_name)

  return sorted(out)


def has_cc_imported_symbols(name):
  # Check for `module_py`.

  if name + "._module_py" in sys.modules:
    return True
  pieces = name.split(".")

  if len(pieces) > 1:
    sub = pieces[-1]
    test = ".".join(pieces[:-1] + ["_{}_py".format(sub)])

    if test in sys.modules:
      raise RuntimeError(
        (
          "The module `{}` should not exist; instead, only `{}` should " "exist"
        ).format(test, name)
      )

  return False


def write_module(f_name, name, version, verbose):
  if verbose:
    print("Write: {}".format(name))
  with open(f_name, "w") as f:
    f.write(".. GENERATED FILE DO NOT EDIT\n")
    f.write("\n")
    # rst_name = name.replace("_", "\\_") + " " + version + " Documentation"
    rst_name = "API Documentation"
    f.write("=" * len(rst_name) + "\n")
    f.write("{}\n".format(rst_name))
    f.write("=" * len(rst_name) + "\n")
    f.write("\n")

    f.write("\nClasses\n")
    f.write("========\n\n")

    f.write(".. toctree::\n")
    f.write("    :maxdepth: 1\n")
    f.write("\n")

    f.write(".. automodule:: {}\n".format(name))
    f.write(".. autosummary:: \n")
    f.write("    :toctree: _generate\n\n")

    for i in CLASSES:
      f.write("    {}\n".format(i))

    f.write("\nFunctions\n")
    f.write("=========\n\n")

    for i in FUNCTIONS:
      f.write(".. autofunction:: {}\n".format(i))


def write_examples(output_dir, verbose=False):
  """
  Generate an examples.rst file that links to all .ipynb files in the input_dir.

  Parameters:
      output_dir (str): Path to the directory where examples.rst will be written.
  """
  if not isabs(output_dir):
    raise RuntimeError("Please provide an absolute path: {}".format(output_dir))
  if not os.path.isdir(output_dir):
    raise ValueError(
      f"The output directory '{output_dir}' does not exist or is not a directory."
    )

  # Find all .ipynb files in the examples/ directory
  examples_dir = join(output_dir, "examples")
  notebooks = [f for f in os.listdir(examples_dir) if f.endswith(".ipynb")]

  if not notebooks:
    raise ValueError(
      f"No .ipynb files were found in the input directory '{examples_dir}'."
    )

  # Sort notebooks alphabetically
  notebooks.sort()

  # Path to the output file
  output_file = os.path.join(output_dir, "examples.rst")

  # Write the examples.rst content
  with open(output_file, "w") as rst_file:
    # Header
    rst_file.write("Examples\n========\n\n")
    rst_file.write(
      "The following example notebooks are included in this documentation:\n\n"
    )

    # Toctree with :titlesonly:
    rst_file.write(".. toctree::\n")
    rst_file.write("   :maxdepth: 1\n")
    rst_file.write("   :titlesonly:\n\n")

    # Add each notebook to the toctree
    for notebook in notebooks:
      notebook_base = os.path.splitext(notebook)[0]  # Remove .ipynb extension
      rst_file.write(f"   examples/{notebook_base}\n")

  if verbose:
    print(f"examples.rst has been written to {output_file}")


def write_doc_modules(output_dir, verbose=False):
  if not isabs(output_dir):
    raise RuntimeError("Please provide an absolute path: {}".format(output_dir))
  features_file = join(output_dir, "features.rst")
  if verbose:
    print("Writing features to: {}".format(features_file))
  import pyvinecopulib as pv

  version = pv.__version__
  write_module(features_file, "pyvinecopulib", version, verbose)


def _die(s):
  print(s, file=sys.stderr)
  exit(1)


def gen_main(input_dir, strict, src_funcs=None, extra_files=None):
  """Main entry point for generation.
  Args:
      input_dir: Directory which contains initial input files.
      strict: Determines if Sphinx warnings should be interpreted as errors.
      src_funcs: (optional) Callable of form `f(src_dir)` which will introduce
          additional source files to `src_dir`.
  """
  parser = argparse.ArgumentParser()
  parser.add_argument(
    "--out_dir",
    type=str,
    required=True,
    help="Output directory. Does not have to exist beforehand.",
  )
  parser.add_argument(
    "--debug",
    action="store_true",
    help="If enabled, leaves intermediate files that are otherwise " "deleted.",
  )
  parser.add_argument(
    "--verbose",
    action="store_true",
    help="If enabled, print more information.",
  )
  args = parser.parse_args()
  out_dir = args.out_dir

  if out_dir == "<test>":
    out_dir = join(os.environ["TEST_TMPDIR"], "doc")

  if not isabs(out_dir):
    _die("--out_dir must be absolute path: {}".format(out_dir))
  # Use default temp directory, handling both `bazel run` and `bazel test`.
  tmp_dir = tempfile.mkdtemp(dir=os.environ.get("TEST_TMPDIR"))
  doctree_dir = join(tmp_dir, "doctrees")
  src_dir = join(tmp_dir, "src")
  # Symlink inputs to src dir (so that we can also generate doc modules).
  mkdir(src_dir)

  for f in listdir(input_dir):
    src_f = join(src_dir, f)
    src_path = join(input_dir, f)
    if args.verbose:
      print("Symlink: {} -> {}".format(src_path, src_f))
    symlink(join(input_dir, f), src_f)

  # Symlink additional files or directories
  if extra_files is not None:
    for extra in extra_files:
      extra_path = abspath(extra)
      extra_dest = join(src_dir, os.path.basename(extra))
      if args.verbose:
        print("Symlink: {} -> {}".format(extra_path, extra_dest))
      if os.path.isdir(extra_path):
        symlink(extra_path, extra_dest)
      elif os.path.isfile(extra_path):
        symlink(extra_path, extra_dest)

  # Optionally generate additional input files as source.
  if src_funcs is not None:
    for src_func in src_funcs:
      src_func(src_dir)

  if strict:
    # Turn warnings into errors; else be quiet.
    warning_args = ["-W", "-N", "-q"]
  else:
    warning_args = [
      "-N",
      "-Q",  # Be very quiet.
      "-T",  # Traceback (for plugin)
    ]

  os.environ["LANG"] = "en_US.UTF-8"
  sphinx_args = (
    [
      "-b",
      "html",  # HTML output.
      "-a",
      "-E",  # Don't use caching.
      "-d",
      doctree_dir,
    ]
    + warning_args
    + [
      src_dir,  # Source dir.
      out_dir,
    ]
  )
  try:
    from sphinx.cmd.build import main as sphinx_main
  except ImportError:
    from sphinx import main as sphinx_main

  sphinx_main(sphinx_args)

  if not args.debug:
    rmtree(tmp_dir)
  else:
    print("DEBUG: Temporary files: {}".format(tmp_dir))


def main():
  input_dir = dirname(abspath(__file__))
  gen_main(
    input_dir=input_dir,
    strict=False,
    src_funcs=[write_doc_modules, write_examples],
    extra_files=EXTRA_FILES,
  )


if __name__ == "__main__":
  main()
