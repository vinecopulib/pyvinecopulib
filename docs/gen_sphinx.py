"""
Generates documentation for `pyvinecopulib`.
"""

import argparse
import os
import sys
import tempfile
from os import listdir, mkdir, symlink
from os.path import abspath, dirname, isabs, join
from shutil import rmtree

CLASSES = [
    "BicopFamily",
    "Bicop",
    "FitControlsBicop",
    "Vinecop",
    "FitControlsVinecop",
    "CVineStructure",
    "DVineStructure",
    "RVineStructure"
]

FUNCTIONS = [
    "to_pseudo_obs",
    "simulate_uniform",
    "ghalton",
    "sobol"
]

EXCLUDE = [
]


def get_submodules(name):
    prefix = name + "."
    out = []

    for s_name in sys.modules.keys():
        if s_name in EXCLUDE:
            continue

        if not s_name.startswith(prefix):
            continue
        sub = s_name[len(prefix):]
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
                ("The module `{}` should not exist; instead, only `{}` should "
                 "exist").format(test, name))

    return False


def write_module(f_name, name, version, verbose):
  if verbose:
    print("Write: {}".format(name))
  with open(f_name, "w") as f:
    f.write(".. GENERATED FILE DO NOT EDIT\n")
    f.write("\n")
    rst_name = name.replace("_", "\\_") + " " + version + " Documentation"
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


def write_doc_modules(output_dir, verbose=False):
  if not isabs(output_dir):
    raise RuntimeError("Please provide an absolute path: {}".format(output_dir))
  index_file = join(output_dir, "index.rst")
  import pyvinecopulib as pv

  version = pv.__version__
  write_module(index_file, "pyvinecopulib", version, verbose)


def _die(s):
  print(s, file=sys.stderr)
  exit(1)


def gen_main(input_dir, strict, src_func=None):
  """Main entry point for generation.
  Args:
      input_dir: Directory which contains initial input files.
      strict: Determines if Sphinx warnings should be interpreted as errors.
      src_func: (optional) Callable of form `f(src_dir)` which will introduce
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
    symlink(join(input_dir, f), src_f)
  # Optionally generate additional input files as source.

  if src_func:
    src_func(src_dir)
  print("Generating documentation...")

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
  # Generate.
  gen_main(input_dir=input_dir, strict=False, src_func=write_doc_modules)


if __name__ == "__main__":
    main()
