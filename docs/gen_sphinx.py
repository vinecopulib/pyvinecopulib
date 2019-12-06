"""
Generates documentation for `pyvinecopulib`.
"""

import sys
from os.path import abspath, dirname, isabs, join

import pyvinecopulib

from sphinx_base import gen_main

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


def write_module(f_name, name, verbose):
    if verbose:
        print("Write: {}".format(name))
    with open(f_name, 'w') as f:
        f.write(".. GENERATED FILE DO NOT EDIT\n")
        f.write("\n")
        rst_name = name.replace("_", "\\_")
        f.write("{}\n".format(rst_name))
        f.write("=" * len(rst_name) + "\n")
        f.write("\n")

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
        raise RuntimeError(
            "Please provide an absolute path: {}".format(output_dir))
    index_file = join(output_dir, "index.rst")
    write_module(index_file, "pyvinecopulib", verbose)


def main():

    input_dir = dirname(abspath(__file__))
    # Generate.
    gen_main(
        input_dir=input_dir, strict=False, src_func=write_doc_modules)
    # TODO(eric.cousineau): Do some simple linting if this is run under
    # `bazel test` (e.g. scan for instances of `TemporaryName`, scan for raw
    # C++ types in type signatures, etc).


if __name__ == "__main__":
    main()
