# -*- coding: utf-8 -*-
#
# pyvinecopulib documentation build configuration file

# Sphinx extension modules

import inspect
import os
import re
from typing import Any

import sphinx.ext.autodoc as autodoc
import sphinx.util.inspect as sphinxinspect
from sphinx.ext.autodoc import AttributeDocumenter, ModuleDocumenter

import pyvinecopulib as pv

# -- Monkey-patch the autodoc module for nanobind compatibility ------------


def isnbfunc(obj: Any) -> bool:
  """Check if the object is nanobind.nb_func."""
  return (
    hasattr(type(obj), "__module__")
    and type(obj).__module__ == "nanobind"
    and type(obj).__name__ == "nb_func"
  )


def isfunction(obj: Any) -> bool:
  """Check if the object is a user-defined function.
  Partial objects are unwrapped before checking them.
  .. seealso:: :external+python:func:`inspect.isfunction`
  """
  return inspect.isfunction(sphinxinspect.unpartial(obj)) or isnbfunc(obj)


def isroutine(obj: Any) -> bool:
  """Check if the object is a kind of function or method.

  Partial objects are unwrapped before checking them.

  .. seealso:: :external+python:func:`inspect.isroutine`
  """
  return inspect.isroutine(sphinxinspect.unpartial(obj)) or isnbfunc(obj)


sphinxinspect.isfunction = isfunction
sphinxinspect.isroutine = isroutine

# show the body of the function
assert autodoc.inspect.isfunction is isfunction
assert autodoc.inspect.isroutine is isroutine


# Define the patched method
@classmethod
def patched_can_document_member(
  cls, member: Any, membername: str, isattr: bool, parent: Any
) -> bool:
  """
  Patched version of AttributeDocumenter's can_document_member.
  """
  if isinstance(parent, ModuleDocumenter):
    return False
  if sphinxinspect.isroutine(member):
    return False  # New behavior: routines are not attributes
  if sphinxinspect.isattributedescriptor(member):
    return True
  return not sphinxinspect.isroutine(member) and not isinstance(member, type)


# Apply the patch
AttributeDocumenter.can_document_member = patched_can_document_member

assert inspect.getsource(
  AttributeDocumenter.can_document_member
) == inspect.getsource(patched_can_document_member)


# -- General configuration ------------------------------------------------

extensions = [
  "sphinx.ext.autodoc",
  "sphinx_autodoc_typehints",
  "sphinx.ext.githubpages",
  "sphinx.ext.mathjax",
  "sphinx_rtd_theme",
  "sphinx.ext.autosummary",
  "sphinx.ext.napoleon",
  "nbsphinx",
  "myst_parser",
]

napoleon_include_init_with_doc = True
napoleon_use_rtype = False
napoleon_custom_sections = [("Usage", "Usage")]
autosummary_generate = True
nbsphinx_execute = "never"

# Warning categories silenced for fresh single-pass builds:
#   - `autosummary`: stub files don't exist when the toctree first reads
#     them; they're generated mid-build, so a clean build always emits this.
#   - `myst.header`: index.rst inlines a slice of README.md (`:start-line: 8`)
#     that starts at an H2, which myst-parser flags. Cosmetic.
suppress_warnings = ["autosummary", "myst.header"]

# The suffix(es) of source filenames.
source_suffix = {
  ".rst": "restructuredtext",
  ".md": "markdown",
}

# For the templates.
templates_path = ["_templates"]

# The master toctree document.
master_doc = "index"

# Don't scan build outputs or jupyter checkpoint dirs as source documents.
# Without `_build` here, sphinx would re-discover its own output (which
# includes converted notebooks) and register every notebook twice. Note
# `_generate` is intentionally NOT excluded — autosummary writes stub
# files there and needs sphinx to read them back as source documents.
exclude_patterns = [
  "_build",
  "**/.ipynb_checkpoints",
]

# General information about the project.
project = "pyvinecopulib"
copyright = "2024, Thomas Nagler and Thibault Vatter"
author = "Thomas Nagler and Thibault Vatter"

# The version info.
release = pv.__version__
version = ".".join(release.split(".")[:3])

# Note: don't set `html_extra_path = ["examples"]` here. nbsphinx processes
# the staged docs/examples/*.ipynb as documents (toctree-targets); putting
# the same directory in html_extra_path would treat it as raw extra files
# and shadow nbsphinx's document registration, leaving the toctree entries
# unresolved.

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"

html_static_path = ["_static"]

html_copy_source = False

html_show_copyright = False

html_show_sphinx = False

add_module_names = False

pygments_style = "sphinx"

html_logo = "_static/pyvinecopulib.png"


def process_cross_references(content: str, is_docstring: bool = True) -> str:
  classes = [
    "BicopFamily",
    "Bicop",
    "Vinecop",
    "CVineStructure",
    "DVineStructure",
    "RVineStructure",
    "FitControlsBicop",
    "FitControlsVinecop",
    "Kde1d",
  ]

  meth_ref = r":meth:`" if is_docstring else r"{py:meth}`"
  cls_ref = r":class:`" if is_docstring else r"{py:class}`"

  # Try convert references to methods or classes to cross-references
  content = re.sub(r"``(\w+)\.(\w+)\(\)``", rf"{meth_ref}\1.\2`", content)
  content = re.sub(r"``(\w+)\.(\w+)``", f"{cls_ref}\1.\2`", content)
  for cls in classes:
    content = re.sub(rf"``{cls}``", rf"{cls_ref}{cls}`", content)

  return content


def autodoc_process_docstring(app, what, name, obj, options, lines):
  # print(f"Processing: {what}, {name}")

  # Join the existing lines and try to reformat the docstring
  docstring = "\n".join(lines)

  # Process cross-references
  docstring = process_cross_references(docstring)

  # Clear lines and replace with the cleaned, structured overloads
  lines.clear()
  lines.extend(docstring.splitlines())


def preprocess_markdown(app, docname, source):
  """
  Preprocess Markdown files before inclusion in the documentation.
  """
  # Apply only to specific Markdown files, e.g., CHANGELOG.md
  if docname == "CHANGELOG":
    # Inject the currentmodule directive at the beginning
    source[0] = (
      "```{eval-rst}\n.. currentmodule:: pyvinecopulib\n```\n" + source[0]
    )
    # Process cross-references without requiring the module name
    source[0] = process_cross_references(source[0], is_docstring=False)


# --- Dynamic source staging --------------------------------------------------
# Equivalent of the (now-retired) docs/gen_sphinx.py driver: stage the repo's
# README, CHANGELOG, and examples/ dir into the docs source tree at build
# time, and write the generated features.rst / examples.rst toctrees. This
# lets plain `sphinx-build docs <out>` work end-to-end. All staged paths are
# gitignored.

# Public API surface, grouped by subpackage. Documentation references the
# canonical location of each symbol (pyvinecopulib.<subpkg>.<name>); the
# top-level aliases still work but are intentionally not documented.
DOCSTRING_SUBPACKAGES = {
  "core": {
    "classes": [
      "Bicop",
      "FitControlsBicop",
      "Vinecop",
      "FitControlsVinecop",
      "CVineStructure",
      "DVineStructure",
      "RVineStructure",
    ],
    "functions": [],
  },
  "families": {
    "classes": ["BicopFamily"],
    "functions": [],
  },
  "utils": {
    "classes": ["Kde1d"],
    "functions": [
      "to_pseudo_obs",
      "simulate_uniform",
      "wdm",
      "ghalton",
      "sobol",
      "pairs_copula_data",
      "benchmark",
    ],
  },
}

# Kept for backward compatibility with downstream tools that consume these
# names. Synthesized from DOCSTRING_SUBPACKAGES.
DOCSTRING_CLASSES = [
  cls for sub in DOCSTRING_SUBPACKAGES.values() for cls in sub["classes"]
]
DOCSTRING_FUNCTIONS = [
  fn for sub in DOCSTRING_SUBPACKAGES.values() for fn in sub["functions"]
]


def _stage_repo_files(docs_dir, repo_root):
  """Symlink (or copy on Windows) repo-root files into the docs source dir
  so the toctree in index.rst can reference them by name."""
  import shutil

  to_stage = ["README.md", "examples", "CHANGELOG.md", "CONTRIBUTING.md"]
  for name in to_stage:
    src = os.path.join(repo_root, name)
    dst = os.path.join(docs_dir, name)
    if os.path.lexists(dst):
      # Already present (a previous build, or hand-placed). Trust it.
      continue
    try:
      os.symlink(src, dst)
    except (OSError, NotImplementedError):
      # Fallback for Windows runners without symlink permission.
      if os.path.isdir(src):
        shutil.copytree(src, dst)
      else:
        shutil.copy2(src, dst)


def _write_features_rst(out_path):
  """Generate the auto-summary RST for the public API.

  One section per subpackage (core / families / utils), each with an
  autosummary toctree of its classes and (if any) explicit autofunctions.
  """
  rst_name = "API Documentation"
  bar = "=" * len(rst_name)
  with open(out_path, "w") as f:
    f.write(".. GENERATED FILE DO NOT EDIT\n\n")
    f.write(f"{bar}\n{rst_name}\n{bar}\n\n")
    for subpkg, contents in DOCSTRING_SUBPACKAGES.items():
      module = f"pyvinecopulib.{subpkg}"
      header = f"``{module}``"
      f.write(f"{header}\n{'-' * len(header)}\n\n")
      f.write(f".. automodule:: {module}\n\n")
      classes = contents.get("classes", [])
      if classes:
        f.write("Classes\n^^^^^^^\n\n")
        f.write(".. autosummary::\n    :toctree: _generate\n\n")
        for cls in classes:
          f.write(f"    {module}.{cls}\n")
        f.write("\n")
      functions = contents.get("functions", [])
      if functions:
        f.write("Functions\n^^^^^^^^^\n\n")
        for fn in functions:
          f.write(f".. autofunction:: {module}.{fn}\n")
        f.write("\n")


def _write_examples_rst(out_path, examples_dir):
  """Generate the toctree linking to all example notebooks."""
  if not os.path.isdir(examples_dir):
    return
  notebooks = sorted(
    f for f in os.listdir(examples_dir) if f.endswith(".ipynb")
  )
  if not notebooks:
    return
  with open(out_path, "w") as f:
    f.write("Examples\n========\n\n")
    f.write(
      "The following example notebooks are included in this documentation:\n\n"
    )
    f.write(".. toctree::\n   :maxdepth: 1\n   :titlesonly:\n\n")
    for nb in notebooks:
      f.write(f"   examples/{os.path.splitext(nb)[0]}\n")


# Stage README/CHANGELOG/CONTRIBUTING/examples + generate features.rst /
# examples.rst eagerly at conf.py load time. Doing this at module level
# (rather than in a `builder-inited` hook) is required because
# sphinx.ext.autosummary's stub-generation also runs at builder-inited
# and races with us — if features.rst doesn't exist before autosummary's
# scan, the `:toctree:` stubs are never written, the per-class pages
# never render, and on a single-pass build the cross-links from the API
# index (features.html) point at non-existent _generate/*.html files.
_DOCS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_DOCS_DIR)
_stage_repo_files(_DOCS_DIR, _REPO_ROOT)
_write_features_rst(os.path.join(_DOCS_DIR, "features.rst"))
_write_examples_rst(
  os.path.join(_DOCS_DIR, "examples.rst"),
  os.path.join(_DOCS_DIR, "examples"),
)


# Register Sphinx setup with recommonmark configuration and autodoc
def setup(app):
  """
  Configure Sphinx to handle autodoc and preprocess Markdown files.
  """
  # Register the autodoc docstring processor
  app.connect("autodoc-process-docstring", autodoc_process_docstring)

  # Register the Markdown preprocessor
  app.connect("source-read", preprocess_markdown)


# Allow .md files to be included
source_suffix = {
  ".rst": "restructuredtext",
  ".md": "markdown",
}
