# pyvinecopulib documentation build configuration file.

import inspect
import os
import re
from typing import Any

import sphinx.ext.autodoc as autodoc
import sphinx.util.inspect as sphinxinspect
from sphinx.ext.autodoc import AttributeDocumenter, ModuleDocumenter

import pyvinecopulib as pv

# Monkey-patch sphinx's autodoc to recognise nanobind's nb_func as a routine.


def isnbfunc(obj: Any) -> bool:
  return (
    hasattr(type(obj), "__module__")
    and type(obj).__module__ == "nanobind"
    and type(obj).__name__ == "nb_func"
  )


def isfunction(obj: Any) -> bool:
  return inspect.isfunction(sphinxinspect.unpartial(obj)) or isnbfunc(obj)


def isroutine(obj: Any) -> bool:
  return inspect.isroutine(sphinxinspect.unpartial(obj)) or isnbfunc(obj)


sphinxinspect.isfunction = isfunction
sphinxinspect.isroutine = isroutine

assert autodoc.inspect.isfunction is isfunction
assert autodoc.inspect.isroutine is isroutine


@classmethod
def patched_can_document_member(
  cls, member: Any, membername: str, isattr: bool, parent: Any
) -> bool:
  if isinstance(parent, ModuleDocumenter):
    return False
  if sphinxinspect.isroutine(member):
    return False
  if sphinxinspect.isattributedescriptor(member):
    return True
  return not sphinxinspect.isroutine(member) and not isinstance(member, type)


AttributeDocumenter.can_document_member = patched_can_document_member

assert inspect.getsource(
  AttributeDocumenter.can_document_member
) == inspect.getsource(patched_can_document_member)


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

# `autosummary`: stub files don't exist on first toctree read (single-pass build).
# `myst.header`: index.rst inlines a slice of README.md starting at H2.
suppress_warnings = ["autosummary", "myst.header"]

source_suffix = {
  ".rst": "restructuredtext",
  ".md": "markdown",
}

templates_path = ["_templates"]
master_doc = "index"

# Don't exclude `_generate` — autosummary writes stubs there and sphinx
# reads them back as source documents.
exclude_patterns = ["_build", "**/.ipynb_checkpoints"]

project = "pyvinecopulib"
copyright = "2024, Thomas Nagler and Thibault Vatter"
author = "Thomas Nagler and Thibault Vatter"

release = pv.__version__
version = ".".join(release.split(".")[:3])

# Don't set `html_extra_path = ["examples"]` — nbsphinx registers the
# staged notebooks as documents, and html_extra_path would shadow that.

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

  content = re.sub(r"``(\w+)\.(\w+)\(\)``", rf"{meth_ref}\1.\2`", content)
  content = re.sub(r"``(\w+)\.(\w+)``", f"{cls_ref}\1.\2`", content)
  for cls in classes:
    content = re.sub(rf"``{cls}``", rf"{cls_ref}{cls}`", content)

  return content


def autodoc_process_docstring(app, what, name, obj, options, lines):
  docstring = process_cross_references("\n".join(lines))
  lines.clear()
  lines.extend(docstring.splitlines())


def preprocess_markdown(app, docname, source):
  if docname == "CHANGELOG":
    source[0] = (
      "```{eval-rst}\n.. currentmodule:: pyvinecopulib\n```\n" + source[0]
    )
    source[0] = process_cross_references(source[0], is_docstring=False)


# Stage README/CHANGELOG/CONTRIBUTING/examples into the docs source tree and
# write features.rst / examples.rst toctrees so plain `sphinx-build` works.
# Documentation references the canonical pyvinecopulib.<subpkg>.<name>;
# top-level aliases still resolve but aren't documented.
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
  "sklearn": {
    "classes": ["VineDensity", "VineRegressor"],
    "functions": [],
  },
}

DOCSTRING_CLASSES = [
  cls for sub in DOCSTRING_SUBPACKAGES.values() for cls in sub["classes"]
]
DOCSTRING_FUNCTIONS = [
  fn for sub in DOCSTRING_SUBPACKAGES.values() for fn in sub["functions"]
]


def _stage_repo_files(docs_dir, repo_root):
  """Symlink (or copy on Windows) repo-root files into the docs source dir."""
  import shutil

  to_stage = ["README.md", "examples", "CHANGELOG.md", "CONTRIBUTING.md"]
  for name in to_stage:
    src = os.path.join(repo_root, name)
    dst = os.path.join(docs_dir, name)
    if os.path.lexists(dst):
      continue
    try:
      os.symlink(src, dst)
    except (OSError, NotImplementedError):
      if os.path.isdir(src):
        shutil.copytree(src, dst)
      else:
        shutil.copy2(src, dst)


def _write_features_rst(out_path):
  """Generate API documentation RST: one section per subpackage."""
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
        f.write(".. autosummary::\n    :toctree: _generate\n\n")
        for fn in functions:
          f.write(f"    {module}.{fn}\n")
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


# Stage docs and write features.rst / examples.rst at conf.py load time.
# Must be eager (not `builder-inited`), otherwise it races with autosummary's
# stub generation and the per-class pages never render on a single-pass build.
_DOCS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_DOCS_DIR)
_stage_repo_files(_DOCS_DIR, _REPO_ROOT)
_write_features_rst(os.path.join(_DOCS_DIR, "features.rst"))
_write_examples_rst(
  os.path.join(_DOCS_DIR, "examples.rst"),
  os.path.join(_DOCS_DIR, "examples"),
)


def setup(app):
  app.connect("autodoc-process-docstring", autodoc_process_docstring)
  app.connect("source-read", preprocess_markdown)


source_suffix = {
  ".rst": "restructuredtext",
  ".md": "markdown",
}
