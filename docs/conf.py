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
nbsphinx_execute = "always"

# The suffix(es) of source filenames.
source_suffix = {
  ".rst": "restructuredtext",
  ".md": "markdown",
}

# For the templates.
templates_path = ["_templates"]

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "pyvinecopulib"
copyright = "2024, Thomas Nagler and Thibault Vatter"
author = "Thomas Nagler and Thibault Vatter"

# The version info.
release = pv.__version__
version = ".".join(release.split(".")[:3])

# Specify additional files to copy
html_extra_path = [
  os.path.abspath(os.path.join(os.path.dirname(__file__), "../CHANGELOG.md"))
]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"

html_static_path = ["_static"]

html_copy_source = False

html_show_copyright = False

html_show_sphinx = False

add_module_names = False

pygments_style = "sphinx"

html_logo = "_static/pyvinecopulib.png"

html_extra_path = ["../examples/"]


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
