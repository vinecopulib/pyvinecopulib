# -*- coding: utf-8 -*-
#
# pyvinecopulib documentation build configuration file

# Sphinx extension modules

import inspect
import re
from typing import Any

import sphinx.ext.autodoc as autodoc
import sphinx.util.inspect as sphinxinspect
from pkg_resources import get_distribution
from sphinx.ext.autodoc import AttributeDocumenter, ModuleDocumenter

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
]

napoleon_include_init_with_doc = True
napoleon_use_rtype = False
napoleon_custom_sections = [("Usage", "Usage")]
autosummary_generate = True

# The suffix(es) of source filenames.
source_suffix = ".rst"

# For the templates.
templates_path = ["_templates"]

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "pyvinecopulib"
copyright = "2024, Thomas Nagler and Thibault Vatter"
author = "Thomas Nagler and Thibault Vatter"

# The version info.
release = get_distribution("pyvinecopulib").version
version = ".".join(release.split(".")[:3])

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"

html_static_path = ["_static"]

html_copy_source = False

html_show_copyright = False

html_show_sphinx = False

add_module_names = False

pygments_style = "sphinx"

html_logo = "_static/pyvinecopulib.png"


def autodoc_process_docstring(app, what, name, obj, options, lines):
  print(f"Processing: {what}, {name}")

  # Join the existing lines and try to reformat the docstring
  docstring = "\n".join(lines)

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

  # Try convert references to methods or classes to cross-references
  docstring = re.sub(r"``(\w+)\.(\w+)\(\)``", r":meth:`\1.\2`", docstring)
  docstring = re.sub(r"``(\w+)\.(\w+)``", r":class:`\1.\2`", docstring)
  for cls in classes:
    docstring = re.sub(r"``" + cls + "``", r":class:`" + cls + "`", docstring)

  # overloaded_classes_constructors = {
  #   "pyvinecopulib." + cls + ".__init__"
  #   for cls in [
  #     "Bicop",
  #     "Vinecop",
  #     "CVineStructure",
  #     "DVineStructure",
  #     "RVineStructure",
  #   ]
  # }
  # if name in overloaded_classes_constructors:
  #   docstring = process_overloaded(docstring)
  # Clear lines and replace with the cleaned, structured overloads
  lines.clear()
  lines.extend(docstring.splitlines())


# Register the event handler with Sphinx
def setup(app):
  app.connect("autodoc-process-docstring", autodoc_process_docstring)