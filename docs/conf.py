# -*- coding: utf-8 -*-
#
# pyvinecopulib documentation build configuration file

# Sphinx extension modules
from pkg_resources import get_distribution

# -- General configuration ------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.githubpages',
    'sphinx.ext.mathjax',
    'sphinx_rtd_theme',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
]

napoleon_include_init_with_doc = True
autosummary_generate = True

# The suffix(es) of source filenames.
source_suffix = '.rst'

# For the templates.
templates_path = ['_templates']

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'pyvinecopulib'
copyright = u'2019, Thomas Nagler and Thibault Vatter'
author = u'Thomas Nagler and Thibault Vatter'

# The version info.
release = get_distribution('pyvinecopulib').version
version = '.'.join(release.split('.')[:3])

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']

html_copy_source = False

html_show_copyright = False

html_show_sphinx = False

add_module_names = False

pygments_style = 'sphinx'

html_logo = '_static/pyvinecopulib.png'
