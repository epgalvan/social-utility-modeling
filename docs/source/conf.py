# Configuration file for the Sphinx documentation builder.

import os, sys

# -- Project information

project = 'Social Utility Modeling'
copyright = '2023, Elijah Galvan'
author = 'Elijah Galvan'

release = '0.1'
version = '0.1.0'

# -- General configuration

sys.path.append(os.path.abspath('sphinx_design'))

extensions = ['dropdown', 'tabs', 'grids']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
