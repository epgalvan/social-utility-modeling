# Configuration file for the Sphinx documentation builder.

import os

# -- Project information

project = 'Social Utility Modeling'
copyright = '2023, Elijah Galvan'
author = 'Elijah Galvan'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = ["sphinx_design"]

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
