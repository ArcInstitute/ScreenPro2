import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

from screenpro import __version__


# -- General configuration ---------------------------------------------

# https://brendanhasz.github.io/2019/01/05/sphinx.html

# General information about the project.
project = 'ScreenPro2'
author = "ScreenPro2 Development Team"
copyright = "2022-2024 – ScreenPro2 Developers Team"
# Gilbart Lab, UCSF / Arc Institute.
# Multi-Omics Tech Center, Arc Insititue.

repository_url = "https://github.com/ArcInstitute/ScreenPro2"

version = __version__
release = __version__

language = 'en'

# improve class docs
# By default, Sphinx ignores docstrings for __init__.py and shows only top-level docstrings for the class itself. 

autoclass_content = 'both'

autodoc_mock_imports = [
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "seaborn",
    "scanpy",
    "anndata",
    "screenpro",
    "pyarrow",
    "biobear",
    "click",
    "polars",
    "biobear",
    "numba",
    "bokeh",
    "pydeseq2",
    "watermark"
]

# Bibliography settings
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"
bibtex_default_style = 'unsrt'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
source_suffix = ['.rst', '.md']

# The master toctree document.
master_doc = 'index'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# Enabling napoleon
napoleon_google_docstring = True
napoleon_numpy_docstring = False

extensions = [
    'sphinx.ext.autodoc', 
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    "sphinx.ext.extlinks",
    'sphinx.ext.viewcode', 
    "sphinxcontrib.bibtex",
    'myst_parser',
]

# -- Options for HTML output -------------------------------------------
# Activate the theme.
html_theme = 'sphinx_rtd_theme'
# html_theme = "furo"

# Theme options are theme-specific and customize the look and feel of a
# theme further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']


# -- Options for HTMLHelp output ---------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'screenprodoc'


# -- Options for LaTeX output ------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',

    # Latex figure (float) alignment
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass
# [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'screenpro.tex',
     'ScreenPro2 Documentation',
     'ScreenPro2 Development Team', 
     'manual'),
]


# -- Options for manual page output ------------------------------------

htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", doc_title, author, "manual")]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (
        master_doc,
        project,
        doc_title,
        author,
        project,
        "One line description of project.",
        "Miscellaneous",
    )
]

intersphinx_mapping = {
    'anndata': ('https://anndata.readthedocs.io/en/latest/', None)
}

# extlinks config
extlinks = {
    "issue": ("https://github.com/scverse/scanpy/issues/%s", "issue%s"),
    "pr": ("https://github.com/scverse/scanpy/pull/%s", "pr%s"),
}