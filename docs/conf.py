#!/usr/bin/env python
#
# HL-gaps-pub documentation build configuration file, created by
# sphinx-quickstart on Fri Jun  9 13:47:02 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import os
import sys

# If extensions (or modules to document with autodoc) are in another
# directory, add these directories to sys.path here. If the directory is
# relative to the documentation root, use os.path.abspath to make it
# absolute, like shown here.
#
sys.path.insert(0, os.path.abspath("../"))


import hl_gaps_pub

# -- General configuration ---------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = "1.0"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "sphinx_tabs.tabs",
    "sphinx_click",
    'nbsphinx',  # nbsphinx extension
]

# -- Options for nbsphinx extension ---------------------------------------

# Allow nbgallery directive
nbsphinx_allow_directives = {
    'gallery',  # Allows nbgallery (which is an alias for gallery)
    'code-block', # Allows code blocks inside notebooks
}

# Example of other nbsphinx options (you might already have some of these)
nbsphinx_execute = 'never'  # Don't execute notebooks by default (use :timeout:)
# nbsphinx_timeout = 600  # Global timeout (you can override per-notebook)

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = [".rst", ".md"]
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "HL-gaps-pub"
copyright = "2025, Sascha Thinius"
author = "Sascha Thinius"

# The version info for the project you're documenting, acts as replacement
# for |version| and |release|, also used in various other places throughout
# the built documents.
#
# The short X.Y version.
version = hl_gaps_pub.__version__
# The full version, including alpha/beta/rc tags.
release = hl_gaps_pub.__version__

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = "en"
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']  

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output -------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

# Theme options are theme-specific and customize the look and feel of a
# theme further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as
# html_title.
# html_short_title = None

# The name of an image file (relative to this directory) to place at the
# top of the sidebar.
# html_logo = None

# The name of an image file (within the static path) to use as favicon
# of the docs.  This file should be a Windows icon file (.ico) being
# 16x16 or 32x32 pixels large.
# html_favicon = None

# Add any paths that contain custom static files (such as style sheets)
# here, relative to this directory. They are copied after the builtin
# static files, so a file named "default.css" will overwrite the builtin
# "default.css".
html_static_path = ["_static"]
html_logo = "_static/flask.png"
html_theme_options = {
    #"light_logo": "logo-light-mode.png",
    #"dark_logo": "logo-dark-mode.png",
    "footer_icons": [
        {
            "name": "GitLab",
            "url": "https://github.com/sthinius87/HL-gaps-pub",
            "html": """
                <svg stroke="currentColor" fill="none" stroke-width="2" viewBox="0 0 24 24" stroke-linecap="round" stroke-linejoin="round" height="1em" width="1em" xmlns="http://www.w3.org/2000/svg"><path d="M22.65 14.39L12 22.13 1.35 14.39a.84.84 0 0 1-.3-.94l1.22-3.78 2.44-7.51A.42.42 0 0 1 4.82 2a.43.43 0 0 1 .58 0 .42.42 0 0 1 .11.18l2.44 7.49h8.1l2.44-7.51A.42.42 0 0 1 18.6 2a.43.43 0 0 1 .58 0 .42.42 0 0 1 .11.18l2.44 7.51L23 13.45a.84.84 0 0 1-.35.94z"></path></svg>
            """,
            "class": "",
        },
    ],
    "source_edit_link": "https://github.com/sthinius87/HL-gaps-pub/docs/{filename}",
}

# If not "", a "Last updated on:" timestamp is inserted at every page
# bottom, using the given strftime format.
# html_last_updated_fmt = "%b %d, %Y"

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {
#     "**": ["globaltoc.html", "relations.html", "sourcelink.html", "searchbox.html"]
# }

# Additional templates that should be rendered to pages, maps page names
# to template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer.
# Default is True.
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer.
# Default is True.
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages
# will contain a <link> tag referring to it.  The value of this option
# must be the base URL from which the finished HTML is served.
# html_use_opensearch = ""

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

html_context = {
    "display_gitlab": True,  # Integrate Gitlab
    "gitlab_host": "HL_gaps_pub",
    "gitlab_user": "ifam418",  # Username
    "gitlab_repo": "HL-gaps-pub",  # Repo name
    "gitlab_version": "master",  # Version
    "conf_py_path": "/docs/",  # Path in the checkout to the docs root
}

# -- Options for HTMLHelp output ---------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "HL-gaps-pubdoc"


# -- Options for LaTeX output ------------------------------------------

latex_elements = {
    # The paper size ("letterpaper" or "a4paper").
    #
    # "papersize": "letterpaper",
    # The font size ("10pt", "11pt" or "12pt").
    #
    # "pointsize": "10pt",
    # Additional stuff for the LaTeX preamble.
    #
    "preamble": r"\setcounter{tocdepth}{8}"
    # Latex figure (float) alignment
    #
    # "figure_align": "htbp",
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass
# [howto, manual, or own class]).
latex_documents = [
    (
        master_doc,
        "HL-gaps-pub.tex",
        r"HL-gaps-pub Documentation",
        "Sascha Thinius",
        "manual",
     ),
]


# -- Options for manual page output ------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (
        master_doc,
        "HL-gaps-pub",
        "HL-gaps-pub Documentation",
        [author],
        1,
    )
]


# -- Options for Texinfo output ----------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "HL-gaps-pub",
        "HL-gaps-pub Documentation",
        author,
        "HL-gaps-pub",
        "One line description of project.",
        "Miscellaneous",
    ),
]

# -- Options for matplotlib output -------------------------------------

# Default value for the include-source option
# plot_include_source=False

# Whether to show a link to the source in HTML.
# plot_html_show_source_link=False

# Code that should be executed before each plot.
# plot_pre_code=

# Base directory, to which plot:: file names are relative to.
# (If None or empty, file names are relative to the directory where the file containing the directive is.)
# plot_basedir

# File formats to generate. List of tuples or strings:
# [(suffix, dpi), suffix, ...]
# that determine the file format and the DPI. For entries whose DPI was omitted,
# sensible defaults are chosen. When passing from the command line through sphinx_build
# the list should be passed as suffix:dpi,suffix:dpi, ....
# plot_formats=

# Whether to show links to the files in HTML.
# plot_html_show_formats=False

# A dictionary containing any non-standard rcParams that should be applied before each plot.
# plot_rcparams=

# By default, rcParams are applied when context option is not used in a plot directive.
# This configuration option overrides this behavior and applies rcParams before each plot.
# plot_apply_rcparams=

# By default, the working directory will be changed to the directory of the example,
# so the code can get at its data files, if any. Also its path will be added to sys.path
# so it can import any helper modules sitting beside it. This configuration option
# can be used to specify a central directory (also added to sys.path) where data files
# and helper modules for all code are located.
# plot_working_directory=

# Provide a customized template for preparing restructured text.
# plot_template=
