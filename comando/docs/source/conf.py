"""Configuration file for building COMANDO's documentation."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHOR: Marco Langiu

# ececute with:
# export SPHINX_APIDOC_OPTIONS="members,undoc-members,show-inheritance"
# sphinx-apidoc --separate --force --implicit-namespaces --module-first \
#   -o source/api ../comando && make clean html latexpdf
import os
import sys
sys.path.insert(0, os.path.abspath(''))
sys.path.insert(0, os.path.abspath('../..'))

# import recommonmark
# from recommonmark.parser import CommonMarkParser
# from recommonmark.transform import AutoStructify
from pypandoc import convert_file

from inject_dynamic_member_doc import inject


with open('interfaces.rst', 'w') as f:
    f.write('.. _interfaces:\n\n')
    f.write(convert_file('../../comando/interfaces/README.md', 'rst'))

# -- Project information -----------------------------------------------------

project = 'COMANDO'
copyright = '2020, Marco Langiu'
author = 'Marco Langiu'

# The full version, including alpha/beta/rc tags
import comando
release = comando.__version__

# -- General configuration ---------------------------------------------------
# The master toctree document.
master_doc = 'index'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',  # Render Latex in html
    'sphinx_math_dollar',
    'sphinx.ext.napoleon',  # Teaches Sphinx to understand NumPy docstrings
    'sphinx.ext.intersphinx',  # Map to other project's sphinx dicumentation
    # 'recommonmark'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['*debug*']

# INCLUDING MARKDOWN (NEVER MANAGED TO GET THIS TO WORK!)
# source_parsers = {
#     '.md': CommonMarkParser,
# }
# source_suffix = ['.rst', '.md']
source_suffix = ['.rst']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
# html_theme = 'sphinx_pdj_theme'
# html_theme = 'groundwork'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# from https://sphinxguide.readthedocs.io/en/latest/sphinx_basics/settings.html
# -- Options for LaTeX output ---------------------------------------------
with open('preamble.tex', 'r') as f:
    PREAMBLE = f.read()
latex_engine = 'pdflatex'
# latex_engine = 'xelatex'
latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    'papersize': 'a4paper',
    'releasename': "COMANDO",
    # Sonny, Lenny, Glenn, Conny, Rejne, Bjarne and Bjornstrup
    # 'fncychap': '\\usepackage[Lenny]{fncychap}',
    # 'fncychap': '\\usepackage{fncychap}',
    'fontpkg': '\\usepackage{amsmath,amsfonts,amssymb,amsthm}',

    'figure_align': 'htbp',
    # The font size ('10pt', '11pt' or '12pt').
    #
    'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    'preamble': PREAMBLE,

    'maketitle': r"""
\pagenumbering{Roman} %%% to avoid page 1 conflict with actual page 1

\begin{titlepage}
    \centering

    \vspace*{40mm} %%% * is used to give space from top
    \textbf{\Huge {The COMANDO Framework}}

    \Large \textbf{{Marco Langiu}}

\end{titlepage}

\clearpage
\pagenumbering{roman}
\tableofcontents
\listoffigures
\listoftables
\clearpage
\pagenumbering{arabic}
""",
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
    'sphinxsetup': r"""
    hmargin={0.7in,0.7in},
    vmargin={1in,1in},
    verbatimwithframe=true,
    TitleColor={rgb}{0,0,0},
    HeaderFamily=\rmfamily\bfseries,
    InnerLinkColor={rgb}{0,0,1},
    OuterLinkColor={rgb}{0,0,1}
""",
    'tableofcontents': ' ',
}

# latex_logo = '_static/logo.jpg'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'main.tex', 'COMANDO Documentation', 'Marco Langiu', 'manual')
]


# -- Extension configuration -------------------------------------------------
# Display todos by setting to True
todo_include_todos = True


mathjax_config = {
    'tex2jax': {
        'inlineMath': [["\\(", "\\)"]],
        'displayMath': [["\\[", "\\]"]],
        'packages': ['base', 'require'],
    },
}

# https://www.sphinx-doc.org/en/stable/ext/napoleon.html#configuration
napoleon_google_docstring = False
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = True
napoleon_use_param = False
napoleon_custom_sections = [('Arguments', 'Parameters'),
                            ('Options', 'Parameters')]


intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'sympy': ('https://docs.sympy.org/latest', None),
    'pandas': ('http://pandas.pydata.org/pandas-docs/dev', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
    'matplotlib': ('https://matplotlib.org', None)
}


def run_apidoc(_):
    """Run sphinx-apidoc."""
    import subprocess
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    output_path = os.path.join(cur_dir, 'api')
    cmd_path = 'sphinx-apidoc'
    subprocess.check_call([cmd_path, '-e', '-o', output_path, '../../comando', '--force'])
    inject()


# This is the expected signature of the handler for this event, cf doc
# blacklist = ['sympy']
def autodoc_skip_member_handler(app, what, name, obj, skip, options):
    """Skip sympy modules included in COMANDO's namespace."""
    try:
        from_sympy = obj.__module__.startswith('sympy')
    except AttributeError:
        from_sympy = False
    return skip or from_sympy


# Automatically called by sphinx at startup
def setup(app):
    """Connect the autodoc-skip-member event from apidoc to the callback."""
    app.connect('builder-inited', run_apidoc)
    app.connect('autodoc-skip-member', autodoc_skip_member_handler)

    # # Mention this on the reconmark doc
    # app.add_config_value('recommonmark_config', {
    #     # 'url_resolver': lambda url: github_doc_root + url,
    #     'auto_toc_tree_section': 'Contents',
    #     'enable_math': False,
    #     'enable_inline_math': False,
    #     'enable_eval_rst': True,
    #     # 'enable_auto_doc_ref': True,  DEPRECATED
    # }, True)
    # app.add_transform(AutoStructify)