# -*- coding: utf-8 -*-
import datetime
import os
import sys
from pathlib import Path

import git
import sphinx_rtd_theme
from sphinx.ext.autodoc import (
    ClassLevelDocumenter,
    DataDocumenter,
    InstanceAttributeDocumenter,
)
from sphinx.ext.napoleon.docstring import GoogleDocstring

import qalgebra


DOCS_SOURCES = Path(__file__).parent
ROOT = DOCS_SOURCES / '..' / '..'  # project root

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('_extensions'))

# -- Generate API documentation ------------------------------------------------


def run_apidoc(app):
    """Generage API documentation"""
    import better_apidoc

    better_apidoc.APP = app
    better_apidoc.main(
        [
            'better-apidoc',
            '-t',
            str(DOCS_SOURCES / '_templates'),
            '--force',
            '--no-toc',
            '--separate',
            '-o',
            str(DOCS_SOURCES / 'API'),
            os.path.join(ROOT / 'src' / 'qalgebra'),
        ]
    )


# -- General configuration -----------------------------------------------------

# Report broken links as warnings
nitpicky = True
nitpick_ignore = [('py:class', 'callable')]

extensions = [
    'graphviz_ext',
    'inheritance_diagram',
    'doctr_versions_menu',
    'nbsphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    'sphinx_copybutton',
    'dollarmath',
]


if os.getenv('SPELLCHECK'):
    extensions.append('sphinxcontrib.spelling')
    spelling_show_suggestions = True
    spelling_lang = os.getenv('SPELLCHECK')
    spelling_word_list_filename = 'spelling_wordlist.txt'
    spelling_ignore_pypi_package_names = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3.8', None),
    'sympy': ('https://docs.sympy.org/latest/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
    'matplotlib': ('https://matplotlib.org/', None),
    'qutip': ('http://qutip.org/docs/latest/', None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

source_suffix = '.rst'
master_doc = 'index'
project = 'QAlgebra'
year = str(datetime.datetime.now().year)
author = 'Michael Goerz'
copyright = '{0}, {1}'.format(year, author)
version = qalgebra.__version__
release = version
git_tag = "v%s" % version
if version.endswith('dev'):
    try:
        last_commit = str(git.Repo(ROOT).head.commit)[:7]
        release = "%s (%s)" % (version, last_commit)
        git_tag = str(git.Repo(ROOT).head.commit)
    except git.exc.InvalidGitRepositoryError:
        git_tag = "master"
numfig = True

pygments_style = 'friendly'
extlinks = {
    'issue': ('https://github.com/qalgebra/qalgebra/issues/%s', '#'),
    'pr': ('https://github.com/qalgebra/qalgebra/pull/%s', 'PR #'),
}

# autodoc settings
autoclass_content = 'class'
autodoc_member_order = 'bysource'
autodoc_mock_imports = []  # e.g.: 'numpy', 'scipy', ...


html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {'**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html']}
html_short_title = '%s-%s' % (project, version)


# Mathjax settings
mathjax_path = (
    'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.6/MathJax.js'
)
mathjax_config = {
    'extensions': ['tex2jax.js'],
    'jax': ['input/TeX', 'output/SVG'],
    'TeX': {
        'extensions': ["AMSmath.js", "AMSsymbols.js"],
        'Macros': {
            'tr': ['{\\operatorname{tr}}', 0],
            'Tr': ['{\\operatorname{tr}}', 0],
            'diag': ['{\\operatorname{diag}}', 0],
            'abs': ['{\\operatorname{abs}}', 0],
            'pop': ['{\\operatorname{pop}}', 0],
            'SLH': ['{\\operatorname{SLH}}', 0],
            'aux': ['{\\text{aux}}', 0],
            'opt': ['{\\text{opt}}', 0],
            'tgt': ['{\\text{tgt}}', 0],
            'init': ['{\\text{init}}', 0],
            'lab': ['{\\text{lab}}', 0],
            'rwa': ['{\\text{rwa}}', 0],
            'fwhm': ['{\\text{fwhm}}', 0],
            'bra': ['{\\langle#1\\vert}', 1],
            'ket': ['{\\vert#1\\rangle}', 1],
            'Bra': ['{\\left\\langle#1\\right\\vert}', 1],
            'Braket': [
                '{\\left\\langle #1\\vphantom{#2} \\mid #2\\vphantom{#1}\\right\\rangle}',
                2,
            ],
            'Ket': ['{\\left\\vert#1\\right\\rangle}', 1],
            'mat': ['{\\mathbf{#1}}', 1],
            'op': ['{\\hat{#1}}', 1],
            'Op': ['{\\hat{#1}}', 1],
            'dd': ['{\\,\\text{d}}', 0],
            'daggered': ['{^{\\dagger}}', 0],
            'transposed': ['{^{\\text{T}}}', 0],
            'Liouville': ['{\\mathcal{L}}', 0],
            'DynMap': ['{\\mathcal{E}}', 0],
            'identity': ['{\\mathbf{1}}', 0],
            'Norm': ['{\\lVert#1\\rVert}', 1],
            'Abs': ['{\\left\\vert#1\\right\\vert}', 1],
            'avg': ['{\\langle#1\\rangle}', 1],
            'Avg': ['{\\left\langle#1\\right\\rangle}', 1],
            'AbsSq': ['{\\left\\vert#1\\right\\vert^2}', 1],
            'Re': ['{\\operatorname{Re}}', 0],
            'Im': ['{\\operatorname{Im}}', 0],
            'Real': ['{\\mathbb{R}}', 0],
            'Complex': ['{\\mathbb{C}}', 0],
            'Integer': ['{\\mathbb{N}}', 0],
        },
    },
}


# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# -- Extensions to the  Napoleon GoogleDocstring class ---------------------


# first, we define new methods for any new sections and add them to the class
def parse_keys_section(self, section):
    return self._format_fields('Keys', self._consume_fields())


GoogleDocstring._parse_keys_section = parse_keys_section


def parse_attributes_section(self, section):
    return self._format_fields('Attributes', self._consume_fields())


GoogleDocstring._parse_attributes_section = parse_attributes_section


def parse_class_attributes_section(self, section):
    return self._format_fields('Class Attributes', self._consume_fields())


GoogleDocstring._parse_class_attributes_section = (
    parse_class_attributes_section
)

# we now patch the parse method to guarantee that the the above methods are
# assigned to the _section dict
def patched_parse(self):
    self._sections['keys'] = self._parse_keys_section
    self._sections['class attributes'] = self._parse_class_attributes_section
    self._unpatched_parse()


GoogleDocstring._unpatched_parse = GoogleDocstring._parse
GoogleDocstring._parse = patched_parse


# -- Monkeypatch for instance attribs (sphinx bug #2044) -----------------------


def iad_add_directive_header(self, sig):
    ClassLevelDocumenter.add_directive_header(self, sig)


InstanceAttributeDocumenter.add_directive_header = iad_add_directive_header


# -- Documenter for Singletons -------------------------------------------------


class SingletonDocumenter(DataDocumenter):
    directivetype = 'data'
    objtype = 'singleton'
    priority = 20

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        return isinstance(member, qalgebra.utils.singleton.SingletonType)


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    'collapse_navigation': True,
    'display_version': True,
}

inheritance_graph_attrs = dict(size='""')
graphviz_output_format = 'svg'

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = 'favicon.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None
nbsphinx_prolog = r"""
{% set docname = env.doc2path(env.docname, base='docs') %}

.. only:: html

    .. role:: raw-html(raw)
        :format: html

    :raw-html:`<a href="http://nbviewer.jupyter.org/github/qalgebra/qalgebra/blob/<<GIT_TAG>>/{{ docname }}" target="_blank"><img alt="Render on nbviewer" src="https://img.shields.io/badge/render%20on-nbviewer-orange.svg" style="vertical-align:text-bottom"></a>&nbsp;<a href="https://mybinder.org/v2/gh/qalgebra/qalgebra/<<GIT_TAG>>?filepath={{ docname }}}" target="_blank"><img alt="Launch Binder" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>`
""".replace(
    '<<GIT_TAG>>', git_tag
)

# -- Options for LaTeX output -------------------------------------------------

# latex_engine = 'lualatex'
latex_elements = {
    'preamble': r'''
\usepackage[titles]{tocloft}
\cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
\usepackage{emptypage}
''',
    'fncychap': r'\usepackage[Bjornstrup]{fncychap}',
    'printindex': r'\footnotesize\raggedright\printindex',
    'babel': '',
}
latex_show_urls = 'no'

# -----------------------------------------------------------------------------


def setup(app):
    app.add_autodocumenter(SingletonDocumenter)
    app.connect('builder-inited', run_apidoc)
