# -*- coding: utf-8 -*-
#
# Sphinx documentation build configuration file, created by
# sphinx-quickstart.py on Sat Mar  8 21:47:50 2008.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default value; values that are commented out
# serve to show the default value.
from __future__ import print_function

import sys, os, re
import contextlib

@contextlib.contextmanager
def cd(newpath):
    """
    Change the current working directory to `newpath`, temporarily.

    If the old current working directory no longer exists, do not return back.
    """
    oldpath = os.getcwd()
    os.chdir(newpath)
    try:
        yield
    finally:
        try:
            os.chdir(oldpath)
        except OSError:
            # If oldpath no longer exists, stay where we are.
            pass

# Check Sphinx version
import sphinx
if sphinx.__version__ < "1.3":
    raise RuntimeError("Sphinx 1.3 or newer required")

# Environment variable to know if the docs are being built on rtd.
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
print
print("Building on ReadTheDocs: {}".format(on_rtd))
print
print("Current working directory: {}".format(os.path.abspath(os.curdir)))
print("Python: {}".format(sys.executable))

if on_rtd:
    # Build is not via Makefile (yet).
    # So we manually build the examples and gallery.
    import subprocess
    with cd('..'):
        # The Makefile is run from kgen/doc, so we need to move there
        # from kgen/doc/source (which holds conf.py).
        py = sys.executable
        #subprocess.call([py, 'make_gallery.py'])
        subprocess.call([py, 'make_examples_rst.py', '../examples', 'source'])

# If your extensions are in another directory, add it here.
# These locations are relative to conf.py
sys.path.append(os.path.abspath('../sphinxext'))

# General configuration
# ---------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    #'sphinxcontrib.bibtex',
    #'IPython.sphinxext.ipython_console_highlighting',
    #'IPython.sphinxext.ipython_directive',
]


# generate autosummary pages
autosummary_generate=True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates','../rst_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
source_encoding = 'utf-8'

# The master toctree document.
master_doc = 'index'

# General substitutions.
project = 'KGen'
copyright = '2016, KGen Developers'

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
#
# The short X.Y version.
#import kgen
# TODO: need to read version from src
version = '0.7.1'
# The full version, including dev info
release = '0.7.1'

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of documents that shouldn't be included in the build.
# unused_docs = ['reference/pdf_reference']

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).

add_module_names = False

# show_authors = True

# The name of the Pygments (syntax highlighting) style to use.
#pygments_style = 'friendly'
pygments_style = 'sphinx'

# A list of prefixs that are ignored when creating the module index. (new in Sphinx 0.6)
modindex_common_prefix=['kgen.']

doctest_global_setup="import kgen"

# treat ``x, y : type`` as vars x and y instead of default ``y(x,) : type``
napoleon_use_param = False

# Options for HTML output
# -----------------------

if not on_rtd:
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

#html_theme_options = {
#    "rightsidebar": "true",
#    "relbarbgcolor: "black"
#}


# The style sheet to use for HTML and HTML Help pages. A file of that name
# must exist either in Sphinx' static/ path, or in one of the custom paths
# given in html_static_path.
html_style = 'networkx.css'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Content template for the index page.
#html_index = 'index.html'
html_index = 'contents.html'

# Custom sidebar templates, maps page names to templates.
#html_sidebars = {'index': 'indexsidebar.html'}

# Additional templates that should be rendered to pages, maps page names to
# templates.
#html_additional_pages = {'index': 'index.html','gallery':'gallery.html'}
#html_additional_pages = {'gallery':'gallery.html'}

# If true, the reST sources are included in the HTML build as _sources/<name>.
html_copy_source = False

html_use_opensearch = 'http://kgen.github.io'

# Output file base name for HTML help builder.
htmlhelp_basename = 'KGen'

pngmath_use_preview = True

# Options for LaTeX output
# ------------------------

# The paper size ('letter' or 'a4').
latex_paper_size = 'letter'

# The font size ('10pt', '11pt' or '12pt').
#latex_font_size = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [('tutorial/index', 'kgen_tutorial.tex',
                    'KGen Tutorial',
                    'Aric Hagberg, Dan Schult, Pieter Swart', 'howto', 1),
                   ('reference/pdf_reference', 'kgen_reference.tex',
                    'KGen Reference',
                    'Aric Hagberg, Dan Schult, Pieter Swart', 'manual', 1)]

#latex_appendices = ['installing']#,'legal'],'citing','credits','history']

#latex_appendices = ['credits']

# Intersphinx mapping
intersphinx_mapping = {'http://docs.python.org/': None,
                       'http://docs.scipy.org/doc/numpy/': None,
                      }

# For trac custom roles

default_role = 'obj'
trac_url = 'https://notdefined/trac/'
mathjax_path = 'https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML'

numpydoc_show_class_members = False

# add path for source code
sys.path.insert(0, os.path.abspath('../../../base'))
