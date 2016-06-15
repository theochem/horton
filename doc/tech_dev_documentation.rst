..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2016 The HORTON Development Team
    :
    : This file is part of HORTON.
    :
    : HORTON is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : HORTON is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

Writing documentation
#####################

Introduction
============

All the documentation is located in the ``doc`` directory. We use the `Sphinx
<http://sphinx.pocoo.org/>`_ formatting engine to compile the `documentation
source code` into fancy formatted HTML or PDF.

The source files have the extension ``.rst``, and are written in the
`ReStructuredText <http://docutils.sourceforge.net/rst.html>`_ (RST) format.
RST is in some sense comparable to latex, but more intuitive to use.
It also has some specific advantages for documenting software.

All ``.rst`` files are part of the source tree, just like the actual source
code. Git is also used to keep track of changes in the documentation. Whenever
you add a new feature, please add the corresponding documentation to explain how
your new feature can be used effectively. When you add a significant feature,
also update the file ``ref_features.rst`` and ``ref_releases.rst``.


Building documentation
======================

There is a makefile to generate the documentation based in the source code:

.. code-block:: bash

    cd doc; make html


If you don't want to rebuild the documentations with sphinx every time you make
a change, you can use the `sphinx-autobuild` tool available through PyPI.
Installation is pretty much like any other PyPI package:

.. code-block:: bash

    pip install --user sphinx-autobuild


If you are using `sphinx-autobuild`  the command is as follows:

.. code-block:: bash

    (cd doc; firefox http://localhost:8000 &; make livehtml)

This sets up a server at `localhost:8000` and rebuilds the website whenever you
make a change to the source files. Just like any other process, you can stop it
with `Ctrl-C`

Common issues
=============

When sphinx reports errors or warnings, please fix these in the ``.rst`` files
and the doc strings in the source code. Keep in mind that only the errors and
warnings are shown for files that changed since the last call to ``make html``.
If you want to see all errors and warnings, then run ``make clean; make html``.

The following problems are often encountered:

**Duplicate labels**
    When the same label is defined in two places, they become useless. To avoid
    such name clashes, add a unique prefix to all labels in one rst file. This
    is a bad example (once found in the file ``download_and_install_mac.rst``)::

        .. _compile_install::

    It is safer to use instead::

        .. _mac_compile_install::

**Indentation and empty-line errors**
    RestructeredText is sensitive to indentation and blank lines. For example
    when making bullet points, the following formatting must be used:

    .. code-block:: rst

        Some paragraph before.

        * This is a bullet point with a lot of text that spans several lines.
          Blah blah blah.
        * This is the next

        Some paragraph after

    This won't work:

    .. code-block:: rst

        * This is a bullet point with a lot of text that spans several lines.
        Blah blah blah.
        * This is the next

    This won't work either:

    .. code-block:: rst

        Some paragraph before.
        * This is a bullet point with a lot of text that spans several lines.
          Blah blah blah.
        * This is the next
        Some paragraph after

**Overriden methods in subclasses do not get inherited docstrings**
    Please use the :py:func:`horton.utils.doc_inherit` decorator.

**A signature should be added manually to the `__init__` method in Cython**
    The current behavior of Cython is not compatible with PEP257. The signature gets
    assigned to the class docstring instead of the __init__ docstring. This is a bug in
    Cython that should eventually be fixed. See...
    See https://groups.google.com/forum/#!searchin/cython-users/embedsignature/cython-users/sXtBF2nh5RI/GjcirKOWqBcJ
