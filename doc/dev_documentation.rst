Writing documentation
#####################

All the documentation is located in the ``doc`` directory. We use the `Sphinx
<http://sphinx.pocoo.org/>`_ formatting engine to compile the `documentation
source code` into fancy formatted HTML or PDF.

The source files have the extension ``.rst``, and are written in the
`ReStructuredText <http://docutils.sourceforge.net/rst.html>`_ (RST) format.
RST is in some sense comparable to latex, but more intuitive to use.
It also has some specific advantages for documenting software.

All ``.rst``-files are part of the source tree, just like the actual source
code. Git is also used to keep track of changes in the documentation.

There is a makefile to generate the documentation based in the source code::

    toony@poony ~/.../horton:master> cd doc
    toony@poony ~/.../horton/doc:master> make html
    toony@poony ~/.../horton/doc:master> make pdf

Whenever you add a new feature, make sure that at least files
``lib_horton*.rst`` are up to date. With more serious work, please also write
a tutorial, e.g. like this one, to explain how your new feature can be used
effectively. If you added a significant feature, also update the file
``ref_features.rst``.
