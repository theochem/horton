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

Whenever you add a new feature, please add the corresponding documentation to
explain how your new feature can be used effectively. When you add a significant
feature, also update the file ``ref_features.rst`` and ``ref_releases.rst``.
