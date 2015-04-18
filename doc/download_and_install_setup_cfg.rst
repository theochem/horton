.. _setup_cfg:

Controlling dynamic/static linking against LibXC, LibInt2 and BLAS
##################################################################

Introduction
============

The script ``setup.py`` also compiles C++ code into
Python extensions. Some of these extension are linked against LibXC, LibInt2,
a BLAS library or a combination of these. The ``setup.py`` script attempts
to detect all compiler and linker settings for these libraries automatically.

For every library, the following attempts are made (in given order) to detect
the compiler and linker flags. If a given step results in such flags, further
steps are not considered.

1. If a file ``setup.cfg`` is present in the root of the source tree, and it
   contains compile and link options for a given library, as documented below,
   the automatic configuration for that library is skipped. One can override the
   settings in ``setup.cfg`` with equivalent environment variables.

2. For LibXC and LibInt2, ``setup.py`` will try to use the static link libraries
   in the ``depends`` directory.

3. An attempt is made to get all compiler and linker flags with the program
   ``pkg-config``. See http://www.freedesktop.org/wiki/Software/pkg-config/ for
   more details.

4. Based on the operating systems and CPU architecture, a ``setup.*.cfg`` file
   is taken from the directory ``data/setup_cfgs``. If the right file is
   present, it is used.

5. A default library name is used for dynamic linking: ``xc``, ``int2`` and
   ``atlas``, for LibXC, LibInt2 and BLAS, respectively.

However, if for some library, these steps do not result in decent
compiler/linker options or the dependencies are not installed, the compilation
will fail and you will get a compiler error message.

The following sections explain how one can override the default guesses of the
``setup.py`` script.


``setup.cfg`` and environment variables
=======================================

For each library, a section can be added to the file ``setup.cfg`` to configure
compiler and linker options. By default, no such file is present in the root
of the source tree, so you have to create a new one. Several examples can be
found in ``data/setup_cfgs``. For example, this is the file to compile the
extensions properly on a 64 bit version of Fedora 21:

.. literalinclude:: ../data/setup_cfgs/setup.Linux-Fedora-21-x86_64.cfg

In each section one can define the following variables: ``include_dirs``,
``library_dirs``, ``libraries``, ``extra_objects``, ``extra_compile_args`` and
``extra_link_args``. They correspond to the optional arguments of the
``Extension`` class of the Distutils package, as described here:
https://docs.python.org/2/distutils/setupscript.html#describing-extension-modules
The purpose of each keyword is summarized below:

``include_dirs``
    A list of non-standard directories containing C/C++ header files

``library_dirs``
    A list of non-standard directories with shared objects (dynamic link
    libraries). When any of the ``*_LIBRARY_DIRS`` variables or the ``-L``
    option of ``setup.py build_ext`` are used, keep in mind that at runtime, the
    dynamic link loader must be informed of these directories. For example, on
    Linux, the variable ``LD_LIBRARY_PATH`` must be set accordingly.

``libraries``
    A list of shard objects to link with.

``extra_objects``
    Extra object files for static linking.

``extra_compile_args``
    Extra arguments given to the compiler when compiling source code.

``extra_link_args``
    Extra arguments given to the compiler when linking object files.


Multiple paths or names for one keyword are separated by a colon (``:``).

Instead of setting the library configuration in the file ``setup.cfg``, one may
also set the following environment variables:

* ``LIBXC_LIBRARY_DIRS``
* ``LIBXC_LIBRARIES``
* ``LIBXC_EXTRA_OBJECTS``
* ``LIBXC_EXTRA_COMPILE_ARGS``
* ``LIBXC_EXTRA_LINK_ARGS``
* ``LIBINT2_LIBRARY_DIRS``
* ``LIBINT2_LIBRARIES``
* ``LIBINT2_EXTRA_OBJECTS``
* ``LIBINT2_EXTRA_COMPILE_ARGS``
* ``LIBINT2_EXTRA_LINK_ARGS``
* ``BLAS_LIBRARY_DIRS``
* ``BLAS_LIBRARIES``
* ``BLAS_EXTRA_OBJECTS``
* ``BLAS_EXTRA_COMPILE_ARGS``
* ``BLAS_EXTRA_LINK_ARGS``


Other ways of controlling compilation and linker flags
======================================================

Instead of the library-specific variables above, there are also general methods
to configure the compiler and linker.

One may specify a colon-separated list of non-standard directories with
include files in the ``CPATH`` environment variable. For example::

    export CPATH=/usr/include/libint2

In addition, general compiler and linker flags may be set in the ``build_ext``
stage of the installation process. To do so, the installation must be done
in two steps::

    ./setup.py build_ext EXTRA OPTIONS HERE
    ./setup.py install --user

Run ``./setup.py build_ext --help`` for a complete list of options. These
options apply to all extension, so avoid them for static linking.
