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

.. _setup_cfg:

Controlling dynamic/static linking against LibXC, LibInt2 and BLAS
##################################################################

Introduction
============

The script ``setup.py`` also compiles C++ code into Python extensions. Some of
these extension are linked against LibXC, LibInt2, BLAS library or a combination
of these. The ``setup.py`` script attempts to detect all of the compiler and
linker settings for these libraries automatically.

For **each** library, following attempts are made (in the given order) to detect
the compiler and linker flags. If the attempt succeeds, i.e. gives a set of
satisfactory flags, further steps are not considered for that library.\

1. Read ``setup.cfg`` in the root of the source tree for compilation and link
   options for that library, as documented in :ref:`cfgfile`. Read in
   environment variables for compilation and link options for that library.

2. Statically link libraries in the ``depends`` directory (for LibXC and LibInt2).

3. Use program ``pkg-config`` to get all of the compiler and linker flags
   for each library. See http://www.freedesktop.org/wiki/Software/pkg-config/ for
   more details.

4. Read ``setup.*.cfg`` file in the directory ``data/setup_cfgs`` that
   corresponds to your operating system and CPU architecture for compilation and
   link options for that library.

5. Dynamically link libraries using default library names: ``xc``, ``int2`` and
   ``atlas``, for LibXC, LibInt2 and BLAS, respectively.

However, if these steps do not result in decent compiler/linker options for some
library or the dependencies are not installed, then the compilation will fail
and you will get a compiler error message.

The following sections explain how you can override the default guesses of the
``setup.py`` script.

.. _cfgfile:

``setup.cfg`` and environment variables
=======================================

For each library, a section can be added to the file ``setup.cfg`` to configure
the compiler and linker options. By default, no such file is present in the root
of the source tree, so you have to create a new one. Several examples can be
found in ``data/setup_cfgs``. For example, this file would compile and link the
extensions on a 64 bit version of Fedora 21:

.. literalinclude:: ../data/setup_cfgs/setup.Linux-Fedora-21-x86_64.cfg
    :caption: data/setup_cfgs/setup.Linux-Fedora-21-x86_64.cfg


This is a setup.cfg file for a compilation with the Intel MKL libraries on
Fedora:

.. literalinclude:: ../data/setup_cfgs/setup.Linux-Fedora-MKL.cfg
    :caption: data/setup_cfgs/setup.Linux-Fedora-MKL.cfg

In each section you can define the following variables: ``include_dirs``,
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
    A list of shared objects to link with.

``extra_objects``
    Extra object files for static linking.

``extra_compile_args``
    Extra arguments given to the compiler when compiling source code.

``extra_link_args``
    Extra arguments given to the compiler when linking object files.


Multiple paths or names for one keyword are separated by a colon (``:``).

Instead of setting the library configuration in the file ``setup.cfg``, you may
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

You may specify a colon-separated list of non-standard directories with
include files in the ``CPATH`` environment variable. For example::

    export CPATH=/usr/include/libint2

In addition, general compiler and linker flags may be set in the ``build_ext``
stage of the installation process. To do so, the installation must be done
in two steps::

    ./setup.py build_ext EXTRA OPTIONS HERE
    ./setup.py install --user

Run ``./setup.py build_ext --help`` for a complete list of options. These
options apply to all extensions, so avoid them for static linking.
