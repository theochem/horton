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

Troubleshooting
###############

On a fresh Linux or Mac machine, the instructions given in the "Download and
Install" guides should result in a working installation of HORTON, without using
any of the suggestions below. In reality, however, the Unix system of an average
researcher isn't pristine but rahter ranges from *cleverly customized* to
*completely borked*. Such customizations may interfere with the installation of
HORTON. This section provides some guidance for the novice Unix user on how to
get HORTON working on a not-so-well-maintained Unix system.

If you are still stuck after trying the suggestions in this section, do not
hesitate to contact us on the `the HORTON mailing list
<https://groups.google.com/forum/#!forum/horton-discuss>`_.


Introduction
============

Let us assume you have already built and installed all your dependencies.
However, when you try to install HORTON, i.e.

.. code-block:: bash

    ./setup.py install --user

or when you run ``nosetests``, you get an unexpected error message. The problem
is most likely related to finding and using the dependencies. You have to make
sure ``setup.py`` and the HORTON modules can find the right dependencies and are
able to use them. We have seen problems with three types of dependencies: Python
modules, executables, and libraries.


Python modules
==============

If you have installed a python package (e.g. NumPy, SciPy, Cython, H5Py,
SymPy, MatPlotLib, Nosetests, Sphinx, Breathe, Docutils) and you get an error
saying your system cannot find that package, then you need to check the
directories in which Python searches for package. These are stored
in the attribute ``path`` of the ``sys`` module, which can be accessed by:

.. code-block:: bash

    python -c "import sys; import pprint; pprint.pprint(sys.path)"

A typical output can be as follows (but is probably different in your case):

.. code-block:: python

    ['',
     '/usr/lib64/python27.zip',
     '/usr/lib64/python2.7',
     '/usr/lib64/python2.7/plat-linux2',
     '/usr/lib64/python2.7/lib-tk',
     '/usr/lib64/python2.7/lib-old',
     '/usr/lib64/python2.7/lib-dynload',
     '/home/foo/.local/lib/python2.7/site-packages',
     '/usr/lib64/python2.7/site-packages',
     '/usr/lib64/python2.7/site-packages/gtk-2.0',
     '/usr/lib/python2.7/site-packages']

If you installed a Python package in another directory, Python will not be able
to load it. This can be fixed by adding your directory to the ``PYTHONPATH``
variable in your ``${HOME}/.bashrc`` (Linux) or ``${HOME}/.bash_profile`` (Mac),
e.g.

.. code-block:: bash

    export PYTHONPATH=/some/custom/path/lib/python2.7/site-packages:${PYTHONPATH}

The next time you start Python (or any program implemented in with Python), the
packages you installed in a non-standard location will become importable. If the
same Python module or package is installed in multiple directories, the one
found in the first directory in the ``sys.path`` list takes precedence.

A typical problem is that there are multiple lines like these in ``.bashrc`` or
``.bash_profile`` of which the last is overwriting the former ones, e.g.:

.. code-block:: bash

    export PYTHONPATH=/some/custom/path/lib/python2.7/site-packages:${PYTHONPATH}
    # and several lines further ...
    export PYTHONPATH=/some/other/path/lib/python2.7/site-packages

The second ``export`` line overrides the first one because it does not end with
``:${PYTHONPATH}``.

Some examples are given below. Note that, in principle, none of these should be
necessary but they seem to have helped people with a broken installation of
Python:

* Some Mac users needed to set the ``PYTHONPATH`` after installing modules
  through PIP:

  .. code-block:: bash

      export PYTHONPATH=${HOME}/Library/Python/2.7/lib/python/site-packages:${PYTHONPATH}

  or their system site-packages:

  .. code-block:: bash

      export PYTHONPATH=/Library/Python/2.7/lib/python/site-packages:${PYTHONPATH}

* Similarly, a few Linux users needed to set ``PYTHONPATH`` after installation
  through PIP:

  .. code-block:: bash

      export PYTHONPATH=${HOME}/.local/lib/python2.7/site-packages:${PYTHONPATH}

  or

  .. code-block:: bash

      export PYTHONPATH=/lib/python2.7/site-packages:${PYTHONPATH}

  or

  .. code-block:: bash

      export PYTHONPATH=/lib64/python2.7/site-packages


Excecutables
============

During the installation (or when building the documentation) HORTON will use
some executables, e.g. a compiler, ``sphinx-build``, etc. These executables must
be in one of the directories in the ``PATH`` environment variable. The essential
changes to the ``PATH`` variable were already discussed in the "Download and
install" guides but if your system is somehow broken, more changes may be
needed.

The contents of ``PATH`` can be accessed by:

.. code-block:: bash

    echo $PATH

In unfavorable circumstances, some directories may be missing from the ``PATH``,
e.g because it got carelessly overwritten in ``${HOME}/.bashrc`` (Linux) or
``${HOME}/.bash_profile`` (Mac). For example, the following should be avoided:

.. code-block:: bash

    export PATH=/some/custom/path/bin

Instead, make sure the existing ``PATH`` variable is included as follows:

.. code-block:: bash

    export PATH=/some/custom/path/bin:${PATH}

If the same executable name occurs in several directories in the ``PATH``, the
one in the first directory takes precedence.

The following examples are in principle not needed but they seemed to be helpful
for some:

* Mac users that uses python scripts might do

  .. code-block:: bash

      # Already mentioned in "Download and install" guide:
      export PATH=${HOME}/Library/Python/2.7/bin:${PATH}
      # Should already be in the PATH anyway, unless your system is broken:
      export PATH=/Library/Python/2.7/bin:${PATH}

* Similarly, Linux users may do

  .. code-block:: bash

      # Already mentioned in "Download and install" guide:
      export PATH=${HOME}/.local/bin:${PATH}
      # Should already be in the PATH anyway, unless your system is broken:
      export PATH=/usr/bin:${PATH}

When you forgot where you installed a dependency, the ``find`` command may help
you find the appropriate directory. The following example will search for
location of the ``sphinx-build`` executable:

.. code-block:: bash

    find / | grep sphinx-build


Libraries
=========

You need to make sure ``setup.py`` can find the necessary libraries. You should
consult :ref:`setup_cfg` for a more complete understanding of the library
linking process when installing HORTON. Here, we will show how we solved some
library problems we encountered before.

First, you need to locate the library that can not be found by ``setup.py``. You
can locate libraries in standard directories by using the unix command
``ldconfig``:

.. code-block:: bash

    ldconfig -p | grep libraryname

``ldconfig -p`` prints all cached libraries, and piping to ``grep`` searches through
the results for the library with the ``libraryname``. This only works when a
library is installed in a standard location and the library cache is up-to-date.
If you can not find it with ``ldconfig``, you may try to used the ``find``
command, e.g.:

.. code-block:: bash

    find / | grep libraryname

Here is an example that searches for the Atlas libraries on a cluster:

.. code-block:: bash

    ldconfig -p | grep atlas

which gives

.. code-block:: bash

    libptf77blas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libptf77blas.so.3
    libptf77blas.so (libc6,x86-64) => /usr/lib64/atlas/libptf77blas.so
    libptcblas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libptcblas.so.3
    libptcblas.so (libc6,x86-64) => /usr/lib64/atlas/libptcblas.so
    liblapack.so.3 (libc6,x86-64) => /usr/lib64/atlas/liblapack.so.3
    liblapack.so (libc6,x86-64) => /usr/lib64/atlas/liblapack.so
    libf77blas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libf77blas.so.3
    libf77blas.so (libc6,x86-64) => /usr/lib64/atlas/libf77blas.so
    libclapack.so.3 (libc6,x86-64) => /usr/lib64/atlas/libclapack.so.3
    libclapack.so (libc6,x86-64) => /usr/lib64/atlas/libclapack.so
    libcblas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libcblas.so.3
    libcblas.so (libc6,x86-64) => /usr/lib64/atlas/libcblas.so
    libatlas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libatlas.so.3
    libatlas.so (libc6,x86-64) => /usr/lib64/atlas/libatlas.so

All the libraries are located in ``/usr/lib64/atlas/``. Notice that all the
libraries use the x86-64 instruction set.

Next, we need to find the include directory. You can find this with the ``find``
command function. Usually, the include directory is almost same as the library
directory, except instead of the ``lib`` or ``lib64``, it reads ``include``.
Continuing the above example,

.. code-block:: bash

    ls -d /usr/include/*atlas*

will give the list of directories that includes the word ``atlas``. The output
gives:

.. code-block:: bash

   /usr/include/atlas
   /usr/include/atlas-x86_64-base

Since we used the x86-64 instruction set, we select the directory that would
correspond with that instruction set, i.e. ``/usr/include/atlas-x86_64-base``.
(This should not matter too much as header files are normally indepedent of the
architecture.)

In the above list of libraries associated with atlas, we have ``ptf77blas``,
``ptcblas``, ``lapack``, ``f77blas``, ``clapack``, ``cblas``, and ``atlas``.
Though we can include all these libraries, HORTON only uses ``atlas`` and
``cblas``. Therefore, the resulting ``setup.cfg`` file includes

.. code-block:: bash

  [blas]
  library_dirs=/usr/lib64/atlas
  libraries=atlas:cblas
  include_dirs=/usr/include/atlas-x86_64-base

Similarly, we can repeat the process for the LibXC and Libint2, where the
libraries that are needed are only ``libxc`` and ``libint``, respectively. See
:ref:`setup_cfg` for more details.
