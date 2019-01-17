..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2019 The HORTON Development Team
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
researcher isn't pristine but rather ranges from *cleverly customized* to
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
able to use them. We have seen problems with several types of dependencies: missing Python
packages, executables, libraries and failing tests.


Python packages
===============

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

* LibXC-3

  When you get the following error message upon running ``./setup.py``, you are trying
  to compile HORTON with LibXC-3:

  .. code-block:: bash

      gcc -pthread -fno-strict-aliasing -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv -DNDEBUG -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv -fPIC -I/usr/lib64/python2.7/site-packages/numpy/core/include -I. -I/usr/include/python2.7 -c horton/meanfield/cext.cpp -o build/temp.linux-x86_64-2.7/horton/meanfield/cext.o -std=c++11
      In file included from /usr/lib64/python2.7/site-packages/numpy/core/include/numpy/ndarraytypes.h:1777:0,
                       from /usr/lib64/python2.7/site-packages/numpy/core/include/numpy/ndarrayobject.h:18,
                       from /usr/lib64/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h:4,
                       from horton/meanfield/cext.cpp:449:
      /usr/lib64/python2.7/site-packages/numpy/core/include/numpy/npy_1_7_deprecated_api.h:15:2: warning: #warning "Using deprecated NumPy API, disable it by " "#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION" [-Wcpp]
       #warning "Using deprecated NumPy API, disable it by " \
        ^~~~~~~
      horton/meanfield/cext.cpp: In function ‘PyObject* __pyx_pf_6horton_9meanfield_4cext_12LibXCWrapper_4refs___get__(__pyx_obj_6horton_9meanfield_4cext_LibXCWrapper*)’:
      horton/meanfield/cext.cpp:2110:74: error: cannot convert ‘func_reference_type* const*’ to ‘const char*’ for argument ‘1’ to ‘PyObject* PyString_FromString(const char*)’
         __pyx_t_1 = __Pyx_PyBytes_FromString((__pyx_v_self->_func.info[0]).refs); if (unlikely(!__pyx_t_1)) __PYX_ERR(0, 117, __pyx_L1_error)
                                                                                ^
      error: command 'gcc' failed with exit status 1

  The solution is to install LibXC-2.2.2. This can always be done by running
  ``./tools/qa/install_libxc-2.2.2.sh`` before running ``./setup.py``.


Failing tests
=============

The following failing tests are symptoms of specific problems:

* ``horton.meanfield.test.test_libxc.test_dot_hessian_o3lyp_cs_polynomial``. This is most
  likely caused by linking against a LibXC that has been compiled with too aggressive
  optimization flags. Use the script ``/toos/qa/install_libxc-2.2.2.sh`` to build a more
  modest version of LibXC, which can then be used to compile HORTON.

* When tests fail due to missing files, you may have installed HORTON from an incomplete
  source tar.gz file, e.g. due to a failed download. In that case, it is best to remove
  your old installation, to make sure no broken files remain, and to reinstall from
  scratch.
