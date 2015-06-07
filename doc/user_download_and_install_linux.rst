..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2015 The HORTON Development Team
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

Download and installation on Linux (Fedora and Ubuntu)
######################################################


Disclaimer
==========

HORTON has been tested on Fedora and Ubuntu Linux. If you are running
any other operating system, some of the instructions below may not work.


Download the code
=================

The latest stable source code release of HORTON can be downloaded here:

    https://github.com/theochem/horton/releases/download/2.0.0/horton-2.0.0.tar.gz

Choose a suitable directory, e.g. ``~/build``, download and unpack the archive

.. code-block:: bash

    mkdir -p ~/build
    cd ~/build
    curl -O https://github.com/theochem/horton/releases/download/2.0.0/horton-2.0.0.tar.gz
    tar -xvzf horton-2.0.0.tar.gz
    cd horton-2.0.0


Dependencies for building, installing and testing HORTON
========================================================

In order to compile and test HORTON, you need to install relatively recent
versions of the following programs/libraries:

* Python >= 2.7, < 3.0: http://www.python.org/ (also install `development files`)
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* GCC, G++ and GFortran >= 4.5: http://gcc.gnu.org/ (In principle the Intel compilers or
  any other of your favorite compilers should work. The GNU compilers are used for
  development and testing.)
* Atlas >= 3.10.1: http://math-atlas.sourceforge.net/ (or any other BLAS implementation that you like more)
* Numpy >= 1.9.1: http://www.numpy.org/
* Scipy >= 0.10.0: http://www.scipy.org/
* Cython >= 0.17.1 : http://www.cython.org/
* h5py >= 2.2.1: http://www.h5py.org/
* Sympy >= 0.7.1: http://code.google.com/p/sympy/
* Matplotlib >= 1.0: http://matplotlib.org/
* LibXC >= 2.2.2: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc
* LibInt2 >= 2.0.3: http://sourceforge.net/p/libint/home
* Curl: http://curl.haxx.se/


Installing dependencies with a package manager
----------------------------------------------

With popular Linux distributions, most of these dependencies can be installed
with a package manager:

* **Ubuntu Linux 15.04**:

  .. code-block:: bash

    sudo apt-get install python-dev gcc g++ gfortran python-numpy python-h5py \
                         python-scipy cython python-nose python-sympy \
                         python-matplotlib libxc-dev libatlas-base-dev curl

  Note that Ubuntu 15.04 does not have a LibInt2 package and will install a
  Numpy version that is too old. Go to the section
  :ref:`linux_manual_dependency_install` to resolve these issues.

* **Fedora Linux 22**:

  .. code-block:: bash

    sudo dnf install python-devel gcc gcc-c++ gcc-gfortran numpy h5py scipy \
                     Cython python-sphinx python-nose sympy python-matplotlib \
                     libint2-devel libxc-devel atlas-devel curl

* **Fedora Linux 20 and 21**:

  .. code-block:: bash

    sudo yum install python-devel gcc gcc-c++ gcc-gfortran numpy h5py scipy \
                     Cython python-sphinx python-nose sympy python-matplotlib \
                     libint2-devel libxc-devel atlas-devel curl

  Note that the Numpy version in Fedora 21 (and earlier releases) is too old. Go
  to the section :ref:`linux_manual_dependency_install` to resolve these issues.


.. _linux_manual_dependency_install:

Installing dependencies manually
--------------------------------

If the package manager of your Linux distribution does not have the desired
packages (or the right versions), you have to install them manually, e.g.
download and execute an installer, or download and unpack a binary package. On
HPC environments a compilation from scratch is recommended.

**BLAS**

In principle, any BLAS implementation may be used. In case of a custom build,
some environment variables must be set prior to building HORTON, as discussed
in :ref:`linux_compile_install`.


**LibXC**

The directory ``depends`` of the HORTON source tree contains a make file that
will download and build LibXC, which will work on most systems:

.. code-block:: bash

    (cd depends; make libxc)

This results in a LibXC library suitable for static linking. If this fails,
consult your local Linux guru to build LibXC. For more info about LibXC, check
the website: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc

**LibInt2**

The directory ``depends`` of the HORTON source tree contains a make file that
will download and build LibInt2, which will work on most systems:

.. code-block:: bash

    (cd depends; make libint -j4)

The compilation of LibInt2 takes a few minutes and results in a library for
static linking. If this fails, consult your local Linux guru to build LibInt2.
For more info about LibInt2, check the website:
http://sourceforge.net/p/libint/home

**Python dependencies**

In some cases, PIP, the Python package manager, may be a good choice to install
the most recent versions of the Python packages in the list of dependencies.
Here are some examples on how to use ``pip`` to install newer versions of
dependencies on Linux distributions that have outdated packages:

* **Ubuntu Linux 15.04 and 14.04**:

  .. code-block:: bash

      sudo apt-get install python-pip
      pip install --user --upgrade numpy

* **Ubuntu Linux 12.04**:

  .. code-block:: bash

      sudo apt-get install python-pip
      pip install --user --upgrade numpy h5py

* **Fedora Linux 20 and 21**:

  .. code-block:: bash

      sudo yum install python-pip
      pip install --user --upgrade numpy


.. _linux_compile_install:

Compilation and installation
============================

Build and install
-----------------

The regular build and install is done as follows:

.. code-block:: bash

    ./setup.py install --user

The ``setup.py`` script makes a reasonable attempt at configuring the compiler and
linker settings for the LibXC, LibInt2 and BLAS libraries. However, this does
not work in all environments. In case of a failure, or if a configuration other
than the default is desired, read the following section.


Overriding default compiler/linker settings for LibXC, LibInt2 and BLAS
-----------------------------------------------------------------------

The manual configuration of the compiler and linker settings is described here:
:ref:`setup_cfg`. You should read this section if the default build and install
has failed or if you would like to specify which libraries to use.


Runtime configuration
---------------------

You need to set the following variable in ``~/.bashrc`` to use HORTON:

.. code-block:: bash

    export PATH=${HOME}/.local/bin:${PATH}

    # If you used special link options for LibXC, LibInt2 or BLAS, something along
    # the following lines may also be needed:
    # export LD_LIBRARY_PATH=some_dir/with/shared_objects/${LD_LIBRARY_PATH}

If you run HORTON on a headless node, i.e. without an X server, you need to
configure Matplotlib to use a backend that does not require a graphical user
interface. (See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for
more details on the Matplotlib backends.) This can be done by adding the
following line to your ``matplotlibrc`` file:

.. code-block:: text

    backend: agg

This file is located in either ``${HOME}/.matplotlib`` or
``${HOME}/.config/matplotlib``.


Running the tests
=================

To test if HORTON was installed properly and if it can be accessed from any directory,
you should change to a directory outside of the source tree and call nosetests
as follows:

.. code-block:: bash

    (cd ~; nosetests -v horton)

Building the documentation
==========================

Dependencies
------------

If you are interested in generating the documentation from source, the following
packages are also needed:

* Sphinx >= 1.3.1: http://sphinx.pocoo.org/
* Doxygen >= 1.8.6: http://www.doxygen.org/
* Breathe >= 1.2.0: http://breathe.readthedocs.org/en/latest/
* Docutils >= 0.11: http://docutils.sourceforge.net/


Installing the dependencies with a package manager and PIP
----------------------------------------------------------

* **Ubuntu Linux 15.04**:

  .. code-block:: bash

      sudo apt-get install doxygen python-docutils python-pip

* **Fedora Linux 22**:

  .. code-block:: bash

      sudo dnf install doxygen python-docutils python-pip

* **Fedora Linux 20 and 21**:

  .. code-block:: bash

      sudo yum install doxygen python-docutils python-pip

Since Breathe (>=1.2.0) and Sphinx (>=1.3.1) may not be available through the
Fedora or Ubuntu repositories, they should be installed manually. For example,
They are available through PyPI.

**PIP packages**:

.. code-block:: bash

    pip install --user --upgrade sphinx breathe

You must also build LibXC statically in the ``depends`` directory, as explained
above, to generate the list of DFT functionals in the documentation.


Actual build
------------

The documentation is compiled and viewed as follows:

.. code-block:: bash

    (cd doc; make html; firefox _build/html/index.html)
