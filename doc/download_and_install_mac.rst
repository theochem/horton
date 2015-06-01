..
    : Horton is a development platform for electronic structure methods.
    : Copyright (C) 2011-2015 The Horton Development Team
    :
    : This file is part of Horton.
    :
    : Horton is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : Horton is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

Download and Installation on Mac OS X (10.8 - 10.10)
####################################################

Disclaimer
==========

Horton has been tested on Mac OS X 10.8 - 10.10 using MacPorts. If you
are running any other version of OS X or using other package managers,
some of the instructions below may not work.


MacPorts
=========

We strongly recommend that you install all of the packages required by Horton
through MacPorts. We also advise you to uninstall/remove all python packages
installed through package managers other than MacPorts (e.g., Canopy).

The latest version of MacPorts can be downloaded from the web:
https://www.macports.org/install.php. This guide has been tested using
MacPorts 2.3.3 but should also work with newer versions.


Quick guideline through MacPorts
--------------------------------

Here are some basic MacPort commands:

* updating ports (recommended):

.. code-block:: bash

    sudo port -v selfupdate

* upgrade ports:

.. code-block:: bash

    sudo port upgrade outdated

* finding ports (e.g, port_name = python27):

.. code-block:: bash

    sudo port list python27

* searching ports (e.g, port_name = python27):

.. code-block:: bash

    port search python27

* installing ports (e.g, port_name = python27):

.. code-block:: bash

    sudo port install python27

* selecting ports (e.g, select python27 as python):

.. code-block:: bash

    sudo port select --set python python27


Download the code
=================

Stable release (recommended)
----------------------------

The latest stable source code release of Horton can be downloaded here:

    https://github.com/theochem/horton/releases/download/2.0.0/horton-2.0.0.tar.gz

Choose a suitable directory, e.g. ``~/build``, download and unpack the archive:

.. code-block:: bash

    mkdir -p ~/build
    cd ~/build
    curl -O https://github.com/theochem/horton/releases/download/2.0.0/horton-2.0.0.tar.gz
    tar -xvzf horton-2.0.0.tar.gz
    cd horton-2.0.0


Latest development code (experts only)
--------------------------------------

In order to get the latest development version of the source code, and to upload
your own changes, you need to work with git. Git is a version control system
that makes life easy when a group of people are working on a common source code.
All information about git (including downloads and tutorials) can be found here:
http://git-scm.com/. The official public git URL of Horton is:
``git://github.com/theochem/horton.git``. Git can be installed through MacPorts:

.. code-block:: bash

    port install git

In order to `clone` the public Horton repository, run these commands:

.. code-block:: bash

    mkdir -p ~/build
    cd ~/build
    git clone git://github.com/theochem/horton.git
    cd horton

The version history can be updated with the latest patches with the following
command:

.. code-block:: bash

    git pull

There is also a web interface to Horton's git repository:
https://github.com/theochem/horton


Dependencies for building, installing and testing Horton
========================================================

In order to compile and test Horton, you need to install relatively recent
versions of the following programs/libraries:

* GCC, G++ and GFortran >= 4.5: http://gcc.gnu.org/
* Python >= 2.7, < 3.0: http://www.python.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* Atlas >= 3.10.1: http://math-atlas.sourceforge.net/ (or any other BLAS implementation that you like more)
* Numpy >= 1.7.0: http://www.numpy.org/
* Scipy >= 0.10.0: http://www.scipy.org/
* Cython >= 0.17.1 : http://www.cython.org/
* h5py >= 2.2.1: http://www.h5py.org/
* Sympy >= 0.7.1: http://code.google.com/p/sympy/
* Matplotlib >= 1.0: http://matplotlib.org/
* LibXC >= 2.2.2: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc
* LibInt2 >= 2.0.3: http://sourceforge.net/p/libint/home


Installing the dependencies with MacPorts
-----------------------------------------

All dependencies can be installed with MacPorts. We recommend
the following ports:

* ``gcc49``, https://trac.macports.org/browser/trunk/dports/lang/gcc47/Portfile
* ``python27``, https://trac.macports.org/browser/trunk/dports/lang/python27/Portfile
* ``py27-nose``, https://trac.macports.org/browser/trunk/dports/python/py-nose/Portfile
* ``atlas``, https://trac.macports.org/browser/trunk/dports/math/atlas/Portfile
* ``py27-numpy +atlas`` (Numpy with Atlas support), https://trac.macports.org/browser/trunk/dports/python/py-numpy/Portfile
* ``py27-scipy +atlas`` (SciPy with Atlas support), https://trac.macports.org/browser/trunk/dports/python/py-scipy/Portfile
* ``py27-cython``, https://trac.macports.org/browser/trunk/dports/python/py-cython/Portfile
* ``py27-h5py``, https://trac.macports.org/browser/trunk/dports/python/py-h5py/Portfile
* ``py27-sympy``, https://trac.macports.org/browser/trunk/dports/python/py-sympy/Portfile
* ``py27-matplotlib``, https://trac.macports.org/browser/trunk/dports/python/py-matplotlib/Portfile
* ``libxc``, https://trac.macports.org/browser/trunk/dports/science/libxc/Portfile
* ``libint``, https://trac.macports.org/browser/trunk/dports/science/libint/Portfile

These are installed with the following commands. (When MacPorts is installed in user
space, the ``sudo`` can be omitted.):

.. code-block:: bash

    sudo port install gcc49
    sudo port select --set gcc mp-gcc49
    sudo port install python27
    sudo port select --set python python27
    sudo port install py27-nose
    sudo port select --set nosetests nosetests27
    sudo port install atlas
    sudo port install py27-numpy +atlas
    sudo port install py27-scipy +atlas
    sudo port install py27-cython
    sudo port select --set cython cython27
    sudo port install py27-h5py
    sudo port install py27-sympy
    sudo port select --set py-sympy py27-sympy
    sudo port install py27-matplotlib
    sudo port install libxc
    sudo port install libint

The GNU compilers are only used to compile Fortran code as the default C/C++
compiler on the Mac is ``clang``.

.. _mac_manual_dependency_install:

Installing dependencies manually
--------------------------------

**BLAS**

In principle, any BLAS implementation may be used. In case of a custom build,
some environment variables must be set prior to building Horton, as discussed
in :ref:`mac_compile_install`. Also, Keep in mind that MacPorts only supports Atlas
for building NumPy and SciPy.


**LibXC**

The directory ``depends`` of the Horton source tree contains a make file that
will download and build LibXC, which will work on most systems:

.. code-block:: bash

    (cd depends; make libxc)

This results in a libxc library suitable for static linking. If this fails,
consult your local Mac guru to build LibXC. For more info about LibXC, check
the website: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc

**LibInt2**

The directory ``depends`` of the Horton source tree contains a make file that
will download and build LibInt2, which will work on most systems:

.. code-block:: bash

    (cd depends; make libint -j4)

The compilation of libint takes a few minutes and results in a library for
static linking. If this fails, consult your local Mac guru to build LibInt2.
For more info about LibInt2, check the website:
http://sourceforge.net/p/libint/home


Reference atoms
===============

This step can be skipped when compiling the stable release because each stable
release already contains reference atoms.

Several parts of Horton make use of reference atomic computations. These files
are too large to be included in the git revision system. Therefore, they must be
downloaded separately when compiling a development version of Horton:

.. code-block:: bash

    (cd data/refatoms; make all)

.. _mac_compile_install:

Compilation and installation
============================

Build and install
-----------------

The regular build and install is done as follows:

.. code-block:: bash

    ./setup.py install --user

The ``setup.py`` script makes a reasonable attemp configuring the compiler and
linker settings for the LibXC, LibInt2 and BLAS libraries. However, this does
not work in all environments. In case of a faillure, or if a configuration other
than the default is desired, read the following section.


Overriding default compiler/linker settings for LibXC, LibInt2 and BLAS
-----------------------------------------------------------------------

The manual configuration of the compiler and linker settings is described here:
:ref:`setup_cfg`. You should read this section if the default build and install
failed or if you would like to specify which libraries to use.


Runtime Configuration
---------------------

You need to set some environment variables to use Horton. Add the following to
``~/.bash_profile`` if it exists, otherwise add them to ``~/.profile``:


.. code-block:: bash

    export PATH=${HOME}/Library/Python/2.7/bin:${PATH}
    # I did not have to set the following two.
    # The --user option of the setup.py script normally installs stuff in a place
    # where Python will find it without setting environment variables. ~Toon
    export PYTHONPATH=${PYTHONPATH}:${HOME}/path-to-horton-installation/
    export HORTONDATA=${HOME}/path-to-horton-installation/data/

If you run Horton on a headless node, i.e. without an X server, you need to
configure Matplotlib to use a backend that does not require a graphical user
interface. (See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for
more details on the Matplotlib backends.) This can be done by adding the
following line to your ``matplotlibrc`` file:

.. code-block:: text

    backend: agg

This file is located either in ``${HOME}/.matplotlib`` or
``${HOME}/.config/matplotlib``.


Running the tests
=================

To test that Horton was installed properly and that you can can access it from
other directories, you should change to a directory outside of the source tree
and call nosetests as follows:

.. code-block:: bash

    (cd ~; nosetests -v horton)


Building the documentation
==========================

Dependencies
------------

If you are interested in generating the documentation from source, the following
packages are also needed:

* PIP >= 6.1.1: https://pypi.python.org/pypi/pip
* Sphinx >= 1.3.1: http://sphinx.pocoo.org/
* Doxygen >= 1.8.6: http://www.doxygen.org/
* Breathe >= 1.2.0: http://breathe.readthedocs.org/en/latest/
* Docutils >= 0.11: http://docutils.sourceforge.net/


Installing the dependencies with MacPorts and PIP
-------------------------------------------------

Most can be installed directly with MacPorts. The following list of ports is recommended:

* ``doxygen``: https://trac.macports.org/browser/trunk/dports/textproc/doxygen/Portfile
* ``py27-pip``: https://trac.macports.org/browser/trunk/dports/python/py-pip/Portfile

The following commands will install the ports:

.. code-block:: bash

    sudo port install doxygen
    sudo port install py27-pip
    sudo port select --set pip pip27

Since Breathe (>=1.2.0) and Sphinx (>=1.3.1) may not be available through
MacPort, they should be installed through PIP:

.. code-block:: bash

    pip install --user --upgrade sphinx breathe

You must also build LibXC statically in the ``depends`` directory, as explained
above, to generate the list of DFT functionals in the documentation.


Actual build
------------

The documentation is compiled and viewed as follows:

.. code-block:: bash

    (cd doc; make html; open _build/html/index.html)


Common Problems
===============

* If you get errors saying that you do not have a file or script when you have
  clearly installed it beforehand, it may not be named appropriately. You can fix
  this by symbolically linking that file to the appropriate name. E.g.::

     ln -s something something

* If you get an error, saying you have not installed xcode, install xcode.
  E.g.::

    some example
