Getting started
###############


.. contents::


Disclaimer
==========

Horton is mainly developed and tested on Linux systems. If you run any other
operating system, some of the instructions below may not work.


Download the code
=================

Stable release (recommended)
----------------------------

The latest stable source code release of Horton can be downloaded here:
http://users.ugent.be/~tovrstra/horton/horton-1.0.tar.gz. Choose a suitable
directory, e.g. ``~/build``, download and unpack the archive::

    mkdir -p ~/build
    cd ~/build
    wget http://users.ugent.be/~tovrstra/horton/horton-1.0.tar.gz
    tar -xvzf horton-1.0.tar.gz
    cd horton-1.0


Latest development code (experts only)
--------------------------------------

In order to get the latest development version of the source code, and to upload
your own changes, you need to work with git. Git is a version control system
that makes life easy when a group of people are working on a common source code.
All information about git (including downloads and tutorials) can be found here:
http://git-scm.com/. The official git URL of Horton is:
git://github.com/theochem/horton.git. In order to `clone` the public Horton
repository, run this command::

    git clone git://github.com/theochem/horton.git
    cd horton

The version history can be updated with the latest patches with the following
command::

    git pull

There is also a web interface to Horton's git repository:
https://github.com/theochem/horton


Common dependencies
===================

In order to compile and test Horton, one needs to
install relatively recent versions of the following programs/libraries:

* Python >= 2.7, < 3.0: http://www.python.org/ (also install `development files`)
* GCC, G++ and GFortran >= 4.5: http://gcc.gnu.org/ (In principle the Intel compilers or any other of your favorite compilers should work. The GNU compilers are used for development and testing.)
* Numpy > 1.0: http://www.scipy.org/
* h5py >= 2.1.1: http://code.google.com/p/h5py/
* Scipy >= 0.10.0: http://www.scipy.org/
* Cython >= 0.17.1 : http://www.cython.org/
* Sphinx > 1.0: http://sphinx.pocoo.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* Sympy >= 0.7.1: http://code.google.com/p/sympy/
* The ``patch`` program >= 2.0: http://savannah.gnu.org/projects/patch/ (or any of its equivalents)

Optionally, one may also install Matplotlib:

* Matplotlib >= 1.0: http://matplotlib.org/

On a decent operating system, these programs/libraries can be easily installed
with a package manager. First check this possibility before manually installing
the dependencies.

On Ubuntu Linux::

    sudo apt-get install python-dev gcc g++ gfortran python-numpy python-h5py python-scipy cython python-sphinx python-nose python-sympy patch python-matplotlib

On Fedora Linux::

    sudo yum install python-devel gcc gcc-c++ gcc-gfortran numpy h5py scipy Cython python-sphinx python-nose sympy patch python-matplotlib

If the package manager of your operating system does not have the desired
packages (or the right versions), one has to install them manually, e.g.
download and execute an installer, or download and unpack a binary package. On
HPC environments a compilation from scratch is recommended. In some cases, Pip,
the Python package manager, may be a good choice to install the most recent
versions of the Python packages in the list of dependencies. Assuming that you
have installed some compilers, the Python development files and HDF5 development
files, the following command installs the remaining dependencies in your home
directory::

    pip install --user numpy scipy cython h5py sphinx nose sympy


Specific dependencies
=====================

The directory ``depends`` of the Horton source tree is used to build specific
dependencies from source. For the moment, there are two such dependencies,
namely `libint2 <http://sourceforge.net/p/libint/>`_ and `libxc
<http://www.tddft.org/programs/octopus/wiki/index.php/Libxc>`_
[marques2012]_. The directory ``depends``
contains a ``Makefile`` that takes care of downloading the right version and
compiling it. The following should get you both libint and libxc::

    cd depends
    make libint
    make libxc
    cd ..

The compilation of libint takes a few minutes.


Reference atoms
===============

Several parts of Horton make use of reference atomic computations. These files
are too large to be included in the revision system. Therefore they must be
downloaded separately as follows::

    cd data/refatoms
    wget http://users.ugent.be/~tovrstra/horton/refatoms.tar.bz2
    tar -xvjf refatoms.tar.bz2
    cd ../..


Compilation and installation
============================

One may either install Horton in the home directory, or perform an in-pace build
in the source tree.

The **regular build+install** is as done follows::

    ./setup.py install --user

The ``horton-*.py`` scripts are installed in ``~/.local/bin`` and you have to
add this directory to your ``PATH`` environment variable to make them accessible
from any directory.

The **in-place build** is useful for testing purposes, and is done as follows::

    ./setup.py build_ext -i

The documentation is compiled and viewed as follows::

    cd doc
    make html
    firefox _build/html/index.html
    cd ..


Testing
=======

A bunch of validation routines are included in Horton. To run the tests, first
perform an **in-place build** and run ``nosetests`` afterwards::

    ./setup.py build_ext -i
    nosetests -v

If all tests pass, the screen output should end with ``OK``. If at some point,
something during the build process fails, clean up the source tree with the
``cleanfiles.sh`` script and try again. If you then still have problems,
please report them on the `Horton mailing list
<https://groups.google.com/d/forum/horton-discuss>`_.
