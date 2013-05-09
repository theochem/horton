Getting started
###############

Downloading the code
====================

In order to get the latest version of the source code, and to upload your own
changes, you need to work with git. Git is a version control system that
makes life easy when a group of people are working on a common source code. It
also makes sense to use it for personal projects. All information about git
(including downloads and tutorials) can be found here: http://git-scm.com/. The
official git URL of Horton is::

    git://github.com/theochem/horton.git

In order to `clone` the public Horton repository, run this command::

    git clone git://github.com/theochem/horton.git

The version history can be updated with the latest patches with the following
command::

    git pull

There is also a web interface to Horton's git repository:
https://github.com/theochem/horton

Common dependencies
===================

In order to compile and test Horton and its documentation, one needs to
install relatively recent versions of the following programs/libraries:

* Python >= 2.7, < 3.0: http://www.python.org/ (also install `development files`)
* Numpy > 1.0: http://www.scipy.org/
* h5py >= 2.1.1: http://code.google.com/p/h5py/
* Scipy >= 0.10.0: http://www.scipy.org/
* Cython >= 0.17.1 : http://www.cython.org/
* Sphinx > 1.0: http://sphinx.pocoo.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* Sympy >= 0.7.1: http://code.google.com/p/sympy/
* Matplotlib >= 1.0: http://matplotlib.org/

On a decent operating system, these programs/libraries can be easily installed
with a package manager. First check that option before manually installing this
software.

On Ubuntu Linux::

    sudo apt-get install python-dev python-numpy python-h5py python-scipy cython python-sphinx python-nose python-sympy python-matplotlib

On Fedora Linux::

    sudo apt-get install python-devel numpy h5py scipy Cython python-sphinx python-nose sympy python-matplotlib


Specific dependencies
=====================

The directory ``depends`` of the Horton source tree is used to build specific
dependencies from source. For the moment, there are two such dependencies,
namely `libint2 <http://sourceforge.net/p/libint/>`_ and libxc. The directory ``depends``
contains a ``Makefile`` that can take care of downloading the right version and
compiling. The following should take care of everything, assuming that you have
installed all the libint2 dependencies::

    cd depends
    make -j4 libint
    make -j4 libxc
    cd ..

The option ``-j4`` instructs make to build code in parallel on four cores.


Reference atoms
===============

Several parts of Horton make use of reference atomic computations. These files
are too large to be included in the revision system. Therefore they must be
downloaded separately as follows::

    cd data/refatoms
    wget http://users.ugent.be/~tovrstra/horton_refatoms.tar.bz2
    tar -xvjf horton_refatoms.tar.bz2
    cd ..


Compilation and installation
============================

One may either install Horton in the home directory, or perform an in-pace build
in the source tree.

The **regular build+install** is as done follows::

    ./setup.py install --user

This will put a `compiled` version of Horton in the ``build`` directory. Then a
copy of this directory is placed under
``~/.local/lib/python2.*/site-packages/horton``.

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
``cleanfiles.sh`` script and try again.


Basic example
=============

This is a basic example computation in Horton. The input file is just
a small Python main program that uses the Horton library. The script
``examples/001_hf_water/run.py`` performs a HF/3-21G computation on water and
partitions the density with the Becke scheme:

.. literalinclude:: ../examples/001_hf_water/run.py

The molecular geometry is loaded from an XYZ file. It is also possible (yet less
convenient) to directly type the atomic coordinates in the python script.

In addition to this Python interface of Horton, several Python scripts were also
developed (``horton-*.py``) that are easier to use but only have a limited
functionality compared to the Python interface
