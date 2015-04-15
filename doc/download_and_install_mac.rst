Getting started on Mac OS X (10.8 - 10.10)
######################################################


.. contents::


Disclaimer
==========

Horton has been tested on Mac OS X 10.8 - 10.10 using MacPorts. If you run any other version of OS X, 
some of the instructions below may not work. 

MacPorts
=========

We strongly recommend to install all the packages required by Horton through MacPorts. We also advise to 
uninstall/remove all python packages installed outside MacPorts (e.g., Canopy). 

The latest version of MacPorts can be downloaded from the web:

https://www.macports.org/install.php

Quick guideline through MacPorts
--------------------------------

Here are some basic MacPort commands:

* updating ports (recommended)::

    sudo port -v selfupdate 

* finding ports (e.g, port_name = python27)::
   
    sudo port list python27

* searching ports (e.g, port_name = python27)::
   
    port search python27

* installing ports (e.g, port_name = python27)::
   
    sudo port install python27

* selecting ports (e.g, select python27 as python)::
   
    sudo port select --set python python27 

Download the code
=================

Stable release (recommended)
----------------------------

The latest stable source code release of Horton can be downloaded here:

    https://github.com/theochem/horton/releases/download/1.2.1/horton-1.2.1.tar.gz

Choose a suitable directory, e.g. ``~/build``, download and unpack the archive::

    mkdir -p ~/build
    cd ~/build
    wget https://github.com/theochem/horton/releases/download/1.2.1/horton-1.2.1.tar.gz
    tar -xvzf horton-1.2.1.tar.gz
    cd horton-1.2.1

Note: If you don't have wget installed on your computer, you can easily get it from MacPorts.

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

Note: If you don't have git installed on your computer, you can easily get it from MacPorts.

Common dependencies
===================

In order to compile and test Horton, one needs to
install relatively recent versions of the following programs/libraries:

* GCC, G++ and GFortran >= 4.5: http://gcc.gnu.org/ (suggested port: ``gcc47``, https://www.macports.org/ports.php?by=name&substr=gcc)
* Python >= 2.7, < 3.0: http://www.python.org/ (suggested port: ``python27``, https://www.macports.org/ports.php?by=name&substr=python27)
* Atlas  https://www.macports.org/ports.php?by=library&substr=atlas (suggested port: ``atlas``, https://trac.macports.org/browser/trunk/dports/math/atlas/Portfile)
* Numpy > 1.0: http://www.scipy.org/ (install numpy with atlas support:``sodp port install py27-numpy+atlas``, https://www.macports.org/ports.php?by=name&substr=numpy)
* h5py >= 2.2.1: http://www.h5py.org/ (suggested port: ``py27-h5py``, https://www.macports.org/ports.php?by=name&substr=h5py)
* Scipy >= 0.10.0: http://www.scipy.org/ (suggested port: ``py27-scipy``, https://trac.macports.org/browser/trunk/dports/python/py-scipy/Portfile)
* Cython >= 0.17.1 : http://www.cython.org/ (suggested port: ``py27-cython``, https://trac.macports.org/browser/trunk/dports/python/py-cython/Portfile) 
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/ (suggested port: ``py27-nose``, https://www.macports.org/ports.php?by=library&substr=py27-nose)
* Sympy >= 0.7.1: http://code.google.com/p/sympy/ (suggested port: ``py27-sympy``, https://trac.macports.org/browser/trunk/dports/python/py-sympy/Portfile)
* The ``patch`` program >= 2.0: http://savannah.gnu.org/projects/patch/ (or any of its equivalents)

Optionally, one may also install Matplotlib:

* Matplotlib >= 1.0: http://matplotlib.org/

If one is interested in generating the documentation from source, the following
packages are also needed:

* Sphinx > 1.0: http://sphinx.pocoo.org/
* Doxygen >= 1.8.6: http://www.doxygen.org/
* Breathe >= 1.2.0: http://breathe.readthedocs.org/en/latest/
* Docutils >= 0.11: http://docutils.sourceforge.net/

Most of these programs/libraries can be installed through MacPorts. 
First check this possibility before manually installing
the dependencies.

Specific dependencies
=====================

Before you start your Horton installation, please make sure that your ``setup.py`` has a correct path 
to your blas library. Most probably it should look like::

    blas_include_dirs = ["/opt/local/include"]                                                                
    blas_lib_dirs = ['/opt/local/include/atlas']  

The directory ``depends`` of the Horton source tree is used to build specific
dependencies from source. For the moment, there are two such dependencies,
namely `libint2 <http://sourceforge.net/p/libint/>`_ and `libxc
<http://www.tddft.org/programs/octopus/wiki/index.php/Libxc>`_
[marques2012]_. The directory ``depends``
contains a ``Makefile`` that takes care of downloading the right version and
compiling it. The following should get you the proper versions of libint and
libxc::

    cd depends
    make libint
    make libxc
    cd ..

The compilation of libint takes a few minutes. These commands will build
libraries suitable for static linking.

.. note::

    Alternatively, it is also possible to link libint and libxc dynamically. This
    requires some familiarity with software compilation on Unix systems. Make
    sure you have the following versions installed:

    * libint (for mpqc) >= 2.0.3-stable
    * libxc >= 2.1.0


Reference atoms
===============

This step can be skipped when compiling a stable release because each stable
release already contains reference atoms.

Several parts of Horton make use of reference atomic computations. These files
are too large to be included in the git revision system. Therefore, they must be
downloaded separately when compiling a development version of Horton::

    cd data/refatoms
    make all
    cd ../..


Compilation and installation
============================

The regular build and install is as done follows::

    ./setup.py install --user

The ``horton-*.py`` scripts are installed in ``~/.local/bin`` and you have to
add this directory to your ``PATH`` environment variable to make them accessible
from any directory.

.. note::

    When libint and libxc are compiled for static linking (as explained above),
    these libraries are found automatically. In case of dynamic linking,
    it may be necessary to specify explicitly the location of the shared objects
    and the header files with the options ``-I`` and ``-L`` of the setup script.

The documentation is compiled and viewed as follows::

    cd doc
    make html
    open _build/html/index.html
    cd ..

Environmental variables
===============================

We need to set up ``PYTHONPATH``::
    
    export PYTHONPATH="$PYTHONPATH:/Users/your-username/path-to-horton-installation/"

and ``HORTONDATA``::
   
    export HORTONDATA=/Users/your-username/path-to-horton-installation/data/

Note: If you are using Breathe for documentation purposes you also need to add the path-to-Breathe to your ``PYTHONPATH``

Running the tests
=================

Move to a directory outside the source tree and call nosetests as follows::

    cd ~
    nosetests-2.7  -v horton

In case one is testing horton on a system without an X Server, one has to
configure matplotlib to use a backend that does not rely on an X Server. This
can be done by adding a line ``backend: agg`` to the file
``~/.matplotlib/matplotlibrc``.
