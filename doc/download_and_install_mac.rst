Download and installation on Mac OS X (10.8 - 10.10)
####################################################

Disclaimer
==========

Horton has been tested on Mac OS X 10.8 - 10.10 using MacPorts. If you
run any other version of OS X, some of the instructions below may not
work.


MacPorts
=========

We strongly recommend to install all the packages required by Horton
through MacPorts. We also advise to uninstall/remove all python packages
installed outside MacPorts (e.g., Canopy).

The latest version of MacPorts can be downloaded from the web:
https://www.macports.org/install.php. This guide was tested with MacPorts 2.3.3 but
it should also work with newer version.


Some hints for installing MacPorts
----------------------------------

It is recommended to place the MacPort ``bin`` and ``sbin`` directories in the beginning
of the ``PATH`` variable, such that they take precedence over other directories in the
``PATH``. This reduces the risk of conflicts with other software. For example, put the
following in your ``~/.bash_profile``::

    MACPORTSHOME=/opt/local
    export PATH=${MACPORTSHOME}/bin:${MACPORTSHOME}/sbin:${PATH}
    export CPATH=${MACPORTSHOME}/include:${CPATH}

The ``CPATH`` variable is required for the compilation of Horton. It tells the compiler
where to look for C/C++ header files (besides the standard locations).

It is possible to install MacPorts in user space. This offers some advantages, e.g. no
need for ``sudo`` and root permissions, easier to clean up in case of epic failures, etc.
Only do this if you feel comfortable with non-standard installs. Otherwise, just use the
default installation of MacPorts.

To do a user-space install, download the MacPorts source code, unpack it, and install
as follows::

    export MACPORTSHOME=${HOME}/macports
    ./configure --prefix=${MACPORTSHOME} --with-applications-dir=${MACPORTSHOME}/Applications --with-install-user=${USER} --with-install-group=$(id -gn) --with-no-root-privileges
    make install

Finally, set the ``MACPORTSHOME`` directory in the ``.bash_profile`` and open a fresh
terminal such that your new ``.bash_profile`` is in effect.


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
    curl -O https://github.com/theochem/horton/releases/download/1.2.1/horton-1.2.1.tar.gz
    tar -xvzf horton-1.2.1.tar.gz
    cd horton-1.2.1


Latest development code (experts only)
--------------------------------------

In order to get the latest development version of the source code, and to upload
your own changes, you need to work with git. Git is a version control system
that makes life easy when a group of people are working on a common source code.
All information about git (including downloads and tutorials) can be found here:
http://git-scm.com/. The official public git URL of Horton is:
``git://github.com/theochem/horton.git.`` If you don't have git installed on
your computer, you can easily get it from MacPorts::

    port install git

In order to `clone` the public Horton repository, run this command::

    git clone git://github.com/theochem/horton.git
    cd horton

The version history can be updated with the latest patches with the following
command::

    git pull

There is also a web interface to Horton's git repository:
https://github.com/theochem/horton


Dependencies for building, installing and testing Horton
========================================================

In order to compile and test Horton, one needs to
install relatively recent versions of the following programs/libraries:

* GCC, G++ and GFortran >= 4.5: http://gcc.gnu.org/
* Python >= 2.7, < 3.0: http://www.python.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* Atlas >= 3.10.1: http://math-atlas.sourceforge.net/ (or any other BLAS implementation that you like more)
* Numpy > 1.0: http://www.scipy.org/
* Scipy >= 0.10.0: http://www.scipy.org/
* Cython >= 0.17.1 : http://www.cython.org/
* h5py >= 2.2.1: http://www.h5py.org/
* Sympy >= 0.7.1: http://code.google.com/p/sympy/
* Matplotlib >= 1.0: http://matplotlib.org/
* LibXC >= 2.1.0: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc
* LibInt2 >= 2.0.3: http://sourceforge.net/p/libint/home


Installing the dependencies with MacPorts
-----------------------------------------

All dependencies can be installed with MacPorts, except for LibInt2. We recommend the
following ports:

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

These are installed with the following commands. (When MacPorts is installed in user
space, the ``sudo`` can be omitted.)::

    sudo port install gcc49
    sudo port select --set gcc mp-gcc49
    sudo port install python27
    sudo port select --set python python27
    sudo port install py27-nose
    sudo port select --set nosetests nosetests27
    sudo port install atlas
    sudo port install py27-numpy +atlas
    sudo port install py27-scipy +atlas
    sudo port install cython
    sudo port select --set cython cython27
    sudo port install py27-h5py
    sudo port install py27-sympy
    sudo port select --set py-sympy py27-sympy
    sudo port install py27-matplotlib
    sudo port install libxc

LibInt2 cannot be installed with MacPorts yet and must be installed manually, as
explained in the next section. The GNU compilers are in fact only used to compile
Fortran code as the default C/C++ compiler on the Mac is ``clang``.


Installing dependencies manually
--------------------------------

**BLAS**

In principle, any BLAS implementation may be used. In case of a custom build,
some environment variables must be set prior to building Horton, as discussed
below. Also keep in mind that MacPorts only supports Atlas for building NumPy and SciPy.


**LibXC**

The directory ``depends`` of the Horton source tree contains a make file that
will download and LibXC, which will work on most systems::

    (cd depends; make libxc)

This results in a libxc library suitable for static linking. If this fails,
consult your local Mac guru to build LibXC. For more info about LibXC, check
the website: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc

**LibInt2**

The directory ``depends`` of the Horton source tree contains a make file that
will download and LibInt2, which will work on most systems::

    (cd depends; make libint -j4)

The compilation of libint takes a few minutes and results in a library for
static linking. If this fails, consult your local Mac guru to build LibInt2.
For more info about LibInt2, check the website:
http://sourceforge.net/p/libint/home


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

Build and install
-----------------

The regular build and install is done as follows::

    ./setup.py install --user

The ``setup.py`` script does a reasonable attempt to configure the compiler and
linker settings for the LibXC, LibInt2 and BLAS libraries. However, this does
not work in all environments. In case of a faillure, or if another configuration
than the default is desired, read the following section.


Overriding default compiler/linker settings for LibXC, LibInt2 and BLAS
-----------------------------------------------------------------------

The manual configuration of the compiler and linker settings is described here:
:ref:`setup_cfg`. Only read this section if the default build and install did
not work.


Runtime Environmental variables
-------------------------------

We need to set the following variables in ``~/.bash_profile to use Horton::

    export PATH=${HOME}/Library/Python/2.7/bin:${PATH}
    # I did not have to set the following two.
    # The --user option of the setup.py script normally installs stuff in a place
    # where Python will find it without setting environment variables. ~Toon
    export PYTHONPATH=${PYTHONPATH}:${HOME}/path-to-horton-installation/
    export HORTONDATA=${HOME}/path-to-horton-installation/data/


Running the tests
=================

Change to a directory outside the source tree and call nosetests as follows::

    (cd ~; nosetests -v horton)

In case one is testing horton on a system without an X Server, one has to
configure matplotlib to use a backend that does not rely on an X Server. This
can be done by adding a line ``backend: agg`` to the ``matplotlibrc`` file.
This file is located in ``~/.matplotlib`` or ``~/.config/matplotlib``.


Building the documentation
==========================

Dependencies
------------

If one is interested in generating the documentation from source, the following
packages are also needed:

* Sphinx > 1.0: http://sphinx.pocoo.org/
* Doxygen >= 1.8.6: http://www.doxygen.org/
* Breathe >= 1.2.0: http://breathe.readthedocs.org/en/latest/
* Docutils >= 0.11: http://docutils.sourceforge.net/
* A latex distribution (Texlive)
* DVIpng >= 1.14: http://savannah.nongnu.org/projects/dvipng/
* The Preview style for Latex (preview.sty)


Installing the dependencies with MacPorts and PIP
-------------------------------------------------

Most can be installed directly with MacPorts. The following list of ports is recommended:

* ``py27-sphinx`` (also install Docutils): https://trac.macports.org/browser/trunk/dports/python/py-sphinx/Portfile
* ``doxygen``: https://trac.macports.org/browser/trunk/dports/textproc/doxygen/Portfile
* ``dving`` (installs texlive as dependency): https://trac.macports.org/browser/trunk/dports/tex/dvipng/Portfile
* ``texlive-latex-extra`` (contains ``preview.sty``): https://trac.macports.org/browser/trunk/dports/tex/texlive-latex-extra/Portfile
* ``py27-pip`` (neede for Breathe): https://trac.macports.org/browser/trunk/dports/python/py-pip/Portfile

For Breathe, one can use the PIP installer. The following commands will install everything
as suggested::

    sudo port install py27-sphinx
    sudo port select --set sphinx py27-sphinx
    sudo port install doxygen
    sudo port install dviping
    sudo port install texlive-latex-extra
    sudo port install py27-pip
    sudo port select --set pip pip27
    pip install breathe --user

One must also build LibXC statically in the ``depends`` directory, as explained
above, to generate the list of DFT functionals in the documentation.


Actual build
------------

The documentation is compiled and viewed as follows::

    cd doc
    make html
    open _build/html/index.html
    cd ..
