Download and installation on Linux (Fedora and Ubuntu)
######################################################


Disclaimer
==========

Horton has been tested on Fedora and Ubuntu Linux. If you run any
other operating system, some of the instructions below may not work.


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
http://git-scm.com/. The official git URL of Horton is:
``git://github.com/theochem/horton.git.`` On most Lunx distributions, git can be
installed with the package manager:

**Ubuntu Linux**::

    sudo apt-get install git

**Fedora Linux**::

    sudo yum install git

In order to `clone` the public Horton repository, run these commands::

    mkdir -p ~/build
    cd ~/build
    git clone git://github.com/theochem/horton.git
    cd horton

The version history can be updated with the latest patches with the following
command::

    git pull

There is also a web interface to Horton's git repository:
https://github.com/theochem/horton


Dependencies for building, installing and testing Horton
========================================================

In order to compile and test Horton, one needs to install relatively recent
versions of the following programs/libraries:

* Python >= 2.7, < 3.0: http://www.python.org/ (also install `development files`)
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* GCC, G++ and GFortran >= 4.5: http://gcc.gnu.org/ (In principle the Intel compilers or
  any other of your favorite compilers should work. The GNU compilers are used for
  development and testing.)
* Atlas >= 3.10.1: http://math-atlas.sourceforge.net/ (or any other BLAS implementation that you like more)
* Numpy >= 1.7.0: http://www.numpy.org/
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

With the popular Linux distributions, most of these dependencies can be
installed with a package manager:

**Ubuntu Linux** (does not have a libint2 package, see next section)::

    sudo apt-get install python-dev gcc g++ gfortran python-numpy python-h5py python-scipy cython python-nose python-sympy python-matplotlib libxc-dev libatlas-base-dev curl

**Fedora Linux**::

    sudo yum install python-devel gcc gcc-c++ gcc-gfortran numpy h5py scipy Cython python-sphinx python-nose sympy python-matplotlib libint2-devel libxc-devel atlas-devel curl


Installing dependencies manually
--------------------------------

If the package manager of your Linux distribution does not have the desired
packages (or the right versions), one has to install them manually, e.g.
download and execute an installer, or download and unpack a binary package. On
HPC environments a compilation from scratch is recommended.

**BLAS**

In principle, any BLAS implementation may be used. In case of a custom build,
some environment variables must be set prior to building Horton, as discussed
below.


**LibXC**

The directory ``depends`` of the Horton source tree contains a make file that
will download and LibXC, which will work on most systems::

    (cd depends; make libxc)

This results in a libxc library suitable for static linking. If this fails,
consult your local Linux guru to build LibXC. For more info about LibXC, check
the website: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc

**LibInt2**

The directory ``depends`` of the Horton source tree contains a make file that
will download and LibInt2, which will work on most systems::

    (cd depends; make libint -j4)

The compilation of libint takes a few minutes and results in a library for
static linking. If this fails, consult your local Linux guru to build LibInt2.
For more info about LibInt2, check the website:
http://sourceforge.net/p/libint/home

**Python dependencies**

In some cases, Pip, the Python package manager, may be a good choice to install
the most recent versions of the Python packages in the list of dependencies.
Assuming that you have installed some compilers, the Python development files
and HDF5 development files, the following command installs the remaining
dependencies in your home directory::

    pip install --user numpy scipy cython h5py sphinx nose sympy


Reference atoms
===============

This step can be skipped when compiling a stable release because each stable
release already contains reference atoms.

Several parts of Horton make use of reference atomic computations. These files
are too large to be included in the git revision system. Therefore, they must be
downloaded separately when compiling a development version of Horton::

    (cd data/refatoms; make all)


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


Runtime Configuration
---------------------

We need to set the following variable in ``~/.bashrc`` to use Horton::

    export PATH=${HOME}/.local/bin:${PATH}

    # If you used special link options for LibXC, LibInt2 or BLAS, something along
    # the following lines may also be needed:
    # export LD_LIBRARY_PATH=some_dir/with/shared_objects/${LD_LIBRARY_PATH}

If you run Horton on a headless node, i.e. without an X server, you need to
configure Matplotlib to use a backend that does not require a graphical user
interface. (See http://matplotlib.org/faq/usage_faq.html#what-is-a-backend for
more details on the Matplotlib backends.) This can be done by adding the
following line to your ``matplotlibrc`` file::

    backend: agg

This file is located either in ``${HOME}/.matplotlib`` or
``${HOME}/.config/matplotlib``.


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


Installing the dependencies with a package manager and PIP
----------------------------------------------------------

**Ubuntu Linux**::

    sudo apt-get install python-sphinx doxygen preview-latex-style python-docutils python-pip dvipng preview-latex-style

**Fedora Linux**::

    sudo yum install python-sphinx doxygen tex-preview python-docutils python-pip dvipng tex-preview

Since Breathe is relatively new, it must be installed manually. For example, it
is available through PyPI and can be installed as follows::

    pip install --user breathe

One must also build LibXC statically in the ``depends`` directory, as explained
above, to generate the list of DFT functionals in the documentation.



Actual build
------------

The documentation is compiled and viewed as follows::

    (cd doc; make html; firefox _build/html/index.html)
