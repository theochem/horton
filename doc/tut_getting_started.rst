Getting started
###############

Downloading the code
====================

In order to get the latest version of the source code, and to upload your own
changes, you need to work with git. Git is a version control system that
makes life easy when a group of people are working on a common source code. It
also makes sense to use it for personal projects. All information about git
(including downloads and tutorials) can be found here: http://git-scm.com/. The
official git URL of Horton is:

    git://github.com/theochem/horton.git

In order to `clone` the official Horton repository, run this command::

    git clone git://github.com/theochem/horton.git

The version history can be updated with the latest patches with the following
command::

    git pull

There is also a web interface to Horton's git repository:
https://github.com/theochem/horton

Compilation
===========

In order to `compile` the Horton C-extension and documentation, one needs to
install relatively recent versions of the following programs/libraries:

* Python < 3.0: http://www.python.org/ (also install `development files`)
* Numpy > 1.0: http://www.scipy.org/ (scipy is needed as well)
* Cython > (todo): http://www.cython.org/
* Sphinx > 1.0: http://sphinx.pocoo.org/

On a decent operating system, these programs/libraries can be easily installed
with a package manager. First check that option before manually installing this
software. On Ubuntu, this one-liner will take care of it::

    sudo apt-get install python-dev python-numpy cython python-sphinx

One may either do a regular installation in the home directory, or an in-pace
build in the source tree.

The **regular build/install** is as done follows::

    ./setup.py install --home=~

This will put a `compiled` version of Horton in the ``build`` directory. Then a
copy of this directory is placed under ``~/lib/python``. Configure your
environment variables (e.g. modify .bashrc) such that
``$PYTHONPATH=$HOME/lib/python``.

The **in-place build** is useful for testing purposes, and is done as follows::

    ./setup.py build_ext -i

The documentation is compiled and viewed as follows::

    cd doc
    make html
    make pdf
    firefox _build/html/index.html


Testing
=======

A bunch of validation routines are included in Horton. To run these tests, one
must install the Nosetests testing framework. It can be found here:

http://somethingaboutorange.com/mrl/projects/nose/0.11.2/

We recommend that you use a package manager to install Nosetests. On Ubuntu,
this software can be installed as follows::

    sudo apt-get install python-nose

Once this Python package is installed, perform an **in-place build** and run
nosetests afterwards::

    ./setup.py build_ext -i
    nosetests -v

If all runs well, the screen output should end with 'OK'.


Basic example
=============

(For now, this is just an idea. It does not work yet.)

This is a basic example computation in Horton. The input file is just
a small Python main program that uses the Horton library. The following script
performs a HF/STO-3G computation on HF::

    from horton import *
    import numpy as np

    system = System(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])*angstrom, [1, 9], 'STO-3G')
    hamiltonian = Hamiltonian(system, [KineticEnergy(),  Hartree(), Fock()])
    wfn = minimize(hamiltonian)
    print wfn.energy/kjmol


In many cases, the molecule of interest is stored in an external file, e.g.
an XYZ file. To avoid retyping this molecule into the script, on may work as
follows::

    from horton import *

    system = System.from_file('hcl.xyz', 'STO-3G')
    hamiltonian = Hamiltonian(system, [KineticEnergy(),  Hartree(), Fock()])
    wfn = minimize(hamiltonian)
    print wfn.energy/kjmol

The kinetic energy may be omitted. If not present, it will be added
automatically. There is also a shortcut to combine the Hartree and the Fock
potential::

    from horton import *

    system = System.from_file('hcl.xyz', 'STO-3G')
    hamiltonian = Hamiltonian(system, HartreeFock())
    wfn = minimize(hamiltonian)
    print wfn.energy/kjmol
