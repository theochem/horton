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

    git://todo

In order to `clone` the official Horton repository, run this command::

    git clone git://todo

The version history can be updated with the latest patches with the following
command::

    git pull


Compilation
===========

In order to `compile` the Horton C-extension and documentation, one needs to
install relatively recent versions of the following programs/libraries:

* Python < 3.0: http://www.python.org/ (also install `development files`)
* Numpy > 1.0: http://www.scipy.org/
* Cython > (todo): http://www.cython.org/
* Sphinx > 1.0: http://sphinx.pocoo.org/

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
must install the nosetests python package. It can be found here:

http://somethingaboutorange.com/mrl/projects/nose/0.11.2/

Once this python package is installed, perform an **in-place build** and run
the following in the root of the source tree::

    nosetests -v

If all runs well, the screen output should end with 'OK'.
