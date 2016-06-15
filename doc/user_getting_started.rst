..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2016 The HORTON Development Team
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

.. _user_getting_started:

Getting Started
###############

How to use HORTON?
==================

HORTON is essentially a Python library that can be used for performing electronic structure
calculations as well as interpreting these calculations (i.e. post-processing). There are two different
ways to use HORTON. The most versatile approach is to write Python scripts
that use HORTON as a Python library. This gives you full access to all the features available in
HORTON; however, this requires some knowledge of the Python programming language.
Alternatively, part of HORTON's functionality is accessible through built-in
python scripts whose performance can be controlled through command line arguments.
Obviously, this requires less programming knowledge.


Running HORTON as a Python library
----------------------------------

There will be many examples in the following sections demonstrating how
HORTON can be used when writing your own Python scripts. These scripts should all
start with the following lines:

.. code-block:: python

    #!/usr/bin/env python

    # Import the HORTON library
    from horton import *

    # Import some other stuff (optional)
    import numpy as np, h5py as h5, matplotlib.pyplot as pt

    # Actual Python script

This header is then followed by some Python code that does the actual computation
of interest. In such a script, you basically write the `main` program to do your
calculations, exploiting the components that HORTON offers. The HORTON library
is designed such that all its features are as modular as possible, allowing you to
combine them in various ways.

Before running your script, say ``run.py``, we recommend that you to make it
executable (this needs to be done only once for every script):

.. code-block:: bash

    chmod +x run.py

Now, when your script has completed, you can run it as follows:

.. code-block:: bash

    ./run.py

Do not use ``horton.py`` as your script name; this will cause trouble when loading
the ``horton`` library (due to a namespace collision).

.. _using_horton_as_a_script:

Running HORTON as ``horton-*.py`` scripts
-----------------------------------------

The built-in HORTON scripts all have the ``horton-*.py`` filename pattern.
Through command line arguments, one can control the actual calculations performed by these scripts. Basic
information on how to use each built-in script can be obtained by using the ``--help`` flag. For example,

.. code-block:: bash

    horton-convert.py --help


Writing a basic HORTON Python script
====================================

HORTON scripts just run with a regular Python interpreter
(like `ASE <https://wiki.fysik.dtu.dk/ase/>`_ and unlike
`PSI4 <http://www.psicode.org/>`_, which uses a modified Python interpreter).
This means that you need to have a basic knowledge of Python. In addition, it will be helpful to be
familiar with popular Python packages for scientific computing. The links below provide
some resources to broaden your Python knowledge:

* `Python <https://www.python.org/>`_, `Python documentation <https://docs.python.org/2.7>`_, `Getting started with Python <https://www.python.org/about/gettingstarted/>`_
* `NumPy <http://www.numpy.org/>`_ (array manipulation), `Getting started with NumPy <http://docs.scipy.org/doc/numpy/user/>`_
* `SciPy <http://www.scipy.org/>`_ (scientific computing library), `The SciPy tutorial <http://docs.scipy.org/doc/scipy/reference/tutorial/index.html>`_
* `Matplotlib <http://matplotlib.org/>`_ (plotting)
* `H5Py <http://www.h5py.org/>`_ (load/dump arrays from/to binary files)
* `Programming Q&A <http://stackoverflow.com/>`_
* `Python code snippets <http://code.activestate.com/recipes/langs/python/>`_

The following sections go through some basic features that will appear in many
other examples in the documentation.


Atomic Units
------------

Internally, HORTON works exclusively in atomic units. If you want to convert a
value from a different unit to atomic units, multiply it with the appropriate unit
constant, e.g. the following snippet sets ``length`` to 5 Angstrom and prints it in
atomic units:

.. code-block:: python

    length = 5*angstrom    # recording 5 angstrom in atomic units
    print length

Conversely, if you want to print a value in a different unit than atomic units,
divide it by the appropriate constant. For example, the following prints an
energy of 0.001 Hartree in kJ/mol:

.. code-block:: python

    energy = 0.001            # recording energy in Hartree
    print energy/kjmol        # printing energy in kJ/mol

An overview of all units can be found in :py:mod:`horton.units`.

There are two special cases:

1. Angles are in radians, but you can use the ``deg`` unit to work with degrees,
   for example, ``90*deg`` and ``np.pi/2`` are equivalent.
2. Temperatures are in Kelvin.


Array Indexing
--------------

All arrays and list-like objects in Python use zero-based indexing. This means
that the first element of a vector is accessed as follows:

.. code-block:: python

    vector = np.array([1, 2, 3])
    print vector[0]

This convention also applies to all array-like (and list-like) objects in HORTON, e.g. the
first orbital in a Slater determinant has index 0.


Loading/Dumping Data from/to a File
-----------------------------------

All input and output of data in HORTON is managed through the :py:class:`~horton.io.iodata.IOData`
class. To load data, call the :py:meth:`~horton.io.iodata.IOData.from_file` method
of the :py:class:`~horton.io.iodata.IOData` class, e.g.:

.. code-block:: python

    mol = IOData.from_file('water.xyz')

The information read from the file is accessible through attributes of the ``mol`` object. For
example, the following prints the coordinates of the nuclei in Angstrom:

.. code-block:: python

    print mol.coordinates/angstrom

To write data into a file, first create an instance of the :py:class:`~horton.io.iodata.IOData`
class, then set the appropriate attributes, and write the content into a file by
calling the :py:meth:`~horton.io.iodata.IOData.to_file` method of the
:py:class:`~horton.io.iodata.IOData` class. For example, the following snippet creates a ``.xyz``
file for a Neon atom:

.. code-block:: python

    mol = IOData(title='Neon')
    mol.coordinates = np.array([[0.0, 0.0, 0.0]])
    mol.numbers = np.array([10])
    mol.to_file('neon.xyz')

For a complete list of supported input/output file formats please refer to
:ref:`ref_file_formats`; this includes a list of :py:class:`~horton.io.iodata.IOData`
attributes supported by each file format. A definition of all possible
:py:class:`~horton.io.iodata.IOData` attributes can be found in
:py:class:`horton.io.iodata.IOData`.


Periodic Table
--------------

HORTON has a periodic table of elements alongside several atomic properties that may
come in handy in computations. For more details please refer to :py:mod:`horton.periodic`.
The following example prints some information for Carbon atom:

.. code-block:: python

    print periodic[6].mass/amu                # The mass in atomic mass units
    print periodic['C'].cov_radius/angstrom   # The covalent radius in angstrom
    print periodic['c '].c6                   # The C6 coefficient in atomic units

As demonstrated above, you can be relatively sloppy with the index when referring to elements of the periodic table.


A Complete Example
------------------

This first example is kept very simple in order to illustrate the basics of a HORTON
Python script. (It neither performs an electronic calculation nor does post-processing.) This example
loads a ``.xyz`` file and computes the molecular mass. Finally, it writes the
data read from the ``.xyz`` file and the calculated mass into a ``.h5`` file, using HORTON's
internal data format.

.. literalinclude:: ../data/examples/getting_started/first.py
    :caption: data/examples/getting_started/first.py

Note that the ``context.get_fn('test/water.xyz')`` expression is used to look up a data
file from the HORTON data directory. If you want to use your own file, load the
molecule as follows:

.. code-block:: python

    mol = IOData.from_file('your_file.xyz')
