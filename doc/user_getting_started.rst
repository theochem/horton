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

Getting started
###############

How to run HORTON
=================

HORTON is essentially a Python library that can be used to do electronic structure
calculations and the post-processing of such calculations. There are two different
ways of working with HORTON. The most versatile approach is to write Python scripts
that use HORTON as a Python library. This gives you full access to all features in
HORTON but it also requires some understanding of the Python programming language.
Alternatively, some of the functionality is also accessible through built-in
python scripts whose behaviour can be controlled through command line arguments.
This requires less technical background from the user.


Running HORTON as a Python library
----------------------------------

There will be many examples in the following sections that explain what can be
done with HORTON by writing your own Python scripts. These scripts should all
start with the following lines:

.. code-block:: python

    #!/usr/bin/env python

    # Import the HORTON library
    from horton import *

    # Optionally import some other stuff
    import numpy as np, h5py as h5, matplotlib.pyplot as pt

    # Actual script

This header is then followed by some Python code that does the actual computation
of interest. In such a script, you basically write the `main` program of the
calculation of interest, using the components that HORTON offers. The Horton library
is designed such that all features are as modular as possible, allowing you to
combine them in various ways.

Before you run your script, say ``run.py``, we recommend that you make it executable
(once):

.. code-block:: bash

    $ chmod +x run,py

When your script it completed, run it as follows:

.. code-block:: bash

    $ ./run.py

Do not use ``horton.py`` as a script name because will cause trouble when loading
the ``horton`` library (due to a namespace collision).


Running HORTON with provided ``horton-*.py`` scripts
----------------------------------------------------

The builtin HORTON scripts all have the following filename pattern: ``horton-*.py``.
Through command line arguments, one can control the actual calculation. Basic
usage information is obtained with the ``--help`` flag, e.g.:

.. code-block:: bash

    $ horton-convert.py --help


Writing your first HORTON Python script
=======================================

HORTON scripts just run with a regular Python interpreter
(like `ASE <https://wiki.fysik.dtu.dk/ase/>`_ and unlike
`PSI4 <http://www.psicode.org/>`_, which uses a modified Python interpreter).
This means you first need to learn basic Python. It may also be of interest to become
familiar with the Python packages that are popular in scientific computing. Here are
some links to broaden your background knowledge:

* `Python <https://www.python.org/>`_, `Python documentation <https://docs.python.org/2.7>`_, `Getting started with Python <https://www.python.org/about/gettingstarted/>`_
* `NumPy <http://www.numpy.org/>`_ (array manipulation), `Getting started with NumPy <http://docs.scipy.org/doc/numpy/user/>`_
* `SciPy <http://www.scipy.org/>`_ (scientific computing library), `The SciPy tutorial <http://docs.scipy.org/doc/scipy/reference/tutorial/index.html>`_
* `Matplotlib <http://matplotlib.org/>`_ (plotting)
* `H5Py <http://www.h5py.org/>`_ (load/dump arrays from/to binary files)
* `Programming Q&A <http://stackoverflow.com/>`_
* `Python code snippets <http://code.activestate.com/recipes/langs/python/>`_

The following sections go through some basic features that will appear in many
other examples in the documentation.


Atomic units
------------

HORTON internally works exclusively in atomic units. If you want to convert a
value from a different unit to atomic units, multiply it with the appropriate unit
constant, e.g. the following sets ``length`` to 5 Angstrom and prints it out in
atomic units:

.. code-block:: python

    length = 5*angstrom
    print length

Conversely, if you want to print a value in a different unit than atomic units,
divide it by the appropriate constant. For example, the following prints an
energy of 0.001 Hartree in kJ/mol:

.. code-block:: python

    energy = 0.001
    print energy/kjmol

An overview of all units can be found here: :py:mod:`horton.units`.

There are two special cases:

1. Angles are in radians but you can use the ``deg`` unit the work with degrees,
   e.g. ``90*deg`` and ``np.pi/2`` are equivalent.
2. Temperatures are in Kelvin. (There is no atomic unit for temperature.)


Array indexing
--------------

All arrays (and list-like) objects in Python use zero-based indexing. This means
that the first element of a vector is accessed as follows:

.. code-block:: python

    vector = np.array([1, 2, 3])
    print vector[0]

This convention also applies to all array-like objects in HORTON, e.g. the
first orbital in a Slater determinant has index 0.


Loading/Dumping data from/to a file
-----------------------------------

All input and output of data in HORTON is managed through the :py:class:`~horton.io.iodata.IOData`
class. Data can be loaded with a call to the class method :py:meth:`~horton.io.iodata.IOData.from_file`
of the :py:class:`~horton.io.iodata.IOData` class, e.g.:

.. code-block:: python

    mol = IOData.from_file('water.xyz')

The data read from file are accessible as attributes of the ``mol`` object. For
example, the following prints the coordinates of the nuclei in Angstrom:

.. code-block:: python

    print mol.coordinates/angstrom

Writing data to a file is done by first creating a :py:class:`~horton.io.iodata.IOData`
instance, followed by setting the appropriate attributes and finally followed by
a call to the :py:meth:`~horton.io.iodata.IOData.to_file` method of the
:py:class:`~horton.io.iodata.IOData` class. For example, the following creates a ``.xyz``
file with a Neon atom:

.. code-block:: python

    mol = IOData(title='Neon')
    mol.coordinates = np.array([[0.0, 0.0, 0.0]])
    mol.numbers = np.array([10])
    mol.to_file('neon.xyz')

A Complete list of supported file formats for data input or output can be found
here: :ref:`ref_file_formats`, which includes a list of :py:class:`~horton.io.iodata.IOData`
attributes supported by each file format. A definition of all possible
:py:class:`~horton.io.iodata.IOData` attributes can be found here:
:py:class:`horton.io.iodata.IOData`.


The periodic table
------------------

HORTON has a periodic table of elements with several atomic properties that may
be of interest for computations. Details can be found here: :py:mod:`horton.periodic`.
The following example prints some information for Carbon:

.. code-block:: python

    print periodic[6].mass/amu                # The mass in atomic mass units
    print periodic['C'].cov_radius/angstrom   # The covalent radius in angstrom
    print periodic['c '].c6                   # The C6 coefficient in atomic units

Note that you can be relatively sloppy with the index to refer to elements of the periodic table.


A complete example
------------------

This first example is kept very simple, just to illustrate the basics of a HORTON
Python script. (It does not do an electronic calculation yet.) The following example
loads an ``.xyz`` file and computes the molecular mass. Finally, it writes out
data read from the ``.xyz`` file and the mass into a ``.h5`` file, using HORTON's
internal data format.

.. literalinclude:: ../data/examples/getting_started/first.py
    :caption: data/examples/getting_started/first.py

Note that the part ``context.get_fn('test/water.xyz')`` is used to look up a data
file from the HORTON data directory. If you want to use your own file, load the
molecule as follows instead:

.. code-block:: python

    mol = IOData.from_file('your_file.xyz')
