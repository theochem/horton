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

How to use the matrix package?
##############################

Introduction
============

In Horton, we use our own Matrix package (abstraction layer) to implement the
linear algebra code (low-level) into the quantum chemistry code (high-level). There
are two main reasons for such implementation:

1. We can develop new linear algebra operations without breaking the higher-level
   code. For example, if we implemented Cholesky decomposition for the transformation
   of 2-electron integrals, we should be able to use this method in the Geminals
   code without rewriting the Geminals code.
2. Numpy, the conventional linear algbebra package in Python, may not manage
   memory effectively for large structures, such as the 2-electron
   integrals. For example, Numpy allocates temporary memory when evaluating ``a += b``, which
   is costly when  ``a`` or ``b`` is very large. We developed a Matrix package
   (using the Numpy package) geared towards memory management of quantum
   chemistry calculations, and we use that package (Matrix package) instead of Numpy
   directly.

Though such an abstraction layer seems pedantic and ostentatious, and requires
a tedious implementation of all new operations into the Matrix package, the
small penalty in performance and the ease of implementing different concepts
into the low-level areas of Horton make it well worth the effort.

Below conceptualizes the organization of Horton:

.. image:: matrix_concept.png

How to use this abstraction layer
=================================

The Matrix package is organized (in ``horton/matrix``) as follows:

At the top level, the package is split into modules in which the
implementation of the data storage and manipulation (backend). So far, there are
two such modules, a dense Numpy storage and Cholesky decomposition of the
Numpy storage. Because these modules treat data in fundamentally different
ways, each method used by the higher-level module (listed in the `base.py`)
must be implemented.

Then, the code for each module is organized by the type of data that is
manipulated within a class. So far, classes for 1-index tensor, 2-index tensor,
3-index tensor, 4-index tensor, and wavefunction expansion have been
implemented.

To avoid reallocation of memory, we use a :class:`.LinalgFactory` instance to
create these objects (e.g. 2-index tensor). This instance will allocate memory
for these objects, and then, the operations performed on these objects will modify
their own attributes directly. Most of the operations are in-place, i.e.
modifies their own data based on input and returns no output.

For example, to create and modify a dense two index tensor, we first create a
LinalgFactory instance and call it ``lf``:

.. code-block:: python

    lf = DenseLinalgFactory()

We use ``lf`` to allocate a two index tensor object:

.. code-block:: python

    A = lf.create_two_index() #matrix

We can also pass ``lf`` as an argument to other parts of the code, which
uses it to allocate an object.

.. code-block:: python

    er = obasis.compute_electron_repulsion(lf)

We can modify the two-tensor object by some in-place operations:

.. code-block:: python

    #A = A + B
    A.iadd(B)

    #A = A * B
    A.idot(A)

    #A = B + C (NOT POSSIBLE)
    #A = A + B
    #A = A + C
    A.iadd(B)
    A.iadd(C)

Note that more complex operations, such as (`A = B + C`), can be broken up into
a series of in-place operations, that deal with explicitly allocated memory.

We can appreciate the simplicity of implementing different modules by playing
with the different modules available (two). For example, we could have used

.. code-block:: python

    lf = CholeskyLinalgFactory()

in place of the ``DenseLinalgFactory`` above. Making this change will not change
any of the preceeding code, provided that the same objects are implemented into
this module as well.

We can also allocate different objects, if implemented, using ``lf``:

.. code-block:: python

    A4 = lf.create_four_index() #4_rank_tensor
    wfn = lf.create_expansion() #wavefunction expansion

Many functions and classes have been implemented into the Matrix class. It may
help to read over some of the (hopefully) documented module files to see if a
desired function has already been implemented. In the event that a
desired function has not been implemented, please contact the authors to
make a feature request or for more details on implementing it yourself.
