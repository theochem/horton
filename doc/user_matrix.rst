How to use the matrix package
#############################

Introduction
============

Horton chooses to separate high-level quantum chemistry code from low-level
linear algebra code with an abstraction layer known as the Matrix module. The
reason for this is to keep the algorithm separate from the low-level
implementation details and also to share frequent operations between different
parts of code. This allows Horton to implement new linear algebra operations
(ie Cholesky decomposition of the 2-electron integrals) without rewriting any
high-level code (ie Geminals). Another reason for this module is to keep the
allocation of memory controlled. By default, Numpy will create temporary
objects when evaluating a statement like ``a += b``. This is clearly not
acceptable when ``a`` or ``b`` is a 2-electron integral, since the temporary
object will be hundreds of gigabytes in some cases. The downside to this
approach is a different programming style, a small performance penalty and the
need to implement new operations within the Matrix module.

See below for an image representation of this concept.

.. image:: matrix_concept.png

How to use this abstraction layer
=================================

The Matrix module is split into several types.

First, the module is split according to backend, (ie dense Numpy storage,
Cholesky decomposed Numpy storage, etc) and then according to the type of
object being stored (ie wavefunction expansion, 2-index tensor, 3-index
tensor, 4-index tensor, etc). This distinction can be seen in the files of
the ``horton/matrix`` directory. Each file is a different backend and within
each file, there are classes for different objects.

The instances of the matrix classes are created from a :class:`.LinalgFactory`
object. As a user, you will need to instantiate one of these Factories when
you start writing your program. This is because most of the other parts of
Horton will ask for a LinalgFactory to allocate memory.

You can create a LinalgFactory instance.

.. code-block:: python

    lf = DenseLinalgFactory()

Note that other choices of Factories are valid here too. Any child class of
the "Matrix" class is possible.

.. code-block:: python

    lf = CholeskyLinalgFactory()

You can then use the Factory as an argument to other parts of the code.

.. code-block:: python

    lf = DenseLinalgFactory()
    er = obasis.compute_electron_repulsion(lf)

Allocation of new objects can be done by calling one of the functions
starting with ``create_`` within the Factory.

.. code-block:: python

    A = lf.create_two_index() #matrix
    A4 = lf.create_four_index() #4_rank_tensor
    wfn = lf.create_expansion() #wavefunction expansion


Operations are a little bit unusual because we want to avoid allocating
temporary memory. Most of the operations are of the in-place variety.

.. code-block:: python

    #A = A + B
    A.iadd(B)

    #A = A * B
    A.idot(A)

Operations that are not in-place operations are usually not available because
they create temporary objects. Instead, they can be broken up into a series
of operations, starting with one that explicitly allocates memory.

.. code-block:: python

    #A = B + C (NOT POSSIBLE)
    #A = A + B
    #A = A + C
    A.iadd(B)
    A.iadd(C)

There are convenience functions to add/multiply a string of matrices together
as well. Contact the authors for more details.

Many operations have been implemented in the Matrix class. See
:mod:`.Matrix.dense` for details.
