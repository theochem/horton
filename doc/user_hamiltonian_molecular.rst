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

Molecular Hamiltonians
######################

A molecular Hamiltonian is typically set up in three steps. First, you load or
construct a molecular geometry. Then, you generate a Gaussian basis set for
this molecule. Finally, you can calculate all kinds of matrix elements with that
basis to define a Hamiltonian.


.. _setup-molgeometry:

Specifying the molecular geometry
=================================

HORTON can load/dump the molecular geometry from/to different file formats. The
file format is determined automatically using the extension or prefix of the
filename. For more information on supported file formats, refer to
:ref:`ref_file_formats`.


.. _read-molgeometry:

Reading the molecular geometry from file
----------------------------------------

The molecular geometry can be read from file using the method
:py:meth:`~horton.io.iodata.IOData.from_file` of the ``IOData`` class,

.. code-block:: python

    mol = IOData.from_file(filename)


Constructing a molecular geometry from scratch
----------------------------------------------

You can also generate molecular geometries with Python code. The
following example constructs a ring of Hydrogen atoms and writes it to an
XYZ file

.. literalinclude:: ../data/examples/hamiltonian/hydrogen_ring.py
    :lines: 2-
    :caption: data/examples/hamiltonian/hydrogen_ring.py


.. _user_molecularham_basis:

Specifying the basis set
========================

HORTON supports basis sets consisting of generally contracted Cartesian Gaussian
functions up to a maximum angular momentum of seven (I functions).
HORTON is using the same basis set format as NWChem, and the basis sets can be
downloaded from the EMSL webpage (https://bse.pnl.gov/bse/portal).

HORTON is distributed with most of the popular basis sets. A list of currently
supported built-in basis sets can be found here:
:ref:`ref_gaussian_basis_standard_sets`. The basis set for a given molecule is
constructed with the function :py:func:`~horton.gbasis.gobasis.get_gobasis`


Unique basis set for all atoms
------------------------------

Usually, you want to define a unique (the same) basis set for the whole system.
This can be done by a function call.

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')

where ``mol.coordinates`` and ``mol.numbers`` are read from file (see
:ref:`read-molgeometry`), and ``cc-pvdz`` is the cc-pVDZ basis set.


Specifying different basis sets for different atoms
---------------------------------------------------

In some cases, you may want to specify different basis sets for different atoms. For
example, you might like to use the 3-21G basis set for the hydrogen atom, the
6-31G basis set for the carbon atom, and STO-3G for all remaining atoms:

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         element_map={'H':'3-21g', 'C':'6-31g'})

where `mol.coordinates` and `mol.numbers` are read from file (see
:ref:`read-molgeometry`),  and ``sto-3g``, ``3-21g`` and ``6-31g`` are the basis
set names (see :ref:`ref_gaussian_basis_standard_sets`)

Alternatively, the same result can be obtained by substituting the H and C symbols
with their atomic numbers:

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         element_map={1:'3-21g', 6:'6-31g'})

You can also override the default basis for selected atoms based on their index,
i.e. position in the list of atoms that specify the molecule:

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         index_map={0:'3-21g', 2:'6-31g'})

The above example uses the ``3-21g`` basis for the first atom, the ``6-31g``
basis for the third atom and the ``sto-3g`` basis for all other atoms.


Loading custom basis sets from file
-----------------------------------

You can also use other basis sets besides the ones that are shipped with HORTON.
It is assumed that the basis is available in NWChem format:

.. code-block:: python

    mybasis = GOBasisFamily('myname', filename='mybasis.nwchem'),
    obasis = get_gobasis(mol.coordinates, mol.numbers, mybasis)

Anywhere you can specify a built-in basis set with a string, you can also use
instance of the ``GOBasisFamily`` class (``mybasis`` in the example above), e.g.
in the arguments ``default``, ``element_map`` and ``index_map`` of
``get_gobasis``.


Defining basis sets with Python code
------------------------------------

In some circumstances, it may be useful to generate the basis set with some
Python code. For example, the following code generates an even tempered basis
for Lithium (without polarization functions):

.. literalinclude:: ../data/examples/hamiltonian/even_tempered_li.py
    :lines: 2-
    :caption: data/examples/hamiltonian/even_tempered_li.py

All basis functions in this example are just single s-type primitives, i.e. no
contactions are used. At the end of the example, the basis set is constructed
for a single Li atom in the origin.

Note that ``get_gobasis`` also accepts instances of GOBasisAtom for the
arguments ``default``, ``element_map`` and ``index_map``.


.. _user_molecularham_geom_and_basis:

Loading geometry and basis set info from one file
-------------------------------------------------

When post-processing results from other programs, it may be desirable to use
exactly the same basis set that was used in the other program (Gaussian, ORCA,
PSI4, etc). This can be achieved by loading the geometry, basis set and
wavefunction from one of the following formats: ``.mkl``, ``.molden``,
``.fchk``, or HORTON's internal ``.h5`` format. In principle, it is also
possible to load the basis set from a ``.wfn`` file, but keep in mind that this
format does not support contracted basis functions, so HORTON will then use a
decontracted basis set, which is less efficient.

You simply use the ``IOData.from_file`` method to load the file. The orbital
basis object is then available in the ``obasis`` attribute if the return value.
For example:

.. code-block:: python

    # Load the geometry, orbital basis and wavefunction from a Gaussian
    # formatted checkpoint file:
    mol = IOData.from_file('water.fchk')

    # Print the number of basis functions
    print mol.obasis.nbasis


.. _user_molecularham_matrix_elements:

Computing (Hamiltonian) matrix elements
=======================================

Given a ``GOBasis`` instance (the ``obasis`` object from the
examples in the previous section), you can generate the two-index and four-index
objects for the molecular electronic Hamiltonian. It may also be useful to
construct the overlap operator as the Gaussian basis sets are not orthonormal.

Before computing the matrix elements, you first have to specify how the two- and four-index
objects will be represented. By default, HORTON uses a dense
matrix representation, which is implemented in the ``DenseLinalgFactory`` class.
An instance of this class must be passed to the methods that compute the matrix
elements. Alternatively, you may also use the ``CholeskyLinalgFactory``, which
represents all four-index objects with a Cholesky decomposition. Note that
the four-center electron-repulsion integrals are computed with LibInt
[valeev2014]_.

This is a basic example, assuming some of the preceding code has created the
``obasis`` object:

.. literalinclude:: ../data/examples/hf_dft/rhf_water_dense.py
    :lines: 22-29
    :caption: data/examples/hf_dft/rhf_water_dense.py, lines 22--29

For the nuclear attraction integrals, you also have to specify arrays with atomic
coordinates and nuclear charges.
