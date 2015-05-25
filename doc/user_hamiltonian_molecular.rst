Defining a molecular Hamiltonian
################################

A molecular Hamiltonian is typically set up in three steps. First one loads or
constructs a molecular geometry. Then one generates a Gaussian basis set for
this molecule. Finally, all kinds of matrix elements can be computed with that
basis to define a Hamiltonian.


.. _setup-molgeometry:

Specifying the molecular geometry
=================================

Supported file formats
----------------------

Horton can load/dump the molecular geometry from/to different file formats. The
file format is determined automatically using the extension or prefix of the
filename.

================================ ========================= =============== ======================= ======================
 File format                      Extension                 Prefix           Read (``from_file``)    Write (``to_file``)
================================ ========================= =============== ======================= ======================
 HDF5 file or h5py.Group object   .h5                                        X                       X
 xyz                              .xyz                                       X                       X
 Gaussian                         .fchk                                      X
 Gaussian                         .log                                       X
 Gaussian cube                    .cube                                      X                       X
 Gaussian wfn                     .wfn                                       X
 Molekel                          .mkl                                       X
 Molden                           .molen.input or .molden                    X                       X
 VASP                                                       POSCAR           X                       X
 VASP                                                       CHGCAR, AECCAR   X
 VASP                                                       LOCPOT           X
 CP2K (atom)                      .cp2k                                      X
 CIF (partially working)          .cif                                       X
================================ ========================= =============== ======================= ======================


.. _read-molgeometry:

Reading the molecular geometry from file
----------------------------------------

The molecular geometry can be read from file using the method
:py:meth:`horton.io.molecule.Molecule.from_file` of the
``Molecule`` class,

.. code-block:: python

    mol = Molecule.from_file(filename)


Constructing a molecular geometry from scratch
----------------------------------------------

One may also generate molecular geometries with Python code. The
following example constructs a ring of Hydrogen atoms and writes it to an
XYZ file

.. literalinclude:: ../data/examples/hamiltonian/hydrogen_ring.py
    :lines: 2-
    :caption: data/examples/hamiltonian/hydrogen_ring.py


.. _user_molecularham_basis:

Specifying the basis set
========================

Basis set format
----------------
Horton supports basis sets consisting of generally contracted Cartesian Gaussian functions up to a maximum angular momentum of seven (I functions).
Horton is using the same basis set format as NWChem, and the basis sets can be downloaded from the EMSL webpage (https://bse.pnl.gov/bse/portal).

.. _horton-basis-library:

Basis set library
-----------------
Most of the popular basis sets can be loaded automatically. A list of currently supported built-in sets is tabulated below.

.. cssclass:: table-striped

.. include:: basis.rst.inc

Note that the basis set names are case-insensitive in Horton. These basis sets
were taken from the EMSL library (https://bse.pnl.gov/bse/portal). When
publishing results obtained with these basis sets, please cite the following
references [feller1996]_ and [Didier2007]_.

The basis set for a given molecule is constructed with the function :py:func:`horton.gbasis.gobasis.get_gobasis`


Unique basis set for all atoms
------------------------------

Usually one wants to define a unique (the same) basis set for the whole system, this can be done by a function call.

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')

where ``mol.coordinates`` and ``mol.numbers`` are read from file (see :ref:`read-molgeometry`), and ``cc-pvdz`` is the cc-pVDZ basis set.


Specifying different basis sets for different atoms
---------------------------------------------------

In some cases one wants to specify different basis sets for different atoms. For example, setting the 3-21G basis set for the hydrogen atom and the 6-31G basis set for the carbon atom, and STO-3G for all remainig atoms can be done as follows.

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         element_map={'H':'3-21g', 'C':'6-31g'})

where `mol.coordinates` and `mol.numbers` are read from file (see :ref:`read-molgeometry`),  and ``sto-3g``, ``3-21g`` and ``6-31g`` are the basis set names (see :ref:`horton-basis-library`)

Alternatively, the same result can be obtained by substituting the H and C symbols with their atomic numbers:

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         element_map={1:'3-21g', 6:'6-31g'})

One may also override the default basis for selected atoms based on their index,
i.e. position in the list of atoms that specify the molecule:

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         index_map={0:'3-21g', 2:'6-31g'})

The above example uses the ``3-21g`` basis for the first atom, the ``6-31g``
basis for the third atom and the ``sto-3g`` basis for all other atoms.


Loading custom basis sets from file
-----------------------------------

One can also use other basis sets besides the ones that are shipped with Horton.
It is assumed that the basis is available in NWChem format:

.. code-block:: python

    mybasis = GOBasisFamily('myname', filename='mybasis.nwchem'),
    obasis = get_gobasis(mol.coordinates, mol.numbers, mybasis)

Anywhere one can specify a built-in basis set with a string, one my also use
instance of the ``GOBasisFamily`` class (``mybasis`` in the example above), e.g.
in the arguments ``default``, ``element_map`` and ``index_map`` of
``get_gobasis``.


Defining basis sets with Python code
------------------------------------

In some circumstances it may be useful to generate the basis set with some
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
exactly the same basis set as was used in the other program (Gaussian, ORCA,
PSI4, etc.) This can be acchieved by loading the geometry, basis set and
wavefunction from one of the following formats: ``.mkl``, ``.molden``,
``.fchk``, or Horton's internal ``.h5`` format. In principle, it is also
possible to load the basis set from a ``.wfn`` file, but keep in mind that this
format does not support contracted basis functions, so Horton will then use a
decontracted basis set, which is less efficient.

One simply uses the ``Molecule.from_file`` method to load the file. The orbital
basis object is then available in the ``obasis`` attribute if the return value.
For example:

.. code-block:: python

    # Load the geometry, orbital basis and wavefunction from a Gaussian
    # formatted checkpoint file:
    mol = Molecule.from_file('water.fchk')

    # Print the number of basis functions
    print mol.obasis.nbasis


.. _user_molecularham_matrix_elements:

Computing (Hamiltonian) Matrix elements
=======================================

Given a ``GOBasis`` instance (the ``obasis`` object from the
examples in the previous section), one can generate the two-index and four-index
objects for the molecular electronic Hamiltonian. It may be useful to construct
also the overlap operator as the Gaussian basis sets are not orthonormal.

Before computing the matrix elements, one first has to specify how the two- and four-index
objects will be represented. The default in Horton is to use a dense
matrix representation, which is implemented in the ``DenseLinalgFactory`` class.
An instance of this class must be passed to the methods that compute the matrix
elements. Alternatively, one may also use the ``CholeskyLinalgFactory``, which
represents all four-index objects with a Cholesky decomposition.

This is a basic example, assuming some of the preceding code has created the
``obasis`` object:

.. literalinclude:: ../data/examples/hf_dft/rhf_water_dense.py
    :lines: 18-25
    :caption: data/examples/hf_dft/rhf_water_dense.py, lines 18--25

For the nuclear attraction integrals, one also has to specify arrays with atomic
coordinates and nuclear charges.
