Molecular Hamiltonians
######################

A molecular Hamiltonian is typically set up in three steps. First one loads
a molecular geometry, or constructs the arrays with coordinates from scratch.
Then one generates a Gaussian basis set for this molecule. Finally, all kinds
of matrix elements can be computed with that basis to define a Hamiltonian.


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

The molecular geometry can be read from file using the ``Molecule`` class,

.. code-block:: python

    mol = Molecule.from_file(filename)

with arguments

    :filename: (str) the filename with proper extension or prefix


Constructing a molecular geometry from scratch
----------------------------------------------

One may also generate molecular geometries with Python code. The
following example constructs a ring of Hydrogen atoms and writes it to an
XYZ file

.. code-block:: python

    import numpy as np
    from horton import *

    # define the ring
    natom = 11
    spacing = 1.3 # distance between two neighboring atoms in bohr
    radius = spacing/(2*np.sin(np.pi/natom))

    # define the coordinates and elements
    coordinates = np.zeros((natom, 3))
    numbers = np.ones(natom, dtype=int) # must be integers
    for iatom in xrange(natom):
        angle = (2*np.pi/natom)*iatom
        coordinates[iatom, 0] = radius*np.cos(angle)
        coordinates[iatom, 1] = radius*np.sin(angle)

    # write the molecule to an XYZ file (optional)
    mol = Molecule(coordinates=coordinates, numbers=numbers)
    mol.to_file('ring.xyz')


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

======================== =========================== ===========================
Common basis name           Horton library name         Supported elements
======================== =========================== ===========================
STO-3G                     ``sto-3g``                         H--Kr
STO-6G                   ``sto-6g``                      H--Kr
3-21G                    ``3-21g``                       H--Kr
3-21G*                   ``3-21g*``                      Na--Ar
3-21+G*                  ``3-21+g*``                     Na--Ar
4-31G                    ``4-31g``                       H--Ne, P, S, Cl
6-31G                    ``6-31g``                       H--Zn
6-31G*                   ``6-31g*``                      H--Zn
6-31G**                  ``6-31g**``                     H--Zn
6-31+G                   ``6-31+g``                      H--Ca
6-31+G*                  ``6-31+g*``                     H--Ca
6-31++G*                 ``6-31++g*``                    H--Ca
6-31++G**                ``6-31++g**``                   H, Li--Ca
cc-pVDZ                  ``cc-pvdz``                     H--Ar, Ca-Kr
cc-pVTZ                  ``cc-pvtz``                     H--Ar, Ca-Kr
cc-pVQZ                  ``cc-pvqz``                     H--Ar, Ca-Kr
cc-pCVDZ                 ``cc-cpvdz``                    Li--Ar
cc-pCVTZ                 ``cc-cpvtz``                    Li--Ar
cc-pCVQZ                 ``cc-cpvqz``                    Li--Ar
aug-cc-pVDZ              ``aug-cc-pvdz``                 H--Ar, Sc--Kr
aug-cc-pVTZ              ``aug-cc-pvtz``                 H--Ar, Sc--Kr
aug-cc-pVQZ              ``aug-cc-pvqz``                 H--Ar, Sc--Kr
aug-cc-pCVDZ             ``aug-cc-cpvdz``                Li--Ar
aug-cc-pCVTZ             ``aug-cc-cpvtz``                Li--Ar
aug-cc-pCVQZ             ``aug-cc-cpvqz``                Li--Ar
def2-svpd                ``def2-svpd``                   H--Kr
def2-tzvp                ``def2-tzvp``                   H--Kr
def2-tzvpd               ``def2-tzvpd``                  H--Kr
def2-qzvp                ``def2-qzvp``                   H--Kr
def2-qzvpd               ``def2-qzvpd``                  H--Kr
======================== =========================== ===========================

Note that the basis set names are case-insensitive in Horton. These basis sets
were taken from the EMSL library (https://bse.pnl.gov/bse/portal). When
publishing results obtained with these basis sets, please cite the following
references [feller1996]_ and [Didier2007]_.

The basis sets are loaded using the following function call.

.. code-block:: python

   get_gobasis(coordinates, numbers, default, element_map=None, index_map=None, pure=True)

with arguments

    :coordinates: (float)  A (N, 3) numpy array with Cartesian coordinates of the atoms (see :ref:`read-molgeometry`)
    :numbers: (int) A (N,) numpy vector with the atomic numbers (see :ref:`read-molgeometry`)
    :default: (str) The basis set name applied to each atom.

with optional arguments

    :element_map: (str) A dictionary with element names or numbers as keys, and basis sets as values. These specs override the default basis
    :index_map: (str) A dictionary with atomic indexes (based on the order of the atoms) as keys and basis sets as values
    :pure: (bool) By default pure (spherical) basis functions are used. Set this to false to switch to Cartesian basis functions

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

.. code-block:: python

    import numpy as np
    from horton import *

    # specify the even tempered basis set
    alpha_low = 5e-3
    alpha_high = 5e2
    nbasis = 30
    lnratio = (np.log(alpha_high) - np.log(alpha_low))/(nbasis-1)

    # build a list of "contractions". These aren't real contractions as every
    # constraction only contains one basis function.
    bcs = []
    for ibasis in xrange(nbasis):
        alpha = alpha_low**lnratio
        # arguments of GOBasisContraction:
        #     shell_type, list of exponents, list of contraction coefficients
        bcs.append(GOBasisContraction(0, [alpha], [1.0]))

    # Finish setting up the basis set:
    ba = GOBasisAtom(bcs)
    obasis = get_gobasis(np.array([[0.0, 0.0, 0.0]]), np.array([3]), default=ba)

All basis functions are just single s-type primitives, i.e. no contactions are
used. In the end of the example, the basis set is constructed for a single Li
atom in the origin.

Note that ``get_gobasis`` also accepts instances of GOBasisAtom for the
arguments ``default``, ``element_map`` and ``index_map``, as shown in the above
example.


Computing (Hamiltonian) Matrix elements
=======================================

TODO
