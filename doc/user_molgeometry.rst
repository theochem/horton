Specifying the molecular geometry
#################################

.. contents::
.. _setup-molgeometry:

Supported file formats
======================

Horton can load the molecular geometry from different file formats. The file format is determined automatically using the extension or prefix of the filename.


======================== =========================== ===========================
File format              Extension                   Prefix
======================== =========================== ===========================
h5.Group                 .h5
xyz                      .xyz
Gaussian                 .fchk
Gaussian                 .log
Gaussian cube            .cube
Gaussian wfn             .wfn
Molekel                  .mkl
Molden                   .molen.input or .molden
x                                                    POSCAR
x                                                    CHGCAR, AECCAR
x                                                    LOCPOT
CP2K                     .cp2k
CIF                      .cif
======================== =========================== ===========================

.. _read-molgeometry:

Reading the molecular geometry from file
========================================

The molecular geometry can be read from file using the ``Molecule`` class,

.. code-block:: python

    mol = Molecule.from_file(filename)

with arguments

    :filename: (str) the filename with proper extension or prefix
