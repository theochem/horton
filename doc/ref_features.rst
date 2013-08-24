.. _ref_features:

Features
########

The Horton library
==================

Each package below has a ``test`` subpackage that contains unit tests.


Package ``horton``
------------------

* Datastructures for molecular systems, unit cell and symmetry.
* Checkpointing.
* Unit conversion and physical constants.
* Periodic table.
* Linear algebra abstraction layer.
* Caching tools to avoid recomputation of earlier results.


Package ``horton.espfit``
-------------------------

* Various ESP cost functions to fit charges for 3D periodic systems.
* Evaluation of ESP on grids for 3D periodic systems.


Package ``horton.gbasis``
-------------------------

* Datastructures for Gaussian basis sets (Cartesian and pure).
* Generator for standard basis sets.
* Integrals in terms of Gaussian basis functions.
* Evaluation of density, potential and density gradient on grids.


Package ``horton.grid``
-----------------------

* Molecular integration grids (Becke partitioning with Lebedev Laikov atom
  centered grids).
* Pruned grids for elements H-La, Hf-Rn.
* 1D integration grids.
* 1D Cubic splines.
* Uniform grids (cube files).


Package ``horton.meanfield``
------------------------------

* Definition of HF and DFT Hamiltonians.
* `LibXC <http://www.tddft.org/programs/octopus/wiki/index.php/Libxc>`_ is used
  to provide a large number of exchange and correlation functionals.
  [marques2012]_
* Custom terms that can be added as external potential.
* Wavefunction classes, initial guess and projection.
* SCF algorithms (regular and optimal damping).


Package ``horton.io``
---------------------

The following file formats are supported.

================================================== ================== ==== =====
Format                                             Pattern            Read Write
================================================== ================== ==== =====
Internal checkpointing format based on HDF5        ``*.chk`` ``*.h5`` X    X
CIF (crystalographic information file)             ``*.cif``          X    X
Atomic wavefunctions from CP2K (2.4-r12857)        ``*.cp2k.out``     X
Gaussian cube files                                ``*.cube``         X    X
Gaussian log files (operators)                     ``*.log``          X
Gaussian format checkpoint files                   ``*.fchk``         X
Molden input files (wavefunction)                  ``*.molden.input`` X    X
Molekel file                                       ``*.mkl``          X
Gaussian / GAMESS WFN file                         ``*.wfn``          X
VASP CHGCAR and LOCPOT files                                          X
VASP POSCAR                                                           X    X
XYZ format                                         ``*.xyz``          X    X
================================================== ================== ==== =====


Package ``horton.part``
-----------------------

The following density partitioning schemes are implemented:

* ``becke``: Becke partitioning [becke1988_multicenter]_
* ``h``: Hirshfeld partitioning [hirshfeld1977]_
* ``hi``: Iterative Hirshfeld partitioning [bultinck2007]_
* ``is``: Iterative Stockholder partitioning [lillestolen2008]_
* ``he``: Extended Hirshfeld partitioning [verstraelen2013]_

Partitioning can be carried on isolated and periodic systems. Both all-electron (AE)
and pseudo-potential (PP) densities are supported. The following AIM properties can
be computed:

* Atomic charges, population, pseudo_population, spin charges.
* Cartesian, pure (harmonic) and radial atomic moments.
* Tkatchenko-Scheffler dispersion coefficients. [tkatchenko2009]_
* Decomposition of each AIM in spherical harmonics and the derivation of the
  hartree potential for each component. [becke1988_poisson]_

Note that not all AIM properties work for any combination partitioning scheme
and boundary condition for technical reasons.



Horton command-line scripts
===========================


Density partitioning
--------------------

* ``horton-atomdb.py``: generate database of isolated atoms
* ``horton-cpart.py``: partition densities on uniform grids (cube and LOCPOT files)
* ``horton-wpart.py``: partition densities derived from wavefunctions in a Gaussian basis (fchk, mkl, molden.input)



ESP fitting
-----------

* ``horton-esp-cost.py``: set up a cost function for ESP fitting or testing
* ``horton-esp-fit.py``: derive atomic charges from an ESP cost function
* ``horton-esp-gen.py``: generate the electrostatic potential due to a set of charges
* ``horton-esp-test.py``: test how well charges perform for a given cost function


File conversion
---------------

* ``horton-convert.py``: convert between file formats supported by Horton
* ``horton-hdf2csv.py``: convert the contents of an HDF5 file to a csv file
