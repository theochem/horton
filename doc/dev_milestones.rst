Milestones
##########

This is a formal version of our wish list, grouped into reasonable chunks that
belong together. I probably forgot a bunch of items, so feel free to add
stuff. Ideally we'd implement the most fancy features first, as they will
lead to good papers. Anyway, the minimal program is a milestone that needs to
be done first.

For every item, there are a few possible states:

* TODO (or no label): not started yet

* WIP: work in (slow) progress

* DONE: implementation and test routines are ready

* DOCUMENTED: the implementation and its usage are documented somewhere on this website.


Minimal program
===============

Basic HF single-point computation on small molecules. We will need the
following components:

1. [DONE] A basic unit conversion module, periodic table information, constants.

2. [DONE] Basis set parser (NWChem format)

3. [DONE] Compact array representation of the basis set

4. [DONE] Abstraction layer for large matrices and operations on these objects.
   This is crucial to make the program fast. We have to keep the following in
   mind:

   * Allocation of a big matrix can never happen accidentally. One must call
     allocate or copy routines to make new matrices.

   * Support all kinds of in-place operations.

   * Allocation should only happen prior to an actual computation,
     deallocation only in the end. This may seen trivial, but in Python
     allocation/deallocation is done automatically. This is fine for small
     stuff, but not for big arrays.

   * In the long run, this part is going to be extended with alternative
     implementations, e.g. sparse-matrices, parallel operations, .... The rest
     of the code should not be affected by switching to different matrix
     implementations.

   * The first version does not have to be perfect, but it is should at least
     have an API that allows us to replace this initial version with faster
     alternatives.

5. [DONE] Gaussian integral routines: overlap, kinetic energy, nuclear
   attraction and two-electron. Initially, these are borrowed from `PyQuante
   <http://pyquante.sourceforge.net/>`__ to get started.

6. [DOCUMENTED] transformation between Cartesian and pure basis functions.

7. [WIP] Generic one- and two-body contraction routines to generate the dense
   matrix representation of the operators. In order to test this properly, we
   will have to add some IO routines to read validation results generated with
   other codes

8. [WIP] A system class that can keep track of molecular geometry, wavefunction,
   operators, ...

9. A simple class hierarchy to define the Hamiltonian.

10. Hamiltonian-core guess and basic SCF routine.

The Python/C wrapping will be done with Cython. The Python/Fortran wrapping
will be done with Fwrap.


Less urgent changes in basic version
====================================


1. Our own implementation of the one-electron integrals. Link with `libint2
   <http://sourceforge.net/projects/libint/>`__ for the two-electron integrals.

2. Compact storage and somewhat clever linear aglebra routines for the dense
   matrix objects.

3. `HDF5 <http://en.wikipedia.org/wiki/Hierarchical_Data_Format>`__ storage of
   system class, cfr Gaussian checkpoint file. This will be based on the `h5py
   module <http://code.google.com/p/h5py/>`__.

4. Transformation of wavefunctions from one basis set to the other.


The fun part (students)
=======================

Add more features to make the program more fun:

1. Support a bunch of basis sets. (Figure out to what extent we can redistribute
   part of the `EMSL <https://bse.pnl.gov/bse/portal>`__ data.)

2. Add support for LDA. For this we need to implement Becke's integration
   scheme and good spherical integration routines. We can initially take these
   grid routines from HiPart, but it would be nice to have the Lebedev grids
   implemented in `SymPy <http://code.google.com/p/sympy/>`__, which uses
   `MPMath <http://code.google.com/p/mpmath/>`__ for arbitrary precision
   evaluation.

3. Implement some more functionals. We should make use of `Libxc
   <http://www.tddft.org/programs/octopus/wiki/index.php/Libxc>`__.


Thumper and friends
===================

Computation of nuclear forces. Include Sandra's optimizer.


Constrain the wavefunction
==========================

Arbitrary constraints (atomic charge, bond orders, ...) during the wavefunction
optimization.


Tickle the wavefunction
=======================

Response of the wavefunction to perturbations. This has several applications:

1. Vibrational frequencies.
2. Polarizability.
3. Condensed linear response matrix.
4. ...


Polymers, membranes and crystals
================================

1D, 2D and 3D periodic boundary conditions.


Pseudopotentials
================


Quantum dots
============

Non-uniform Gaussian basis functions.


Cut the crap
============

Decent partitioning schemes.


QM/MM interface
===============

An onion-like scheme, with long-range polarization.
