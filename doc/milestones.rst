Milestones
##########

This is a formal version of our wish list, grouped into reasonable chunks that
belong together. I probably forgot a bunch of items, so feel free to add
stuff. Ideally we'd implement the most fancy features first, as they will
lead to good papers. Anyway, the minimal program is a milestone that needs to
be done first.


Minimal program
===============

Basic HF/STO-3G single-point computation on small molecules. We will need the
following components for this to work:

1. Basis set parser (e.g. Gaussian 94 format)

2. Compact array representation of the basis set

3. Abstraction layer for large matrices and operations on these objects. This
   is crucial to make the program fast. We have to keep the following in mind.

   * Allocation of a big matrix can never happen accidentally. One must call
     allocate or copy routines to make new matrices.

   * Support all kinds of in-place operations.

   * Allocation should only happen prior to an actual computations,
     deallocation only in the end. This may seen trivial, but in Python
     allocation/deallocation is done automatically. This is fine for small
     stuff, but not for big arrays.

   * Most matrices are triangular, so we can use that to reduce the memory
     footprint.

   * In the long run, this part is going to be replaced with a
     sparse-matrix module. The rest of the code should not be affected by
     switching to sparse-matrix code.

4. Gaussian integral routines. We'll probably borrow some from PyQuante to get
   started.

5. Matrix element routines: use basis set structure and Gaussian integral
   routines. We need these to compute the derivatives of the total Lagrangian
   towards the density matrix elements. Expectation values of the
   corresponding operators can be evaluated simultaneously. The following
   operators are needed.

   * Overlap
   * Kinetic energy
   * Hartree & Exchange

6. Hamiltonian-core guess.

7. A Constrained minimizer.

The Python/C wrapping will be done with Cython. The Python/Fortran wrapping
will be done with Fwrap.


The fun part (students)
=======================

Add more features to make the program more fun:

1. Support a bunch of basis sets

2. Add support for LDA. For this we need to implement Becke's integration
   scheme and good spherical integration routines.

3. Implement some more functionals: XXX.

We will need to add a lot of testing here, just to make sure the program works
like it should.


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
3. Softness matrix.
4. ...

Polymers, membranes and crystals
================================

1D, 2D and 3D periodic boundary conditions.


Quantum dots
============

Non-uniform Gaussian basis functions.


Cut the crap
============

Decent partitioning schemes.
