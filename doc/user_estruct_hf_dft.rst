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

.. _user_hf_dft:

How to use HORTON as a Hartree-Fock/DFT program
###############################################

This section explains how the relevant Hartree-Fock and DFT features in HORTON
are used with Python scripts. This is different from most regular quantum
chemistry codes, where you just prepare an input file with a molecular geometry
and some control options. Such an input-file interface for HORTON is under
construction. The script interface, explained here, is more expressive and
laborsome than a regular input file, because you basically have to write the
main program. However, it has the advantage of flexibility, since you can tinker
more with the individual steps in the calculation.

HORTON can do restricted and unrestricted HF/DFT SCF calculations in various
different ways. The following sections cover the typical steps of a calculation:

1) :ref:`user_hf_dft_setup_hamiltonian`
2) :ref:`user_hf_dft_initial_guess`
3) :ref:`user_hf_dft_effective_ham`
4) :ref:`user_hf_dft_occupation`
5) :ref:`user_hf_dft_scf`
6) Optional (if needed for 7 or 8): :ref:`user_hf_dft_get_orbitals`
7) Optional: :ref:`user_hf_dft_to_file`
8) Optional: :ref:`user_hf_dft_preparing_posthf`

The last section contains an overview of :ref:`hf_dft_complete_examples`
that are shipped (and tested) with HORTON. These may be convenient as a starting
point when preparing your own scripts. You may also just dive straight
into these examples and consult the first five sections to better understand how
each example works.

The code snippets below use the following conventions for variable names. The
following are defined by setting up the Hamiltonian:

* ``olp``: two-index object with the overlap integrals.
* ``kin``: two-index object with the kinetic energy integrals.
* ``na``: two-index object with the nuclear attraction integrals.
* ``er``: four-index object with the electron repulsion integrals.

The following names are also used systematically:

* ``exp_alpha``, ``exp_beta``: the expansion coefficient of the alpha and beta
  orbitals in a (local) basis.


.. _user_hf_dft_setup_hamiltonian:

Setting up the (molecular electronic) Hamiltonian
=================================================

See :ref:`user_hamiltonian`.

You first have to compute/load the two- and four-index objects for the one- and
two-body terms in the Hamiltonian. Some DFT implementations in HORTON do not
require pre-computed four-center integrals.


.. _user_hf_dft_initial_guess:

Generating an initial guess of the wavefunction
===============================================

.. _user_hf_dft_core_guess:

Core Hamiltonian guess
----------------------

This guess involves only the one-body terms of the electronic Hamiltonian and
completely neglects electron-electron interactions. The orbitals for such a
one-body Hamiltonian can be computed without any prior guess.

The function :py:func:`~horton.meanfield.guess.guess_core_hamiltonian` can be
used to compute the core hamiltonian guess. This type of guess only relies on a
simple one-body Hamiltonian, usually consisting of the kinetic energy and the
interaction with the external field. In case of an unrestricted calculation,
this means that the same initial orbital guess is made for the alpha and beta
spins.

The guess for a restricted wavefunction is done as follows:

.. literalinclude:: ../data/examples/hf_dft/rhf_water_dense.py
    :lines: 31-35
    :caption: data/examples/hf_dft/rhf_water_dense.py, lines 31--35

For an unrestricted wavefunction, the procedure is very similar:

.. literalinclude:: ../data/examples/hf_dft/uhf_methyl_dense.py
    :lines: 27-32
    :caption: data/examples/hf_dft/uhf_methyl_dense.py, lines 27--32

The arguments ``exp_alpha`` and ``exp_beta`` are treated as output arguments.
Instead of ``kin`` and ``na``, you may provide just any set of one-body
operators to construct different types of guesses.


Randomizing an initial guess
----------------------------

The method :py:meth:`~horton.matrix.dense.DenseExpansion.rotate_random` can be
used to apply a random unitary rotation to all the orbitals (occupied and virtual).
The orbital energies and occupation numbers are not affected. For example, after
constructing a :ref:`user_hf_dft_core_guess`, the following line randomizes the
orbitals:

.. code-block:: python

    # randomly rotate the orbitals (irrespective of occupied or virtual)
    exp_alpha.rotate_random()


Modifying the initial guess
---------------------------

If needed, you may fine-tune the initial guess by making fine-grained
modifications to the orbitals. (These may also be useful for fixing the orbitals
that come out of a failed SCF.)

* The method :py:meth:`~horton.matrix.dense.DenseExpansion.rotate_2orbitals`
  allows you to mix two orbitals. By default, it rotates the HOMO and LUMO
  orbitals by 45 degrees:

  .. code-block:: python

      # Mix HOMO and LUMO orbitals
      exp_alpha.rotate_2orbitals()

      # Rotate 1st and 6th orbital by 30 deg
      exp._alpha.rotate_2orbitals(30*deg, 0, 5)

  Note that HORTON uses radians as unit for angles, i.e. ``30*deg == np.pi/6``.
  Also, HORTON uses zero-based indices, so the arguments ``0, 5``
  refer to the first and the sixth orbital.

* The method :py:meth:`~horton.matrix.dense.DenseExpansion.swap_orbitals` allows
  you to swap several orbitals. It takes as an argument an array where each row
  is a pair of orbitals to swap. For example, the following swaps 1st and 3rd,
  followed by a swap of 2nd and 4th:

  .. code-block:: python

      # Swap some orbitals
      swaps = np.array([[0, 2], [1, 3]])
      exp_alpha.swap_orbitals(swaps)


Reading a guess from a file
---------------------------

You may also load orbitals from an external file. The file formats ``.mkl``,
``.molden``, ``.fchk``, or HORTON's internal ``.h5`` are all viable
sources of orbitals. For example, the orbitals from a Gaussian formatted
checkpoint file, ``*.fchk``, may be loaded as follows:

.. code-block:: python

    # Load fchk file
    mol = IOData.from_file('water.fchk')

    # Print the number of alpha orbitals (occupied and virtual)
    print mol.exp_alpha.nfn

Obviously, if you would like to use these orbitals without projecting them onto
a new basis set (as explained in :ref:`user_hf_dft_project_basis`), you are
forced to continue working in the same basis set, which can be accessed
in this example as ``mol.obasis``. See :ref:`user_molecularham_geom_and_basis`
for more details.


.. _user_hf_dft_project_basis:

Projecting orbitals from a smaller basis onto a larger one
----------------------------------------------------------

Assuming you have obtained (converged) orbitals in a smaller basis set, you can
try to use these as initial guess after projecting the orbitals onto the
larger basis set. This is exactly what the function
:py:func:`~horton.meanfield.project.project_orbitals_mgs` does. The following
snippet assumes that the ``obasis0`` and ``exp_alpha0`` are the small basis set
and a set of orbitals in that basis for the ``IOData`` instance ``mol``.

.. code-block:: python

    # Definition of the bigger basis set
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'aug-cc-pvtz'):

    # Linalg factory for the bigger basis set
    lf1 = DenseLinalgFactory(obasis1.nbasis)

    # Create a expansion object for the alpha orbitals in the large basis
    exp_alpha1 = lf1.create_expansion()

    # The actual projection
    project_orbitals_msg(obasis0, obasis1, exp_alpha0, exp_alpha1)


.. _user_hf_dft_effective_ham:

Defining the effective Hamiltonian
==================================

HORTON implements spin-restricted and spin-unrestricted effective Hamiltonians.
Mathematically speaking, these are models for the energy as function of a set of
density matrices. The implementation also provides an API to compute the
corresponding Fock matrix for every density matrix, i.e. the derivative of the
energy toward the density matrix elements.

* For the restricted case, the alpha and beta density matrices are assumed
  to be identical. Hence the energy is only a function of the alpha density
  matrix. When constructing the Fock matrix, the derivative is divided by two
  to obtain such that the Fock matrix has the conventional orbital energies as
  eigenvalues.

  .. math::
      D^\alpha &\rightarrow E(D^\alpha) \\
               &\rightarrow F^\alpha_{\mu\nu} = \frac{1}{2}\frac{\partial E}{\partial D^\alpha_{\nu\mu}}

* For the unrestricted case, the alpha and beta density matrices are allowed to
  be different. Hence, there are also alpha and beta Fock matrices.

  .. math::
      D^\alpha, D^\beta &\rightarrow E(D^\alpha, D^\beta) \\
                        &\rightarrow F^\alpha_{\mu\nu} = \frac{\partial E}{\partial D^\alpha_{\nu\mu}} \\
                        &\rightarrow F^\beta_{\mu\nu} = \frac{\partial E}{\partial D^\beta_{\nu\mu}}

A generic `API
<http://en.wikipedia.org/wiki/Application_programming_interface>`_ for such
density matrix functionals is implemented in the classes
:py:class:`~horton.meanfield.hamiltonian.REffHam` and
:py:class:`~horton.meanfield.hamiltonian.UEffHam`. The prefixes ``R`` and ``U``
are used (also below) to differentiate between restricted and unrestricted
implementations. A Hartree-Fock or DFT effective Hamiltonian is defined by
constructing an instance of the ``REffHam`` or ``UEffHam`` classes and providing
the necessary energy terms to the constructor.


Supported energy terms for the effective Hamiltonians
-----------------------------------------------------

All classes below take a ``label`` argument to give each term in the Hamiltonian
a name, e.g. used for storing/displaying results. For each class listed below,
follow the hyperlinks to the corresponding documentation for a description of
the constructor arguments.

* Simple one-body terms are specified with
  :py:class:`~horton.meanfield.observable.RTwoIndexTerm`, or
  :py:class:`~horton.meanfield.observable.UTwoIndexTerm`.

* The direct term of a two-body interaction is specified with
  :py:class:`~horton.meanfield.observable.RDirectTerm`, or
  :py:class:`~horton.meanfield.observable.UDirectTerm`.

* The exchange term of a two-body interaction is specified with
  :py:class:`~horton.meanfield.observable.RExchangeTerm`, or
  :py:class:`~horton.meanfield.observable.UExchangeTerm`.

* Functionals of the density (or its derivatives) that require numerical
  integration are all grouped into on term using
  :py:class:`~horton.meanfield.gridgroup.RGridGroup`, or
  :py:class:`~horton.meanfield.gridgroup.UGridGroup`. This makes it possible
  to compute at every SCF iteration the density (and its gradients) only once
  for all terms that depend on the density. This also allows for a similar gain
  in efficiency when building the Fock matrices. The constructor of a
  ``GridGroup`` class takes a numerical integration grid and a list of instances
  of the following classes as arguments:

    * An LDA functional from LibXC can be specified with
      :py:class:`~horton.meanfield.libxc.RLibXCLDA` or
      :py:class:`~horton.meanfield.libxc.ULibXCLDA`.

    * A GGA functional from LibXC can be specified with
      :py:class:`~horton.meanfield.libxc.RLibXCGGA` or
      :py:class:`~horton.meanfield.libxc.ULibXCGGA`.

    * A Hybrid GGA functional from LibXC can be specified with
      :py:class:`~horton.meanfield.libxc.RLibXCHybridGGA` or
      :py:class:`~horton.meanfield.libxc.ULibXCHybridGGA`.

    * A numerical implementation of the Hartree term (using an improved version
      of Becke's Poisson solver) can be used instead of the ``RDirectTerm`` or
      ``UDirectTerm`` classes, which require four-center integrals. The relevant
      classes are
      :py:class:`~horton.meanfield.builtin.RBeckeHartree` or
      :py:class:`~horton.meanfield.builtin.UBeckeHartree`.

  Integration grids are discussed in more detail in the section
  :ref:`user_other_numerical_integration_specify`. A list of the supported LibXC functionals can
  be found in :ref:`ref_functionals`. Note that HORTON does not support
  meta-GGAs (MGGAs) yet.

Using these classes, you can construct the Hartree-Fock or a DFT effective
Hamiltonian.


A few typical examples
----------------------

The examples below assume that some or all of the following variables are
already defined:

* ``obasis``: An orbital basis set.
* ``olp``: two-index object with the overlap integrals.
* ``kin``: two-index object with the kinetic energy integrals.
* ``na``: two-index object with the nuclear attraction integrals.
* ``er``: four-index object with the electron repulsion integrals.
* ``grid``: a numerical integration grid.

If you are unfamiliar with any of these, please read the sections
:ref:`user_hamiltonian` and :ref:`user_other_numerical_integration`. The examples below
also make use of the external argument of
:py:class:`~horton.meanfield.hamiltonian.REffHam` or
:py:class:`~horton.meanfield.hamiltonian.UEffHam` to add the nuclear-nuclear
repulsion energy to the total energy reported by the effective Hamiltonian.

* Restricted Hartree-Fock:

  .. literalinclude:: ../data/examples/hf_dft/rhf_water_dense.py
      :lines: 37-45
      :caption: data/examples/hf_dft/rhf_water_dense.py, lines 37--45

* Unrestricted Hartree-Fock:

  .. literalinclude:: ../data/examples/hf_dft/uhf_methyl_dense.py
      :lines: 35-42
      :caption: data/examples/hf_dft/uhf_methyl_dense.py, lines 35--42

* Restricted Kohn-Sham DFT with the Dirac exchange and the VWN correlation
  functionals:

  .. literalinclude:: ../data/examples/hf_dft/rks_water_lda.py
      :lines: 37-47
      :caption: data/examples/hf_dft/rks_water_lda.py, lines 37--47

* Unrestricted Kohn-Sham DFT with the PBE GGA exchange and correlation
  functionals:

  .. literalinclude:: ../data/examples/hf_dft/uks_methyl_gga.py
      :lines: 38-48
      :caption: data/examples/hf_dft/uks_methyl_gga.py, lines 38--48

* Restricted Kohn-Sham DFT with the Hybrid GGA functional B3LYP:

  .. literalinclude:: ../data/examples/hf_dft/rks_water_hybgga.py
      :lines: 37-46
      :caption: data/examples/hf_dft/rks_water_hybgga.py, lines 37--46

* Unrestricted Kohn-Sham DFT with LDA exchange and correlation
  functionals and with a numerical integration of the Hartree term:

  .. literalinclude:: ../data/examples/hf_dft/uks_methyl_numlda.py
      :lines: 38-48
      :caption: data/examples/hf_dft/uks_methyl_numlda.py, lines 38--48


.. _user_hf_dft_occupation:

Choose an orbital occupation scheme
===================================

Before calling an SCF solver, you have to select a scheme to set the orbital
occupations after each SCF iteration, even when the occupation numbers are to
remain fixed throughout the calculation. You can use any of the following
three options:

* :py:class:`~horton.meanfield.occ.FixedOccModel`. Keep all occupation numbers
  fixed at preset values. For example,

  .. code-block:: python

      # Restricted case
      occ_model = FixedOccModel(np.array([1.0, 1.0, 0.5, 0.5, 0.0]))
      # Unrestricted case
      occ_model = FixedOccModel(np.array([1.0, 1.0, 0.5, 0.5, 0.0]), np.array([1.0, 0.7, 1.0, 0.0, 0.3]))


* :py:class:`~horton.meanfield.occ.AufbauOccModel`. Fill all orbitals according
  to the `Aufbau principle <http://en.wikipedia.org/wiki/Aufbau_principle>`_.
  For example,

  .. code-block:: python

      # Restricted case (three alpha and three beta electrons)
      occ_model = AufbauOccModel(3.0)
      # Unrestricted case (two alpha and three beta electrons)
      occ_model = AufbauOccModel(2.0, 3.0)


* :py:class:`~horton.meanfield.occ.FermiOccModel`. Use the Fermi-smearing method
  to fill up the orbitals. [rabuck1999]_ Only part of the methodology presented
  Rabuck is implemented. See :py:class:`~horton.meanfield.occ.FermiOccModel` for
  details. For example,

  .. code-block:: python

      # Restricted case (three alpha and three beta electrons, 300K)
      occ_model = FermiOccModel(3.0, temperature=300)
      # Unrestricted case (two alpha and three beta electrons, 500K)
      occ_model = AufbauOccModel(2.0, 3.0, temperature=500)


.. _user_hf_dft_scf:

Optimize the wavefunction with a self-consistent field algorithm
================================================================

HORTON supports the following SCF algorithms:

* :py:class:`~horton.meanfield.scf.PlainSCFSolver`: the ordinary SCF solver.
  This method just builds and diagonalizes the Fock matrices at every iteration.

* :py:class:`~horton.meanfield.scf_oda.ODASCFSolver`: the optimal damping SCF
  solver. [cances2001]_ It uses a cubic interpolation to estimate the optimal
  mixing between the old and the new density matrices. This is relatively robust
  but slow.

* :py:class:`~horton.meanfield.scf_cdiis.CDIISSCFSolver`: the (traditional)
  commutator direct inversion of the iterative subspace (CDIIS) algorithm, also
  know as Pulay mixing. [pulay1980]_ This is usually very efficient but
  sensitive to the initial guess.

* :py:class:`~horton.meanfield.scf_ediis.EDIISSCFSolver`: the energy direct
  inversion of the iterative subspace (EDIIS) method. [kudin2002]_ This method
  works well for the initial iterations but becomes numerically unstable close
  to the solution. It typically works better with a relativey poor initial
  guess.

* :py:class:`~horton.meanfield.scf_ediis2.EDIIS2SCFSolver`: a combination of
  CDIIS and EDIIS. [kudin2002]_ This method tries to combine the benefits of
  both approaches.

The plain SCF solver starts from an initial guess of the orbitals and updates
these in-place.

* Usage in the restricted case:

  .. literalinclude:: ../data/examples/hf_dft/rhf_water_dense.py
      :lines: 50-52
      :caption: data/examples/hf_dft/rhf_water_dense.py, lines 50--52

* Usage in the unrestricted case:

  .. literalinclude:: ../data/examples/hf_dft/uhf_methyl_dense.py
      :lines: 47-49
      :caption: data/examples/hf_dft/uhf_methyl_dense.py, lines 47--49

All other solvers start from an initial guess of the density matrix and update
that quantity in-place. The usage pattern is as follow:

* Usage in the restricted case:

  .. literalinclude:: ../data/examples/hf_dft/rks_water_lda.py
      :lines: 52-58
      :caption: data/examples/hf_dft/rks_water_lda.py, lines 52--58

* Usage in the unrestricted case:

  .. literalinclude:: ../data/examples/hf_dft/uks_methyl_lda.py
      :lines: 53-60
      :caption: data/examples/hf_dft/uks_methyl_lda.py, lines 53--60

All SCF solvers support the following two options:

  :threshold: The convergence threshold for the SCF solver. The exact meaning
              of this threshold is algorithm-dependent. In case of the
              ``PlainSCFSolver``, convergence is tested with the RSMD error on
              the Roothaan equations. In all other cases, the square root of
              the commutator of the density and the Fock matrix is tested.

  :maxiter: The maximum allowed number of SCF iterations. If convergence is not
            reached within this number of iterations, a ``NoSCFConvergence``
            error is raised.

Click on the links of the SCF solver classes above for more optional arguments
and their default values.


.. _user_hf_dft_get_orbitals:

Conversion of density and Fock matrix to orbitals
=================================================

The :py:class:`~horton.meanfield.scf.PlainSCFSolver` iteratively updates the
orbitals. Hence, for this SCF solver, there is no need to reconstruct the orbitals
after SCF convergence. For all other SCF algorithms, however, you need to iteratively
update the density matrix and explicitly reconstruct the orbitals after SCF
convergence.

.. note::

    SCF algorithms that update the density matrix may produce a converged
    density matrix with fractional occupations for degenerate orbitals at the
    Fermi level. The code snippets below properly handle such cases as well. For
    more details, refer to
    :py:meth:`~horton.matrix.dense.DenseExpansion.from_fock_and_dm`.

* Usage in the restricted case:

  .. literalinclude:: ../data/examples/hf_dft/rks_water_lda.py
      :lines: 60-66
      :caption: data/examples/hf_dft/rks_water_lda.py, lines 60--67

* Usage in the unrestricted case:

  .. literalinclude:: ../data/examples/hf_dft/uks_methyl_lda.py
      :lines: 62-70
      :caption: data/examples/hf_dft/uks_methyl_lda.py, lines 62--71


.. _user_hf_dft_to_file:

Writing SCF results to a file
=============================

Once you have obtained a converged set of molecular orbitals, they can be stored
in a file. Two file formats can be used: (i) HORTON's internal format based on
HDF5 and (ii) the Molden format. The internal format has the advantage that all
data is stored in full precision and that you can include as much additional
information as you like. The Molden format is useful for visualization purposes
or when you would like to use the orbitals in another quantum chemistry code.
Both formats can be used to store the wavefunction and to use it as an initial
guess afterwards. For more background information on data formats in HORTON,
refer to :ref:`ref_file_formats`.

The two code snippets below show how a wavefunction can be written to a file.
In both examples, the object ``mol`` is an instance of ``IOData`` that was
created earlier in the script by loading an ``.xyz`` file as follows:

.. code-block:: python

    mol = IOData.from_file('some_molecule.xyz')

When the ``.xyz`` file is loaded, the attributes ``numbers`` and ``coordinates``
of the object ``mol`` are already set as required for the ``.molden`` format.

* Usage in the restricted case:

  .. literalinclude:: ../data/examples/hf_dft/rks_water_lda.py
      :lines: 68-80
      :caption: data/examples/hf_dft/rks_water_lda.py, lines 68--80

* Usage in the unrestricted case:

  .. literalinclude:: ../data/examples/hf_dft/uks_methyl_lda.py
      :lines: 72-86
      :caption: data/examples/hf_dft/uks_methyl_lda.py, lines 72--86


.. _user_hf_dft_preparing_posthf:

Preparing for a Post-Hartree-Fock calculation
=============================================

Once the SCF has converged and you have obtained a set of orbitals, you can use
these orbitals to convert the integrals from the atomic-orbital (AO) basis to
the molecular-orbital (MO) basis. Two implementations are available in HORTON:
(i) using all molecular orbitals or (ii) by specifying a frozen core and active
set of orbitals. A full example, ``rhf_n2_dense.py``, which covers both options
and which includes dumping the transformed integrals to a file, is given in the
section :ref:`hf_dft_complete_examples` below.

The conversion to an MO basis is useful for post-HF calculations. For such
purposes, it is also of interest to sum all one-body operators into a single
term:

.. code-block:: python

    lf = DenseLinalgFactory(obasis.nbasis)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    one = kin.copy()
    one.iadd(na)


Transforming the Hamiltonian to the molecular-orbital (MO) basis
----------------------------------------------------------------

The function :py:func:`~horton.orbital_utils.transform_integrals` can be used for
this purpose. It works differently for restricted and unrestricted orbitals:

1. In the case of restricted (Hartree-Fock) orbitals, there is just one
   transformed one-body and one transformed two-body operator:

   .. code-block:: python

       (one_mo,), (two_mo,) = transform_integrals(one, er, 'tensordot', exp_alpha)


2. In the case of unrestricted (Hartree-Fock) orbitals, there are two
   transformed one-body and three transformed two-body operators:

   .. code-block:: python

       one_mo_ops, two_mo_ops = transform_integrals(one, er, 'tensordot', exp_alpha, exp_beta)
       one_mo_alpha, one_mo_beta = one_mo_ops
       two_mo_alpha_alpha, two_mo_alpha_beta, two_mo_beta_beta = two_mo_ops


Reducing the Hamiltonian to an active space
-------------------------------------------

The active-space feature only works for restricted orbitals. The new
one-electron integrals :math:`\tilde{t}_{pq}` become

.. math::

    \tilde{t}_{pq} = t_{pq} + \sum_{i \in \textrm{ncore}} ( 2 \langle pi \vert qi \rangle - \langle pi \vert iq \rangle),

where :math:`t_{pq}` is the element :math:`pq` of the old one-electron integrals
and :math:`\langle pi \vert qi \rangle` is the appropriate two-electron integral
in physicist's notation. The core energy of the active space is calculated as

.. math::

    e_\text{core} = e_\text{nn} + 2\sum_{i \in \textrm{ncore}} t_{ii} + \sum_{i, j \in \textrm{ncore}} (2 \langle ij \vert ij \rangle - \langle ij \vert ji \rangle)

where the two-electron integrals :math:`\langle pq \vert rs \rangle` contain only
the elements with active orbital indices :math:`p,q,r,s`. This type of conversion
is implemented in the function :py:func:`~horton.orbital_utils.split_core_active`.
It is used as follows:

.. code-block:: python

     one_small, two_small, core_energy = split_core_active(one, er,
        external['nn'], exp_alpha, ncore, nactive)


.. _hf_dft_complete_examples:

Complete examples
=================

The following is a basic example of a restricted Hartree-Fock calculation of
water. It contains all the steps discussed in the previous sections.

.. literalinclude:: ../data/examples/hf_dft/rhf_n2_dense.py
    :caption: data/examples/hf_dft/rhf_n2_dense.py
    :lines: 7-

The directory ``data/examples/hf_dft`` contains many more examples that use the
different options discussed above. The following table shows which features are
used in which example.

.. include:: hf_dft_examples.rst.inc

A more elaborate example can be found in ``data/examples/hf_compare``. It
contains a script that systematically computes all elements in the periodic table
(for different charges and multiplicities), and compares the results with
outputs obtained with Gaussian. See the ``README`` file for instructions how to
run this example.
