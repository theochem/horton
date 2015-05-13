.. _user_hf_dft:

How to use Horton as a Hartree-Fock/DFT program
###############################################

**Disclaimer.** The relevant Hartree-Fock and DFT features in Horton are
explained below. This is all very different from a regular quantum chemistry
code where one just prepares an input file with a molecular geometry and some
control options. Instead, Horton is designed to allow one to interfer with most
of the internals, and, in that way, make intersting fundamental research easy
and enjoyable.

Horton can do restricted and unrestricted HF/DFT SCF calculations in various
different ways. The following sections cover the typical steps of a calculation:

1) Definition of the (molecular electronic) Hamiltonian
2) An initial guess of the orbitals
3) A definition of the effective Hamiltonian
4) Choose an orbital occupations scheme
5) Run an SCF algorithm

The last section contains an overview of complete examples
that are shipped (and tested) with Horton. These may be convenient as a starting
point for the preparation of your own scripts. One may also just dive straight
into these examples and consult to first five sections to better understand how
each example works.

The code snippets below use the following conventions for variable names. The
following are defined by setting up the Hamiltonian (see next section):

* ``olp``: two-index object with the overlap integrals.
* ``kin``: two-index object with the kinetic energy integrals.
* ``na``: two-index object with the nuclear attraction integrals.
* ``er``: four-index object with the electron repulsion integrals.

The following names are also used systematically:

* ``exp_alpha``, ``exp_beta``: the expansion coefficient of the alpha and beta
  orbitals in a (local) basis.


Setting up the (molecular electronic) Hamiltonian
=================================================

See :ref:`user_hamiltonian`.

One first has to compute/load the two- and four-index objects for the one- and
two-body terms in the Hamiltonian. Some DFT implementations in Horton do not
require pre-computed four-center integrals.


Generating an initial guess of the wavefunction
===============================================

.. _user_hf_dft_core_guess:

Core Hamiltonian Guess
----------------------

This guess involves only the one-body terms of the electronic Hamiltonian and
completely neglects electron-electron interactions. The orbitals for such a
one-body Hamiltonian can be computed without any prior guess.

The function :py:func:`horton.meanfield.guess.guess_core_hamiltonian` can be used
to compute core hamiltonian guess. The guess does not depend on occupation
numbers and in the case of a unrestricted calculation, the default behavior is
to assign the same orbitals are assigned for the alpha and beta spins.

The guess for a restricted wavefunction is done as follows:

.. code-block:: python

    # Create alpha orbitals
    exp_alpha = lf.create_expansion()

    # Initial guess
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

For a unrestricted wavefunction, the procedure is very similar:

.. code-block:: python

    # Create alpha and beta orbitals
    exp_alpha = lf.create_expansion()
    exp_beta = lf.create_expansion()

    # Initial guess
    guess_core_hamiltonian(olp, kin, na, exp_alpha, exp_beta)

The arguments ``exp_alpha`` and ``exp_beta`` are treated as output arguments.
Instead of ``kin`` and ``na``, one may provide just any set of one-body
operators to construct different types of guesses.


Randomizing an initial guess
----------------------------

The method :py:meth:`horton.matrix.dense.DenseExpansion.rotate_random` can be used to
apply a random unitary rotation all orbitals (occupied and virtual). The orbital
energies and occupation numbers are not affected. For example, after
constructing a :ref:`user_hf_dft_core_guess`, the following line randomizes the
orbitals:

.. code-block:: python

    # randomly rotate the orbitals (irrespective of occupied or virtual)
    exp_alpha.rotate_random()


Modifying the initial guess
---------------------------

If needed one may fine-tune the initial guess by making fine-grained
modifications to the orbitals. (These may also be useful for fixing the orbitals
that come out of a failed SCF.)

* The method :py:meth:`horton.matrix.dense.DenseExpansion.rotate_2orbitals`
  allows one to mix two orbitals. By default it rotates the HOMO and LUMO
  orbitals by 45 degrees:

  .. code-block:: python

      # Mix HOMO and LUMO orbitals
      exp_alpha.rotate_2orbitals()

      # Rotate 1st and 6th orbital by 30 deg
      exp._alpha.rotate_2orbitals(np.pi/6, 0, 5)

* The method :py:meth:`horton.matrix.dense.DenseExpansion.swap_orbitals` allows
  on to swap several orbitals. It takes as an argument an array where each row
  is a pair of orbitals to swap. For example, the following swaps 1st and 3rd,
  followed by a swap of 2nd and 4th:

  .. code-block:: python

      # Swap some orbitals
      swaps = np.array([[0, 2], [1, 3]])
      exp_alpha.swap_orbitals(swaps)


Reading a guess from a file
----------------------------

One may also load orbitals from an external file. The file formats ``.mkl``,
``.molden``, ``.fchk``, or Horton's internal ``.h5`` can are all eligible
sources of orbitals. For example, the orbitals from a Gaussian formatted
checkpoint file may be loaded as follows:

.. code-block:: python

    # Load fchk file
    mol = Molecule.from_file('water.fchk')

    # Print the number of alpha orbitals (occupied and virtual)
    print mol.exp_aplha.nfn

Obviously, if one would like to use these orbitals without projecting them onto
a new basis set (as explained in :ref:`user_hf_dft_project_basis`), one is
forced to continue working in exactly the same basis set, which can be accessed
in this example as ``mol.obasis``. See :ref:`user_molecularham_geom_and_basis`
for more details.


.. _user_hf_dft_project_basis:

Projecting orbitals from a smaller basis onto a larger one
----------------------------------------------------------

Assuming one has obtained (converged) orbitals in a smaller basis set, one can
try to use these as initial guess after projecting the orbitals onto the
larger basis set. This is exactly what the function
:py:func:`horton.meanfield.project.project_orbitals_mgs` does. The following
snippet assumes that the ``obasis0`` and ``exp_alpha0`` are the small basis set
and a set of orbitals in that basis for the ``Molecule`` instance ``mol``.

.. code-block:: python

    # Definition of the bigger basis set
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'aug-cc-pvtz'):

    # Linalg factory for the bigger basis set
    lf1 = DenseLinalgFactory(obasis1.nbasis)

    # Create a expansion object for the alpha orbitals in the large basis
    exp_alpha1 = lf1.create_expansion()

    # The actual projection
    project_orbitals_msg(obasis0, obasis1, exp_alpha0, exp_alpha1)


Effective Hamiltonians
======================

Horton implements spin-restricted and spin-unrestricted effective Hamiltonians.
Mathematically, these are models for the energy as function of a set of
density matrices. The implementation also provides an API to compute for every
density matrix the corresponding Fock matrix, i.e. the derivative of the energy
toward the density matrix elements.

* For the restricted case, the alpha and beta density matrices are assumed
  to be identical. Hence the energy is only a function of the alpha density
  matrix. When constructing the Fock matrix, the derivative is divided by two
  to obtain such that the Fock matrix has conventional orbital energies as
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

This generic API is implemented in the class
:py:class:`horton.meanfield.hamiltonian.REffHam` and
:py:class:`horton.meanfield.hamiltonian.UEffHam`. The prefixes ``R`` and ``U``
are used (also below) to differentiate between restricted and unrestricted
implementations. A Hatree-Fock or DFT effective Hamiltonian is defined by
constructing an instance of the ``REffHam`` or ``UEffHam`` classes and providing
the necessary energy terms to the constructor.


Supported energy terms for the effective Hamiltonians
-----------------------------------------------------

All classes below take a ``label`` argument to give each term in the Hamiltonian
a name, e.g. used for storing/displaying results. For each class listed below,
follow the hyperlinks to he corresponding documentation for a description of the
constructor arguments.

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
  in efficiency when building the Fock matrix/matrices. The constructor of a
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
  :ref:`user_integration_grids_specify`. A list of the supported LibXC functionals can
  be found in :ref:`ref_functionals`. Note that Horton does not support the
  MGGA's yet.

Using these classes, one can construct the Hatree-Fock or a DFT effective
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
:ref:`user_hamiltonian` and :ref:`user_integration_grids`. The examples below
also make use of the external argument of
:py:class:`~horton.meanfield.hamiltonian.REffHam` or
:py:class:`~horton.meanfield.hamiltonian.UEffHam` to add the nuclear-nuclear
repulsion energy to the total energy reported by the effective Hamiltonian.

* Restricted Hartree-Fock:

  .. code-block:: python

      external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
      terms = [
          RTwoIndexTerm(kin, 'kin'),
          RDirectTerm(er, 'hartree'),
          RExchangeTerm(er, 'x_hf'),
          RTwoIndexTerm(na, 'ne'),
      ]
      ham = REffHam(terms, external)


* Unrestricted Hartree-Fock:

  .. code-block:: python

      external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
      terms = [
          UTwoIndexTerm(kin, 'kin'),
          UDirectTerm(er, 'hartree'),
          UExchangeTerm(er, 'x_hf'),
          UTwoIndexTerm(na, 'ne'),
      ]
      ham = UEffHam(terms, external)


* Restricted Kohn-Sham DFT with the Dirac exchange functional:

  .. code-block:: python

      external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
      terms = [
          RTwoIndexTerm(kin, 'kin'),
          RDirectTerm(er, 'hartree'),
          RGridGroup(obasis, grid, [
              RLibXCLDA('x'),
          ]),
          RTwoIndexTerm(na, 'ne'),
      ]
      ham = REffHam(terms, external)


* Unrestricted Kohn-Sham DFT with the PBE GGA exchange and correlation
  functionals:

  .. code-block:: python

      external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
      terms = [
          UTwoIndexTerm(kin, 'kin'),
          UDirectTerm(er, 'hartree'),
          UGridGroup(obasis, grid, [
              ULibXCGGA('x_pbe'),
              ULibXCGGA('c_pbe'),
          ]),
          UTwoIndexTerm(na, 'ne'),
      ]
      ham = UEffHam(terms, external)


* Restricted Kohn-Sham DFT with the Hybrid GGA functional B3LYP:

  .. code-block:: python

      external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
      libxc_term = RLibXCHybridGGA('xc_b3lyp')
      terms = [
          RTwoIndexTerm(kin, 'kin'),
          RDirectTerm(er, 'hartree'),
          RGridGroup(obasis, grid, [libxc_term]),
          RExchangeTerm(er, 'x_hf', libxc_term.get_exx_fraction()),
          RTwoIndexTerm(na, 'ne'),
      ]
      ham = REffHam(terms, external)


* Unrestricted Kohn-Sham DFT with the PBE GGA exchange and correlation
  functionals and with a numerical integration of the Hartree term:

  .. code-block:: python

      external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
      terms = [
          UTwoIndexTerm(kin, 'kin'),
          UGridGroup(obasis, grid, [
              UBeckeHartree(lmax=5, 'hartree_becke'),
              ULibXCGGA('x_pbe'),
              ULibXCGGA('c_pbe'),
          ]),
          UTwoIndexTerm(na, 'ne'),
      ]
      ham = UEffHam(terms, external)


Models for orbital occupations
==============================

Before calling an SCF solver, one has to select a scheme to set the orbital
occupations after each SCF iteration, even when the occuption numbers are to
remain fixed throughout the calculation. One has can use any of the following
three options:

* :py:class:`~horton.meanfield.occ.FixedOccModel`. Keep all occupation numbers
  fixed at preset values. Example usage:

  .. code-block:: python

      # Restricted case
      occ_model = FixedOccModel(np.array([1.0, 1.0, 0.5, 0.5, 0.0]))
      # Unrestricted case
      occ_model = FixedOccModel(np.array([1.0, 1.0, 0.5, 0.5, 0.0]), np.array([1.0, 0.7, 1.0, 0.0, 0.3]))


* :py:class:`~horton.meanfield.occ.AufbauOccModel`. Fill all orbitals according
  to the `Aufbau principle <http://en.wikipedia.org/wiki/Aufbau_principle>`_.
  Example usage:

  .. code-block:: python

      # Restricted case (three alpha and three beta electrons)
      occ_model = AufbauOccModel(3.0)
      # Unrestricted case (two alpha and three beta electrons)
      occ_model = AufbauOccModel(2.0, 3.0)


* :py:class:`~horton.meanfield.occ.FermiOccModel`. Use the Fermi-smearing method
  to fill up the orbitals. [rabuck1999]_
  Example usage:

  .. code-block:: python

      # Restricted case (three alpha and three beta electrons, 300K)
      occ_model = FermiOccModel(3.0, temperature=300)
      # Unrestricted case (two alpha and three beta electrons, 500K)
      occ_model = AufbauOccModel(2.0, 3.0, temperature=500)


Self-consistent field algorithms
================================

Horton supports the following SCF algorithms:

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
  works well for the initial iterations but becomes numerically instable close
  to the solution. It typically works better with a relativey poor initial
  guess.

* :py:class:`~horton.meanfield.scf_ediis2.EDIIS2SCFSolver`: a combination of
  CDIIS and EDIIS. [kudin2002]_ This method tries to combine the benefits of
  both approaches.

TODO: add code snipts of how to call the SCF solver


Complete examples
=================

The following is a basic example of a restricted Hartree-Fock calculation of
water. It contains all the steps discussed in the previous sections.

.. literalinclude:: ../data/examples/hf_dft/rhf_water_dense.py
    :caption: data/examples/hf_dft/rhf_water_dense.py
    :lines: 7-

The directory ``data/examples/hf_dft`` contains many more examples that use the
different options discussed above. The following table shows which features are
used in which example.

.. include:: hf_dft_examples.rst.inc

A more elaborate example can be found in ``data/examples/hf_compare``. It
contains a script that systemtically computes all elements in the periodic table
(for different charges and multiplicities), and compares the results with
outputs obtained with Gaussian. See the ``README`` file for instructions how to
run this example.
