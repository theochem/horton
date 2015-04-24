.. _localization:

The localization of molecular orbitals
######################################

.. contents::


In general, the localization algorithm optimizes the localization functions using an orthogonal transformation between the orbitals :math:`\vert i \rangle` ,

.. math::

    \vert \tilde{i} \rangle = \sum_k \vert k \rangle \exp(-\mathbf{\kappa})_{ki},

where :math:`\mathbf{\kappa}` is the generator of orbital rotations

.. math::

    \mathbf{\kappa} = \sum_{k > l} \kappa_{kl} (a^\dagger_k a_l - a^\dagger_l a_k)

This version of Horton supports

1. :ref:`Pipek-Mezey <pipek-mezey>` localization


.. _pipek-mezey:

Pipek-Mezey localization of restricted Hartree-Fock orbitals
============================================================

In the Pipek-Mezey scheme, the localization function is maximized to obtain a set of orbitals that are on average most local. To evaluate the Pipek-Mezey localization function,

.. math::

    D = \sum_{i} \sum_{A \in \textrm{atoms}} (Q_{ii}^A)^2,

an atomic population matrix :math:`\mathbf{Q}^A` is needed which can be calculated using projectors for atomic basis functions.

If the Mulliken population analysis is used, the Mulliken projectors can be constructed as follow

.. code-block:: python

    mulliken = get_mulliken_operators(obasis, lf)

where

    :obasis: (A ``GOBasis`` instance) the Gaussian orbital basis
    :lf: A ``LinalgFactory`` instance. One of ``DenseLinalgFactory`` or ``CholeskyLinalgFactory``

After defining the projectors, an instance of the ``PipekMezey`` class can be created,

.. code-block:: python

    loc = PipekMezey(lf, occ_model, projector)

where **lf** is again the ``LinalgFactory`` instance and

    :occ_model: (``AufbauOccModel`` instance) an Aufbau occupation model
    :projector: (list of ``TwoIndex`` instances) the projectors for atomic basis functions

To localize the orbitals, use a function call of the ``PipekMezey`` instance

.. code-block:: python

    loc(mocoeff, select[, **kwargs])

with arguments

    :mocoeff: (``Expansion`` instance) the MO coefficient matrix to be localized
    :select: (str) the orbital block to be localized. One of ``occ`` (occupied orbitals), ``virt`` (virtual orbitals)

The **select** argument specifies the orbital block to be localized (either the occupied or the virtual Hartree-Fock orbitals). If both the occupied and virtual Hartree-Fock orbitals are to be localized, two consecutive function calls are required. Note that each function call optimizes the orbitals through orbital rotations within the orbital block of interest (either ``occ`` or ``virt``).

The keyword arguments contain optimization-specific parameters. Specifying keyword arguments is optional. The list of acceptable keyword arguments is as follows

    :maxiter: (int) maximum number of iterations (that is, orbital rotation steps) for localization (default ``2000``)

    :threshold: (float)  localization threshold for objective function (default ``1e-6``)

    :levelshift: (float) level shift of Hessian (default ``1e-8``). Absolute value of elements of the orbital Hessian smaller than **levelshift** are shifted by **levelshift**

    :stepsearch: (dictionary) optimizes an orbital rotation step:

              :method: (str) step search method used. One of ``trust-region`` (default), ``None``,  ``backtracking``
              :optimizer: (str) optimizes step to boundary of trust radius in ``trust-region``. One of ``pcg`` (preconditioned conjugate gradient), ``dogleg`` (Powell's single dogleg step), ``ddl`` (Powell's double-dogleg step) (default ``ddl``)
              :stepa: (float) scaling factor for Newton step. Used in ``backtracking`` and ``None`` method (default ``0.75``)
              :c1: (float) parameter used in the Armijo condition of ``backtracking`` (default ``1e-4``)
              :maxstep: (float) maximum step length/trust radius (default ``0.75``)
              :minstep: (float) minimum step length used in ``backracking`` (default ``1e-6``). If step length falls below **minstep**, the ``backtracking`` line search is terminated and the most recent step is accepted
              :maxiterouter: (int) maximum number of iterations to optimize orbital rotation step  (default ``10``)
              :maxiterinner: (int) maximum number of optimization steps in each step search (used only in ``pcg``, default ``500``)
              :maxeta: (float) upper bound for estimated vs. actual change in ``trust-region`` (default ``0.75``)
              :mineta: (float) lower bound for estimated vs. actual change in ``trust-region`` (default ``0.25``)
              :upscale: (float) scaling factor to increase trust radius in ``trust-region`` (default ``2.0``)
              :downscale: (float) scaling factor to decrease trust radius in ``trust-region`` (default ``0.25``)
              :trustradius: (float) initial trust radius (default ``0.75``)
              :maxtrustradius: (float) maximum trust radius (default ``0.75``)
              :threshold: (float) trust-region optimization threshold, only used in ``pcg`` (default ``1e-8``)

The optimized set of orbitals is stored in **mocoeff** (an ``Expansion`` instance). Note that the initial orbitals **mocoeff** are overwritten.


Example input files
===================

Pipek-Mezey localization of restricted Hartree-Fock orbitals for the water molecule
-----------------------------------------------------------------------------------

This is a basic example on how to perform a Pipek-Mezey localization in Horton. This script performs a Pipek-Mezey localization for the water molecule using the cc-pVDZ basis set and Mulliken projectors.

.. code-block:: python

    from horton import *
    ###############################################################################
    ## Set up molecule, define basis set ##########################################
    ###############################################################################
    mol = Molecule.from_file('mol.xyz')
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
    ###############################################################################
    ## Define Occupation model, expansion coefficients and overlap ################
    ###############################################################################
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(5)
    moceoff = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    ###############################################################################
    ## Construct Hamiltonian ######################################################
    ###############################################################################
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    er = obasis.compute_electron_repulsion(lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    ###############################################################################
    ## Perform initial guess ######################################################
    ###############################################################################
    guess_core_hamiltonian(olp, kin, na, moceoff)
    ###############################################################################
    ## Do a Hartree-Fock calculation ##############################################
    ###############################################################################
    scf_solver = PlainSCFSolver(1e-6)
    scf_solver(ham, lf, olp, occ_model, moceoff)
    ###############################################################################
    ## Define Mulliken projectors #################################################
    ###############################################################################
    mulliken = get_mulliken_operators(obasis, lf)

    ###############################################################################
    ## Pipek-Mezey localizaton ####################################################
    ###############################################################################
    loc = PipekMezey(lf, occ_model, mulliken)
    ###############################################################################
    ## occupied block #############################################################
    ###############################################################################
    loc(mocoeff, 'occ')
    ###############################################################################
    ## virtual block ##############################################################
    ###############################################################################
    loc(mocoeff, 'virt')
