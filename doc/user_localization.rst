.. _localization:

Localization of molecular orbitals
##################################

In general, the localization algorithm optimizes some localization function by
orthogonal transformation of the orbitals. Given orbitals
:math:`\vert i \rangle`, the localized orbital :math:`\vert \tilde{i} \rangle`
can be obtained by some transformation

.. math::

    \vert \tilde{i} \rangle = \sum_k \vert k \rangle \exp(-\mathbf{\kappa})_{ki},

where

.. math::

    \mathbf{\kappa} = \sum_{k > l} \kappa_{kl} (a^\dagger_k a_l - a^\dagger_l a_k)

and :math:`\kappa_{kl}` is determined by the optimization of the localization function.

Many of the localization schemes, and thus
the result of the localization, differ by the localization function. The localization
function somehow measures the localization of the orbitals. So far, Horton only
supports the Pipek-Mezey localization. [pipek1989]_


.. _pipek-mezey:

Pipek-Mezey localization
========================

In the Pipek-Mezey scheme, the Pipek-Mezey localization function, :math:`D`, is maximized.

.. math::

    D = \sum_{i} \sum_{A \in \textrm{atoms}} (Q_{ii}^A)^2,

where :math:`\mathbf{Q}^A` is the atomic population matrix. The atomic population
matrix can be obtained from the overlap of the atomic basis, the expansion of the molecular
orbitals from the atomic basis, the occupation of each molecular orbital, and
some weighted projection of atomic basis function within each atom.

For example, if the Mulliken population analysis is used, the projectors are
obtained through :py:meth:`horton.part.mulliken.get_mulliken_operators`.
Then the Pipek-Mezey localization function and the optimization are obtained through
:py:class:`horton.localization.localization.PipekMezey`. Function call,
:py:meth:`horton.localization.localization.Localization.__call__`, of this
instance results in localization.
Please see documentation in :py:mod:`horton.localization.localization` for more detail.


Example input files
===================

Pipek-Mezey localization of restricted Hartree-Fock orbitals for the water molecule
-----------------------------------------------------------------------------------

This is a basic example on how to perform a Pipek-Mezey localization in Horton. This script performs a Pipek-Mezey localization for the water molecule using the cc-pVDZ basis set and Mulliken projectors.

.. code-block:: python

    #!/usr/bin/env python

    from horton import *

    ###############################################################################
    ## Set up molecule, define basis set ##########################################
    ###############################################################################
    mol = Molecule.from_file(context.get_fn('test/water.xyz'))
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
    ###############################################################################
    ## Define Occupation model, expansion coefficients and overlap ################
    ###############################################################################
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(5)
    orb = lf.create_expansion(obasis.nbasis)
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
    guess_core_hamiltonian(olp, kin, na, orb)
    ###############################################################################
    ## Do a Hartree-Fock calculation ##############################################
    ###############################################################################
    scf_solver = PlainSCFSolver(1e-6)
    scf_solver(ham, lf, olp, occ_model, orb)
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
    loc(orb, 'occ')
    ###############################################################################
    ## virtual block ##############################################################
    ###############################################################################
    loc(orb, 'virt')
