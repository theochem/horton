.. _localization:

The localization of molecular orbitals
######################################

.. contents::

.. _pipek-mezey:

Pipek-Mezey localization of restricted Hartree-Fock orbitals
============================================================

How to set-up a calculation
---------------------------


Example input files
===================

Pipek-Mezey localization of restricted Hartree-Fock orbitals for the water molecule
-----------------------------------------------------------------------------------

This is a basic example on how to perform a Pipek-Mezey localization in Horton. This script performs a Pipek-Mezey localization for the water molecule using the cc-pVDZ basis set.

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
    ## Pipek-Mezey localizaton ####################################################
    ###############################################################################
    mulliken = get_mulliken_operators(obasis, lf)

    loc = PipekMezey(lf, occ_model, mulliken)
    ###############################################################################
    ## occupied block #############################################################
    ###############################################################################
    loc(mocoeff, 'occ')
    ###############################################################################
    ## virtual block ##############################################################
    ###############################################################################
    loc(mocoeff, 'virt')
