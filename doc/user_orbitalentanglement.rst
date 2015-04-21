Orbital entanglement analysis
#############################

.. contents::

.. _orbitalentanglement:

Orbital entanglement and orbital correlation
============================================


Supported features
------------------

.. _orbitalentanglementseniorityzero:

Orbital entanglement and correlation for a seniority zero wavefunction
======================================================================


How to set-up a calculation
---------------------------


How to generate correlation diagrams
------------------------------------


Example input files
===================

Orbital entanglement analysis of an AP1roG wavefunction
-------------------------------------------------------

This is a basic example on how to perform an orbital entanglement analysis in Horton. This script performs an orbital entanglement analysis of the AP1roG wavefunction for the water molecule using the cc-pVDZ basis set.

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
    ## Combine one-electron integrals to single Hamiltonian #######################
    ###############################################################################
    one = kin.copy()
    one.iadd(na)

    ###############################################################################
    ## Do OO-AP1roG optimization ##################################################
    ###############################################################################
    ap1rog = RAp1rog(lf, occ_model)
    energy, g, l = ap1rog(one, er, external['nn'], moceoff, olp, True)

    ###############################################################################
    ## Do orbital entanglement analysis ###########################################
    ###############################################################################
    one_dm = lf.create_one_index()
    one_dm.assign(exp_alpha.occupations)
    twoppqq = lf.create_two_index()
    twopqpq = lf.create_two_index()
    twoppqq.compute_2dm_ap1rog(one_dm, g, l, 'ppqq')
    twopqpq.compute_2dm_ap1rog(one_dm, g, l, 'pqpq')

    entanglement = OrbitalEntanglementAp1rog(lf, one_dm, [twoppqq,twopqpq])
    entanglement()
