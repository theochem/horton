#!/usr/bin/env python

from horton import *
###############################################################################
## Set up molecule, define basis set ##########################################
###############################################################################
# Use the XYZ file from Horton's test data directory.
fn_xyz = context.get_fn('test/water.xyz')
mol = IOData.from_file(fn_xyz)
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
