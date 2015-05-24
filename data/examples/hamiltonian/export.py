#!/usr/bin/env python

from horton import *
###########################################################################################
## Set up molecule, define basis set ######################################################
###########################################################################################
# get the XYZ file from Horton's test data directory
fn_xyz = context.get_fn('test/water.xyz')
mol = Molecule.from_file(fn_xyz)
obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
###########################################################################################
## Define Occupation model, expansion coefficients and overlap ############################
###########################################################################################
lf = DenseLinalgFactory(obasis.nbasis)
occ_model = AufbauOccModel(7)
orb = lf.create_expansion(obasis.nbasis)
olp = obasis.compute_overlap(lf)
###########################################################################################
## Construct Hamiltonian ##################################################################
###########################################################################################
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
###########################################################################################
## Perform initial guess ##################################################################
###########################################################################################
guess_core_hamiltonian(olp, kin, na, orb)
###########################################################################################
## Do a Hartree-Fock calculation ##########################################################
###########################################################################################
scf_solver = PlainSCFSolver(1e-6)
scf_solver(ham, lf, olp, occ_model, orb)
###########################################################################################
## Combine to single one-electron Hamiltonian #############################################
###########################################################################################
one = kin.copy()
one.iadd(na)

###########################################################################################
## Export Hamiltonian in Hartree-Fock molecular orbital basis (all orbitals active) #######
###########################################################################################
dump_fcidump(lf, one, er, external['nn'], orb, 'FCIDUMP')

###########################################################################################
## Export Hamiltonian in Hartree-Fock molecular orbital basis for CAS(8,8) ################
###########################################################################################
dump_fcidump(lf, one, er, external['nn'], orb, 'FCIDUMP8-8',
             **{'nel': 8, 'ncore': 2, 'nactive': 8})
