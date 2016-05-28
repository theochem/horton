#!/usr/bin/env python

from horton import *
from horton.test.common import numpy_seed

###############################################################################
## Set up molecule, define basis set ##########################################
###############################################################################
# get the XYZ file from HORTON's test data directory
fn_xyz = context.get_fn('test/h2.xyz')
mol = IOData.from_file(fn_xyz)
obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
###############################################################################
## Define Occupation model, expansion coefficients and overlap ################
###############################################################################
lf = DenseLinalgFactory(obasis.nbasis)
occ_model = AufbauOccModel(1)
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
## Combine one-electron integrals to single Hamiltonian #######################
###############################################################################
one = kin.copy()
one.iadd(na)

###############################################################################
## Do OO-AP1roG optimization ##################################################
###############################################################################
ap1rog = RAp1rog(lf, occ_model)
with numpy_seed():  # reproducible 'random' numbers to make sure it always works
    energy, g, l = ap1rog(one, er, external['nn'], orb, olp, True)

###############################################################################
## Do PTb calculation #########################################################
###############################################################################
ptb = PTb(lf, occ_model)
with numpy_seed():  # reproducible 'random' numbers to make sure it always works
    energyptb, amplitudes = ptb(one, er, orb, g, **{'eref': energy, 'ecore': external['nn'], 'threshold': 1e-6})
