#! /usr/bin/python
from horton import *

###############################################################################
## Define Occupation model, expansion coefficients and overlap ################
###############################################################################
lf = DenseLinalgFactory(10)
occ_model = AufbauOccModel(5)
modelham = Hubbard(pbc=True)
orb = lf.create_expansion(10)
olp = modelham.compute_overlap(lf)
###############################################################################
# t-param, t = -1 #############################################################
###############################################################################
kin = modelham.compute_kinetic(lf, -1)
###############################################################################
# U-param, U = 2 ##############################################################
###############################################################################
er = modelham.compute_er(lf, 2)
###############################################################################
## Perform initial guess ######################################################
###############################################################################
guess_core_hamiltonian(olp, kin, orb)
terms = [
    RTwoIndexTerm(kin, 'kin'),
    RDirectTerm(er, 'hartree'),
    RExchangeTerm(er, 'x_hf'),
]
ham = REffHam(terms)
###############################################################################
## Do a Hartree-Fock calculation ##############################################
###############################################################################
scf_solver = PlainSCFSolver()
scf_solver(ham, lf, olp, occ_model, orb)
###############################################################################
## Do OO-AP1roG optimization ##################################################
###############################################################################
ap1rog = RAp1rog(lf, occ_model)
energy, g, l = ap1rog(kin, er, 0, orb, olp, True)
