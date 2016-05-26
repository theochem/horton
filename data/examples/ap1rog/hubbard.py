#!/usr/bin/env python

import numpy as np

from horton import *
from horton.test.common import numpy_seed

###############################################################################
## Define Occupation model, expansion coefficients and overlap ################
###############################################################################
lf = DenseLinalgFactory(6)
occ_model = AufbauOccModel(3)
modelham = Hubbard(pbc=True)
orb = lf.create_expansion(6)
olp = modelham.compute_overlap(lf)
###############################################################################
# t-param, t = -1 #############################################################
###############################################################################
kin = modelham.compute_kinetic(lf, -1)
###############################################################################
# U-param, U = 2 ##############################################################
###############################################################################
two = modelham.compute_er(lf, 2)
###############################################################################
## Perform initial guess ######################################################
###############################################################################
guess_core_hamiltonian(olp, kin, orb)
terms = [
    RTwoIndexTerm(kin, 'kin'),
    RDirectTerm(two, 'hartree'),
    RExchangeTerm(two, 'x_hf'),
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
with numpy_seed():  # reproducible 'random' numbers to make sure it always works
    energy, c, l = ap1rog(kin, two, 0, orb, olp, True)
