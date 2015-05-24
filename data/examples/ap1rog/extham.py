#!/usr/bin/env python

from horton import *
###############################################################################
## Define number of occupied orbitals and total number of basis functions #####
###############################################################################
nocc = 5
nbasis = 28
###############################################################################
## Define Occupation model, expansion coefficients and overlap ################
###############################################################################
lf = DenseLinalgFactory(nbasis)
occ_model = AufbauOccModel(nocc)
orb = lf.create_expansion(nbasis)
olp = lf.create_two_index(nbasis)
olp.assign_diagonal(1.0)
orb.assign(olp)
###############################################################################
## Read Hamiltonian from file 'FCIDUMP' #######################################
###############################################################################
one, two, core = load_fcidump(lf, './FCIDUMP')

###############################################################################
## Do OO-AP1roG optimization ##################################################
###############################################################################
ap1rog = RAp1rog(lf, occ_model)
energy, c, l = ap1rog(one, two, core, orb, olp, True)
