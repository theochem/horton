#!/usr/bin/env python

from horton import *

###############################################################################
## Read Hamiltonian from file 'FCIDUMP' #######################################
###############################################################################
mol = Molecule.from_file('./FCIDUMP')

###############################################################################
## Define number of occupied orbitals and total number of basis functions #####
###############################################################################
nocc = 5

###############################################################################
## Define Occupation model, expansion coefficients and overlap ################
###############################################################################
occ_model = AufbauOccModel(nocc)
orb = mol.lf.create_expansion()
olp = mol.lf.create_two_index()
olp.assign_diagonal(1.0)
orb.assign(olp)

###############################################################################
## Do OO-AP1roG optimization ##################################################
###############################################################################
ap1rog = RAp1rog(mol.lf, occ_model)
energy, c, l = ap1rog(mol.one_mo, mol.two_mo, mol.core_energy, orb, olp, True)
