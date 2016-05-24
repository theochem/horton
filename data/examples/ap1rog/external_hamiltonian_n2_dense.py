#!/usr/bin/env python

from horton import *
from horton.test.common import numpy_seed

# Read Hamiltonian from file 'FCIDUMP'
# ------------------------------------
# The required FCIDUMP file can be generated with the script
# data/examples/hf_dft/rhf_n2_dense.py
mol = IOData.from_file('n2.FCIDUMP')

# Define Occupation model, expansion coefficients and overlap
# -----------------------------------------------------------
nocc = 5
occ_model = AufbauOccModel(nocc)
orb = mol.lf.create_expansion()
olp = mol.lf.create_two_index()
olp.assign_diagonal(1.0)
orb.assign(olp)

# Do OO-AP1roG optimization
# -------------------------
ap1rog = RAp1rog(mol.lf, occ_model)
with numpy_seed():  # reproducible 'random' numbers to make sure it always works
    energy, c, l = ap1rog(mol.one_mo, mol.two_mo, mol.core_energy, orb, olp, True)
