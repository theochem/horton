#!/usr/bin/env python

from horton import *
from horton.test.common import numpy_seed

# Set up molecule, define basis set
# ---------------------------------
# get the XYZ file from HORTON's test data directory
fn_xyz = context.get_fn('test/water.xyz')
mol = IOData.from_file(fn_xyz)
# Slighly rescale the coordinates.
mol.coordinates *= 1.01
obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')

# Define Occupation model, expansion coefficients and overlap
# -----------------------------------------------------------
lf = CholeskyLinalgFactory(obasis.nbasis)
occ_model = AufbauOccModel(5)
exp_alpha = lf.create_expansion(obasis.nbasis)
olp = obasis.compute_overlap(lf)

# Construct Hamiltonian
# ---------------------
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
one = kin.copy()
one.iadd(na)
two = obasis.compute_electron_repulsion(lf)
core_energy = compute_nucnuc(mol.coordinates, mol.pseudo_numbers)

# Read an old set of orbitals from a previous OO-AP1roG example
# -------------------------------------------------------------
# In this case, it is assumed that the example is executed after
# water_cholesky_3-21g.py or water_dense_3-21g.py
old = IOData.from_file('checkpoint.h5')

# Re-orthogonalize the orbitals
# -----------------------------
project_orbitals_ortho(old.olp, olp, old.exp_alpha, exp_alpha)

# Do OO-AP1roG optimization
# -------------------------
ap1rog = RAp1rog(lf, occ_model)
with numpy_seed():  # reproducible 'random' numbers to make sure it always works
    energy, g, l = ap1rog(one, two, core_energy, exp_alpha, olp, True)
