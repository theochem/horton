#!/usr/bin/env python

from horton import *

# Set up molecule, define basis set
# ---------------------------------
# get the XYZ file from HORTON's test data directory
fn_xyz = context.get_fn('test/water.xyz')
mol = IOData.from_file(fn_xyz)
obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
lf = CholeskyLinalgFactory(obasis.nbasis)

# Construct Hamiltonian
# ---------------------
mol.lf = lf
mol.kin = obasis.compute_kinetic(lf)
mol.na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
mol.er = obasis.compute_electron_repulsion(lf)
mol.core_energy = compute_nucnuc(mol.coordinates, mol.pseudo_numbers)

# Write to a HDF5 file
# --------------------
mol.to_file('hamiltonian_ao.h5')
