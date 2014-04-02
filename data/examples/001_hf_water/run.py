#!/usr/bin/env python

from horton import *

# Load the coordinates from file
mol = Molecule.from_file('water.xyz')

# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')

# Create a linalg factory
lf = DenseLinalgFactory(obasis.nbasis)

# Initialize the closed-shell wfn
wfn = setup_mean_field_wfn(obasis.nbasis, mol.numbers, lf, charge=0)

# Compute Gaussian integrals
olp = obasis.compute_overlap(lf)
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)

# Initial guess
guess_core_hamiltonian(wfn, olp, kin, na)

# Construct the HF Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
terms = [
    RestrictedOneBodyTerm(kin, 'kin'),
    RestrictedDirectTerm(er, 'hartree'),
    RestrictedExchangeTerm(er, 'x_hf'),
    RestrictedOneBodyTerm(na, 'ne'),
]
ham = RestrictedEffectiveHamiltonian(terms, external)

# Converge WFN with SCF
converged = converge_scf(ham, wfn, lf, olp)

# Compute the energy
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)

# Partition the density with the Becke scheme
grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, mode='keep')
moldens = obasis.compute_grid_density_dm(wfn.dm_full, grid.points)
bp = BeckeWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens)
bp.do_charges()
