#!/usr/bin/env python

from horton import *

# Load the coordinates from file
mol = Molecule.from_file('water.xyz')

# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')

# Create a linalg factory
lf = DenseLinalgFactory(obasis.nbasis)

# Compute Gaussian integrals
olp = obasis.compute_overlap(lf)
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)

# Create alpha orbitals
exp_alpha = lf.create_expansion()

# Initial guess
guess_core_hamiltonian(olp, kin, na, exp_alpha)

# Construct the restricted HF effective Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
terms = [
    RTwoIndexTerm(kin, 'kin'),
    RDirectTerm(er, 'hartree'),
    RExchangeTerm(er, 'x_hf'),
    RTwoIndexTerm(na, 'ne'),
]
ham = REffHam(terms, external)

# Decide how to occupy the orbitals (5 alpha electrons)
occ_model = AufbauOccModel(5)

# Converge WFN with plain SCF
scf_solver = PlainSCFSolver(1e-6)
scf_solver(ham, lf, olp, occ_model, exp_alpha)

# Partition the density with the Becke scheme
# - construct a numerical integration grid
grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, mode='keep')
# - compute the electron density on a grid
dm_alpha = exp_alpha.to_dm()
moldens = 2*obasis.compute_grid_density_dm(dm_alpha, grid.points)
# - initialize the partitioning code
wpart = BeckeWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens)
# - compute and print the charges
wpart.do_charges()
print wpart['charges']
