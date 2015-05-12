#!/usr/bin/env python
#JSON {"lot": "UHF/cc-pVDZ",
#JSON  "scf": "PlainSCFSolver",
#JSON  "linalg": "CholeskyLinalgFactory",
#JSON  "difficulty": 2,
#JSON  "description": "Basic UHF example with Cholesky decomposition of the ERI"}

from horton import *

# Load the coordinates from file
mol = Molecule.from_file('methyl.xyz')

# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pVDZ')

# Create a linalg factory
lf = DenseLinalgFactory(obasis.nbasis)

# Compute Gaussian integrals
olp = obasis.compute_overlap(lf)
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)

# Create alpha orbitals
exp_alpha = lf.create_expansion()
exp_beta = lf.create_expansion()

# Initial guess
guess_core_hamiltonian(olp, kin, na, exp_alpha, exp_beta)

# Construct the restricted HF effective Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
terms = [
    UTwoIndexTerm(kin, 'kin'),
    UDirectTerm(er, 'hartree'),
    UExchangeTerm(er, 'x_hf'),
    UTwoIndexTerm(na, 'ne'),
]
ham = UEffHam(terms, external)

# Decide how to occupy the orbitals (5 alpha electrons, 4 beta electrons)
occ_model = AufbauOccModel(5, 4)

# Converge WFN with plain SCF
scf_solver = PlainSCFSolver(1e-6)
scf_solver(ham, lf, olp, occ_model, exp_alpha, exp_beta)

# Assign results to the molecule object and write it to a file, e.g. for
# later analysis
mol.title = 'UHF computation on water'
mol.energy = ham.cache['energy']
mol.obasis = obasis
mol.exp_alpha = exp_alpha
mol.exp_beta = exp_beta
mol.dm_alpha = ham.cache['dm_alpha']
mol.dm_beta = ham.cache['dm_beta']

# useful for visualization:
mol.to_file('methyl.molden')
# useful for post-processing (results stored in double precision)
mol.to_file('methyl.h5')
