#!/usr/bin/env python
#JSON {"lot": "UKS/6-31G(d)",
#JSON  "scf": "CDIISSCFSolver",
#JSON  "er": "cholesky",
#JSON  "difficulty": 6,
#JSON  "description": "Basic UKS DFT example with MGGA exhange-correlation functional (TPSS)"}

import numpy as np
from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


# Load the coordinates from file.
# Use the XYZ file from HORTON's test data directory.
fn_xyz = context.get_fn('test/methyl.xyz')
mol = IOData.from_file(fn_xyz)

# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, '6-31g(d)')

# Compute Gaussian integrals
olp = obasis.compute_overlap()
kin = obasis.compute_kinetic()
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
er_vecs = obasis.compute_electron_repulsion_cholesky()

# Define a numerical integration grid needed the XC functionals
grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)

# Create alpha orbitals
orb_alpha = Orbitals(obasis.nbasis)
orb_beta = Orbitals(obasis.nbasis)

# Initial guess
guess_core_hamiltonian(olp, kin + na, orb_alpha, orb_beta)

# Construct the restricted HF effective Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
terms = [
    UTwoIndexTerm(kin, 'kin'),
    UDirectTerm(er_vecs, 'hartree'),
    UGridGroup(obasis, grid, [
        ULibXCMGGA('x_tpss'),
        ULibXCMGGA('c_tpss'),
    ]),
    UTwoIndexTerm(na, 'ne'),
]
ham = UEffHam(terms, external)

# Decide how to occupy the orbitals (5 alpha electrons, 4 beta electrons)
occ_model = AufbauOccModel(5, 4)

# Converge WFN with CDIIS SCF
# - Construct the initial density matrix (needed for CDIIS).
occ_model.assign(orb_alpha, orb_beta)
dm_alpha = orb_alpha.to_dm()
dm_beta = orb_beta.to_dm()
# - SCF solver
scf_solver = CDIISSCFSolver(1e-6)
scf_solver(ham, olp, occ_model, dm_alpha, dm_beta)

# Derive orbitals (coeffs, energies and occupations) from the Fock and density
# matrices. The energy is also computed to store it in the output file below.
fock_alpha = np.zeros(olp.shape)
fock_beta = np.zeros(olp.shape)
ham.reset(dm_alpha, dm_beta)
ham.compute_energy()
ham.compute_fock(fock_alpha, fock_beta)
orb_alpha.from_fock_and_dm(fock_alpha, dm_alpha, olp)
orb_beta.from_fock_and_dm(fock_beta, dm_beta, olp)

# Assign results to the molecule object and write it to a file, e.g. for
# later analysis. Note that the CDIIS algorithm can only really construct an
# optimized density matrix and no orbitals.
mol.title = 'UKS computation on methyl'
mol.energy = ham.cache['energy']
mol.obasis = obasis
mol.orb_alpha = orb_alpha
mol.orb_beta = orb_beta
mol.dm_alpha = dm_alpha
mol.dm_beta = dm_beta

# useful for post-processing (results stored in double precision):
mol.to_file('methyl.h5')

# CODE BELOW IS FOR horton-regression-test.py ONLY. IT IS NOT PART OF THE EXAMPLE.
rt_results = {
    'energy': ham.cache['energy'],
    'orb_alpha': orb_alpha.energies,
    'orb_beta': orb_beta.energies,
    'nn': ham.cache["energy_nn"],
    'kin': ham.cache["energy_kin"],
    'ne': ham.cache["energy_ne"],
    'grid': ham.cache["energy_grid_group"],
    'hartree': ham.cache["energy_hartree"],
}
# BEGIN AUTOGENERATED CODE. DO NOT CHANGE MANUALLY.
rt_previous = {
    'energy': -39.836677347925914,
    'grid': -6.447581456959345,
    'hartree': 28.092648155263227,
    'kin': 39.36712413305364,
    'ne': -109.92865317406965,
    'nn': 9.079784942663636,
    'orb_alpha': np.array([
        -10.01932345515595, -0.6135259960982921, -0.36276560205876457,
        -0.3627187398579602, -0.19075935105245162, 0.07150035135891722,
        0.14494635637871764, 0.14498775711923603, 0.510315601262867, 0.5521868550528035,
        0.5522018808304408, 0.6683335215978065, 0.8453792117712211, 0.8454609326443266,
        0.8983531847806548, 1.619269238542316, 1.6192862967691806, 1.9292046668075504,
        2.093858169207324, 2.0946697177558007
    ]),
    'orb_beta': np.array([
        -10.003689426590137, -0.5754877487524174, -0.3559723619217707,
        -0.35584153775415894, -0.06554258134740852, 0.09814780282274696,
        0.1593816751966446, 0.15964525579884276, 0.5688913029872666, 0.5691828001315896,
        0.6008399202646769, 0.7037253351266851, 0.8792123495562243, 0.8812997469246673,
        0.957810069192318, 1.7275803186575, 1.729736365896011, 2.0569028489921184,
        2.1726342557791885, 2.173521534951844
    ]),
}
