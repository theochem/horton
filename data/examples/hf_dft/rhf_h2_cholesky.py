#!/usr/bin/env python
#JSON {"lot": "RHF/6-31G",
#JSON  "scf": "PlainSCFSolver",
#JSON  "linalg": "CholeskyLinalgFactory",
#JSON  "difficulty": 1,
#JSON  "description": "Basic RHF example with Cholesky matrices, includes export of Hamiltonian"}

from horton import *
import numpy as np


# Hartree-Fock calculation
# ------------------------

# Construct a molecule from scratch
mol = IOData.from_file(context.get_fn('test/h2.xyz'))

# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, '6-31G')

# Create a linalg factory
lf = CholeskyLinalgFactory(obasis.nbasis)

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

# Decide how to occupy the orbitals (1 alpha electron)
occ_model = AufbauOccModel(1)

# Converge WFN with plain SCF
scf_solver = PlainSCFSolver(1e-6)
scf_solver(ham, lf, olp, occ_model, exp_alpha)


# Write SCF results to a file
# ---------------------------

# Assign results to the molecule object and write it to a file, e.g. for
# later analysis
mol.title = 'RHF computation on dinitrogen'
mol.energy = ham.cache['energy']
mol.obasis = obasis
mol.exp_alpha = exp_alpha

# useful for visualization:
mol.to_file('h2-scf.molden')
# useful for post-processing (results stored in double precision)
mol.to_file('h2-scf.h5')


# Export Hamiltonian in Hartree-Fock molecular orbital basis (all orbitals active)
# --------------------------------------------------------------------------------

# Transform orbitals
one = kin.copy()
one.iadd(na)
two = er
(one_mo,), (two_mo,) = transform_integrals(one, two, 'tensordot', mol.exp_alpha)

# Prepare an IOData object for writing the Hamiltonian.
mol_all_active = IOData(core_energy=external['nn'], one_mo=one_mo, two_mo=two_mo, lf=lf)
# The Cholesky decomposition can only be stored in the internal format.
mol_all_active.to_file('h2-hamiltonian.h5')

# Assign results to variables for regression testing
# --------------------------------------------------
result_energy = ham.cache['energy']
result_exp_alpha = exp_alpha.energies
result_nn = ham.cache["energy_nn"]
result_kin = ham.cache["energy_kin"]
result_ne = ham.cache["energy_ne"]
result_hartree = ham.cache["energy_hartree"]
result_x_hf = ham.cache["energy_x_hf"]
# --------------------------------------------------
