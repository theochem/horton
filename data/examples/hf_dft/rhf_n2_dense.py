#!/usr/bin/env python
#JSON {"lot": "RHF/cc-pvtz",
#JSON  "scf": "PlainSCFSolver",
#JSON  "linalg": "DenseLinalgFactory",
#JSON  "difficulty": 1,
#JSON  "description": "Basic RHF example with dense matrices, includes export of Hamiltonian"}

from horton import *
import numpy as np


# Hartree-Fock calculation
# ------------------------

# Construct a molecule from scratch
bond_length = 1.098*angstrom
mol = IOData(title='dinitrogen')
mol.coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, bond_length]])
mol.numbers = np.array([7, 7])

# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')

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

# Decide how to occupy the orbitals (7 alpha electrons)
occ_model = AufbauOccModel(7)

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
mol.to_file('n2-scf.molden')
# useful for post-processing (results stored in double precision)
mol.to_file('n2-scf.h5')


# Export Hamiltonian in Hartree-Fock molecular orbital basis (all orbitals active)
# --------------------------------------------------------------------------------

# Transform orbitals
one = kin.copy()
one.iadd(na)
two = er
(one_mo,), (two_mo,) = transform_integrals(one, two, 'tensordot', mol.exp_alpha)

# Write files
mol_all_active = IOData(core_energy=external['nn'], one_mo=one_mo, two_mo=two_mo)
# useful for exchange with other codes
mol_all_active.to_file('n2.FCIDUMP')
# useful for exchange with other HORTON scripts
mol_all_active.to_file('n2.h5')


# Export Hamiltonian in Hartree-Fock molecular orbital basis for CAS(8,8)
# -----------------------------------------------------------------------

# Transform orbitals
one_small, two_small, core_energy = split_core_active(one, er,
    external['nn'], exp_alpha, ncore=2, nactive=8)

# Write files
mol_cas88 = IOData(core_energy=core_energy, one_mo=one_mo, two_mo=two_mo, nelec=8, ms2=0)
# useful for exchange with other codes
mol_cas88.to_file('n2-cas8-8.FCIDUMP')
# useful for exchange with other HORTON scripts
mol_cas88.to_file('n2-cas8-8.h5')
