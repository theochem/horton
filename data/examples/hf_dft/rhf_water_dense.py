#!/usr/bin/env python
#JSON {"lot": "RHF/3-21G",
#JSON  "scf": "PlainSCFSolver",
#JSON  "linalg": "DenseLinalgFactory",
#JSON  "difficulty": 1,
#JSON  "description": "Basic RHF example with dense matrices"}

from horton import *


# Hartree-Fock calculation
# ------------------------

# Load the coordinates from file.
# Use the XYZ file from HORTON's test data directory.
fn_xyz = context.get_fn('test/water.xyz')
mol = IOData.from_file(fn_xyz)

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

# Assign results to the molecule object and write it to a file, e.g. for
# later analysis
mol.title = 'RHF computation on water'
mol.energy = ham.cache['energy']
mol.obasis = obasis
mol.exp_alpha = exp_alpha


# Write SCF results to a file
# ---------------------------

# useful for visualization:
mol.to_file('water-scf.molden')
# useful for post-processing (results stored in double precision)
mol.to_file('water-scf.h5')


# Export Hamiltonian in Hartree-Fock molecular orbital basis (all orbitals active)
# --------------------------------------------------------------------------------

# Transform orbitals
one = kin.copy()
one.iadd(na)
two = er
(one_mo,), (two_mo,) = transform_integrals(one, two, 'tensordot', mol.exp_alpha)

# Prepare an IOData object for writing the Hamiltonian.
mol_all_active = IOData(core_energy=external['nn'], one_mo=one_mo, two_mo=two_mo, lf=lf)
# useful for exchange with other codes
mol_all_active.to_file('water.FCIDUMP')
# Useful for exchange with other HORTON scripts
mol_all_active.to_file('water-hamiltonian.h5')


# Export Hamiltonian in Hartree-Fock molecular orbital basis for CAS(8,8)
# -----------------------------------------------------------------------

# Transform orbitals
one_small, two_small, core_energy = split_core_active(one, er,
    external['nn'], exp_alpha, ncore=2, nactive=8)

# Write files
mol_cas88 = IOData(core_energy=core_energy, one_mo=one_mo, two_mo=two_mo, nelec=8, ms2=0, lf=lf)
# useful for exchange with other codes
mol_cas88.to_file('n2-cas8-8.FCIDUMP')
# useful for exchange with other HORTON scripts
mol_cas88.to_file('n2-hamiltonian-cas8-8.h5')

# Assign results to variables for regression testing
# --------------------------------------------------
result_energy = ham.cache['energy']
result_exp_alpha = exp_alpha.coeffs
# --------------------------------------------------
