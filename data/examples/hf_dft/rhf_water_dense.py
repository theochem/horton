#!/usr/bin/env python
#JSON {"lot": "RHF/3-21G",
#JSON  "scf": "PlainSCFSolver",
#JSON  "er": "dense",
#JSON  "difficulty": 1,
#JSON  "description": "Basic RHF example with dense matrices"}

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


# Hartree-Fock calculation
# ------------------------

# Load the coordinates from file.
# Use the XYZ file from HORTON's test data directory.
fn_xyz = context.get_fn('test/water.xyz')
mol = IOData.from_file(fn_xyz)

# Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')

# Compute Gaussian integrals
olp = obasis.compute_overlap()
kin = obasis.compute_kinetic()
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
er = obasis.compute_electron_repulsion()

# Create alpha orbitals
orb_alpha = Orbitals(obasis.nbasis)

# Initial guess
one = kin + na
guess_core_hamiltonian(olp, one, orb_alpha)

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
scf_solver(ham, olp, occ_model, orb_alpha)

# Assign results to the molecule object and write it to a file, e.g. for
# later analysis
mol.title = 'RHF computation on water'
mol.energy = ham.cache['energy']
mol.obasis = obasis
mol.orb_alpha = orb_alpha


# Write SCF results to a file
# ---------------------------

# useful for visualization:
mol.to_file('water-scf.molden')
# useful for post-processing (results stored in double precision)
mol.to_file('water-scf.h5')


# Export Hamiltonian in Hartree-Fock molecular orbital basis (all orbitals active)
# --------------------------------------------------------------------------------

# Transform orbitals
(one_mo,), (two_mo,) = transform_integrals(one, er, 'tensordot', mol.orb_alpha)

# Prepare an IOData object for writing the Hamiltonian.
mol_all_active = IOData(core_energy=external['nn'], one_mo=one_mo, two_mo=two_mo)
# useful for exchange with other codes
mol_all_active.to_file('water.FCIDUMP')
# Useful for exchange with other HORTON scripts
mol_all_active.to_file('water-hamiltonian.h5')


# Export Hamiltonian in Hartree-Fock molecular orbital basis for CAS(8,8)
# -----------------------------------------------------------------------

# Transform orbitals
one_small, two_small, core_energy = split_core_active(one, er,
    external['nn'], orb_alpha, ncore=2, nactive=8)

# Write files
mol_cas88 = IOData(core_energy=core_energy, one_mo=one_mo, two_mo=two_mo, nelec=8, ms2=0)
# useful for exchange with other codes
mol_cas88.to_file('n2-cas8-8.FCIDUMP')
# useful for exchange with other HORTON scripts
mol_cas88.to_file('n2-hamiltonian-cas8-8.h5')


# CODE BELOW IS FOR horton-regression-test.py ONLY. IT IS NOT PART OF THE EXAMPLE.
rt_results = {
    'energy': ham.cache['energy'],
    'orb_alpha': orb_alpha.energies,
    'nn': ham.cache["energy_nn"],
    'kin': ham.cache["energy_kin"],
    'ne': ham.cache["energy_ne"],
    'hartree': ham.cache["energy_hartree"],
    'x_hf': ham.cache["energy_x_hf"],
}
# BEGIN AUTOGENERATED CODE. DO NOT CHANGE MANUALLY.
import numpy as np  # pylint: disable=wrong-import-position
rt_previous = {
    'energy': -75.585812494951796,
    'orb_alpha': np.array([
        -20.424127702531802, -1.3220952123937071, -0.69192873985781, -0.52710998103304163,
        -0.47641493329615758, 0.26331191570767409, 0.36338901324987316,
        1.2209721590205258, 1.2733593057506916, 1.7850859585681533, 1.8583713969660212,
        2.029058720163976, 3.1043576600512139
    ]),
    'hartree': 46.82886903477147,
    'kin': 75.44720801802484,
    'ne': -198.04980467784898,
    'nn': 9.1571750364299866,
    'x_hf': -8.969259906329114,
}
