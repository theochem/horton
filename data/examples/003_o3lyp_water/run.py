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

# Setup integration grids with default settings
grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)

# Construction of Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
libxc_term = LibXCHybridGGA(wfn, 'xc_o3lyp')
terms = [
    OneBodyTerm(kin, wfn, 'kin'),
    DirectTerm(er, wfn, 'hartree'),
    GridGroup(obasis, grid, wfn, [
        libxc_term,
    ]),
    ExchangeTerm(er, wfn, 'x_hf', libxc_term.get_exx_fraction()),
    OneBodyTerm(na, wfn, 'ne'),
]
ham = Hamiltonian(terms, external)

# Optimal damping SCF cycle
converged = converge_scf_oda(ham, wfn, lf, olp)

# Energy computation
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)
