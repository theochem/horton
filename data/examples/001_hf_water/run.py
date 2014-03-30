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
olp = lf.create_one_body()
kin = lf.create_one_body()
nai = lf.create_one_body()
er = lf.create_two_body()
obasis.compute_overlap(olp)
obasis.compute_kinetic(kin)
obasis.compute_nuclear_attraction(mol.pseudo_numbers, mol.coordinates, nai)
obasis.compute_electron_repulsion(er)

# Initial guess
guess_core_hamiltonian(wfn, olp, kin, nai)

# Construct the HF Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.numbers)}
terms = [
    OneBodyTerm(kin, lf, wfn, 'kin'),
    DirectTerm(er, lf, wfn),
    ExchangeTerm(er, lf, wfn),
    OneBodyTerm(nai, lf, wfn, 'ne'),
]
ham = Hamiltonian(terms, external)

# Converge WFN with SCF
converged = converge_scf(ham, wfn, lf, olp)

# Compute the energy
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)

# Partition the density with the Becke scheme
# TODO: fix this when horton.part rewrite is ready
#grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, mode='keep')
#bp = BeckeWPart(sys, grid)
#bp.do_charges()
