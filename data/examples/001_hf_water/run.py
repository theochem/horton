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
nai = obasis.compute_nuclear_attraction(mol.pseudo_numbers, mol.coordinates, lf)
er = obasis.compute_electron_repulsion(lf)

# Initial guess
guess_core_hamiltonian(wfn, olp, kin, nai)

# Construct the HF Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.numbers)}
terms = [
    OneBodyTerm(kin, wfn, 'kin'),
    DirectTerm(er, wfn),
    ExchangeTerm(er, wfn),
    OneBodyTerm(nai, wfn, 'ne'),
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
