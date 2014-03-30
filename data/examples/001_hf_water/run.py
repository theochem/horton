#!/usr/bin/env python

from horton import *

# Load the coordinates from file
sys = System.from_file('water.xyz', obasis='3-21G')

# Initialize the closed-shell wfn
setup_mean_field_wfn(sys, charge=0)

# Compute Gaussian integrals
olp = sys.get_overlap()
kin = sys.get_kinetic()
nai = sys.get_nuclear_attraction()
er = sys.get_electron_repulsion()

# Initial guess
guess_core_hamiltonian(sys.wfn, olp, kin, nai)

# Construct the HF Hamiltonian
external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
terms = [
    OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
    DirectTerm(er, sys.lf, sys.wfn),
    ExchangeTerm(er, sys.lf, sys.wfn),
    OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
]
ham = Hamiltonian(terms, external)

# Converge WFN with SCF
converged = converge_scf(ham, sys.wfn, sys.lf, olp)

# Compute the energy
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)

# Partition the density with the Becke scheme
grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers, mode='keep')
bp = BeckeWPart(sys, grid)
bp.do_charges()
