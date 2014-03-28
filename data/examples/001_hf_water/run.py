#!/usr/bin/env python

from horton import *

# Load the coordinates from file
sys = System.from_file('water.xyz', obasis='3-21G')

# Initialize the closed-shell wfn
setup_mean_field_wfn(sys, charge=0)

# Initial WFN guess
guess_hamiltonian_core(sys)

# Construct a Hamiltonian
er = sys.get_electron_repulsion()
external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
terms = [
    KineticEnergy(sys.obasis, sys.lf, sys.wfn),
    Hartree(sys.lf, sys.wfn, er),
    HartreeFockExchange(sys.lf, sys.wfn, er),
    ExternalPotential(sys.obasis, sys.lf, sys.wfn, sys.numbers, sys.coordinates),
]
ham = Hamiltonian(sys, terms, external=external)

# Converge WFN with SCF
converged = converge_scf(ham)

# Compute the energy
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)

# Partition the density with the Becke scheme
grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers, mode='keep')
bp = BeckeWPart(sys, grid)
bp.do_charges()
