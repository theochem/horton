#!/usr/bin/env python

from horton import *

# Load the coordinates from file
sys = System.from_file('water.xyz', obasis='3-21G')

# Initialize the closed-shell wfn
setup_mean_field_wfn(sys, charge=0)

# Initial WFN guess
guess_hamiltonian_core(sys)

# Construct a Hamiltonian
ham = Hamiltonian(sys, [HartreeFockExchange()])

# Converge WFN with SCF
converged = converge_scf(ham)

# Compute the energy
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)

# Partition the density with the Becke scheme
grid = BeckeMolGrid(sys, mode='keep')
bp = BeckeWPart(sys, grid)
bp.do_charges()
