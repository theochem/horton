#!/usr/bin/env python

from horton import *

# Load the coordinates from file
sys = System.from_file('water.xyz', obasis='3-21G')

# Initialize the closed-shell wfn
setup_mean_field_wfn(sys, charge=0)

# Initial WFN guess
guess_hamiltonian_core(sys)

# set up cache for Hamiltonian
scf_cache = Cache()

# Construct a Hamiltonian
ham = Hamiltonian(sys, scf_cache, [HartreeFockExchange(scf_cache, sys.lf, sys.wfn,
                                           sys.get_electron_repulsion())])

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
