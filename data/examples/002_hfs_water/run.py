#!/usr/bin/env python


from horton import *

# Loading system
system = System.from_file('water.xyz', obasis='3-21G')

# Allocate wavefuntion
setup_mean_field_wfn(system, charge=0, mult=1)

# Initial guess
guess_hamiltonian_core(system)

# Setup integration grids with default settings
grid = BeckeMolGrid(system.coordinates, system.numbers, system.pseudo_numbers)

#setup cache object for the Hamiltonian
scf_cache = Cache()

# Construction of Hamiltonian
ham = Hamiltonian(system, scf_cache,[Hartree(scf_cache, system.lf, system.wfn,
                                           system.get_electron_repulsion()),
                           DiracExchange(scf_cache, system.lf, system.wfn,
                                           system.get_electron_repulsion())],
                  grid)

# Optimal damping SCF cycle
converged = converge_scf_oda(ham)

# Energy computation
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)
