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

# Construction of Hamiltonian
libxc_term = LibXCHybridGGA(system.lf, system.wfn, 'xc_o3lyp')
xhf_term = HartreeFockExchange(system.lf, system.wfn, system.get_electron_repulsion(),
                               fraction_exchange=libxc_term.get_exx_fraction())
ham = Hamiltonian(system, [xhf_term, libxc_term], grid)

# Optimal damping SCF cycle
converged = converge_scf_oda(ham)

# Energy computation
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)
