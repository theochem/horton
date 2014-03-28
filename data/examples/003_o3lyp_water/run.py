#!/usr/bin/env python


from horton import *

# Loading system
sys = System.from_file('water.xyz', obasis='3-21G')

# Allocate wavefuntion
setup_mean_field_wfn(sys, charge=0, mult=1)

# Initial guess
guess_hamiltonian_core(sys)

# Setup integration grids with default settings
grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers)

# Construction of Hamiltonian
er = sys.get_electron_repulsion()
libxc_term = LibXCHybridGGA(sys.lf, sys.wfn, 'xc_o3lyp')
terms = [
    KineticEnergy(sys.obasis, sys.lf, sys.wfn),
    Hartree(sys.lf, sys.wfn, er),
    libxc_term,
    HartreeFockExchange(sys.lf, sys.wfn, er,
                        fraction_exchange=libxc_term.get_exx_fraction()),
    ExternalPotential(sys.obasis, sys.lf, sys.wfn, sys.numbers, sys.coordinates),
]
ham = Hamiltonian(sys, terms, grid)

# Optimal damping SCF cycle
converged = converge_scf_oda(ham)

# Energy computation
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)
