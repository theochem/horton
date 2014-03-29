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
kin = sys.get_kinetic()
nai = sys.get_nuclear_attraction()
er = sys.get_electron_repulsion()
external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
terms = [
    OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
    DirectTerm(er, sys.lf, sys.wfn),
    GridGroup(sys.obasis, grid, sys.lf, sys.wfn, [
        DiracExchange(sys.wfn),
    ]),
    OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
]
ham = Hamiltonian(sys, terms, external)

# Optimal damping SCF cycle
converged = converge_scf_oda(ham)

# Energy computation
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)
