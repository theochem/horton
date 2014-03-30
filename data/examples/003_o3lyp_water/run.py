#!/usr/bin/env python


from horton import *

# Loading system
sys = System.from_file('water.xyz', obasis='3-21G')

# Allocate wavefuntion
setup_mean_field_wfn(sys, charge=0, mult=1)

# Compute Gaussian integrals
olp = sys.get_overlap()
kin = sys.get_kinetic()
nai = sys.get_nuclear_attraction()
er = sys.get_electron_repulsion()

# Initial guess
guess_core_hamiltonian(sys.wfn, olp, kin, nai)

# Setup integration grids with default settings
grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers)

# Construction of Hamiltonian
external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
libxc_term = LibXCHybridGGA(sys.wfn, 'xc_o3lyp')
terms = [
    OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
    DirectTerm(er, sys.lf, sys.wfn),
    GridGroup(sys.obasis, grid, sys.lf, sys.wfn, [
        libxc_term,
    ]),
    ExchangeTerm(er, sys.lf, sys.wfn,
                 fraction_exchange=libxc_term.get_exx_fraction()),
    OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
]
ham = Hamiltonian(terms, external)

# Optimal damping SCF cycle
converged = converge_scf_oda(ham, sys.wfn, sys.lf, olp)

# Energy computation
log.set_level(log.high)
ham.compute()
log.set_level(log.medium)
