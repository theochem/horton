#!/usr/bin/env python

from horton import *

print 'Loading system'
sys = System.from_file('water.xyz', obasis='3-21G')

print 'Initializing wavefuntion'
sys.init_wfn(charge=0)

print 'Constructing initial guess'
guess_hamiltonian_core(sys)

print 'Constructing Hamiltonian'
ham = Hamiltonian(sys, [HartreeFock()])

print 'SCF cycle'
converged = converge_scf(ham)
print 'Converged:', converged

print 'Computing energies'
ham.compute_energy()
print sys.props

print 'Compute Becke charges'
rtf = LogRTransform(1e-3, 1e1, 100)
int1d = TrapezoidIntegrator1D()
grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=1)
bdp = BeckeDPart(grid)
bdp.do_charges()
print bdp['charges']
