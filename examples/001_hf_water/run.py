#!/usr/bin/env python

from horton import *

print 'Loading system'
sys = System.from_file('water.xyz', obasis='3-21G')

print 'Initializing wavefuntion'
sys.init_wfn(charge=0)

print 'Constructing initial guess'
guess_hamiltionian_core(sys)

print 'Constructing Hamiltonian'
ham = Hamiltonian(sys, [HartreeFock()])

print 'SCF cycle'
converged = converge_scf(ham)
print 'Converged?', converged

print 'Computing energies'
ham.compute_energy()
print sys.props
