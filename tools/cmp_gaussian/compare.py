#!/usr/bin/env python

import os
from glob import glob
from horton import *

for fn_fchk in sorted(glob('*/gaussian.fchk')):
    sys = System.from_file(fn_fchk)
    g09_energy = sys.props['energy']
    guess_hamiltionian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    converge_scf(ham)
    ham.compute_energy()
    error = sys.props['energy'] - g09_energy
    print '%30s % 15.10e % 15.10e % 9.1e' % (fn_fchk, sys.props['energy'], g09_energy, error)
