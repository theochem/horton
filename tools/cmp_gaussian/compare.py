#!/usr/bin/env python

import os, numpy as np
from glob import glob
from horton import *
from horton.hamiltonian.test.common import check_cubic_os_wrapper

#log.set_level(log.high)
log.set_level(log.silent)

#fns_fchk = ['006__C_Q+0_M7/gaussian.fchk']
#fns_fchk = ['006__C_Q+0_M5/gaussian.fchk']
#fns_fchk = ['006__C_Q+2_M1/gaussian.fchk']
fns_fchk = ['005__B_Q-1_M7/gaussian.fchk']
#fns_fchk = ['010_Ne_Q+0_M1/gaussian.fchk']
#fns_fchk = ['006__C_Q-1_M2/gaussian.fchk']
#fns_fchk = ['013_Al_Q-1_M7/gaussian.fchk']
#fns_fchk = sorted(glob('*/gaussian.fchk'))

np.set_printoptions(suppress=True, linewidth=100)
print '                          Case             Horton                G09        H-G  H converged'
print '--------------------------------------------------------------------------------------------'
for fn_fchk in fns_fchk:
    sys = System.from_file(fn_fchk)
    dma0 = sys.wfn.dm_alpha.copy()
    dmb0 = sys.wfn.dm_beta.copy()
    g09_energy = sys.props['energy']
    sys.props.clear()

    guess_hamiltonian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    converged = converge_scf_oda(ham, max_iter=1024, threshold=1e-8, debug=False)
    error = sys.props['energy'] - g09_energy
    print '%30s  % 15.10e  % 15.10e  %+9.4f  %s' % (fn_fchk, sys.props['energy'], g09_energy, error, converged)

    dma1 = sys.wfn.dm_alpha.copy()
    dmb1 = sys.wfn.dm_beta.copy()
    check_cubic_os_wrapper(ham, dma0, dmb0, dma1, dmb1, do_plot=True)
    



