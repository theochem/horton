#!/usr/bin/env python

from horton import *

sys = System.from_file('water.xyz', obasis='3-21G')

sys.init_wfn(charge=0)

guess_hamiltonian_core(sys)

ham = Hamiltonian(sys, [HartreeFock()])

converged = converge_scf(ham)

ham.compute_energy()

rtf = ExpRTransform(1e-3, 1e1, 100)
int1d = TrapezoidIntegrator1D()
grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=1)
bdp = BeckeDPart(grid)
bdp.do_charges()
