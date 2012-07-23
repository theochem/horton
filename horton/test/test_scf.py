# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


import numpy as np
from horton import *


def test_scf_cs():
    fn_fchk = context.get_fn('test/hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    guess_hamiltionian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    converge_scf(ham)
    expected_energies = np.array([
        -2.59083334E+01, -1.44689996E+00, -5.57467136E-01, -4.62288194E-01,
        -4.62288194E-01, 5.39578910E-01,
    ])
    assert abs(sys.wfn.expansion.energies - expected_energies).max() < 1e-5


def test_scf_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    guess_hamiltionian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    converge_scf(ham)
    expected_alpha_energies = np.array([
        -2.76116635E+00, -7.24564188E-01, -1.79148636E-01, -1.28235698E-01,
        -1.28235698E-01, -7.59817520E-02, -1.13855167E-02, 6.52484445E-03,
        6.52484445E-03, 7.52201895E-03, 9.70893294E-01,
    ])
    expected_beta_energies = np.array([
        -2.76031162E+00, -2.08814026E-01, -1.53071066E-01, -1.25264964E-01,
        -1.25264964E-01, -1.24605870E-02, 5.12761388E-03, 7.70499854E-03,
        7.70499854E-03, 2.85176080E-02, 1.13197479E+00,
    ])
    assert abs(sys.wfn.alpha_expansion.energies - expected_alpha_energies).max() < 1e-5
    assert abs(sys.wfn.beta_expansion.energies - expected_beta_energies).max() < 1e-5
