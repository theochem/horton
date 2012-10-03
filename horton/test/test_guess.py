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


def test_guess_hamcore_cs():
    fn_fchk = context.get_fn('test/hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    guess_hamiltonian_core(sys)
    # just a few simple checks
    exp_alpha = sys.wfn.get_exp('alpha')
    assert abs(exp_alpha.energies[0] - (-2.59083334E+01)) > 1e-5 # values from fchk must be overwritten
    assert (exp_alpha.energies.argsort() == np.arange(sys.obasis.nbasis)).all()
    assert sys.wfn._cache.has('dm_alpha')


def test_guess_hamcore_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    guess_hamiltonian_core(sys)
    # just a few simple checks
    exp_alpha = sys.wfn.get_exp('alpha')
    exp_beta = sys.wfn.get_exp('beta')
    assert abs(exp_alpha.energies[0] - (-2.76116635E+00)) > 1e-5 # values from fchk must be overwritten
    assert abs(exp_beta.energies[0] - (-2.76031162E+00)) > 1e-5 # values from fchk must be overwritten
    assert (exp_alpha.energies.argsort() == np.arange(sys.obasis.nbasis)).all()
    assert abs(exp_alpha.energies - exp_beta.energies).max() < 1e-10
    assert abs(exp_alpha.coeffs - exp_beta.coeffs).max() < 1e-10
    assert sys.wfn._cache.has('dm_alpha')
    assert sys.wfn._cache.has('dm_beta')
