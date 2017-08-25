# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


import numpy as np

from horton.meanfield.test.common import load_orbs_alpha, load_orbs_beta, get_obasis, load_na, \
    load_kin, load_olp
from .. import guess_core_hamiltonian


def test_guess_hamcore_cs():
    fname = 'hf_sto3g_fchk'
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    orb_alpha = load_orbs_alpha(fname)
    guess_core_hamiltonian(olp, kin + na, orb_alpha)
    # just a few simple checks
    assert abs(
        orb_alpha.energies[0] - (-2.59083334E+01)) > 1e-5  # values from fchk must be overwritten
    assert (orb_alpha.energies.argsort() == np.arange(olp.shape[0])).all()


def test_guess_hamcore_os():
    fname = 'li_h_3_21G_hf_g09_fchk'
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    orb_alpha = load_orbs_alpha(fname)
    orb_beta = load_orbs_beta(fname)
    guess_core_hamiltonian(olp, kin + na, orb_alpha, orb_beta)
    # just a few simple checks
    assert abs(orb_alpha.energies[0] - (
        -2.76116635E+00)) > 1e-5  # values from fchk must be overwritten
    assert abs(
        orb_beta.energies[0] - (-2.76031162E+00)) > 1e-5  # values from fchk must be overwritten
    assert (orb_alpha.energies.argsort() == np.arange(get_obasis(fname).nbasis)).all()
    assert abs(orb_alpha.energies - orb_beta.energies).max() < 1e-10
    assert abs(orb_alpha.coeffs - orb_beta.coeffs).max() < 1e-10
