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
from nose.plugins.skip import SkipTest

from common import load_olp, load_orbs_alpha, load_orbs_beta
from ..moments import get_ncart_cumul, get_cartesian_powers
from ..utils import get_spin, get_homo_lumo, get_level_shift


def check_spin(fname, sz0, ssq0, eps):
    olp = load_olp(fname)
    orb_alpha = load_orbs_alpha(fname)
    try:
        orb_beta = load_orbs_beta(fname)
        sz, ssq = get_spin(orb_alpha, orb_beta, olp)
        assert abs(sz - sz0) < eps
        assert abs(ssq - ssq0) < eps
        sz, ssq = get_spin(orb_beta, orb_alpha, olp)
        assert abs(sz + sz0) < eps
        assert abs(ssq - ssq0) < eps
    except IOError:
        pass
    sz, ssq = get_spin(orb_alpha, orb_alpha, olp)
    assert abs(sz) < eps
    assert abs(ssq) < eps


def test_spin_li_h():
    check_spin('li_h_3_21G_hf_g09_fchk', 0.5, 0.75, 1e-7)


def test_spin_h3_hfs():
    check_spin('h3_hfs_321g_fchk', 0.5, 0.7530, 1e-4)


def test_spin_h3_pbe():
    check_spin('h3_pbe_321g_fchk', 0.5, 0.7530, 1e-4)


def test_spin_ch3_hf():
    check_spin('ch3_hf_sto3g_fchk', 0.5, 0.7632, 1e-4)


def test_spin_water_hf():
    check_spin('water_sto3g_hf_g03_fchk', 0.0, 0.0, 1e-8)


def check_homo_lumo(fname, homo_energy0, lumo_energy0, eps=1e-8):
    orbs = [load_orbs_alpha(fname)]
    try:
        orbs.append(load_orbs_beta(fname))
    except IOError:
        pass
    homo_energy, lumo_energy = get_homo_lumo(*orbs)
    assert abs(homo_energy - homo_energy0) < eps
    if lumo_energy0 is None:
        assert lumo_energy is None
    else:
        assert abs(lumo_energy - lumo_energy0) < eps


def test_homo_lumo_water_hf():
    check_homo_lumo('water_sto3g_hf_g03_fchk', -3.87671783E-01, 6.03082408E-01)


def test_homo_lumo_ch3_hf():
    check_homo_lumo('ch3_hf_sto3g_fchk', -3.63936540E-01, 3.28562907E-01)


def test_homo_lumo_h():
    check_homo_lumo('atom_001_001_hf_sto3g_fchk', -4.66581850E-01, 3.08024094E-01)


def test_homo_lumo_he():
    check_homo_lumo('helium_hf_sto3g_fchk', -8.76035508E-01, None)


def test_level_shift():
    fname = 'helium_hf_sto3g_fchk'
    overlap = load_olp(fname)
    dm_alpha1 = load_orbs_alpha(fname).to_dm()
    ls_alpha = get_level_shift(dm_alpha1, overlap)
    load_orbs_alpha(fname).from_fock(ls_alpha, overlap)
    dm_alpha2 = load_orbs_alpha(fname).to_dm()
    np.testing.assert_allclose(dm_alpha1, dm_alpha2)


def test_get_ncart_cumul():
    assert get_ncart_cumul(0) == 1
    assert get_ncart_cumul(1) == 4
    assert get_ncart_cumul(2) == 10
    assert get_ncart_cumul(3) == 20


def test_get_cartesian_powers():
    lmax = 4
    cartesian_powers = get_cartesian_powers(lmax)
    assert issubclass(cartesian_powers.dtype.type, int)
    assert cartesian_powers.shape == (get_ncart_cumul(lmax), 3)
    assert (cartesian_powers[0] == [0, 0, 0]).all()
    assert (cartesian_powers[1] == [1, 0, 0]).all()
    assert (cartesian_powers[4] == [2, 0, 0]).all()
    assert (cartesian_powers[5] == [1, 1, 0]).all()
    assert (cartesian_powers[6] == [1, 0, 1]).all()
    assert (cartesian_powers[10] == [3, 0, 0]).all()
    assert (cartesian_powers[11] == [2, 1, 0]).all()
    assert (cartesian_powers[19] == [0, 0, 3]).all()
    assert (cartesian_powers[-1] == [0, 0, 4]).all()

    for lmax in xrange(4):
        tmp = get_cartesian_powers(lmax)
        assert tmp.shape == (get_ncart_cumul(lmax), 3)
        assert (tmp == cartesian_powers[:len(tmp)]).all()


def test_rotate_cartesian_moments():
    raise SkipTest("Need a rotate cartesian moment test")
