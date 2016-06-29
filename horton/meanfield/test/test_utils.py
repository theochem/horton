# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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


from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def check_spin(fn_fchk, sz0, ssq0, eps):
    path_fchk = context.get_fn('test/%s' % fn_fchk)
    mol = IOData.from_file(path_fchk)
    olp = mol.obasis.compute_overlap(mol.lf)
    if hasattr(mol, 'exp_beta'):
        sz, ssq = get_spin(mol.exp_alpha, mol.exp_beta, olp)
        assert abs(sz - sz0) < eps
        assert abs(ssq - ssq0) < eps
        sz, ssq = get_spin(mol.exp_beta, mol.exp_alpha, olp)
        assert abs(sz + sz0) < eps
        assert abs(ssq - ssq0) < eps
    sz, ssq = get_spin(mol.exp_alpha, mol.exp_alpha, olp)
    assert abs(sz) < eps
    assert abs(ssq) < eps


def test_spin_li_h():
    check_spin('li_h_3-21G_hf_g09.fchk', 0.5, 0.75, 1e-7)


def test_spin_h3_hfs():
    check_spin('h3_hfs_321g.fchk', 0.5, 0.7530, 1e-4)


def test_spin_h3_pbe():
    check_spin('h3_pbe_321g.fchk', 0.5, 0.7530, 1e-4)


def test_spin_ch3_hf():
    check_spin('ch3_hf_sto3g.fchk', 0.5, 0.7632, 1e-4)


def test_spin_water_hf():
    check_spin('water_sto3g_hf_g03.fchk', 0.0, 0.0, 1e-8)


def check_homo_lumo(fn_fchk, homo_energy0, lumo_energy0, eps=1e-8):
    mol = IOData.from_file(context.get_fn('test/%s' % fn_fchk))
    exps = [mol.exp_alpha]
    if hasattr(mol, 'exp_beta'):
        exps.append(mol.exp_beta)
    homo_energy, lumo_energy = get_homo_lumo(*exps)
    assert abs(homo_energy - homo_energy0) < eps
    if lumo_energy0 is None:
        assert lumo_energy is None
    else:
        assert abs(lumo_energy - lumo_energy0) < eps


def test_homo_lumo_water_hf():
    check_homo_lumo('water_sto3g_hf_g03.fchk', -3.87671783E-01, 6.03082408E-01)


def test_homo_lumo_ch3_hf():
    check_homo_lumo('ch3_hf_sto3g.fchk', -3.63936540E-01, 3.28562907E-01)


def test_homo_lumo_h():
    check_homo_lumo('atom_001_001_hf_sto3g.fchk', -4.66581850E-01, 3.08024094E-01)


def test_homo_lumo_he():
    check_homo_lumo('helium_hf_sto3g.fchk', -8.76035508E-01, None)


def test_level_shift():
    fn_fchk = context.get_fn('test/helium_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    overlap = mol.obasis.compute_overlap(mol.lf)
    dm_alpha1 = mol.exp_alpha.to_dm()
    ls_alpha = get_level_shift(dm_alpha1, overlap)
    mol.exp_alpha.from_fock(ls_alpha, overlap)
    dm_alpha2 = mol.exp_alpha.to_dm()
    assert abs(dm_alpha1._array - dm_alpha2._array) < 1e-5
