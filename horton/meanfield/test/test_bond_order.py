# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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


def check_bond_orders(fn):
    mol = IOData.from_file(fn)
    operators = get_mulliken_operators(mol.obasis, mol.lf)
    dm_full = mol.get_dm_full()
    if dm_full is not None:
        dm_spin = mol.get_dm_spin()
        if dm_spin is not None:
            dm_alpha = dm_full.copy()
            dm_alpha.iadd(dm_spin)
            dm_alpha.iscale(0.5)
            dm_beta = dm_full.copy()
            dm_beta.iadd(dm_spin, -1)
            dm_beta.iscale(0.5)
            bond_orders, valences, free_valences = compute_bond_orders_os(dm_alpha, dm_beta, operators)
        else:
            dm_alpha = dm_full.copy()
            dm_alpha.iscale(0.5)
            bond_orders, valences, free_valences = compute_bond_orders_cs(dm_alpha, operators)
    else:
        raise NotImplementedError
    assert abs(bond_orders.sum() - 2*mol.numbers.sum()) < 1e-3
    assert bond_orders.shape == (mol.natom, mol.natom)
    assert valences.shape == (mol.natom,)
    assert free_valences.shape == (mol.natom,)
    assert (bond_orders == bond_orders.T).all()
    assert (valences > 0).all()
    return bond_orders, valences, free_valences


def test_bond_order_water_sto3g():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    bond_orders, valences, free_valences = check_bond_orders(fn_fchk)
    assert abs(free_valences).max() < 1e-5


def test_bond_order_h3_321g():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    bond_orders, valences, free_valences = check_bond_orders(fn_fchk)
    assert (free_valences != 0).any()
