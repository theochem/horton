#!/usr/bin/env python
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


from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_base_exceptions():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 110), random_rotate=False, mode='discard')
    dm_full = mol.get_dm_full()
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid.points)
    with assert_raises(ValueError):
        # the default setting is local=true, which is not compatible with mode='discard'.
        dp = WPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 110), random_rotate=False, mode='keep')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 110), random_rotate=False, mode='only')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 110), random_rotate=False, mode='discard')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens, local=False)
