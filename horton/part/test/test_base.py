#!/usr/bin/env python
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


from nose.tools import assert_raises

from .common import load_molecule_npz
from .. base import WPart
from horton.grid import ExpRTransform, RadialGrid, BeckeMolGrid


def test_base_exceptions():
    # load molecule data
    coords, nums, pseudo_nums, dens, points = load_molecule_npz('n2_hfs_sto3g_fchk_exp:1e-3:1e1:100:110.npz')
    # make grid
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(coords, nums, pseudo_nums, (rgrid, 110), random_rotate=False, mode='discard')
    # check the grid points against stored points on which density is evaluated
    assert (abs(points - grid.points) < 1.e-6).all()

    with assert_raises(ValueError):
        # the default setting is local=true, which is not compatible with mode='discard'.
        dp = WPart(coords, nums, pseudo_nums, grid, dens)

    grid = BeckeMolGrid(coords, nums, pseudo_nums, (rgrid, 110), random_rotate=False, mode='keep')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(coords, nums, pseudo_nums, grid, dens)

    grid = BeckeMolGrid(coords, nums, pseudo_nums, (rgrid, 110), random_rotate=False, mode='only')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(coords, nums, pseudo_nums, grid, dens)

    grid = BeckeMolGrid(coords, nums, pseudo_nums, (rgrid, 110), random_rotate=False, mode='discard')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(coords, nums, pseudo_nums, grid, dens, local=False)
