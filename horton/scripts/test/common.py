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


import shutil, os, numpy as np

from horton import context
from horton.grid.cext import UniformGrid
from horton.io.iodata import IOData


__all__ = ['copy_files', 'check_files', 'write_random_lta_cube']


def copy_files(dn, fns):
    for fn in fns:
        shutil.copy(context.get_fn(os.path.join('test', fn)), os.path.join(dn, fn))


def check_files(dn, fns):
    for fn in fns:
        assert os.path.isfile(os.path.join(dn, fn)), "Missing %s" % fn


def write_random_lta_cube(dn, fn_cube):
    '''Write a randomized cube file'''
    # start from an existing cube file
    mol = IOData.from_file(context.get_fn('test/lta_gulp.cif'))
    # Define a uniform grid with only 1000 points, to make the tests fast.
    ugrid = UniformGrid(np.zeros(3, float), mol.cell.rvecs*0.1, np.array([10, 10, 10]), np.array([1, 1, 1]))
    # Write to the file dn/fn_cube
    mol.cube_data = np.random.uniform(0, 1, ugrid.shape)
    mol.grid = ugrid
    mol.to_file(os.path.join(dn, fn_cube))
    return mol
