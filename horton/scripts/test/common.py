# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import shutil, os, numpy as np

from horton import context, UniformGrid, System


__all__ = ['copy_files', 'check_files', 'write_random_lta_cube']


def copy_files(dn, fns):
    for fn in fns:
        shutil.copy(context.get_fn(os.path.join('test', fn)), os.path.join(dn, fn))


def check_files(dn, fns):
    for fn in fns:
        assert os.path.isfile(os.path.join(dn, fn)), "Missing %s" % fn


def write_random_lta_cube(dn, fn_cube):
    sys = System.from_file(context.get_fn('test/lta_gulp.cif'))
    ugrid = UniformGrid(np.zeros(3, float), sys.cell.rvecs*0.1, np.array([10, 10, 10]), np.array([1, 1, 1]))
    cube_data = np.random.uniform(0, 1, ugrid.shape)
    sys.update_grid(ugrid)
    sys.extra['cube_data'] = cube_data
    sys.to_file(os.path.join(dn, fn_cube))
    return sys
