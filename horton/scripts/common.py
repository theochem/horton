# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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



import os, sys, datetime, numpy as np

from horton import UniformIntGrid, angstrom


__all__ = ['reduce_data', 'parse_h5', 'parse_ewald_args', 'parse_pbc',
           'store_args']


def reduce_data(cube_data, ui_grid, factor):
    if (ui_grid.shape % factor != 0).any():
        raise ValueError('The stride is not commensurate with all three grid demsions.')

    new_cube_data = cube_data[::factor, ::factor, ::factor].copy()

    new_shape = ui_grid.shape/factor
    grid_rvecs = ui_grid.grid_cell.rvecs*factor
    new_ui_grid = UniformIntGrid(ui_grid.origin, grid_rvecs, new_shape, ui_grid.pbc)

    return new_cube_data, new_ui_grid


def parse_h5(arg_h5):
    if arg_h5.count(':') != 1:
        raise VallueError('An hdf5 argument must contain one colon.')
    return arg_h5.split(':')


def parse_ewald_args(args):
    rcut = args.rcut*angstrom
    alpha = args.alpha_scale/rcut
    gcut = args.gcut_scale*alpha
    return rcut, alpha, gcut


def parse_pbc(spbc):
    if len(spbc) != 3:
        raise ValueError('The pbc argument must consist of three characters, 0 or 1.')
    result = np.zeros(3, int)
    for i in xrange(3):
        if spbc[i] not in '01':
            raise ValueError('The pbc argument must consist of three characters, 0 or 1.')
        result[i] = int(spbc[i])
    return result


def store_args(args, grp):
    '''Convert the command line arguments to hdf5 attributes'''
    grp.attrs['cmdline'] =  ' '.join(sys.argv)
    grp.attrs['pwd'] = os.getcwd()
    grp.attrs['datetime'] = datetime.datetime.now().isoformat()
    for key, val in vars(args).iteritems():
        if val is not None:
            grp.attrs['arg_%s' % key] = val
