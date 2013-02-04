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


import numpy as np

from horton.cext import Cell
from horton.grid.cext import dot_multi, eval_spline_cube


__all__ = ['UniformIntGrid']


class UniformIntGrid(object):
    def __init__(self, origin, grid_cell, shape, pbc_active=None):
        if grid_cell.nvec != 3:
            raise ValueError('The cell must be a 3D cell.')
        self.origin = origin
        self.grid_cell = grid_cell
        self.shape = shape
        if pbc_active is None:
            self.pbc_active = np.ones(3, int)
        else:
            self.pbc_active = pbc_active

    def eval_spline(self, spline, center, output):
        eval_spline_cube(spline, center, output, self.origin, self.grid_cell, self.shape, self.pbc_active)

    def integrate(self, *args):
        '''Integrate the product of all arguments

           **Arguments:**

           data1, data2, ...
                All arguments must be arrays with the same size as the number
                of grid points. The arrays contain the functions, evaluated
                at the grid points, that must be multiplied and integrated.
        '''
        # This is often convenient for cube grid data:
        args = [arg.ravel() for arg in args]
        # Similar to conventional integration routine:
        return dot_multi(*args)*self.grid_cell.volume
