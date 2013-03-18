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



from horton import UniformIntGrid


def reduce_data(cube_data, ui_grid, factor):
    if (ui_grid.shape % factor != 0).any():
        raise ValueError('The stride is not commensurate with all three grid demsions.')

    new_cube_data = cube_data[::factor, ::factor, ::factor].copy()

    new_shape = ui_grid.shape/factor
    grid_rvecs = ui_grid.grid_cell.rvecs*factor
    new_ui_grid = UniformIntGrid(ui_grid.origin, grid_rvecs, new_shape, ui_grid.pbc)

    return new_cube_data, new_ui_grid
