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


cimport horton.cell

cdef extern from "horton/grid/uniform.h":
    cdef cppclass UniformGrid:
        double origin[3]
        double grid_rvecs[9]
        long shape[3]
        long pbc[3]

        UniformGrid(double* _origin, double* _grid_rvecs, long* _shape, long* _pbc)

        horton.cell.Cell* get_cell()
        horton.cell.Cell* get_grid_cell()

        void set_ranges_rcut(double* center, double rcut, long* ranges_begin, long* ranges_end)
        double dist_grid_point(double* center, long* i)
        void delta_grid_point(double* center, long* i)

    long index_wrap(long i, long high)
