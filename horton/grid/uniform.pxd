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


cimport cell

cdef extern from "uniform.h":
    cdef cppclass UniformIntGrid:
        UniformIntGrid(double* _origin, cell.Cell* _grid_cell, long* _shape, long* _pbc, cell.Cell* _cell) except +

        void copy_origin(double* output)
        void copy_shape(long* output)
        void copy_pbc(long* output)

        void set_ranges_rcut(double* center, double rcut, long* ranges_begin, long* ranges_end)
        double dist_grid_point(double* center, long* i)
        void delta_grid_point(double* center, long* i)

    cdef cppclass UniformIntGridWindow:
        UniformIntGridWindow(UniformIntGrid* ui_grid, long* begin, long* end)

        void copy_begin(long* output)
        void copy_end(long* output)

        void extend(double* small, double* output)

    long index_wrap(long i, long high)

    cdef cppclass Block3Iterator:
        Block3Iterator(const long* begin, const long* end, const long* shape)

        void copy_block_begin(long* output)
        void copy_block_end(long* output)
