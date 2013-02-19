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

cdef extern from "utils.h":
    double dot_multi(long npoint, long nvector, double** data)
    double dot_multi_poly_cube(long npoint, long nvector, double** data,
        double* origin, cell.Cell* grid_cell, long* shape, long* pbc_active,
        cell.Cell* cell, double* center, long mask, double powx, double powy,
        double powz, double powr) except +

    void grid_distances(double *points, double *center, double *distances, long n)
