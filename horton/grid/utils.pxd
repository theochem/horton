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

cimport uniform

cdef extern from "utils.h":
    double dot_multi(long npoint, long nvector, double** data)
    void dot_multi_moments_cube(long nvector, double** data, uniform.UniformGrid* ugrid, double* center, long lmax, long mtype, double* output, long nmoment) except +
    void dot_multi_moments(long npoint, long nvector, double** data, double* points, double* center, long lmax, long mtype, double* output, long nmoment) except +
    void dot_multi_parts(long npoint, long nvector, long noutput, double** data, long* sizes, double* output)

    void grid_distances(double *points, double *center, double *distances, long n)
