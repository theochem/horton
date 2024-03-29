# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2022 The HORTON Development Team
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


cimport horton.grid.uniform as uniform

cdef extern from "horton/grid/utils.h":
    void dot_multi(long npoint, long nvector, double** data, long* segments,
        double* output)
    void dot_multi_moments_cube(long nvector, double** data, uniform.UniformGrid* ugrid,
        double* center, long lmax, long mtype, double* output, long nmoment) except +
    void dot_multi_moments(long npoint, long nvector, double** data, double* points,
        double* center, long lmax, long mtype, long* segments, double* output,
        long nmoment) except +
