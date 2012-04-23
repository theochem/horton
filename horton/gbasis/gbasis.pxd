# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

cdef extern from "gbasis.h":
    cdef cppclass GBasis:
        # Arrays that fully describe the basis set.
        double* centers
        long* shell_map
        long* nprims
        long* shell_types
        double* alphas
        double* con_coeffs
        long ncenter, nshell, nprim_total

        # Auxiliary bits
        long get_nbasis()
        long get_nscales()
        long get_max_shell_nbasis()
        double* get_scales(long iprim)

    cdef cppclass GOBasis:
        GOBasis(double* centers, long* shell_map, long* nprims,
                long* shell_types, double* alphas, double* con_coeffs,
                long ncenter, long nshell, long nprim_total)

        void compute_gobasis_overlap(double* output)
