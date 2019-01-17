# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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


cdef extern from "horton/gbasis/fns.h":
    cdef cppclass GB1DMGridFn:
        long get_nwork()
        long get_max_shell_type()
        long get_max_nbasis()
        long get_dim_work()
        long get_dim_output()
        void reset(long shell_type0, double* r0, double* point) except +
        void add(double coeff, double alpha0, double* scales0)
        void cart_to_pure() except +

        long get_shell_type0()
        double* get_work()

    cdef cppclass GB1DMGridDensityFn:
        GB1DMGridDensityFn(long max_shell_type) except +

    cdef cppclass GB1DMGridGradientFn:
        GB1DMGridGradientFn(long max_shell_type) except +

    cdef cppclass GB1DMGridGGAFn:
        GB1DMGridGGAFn(long max_shell_type) except +

    cdef cppclass GB1DMGridKineticFn:
        GB1DMGridKineticFn(long max_shell_type) except +

    cdef cppclass GB1DMGridHessianFn:
        GB1DMGridHessianFn(long max_shell_type) except +

    cdef cppclass GB1DMGridMGGAFn:
        GB1DMGridMGGAFn(long max_shell_type) except +

    cdef cppclass GB2DMGridFn:
        long get_nwork()
        long get_max_shell_type()
        long get_max_nbasis()
        long get_dim_work()
        long get_dim_output()
        void reset(long shell_type0, double* r0, double* point) except +
        void add(double coeff, double alpha0, double* scales0)
        void cart_to_pure() except +

        long get_shell_type0()
        long get_shell_type1()
        double* get_work()
