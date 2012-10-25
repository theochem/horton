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


# TODO: Derive GB1GridFn and GB2GridFn from GBCalculator (in its own pxd file)
#       Use the same class hierarchy in cext.pyx

cdef extern from "fns.h":
    cdef cppclass GB1GridFn:
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

    cdef cppclass GB1GridDensityFn:
        GB1GridDensityFn(long max_shell_type) except +

    cdef cppclass GB1GridGradientFn:
        GB1GridGradientFn(long max_shell_type) except +

    cdef cppclass GB2GridFn:
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

    cdef cppclass GB2GridHartreeFn:
        GB2GridHartreeFn(long max_shell_type) except +
