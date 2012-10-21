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


cdef extern from "fns.h":
    cdef cppclass GB1GridFn:
        GB1GridFn(long max_shell_type) except +

        long get_nwork()
        long get_max_shell_type()
        long get_max_nbasis()
        void reset(long shell_type0, double* r0, double* point) except +
        void add(double coeff, double alpha0, double* scales0)
        void cart_to_pure() except +

        long get_shell_type0()
        double* get_work()

    cdef cppclass GB1GridDensityFn:
        GB1GridDensityFn(long max_shell_type) except +
