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

cdef extern from "ints.h":
    double gpt_coeff(long k, long n0, long n1, double pa, double pb)
    double gb_overlap_int1d(long n0, long n1, double pa, double pb, double gamma)
    double gob_normalization(double alpha, long* n)

    cdef cppclass GB2Integral:
        long get_max_nbasis()
        void reset(long shell_type0, long shell_type1, double* r0, double* r1)
        void add(double coeff, double alpha0, double alpha1, double* scales0, double* scales1)
        void cart_to_pure()

        long get_shell_type0()
        long get_shell_type1()
        double* get_work()

    cdef cppclass GB2OverlapIntegral:
        GB2OverlapIntegral(long max_nbasis)

    cdef cppclass GB2KineticIntegral:
        GB2KineticIntegral(long max_nbasis)
