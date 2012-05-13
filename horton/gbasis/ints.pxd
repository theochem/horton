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
        long get_nwork()
        long get_max_shell_type()
        long get_max_nbasis()
        void reset(long shell_type0, long shell_type1, double* r0, double* r1)
        void add(double coeff, double alpha0, double alpha1, double* scales0, double* scales1)
        void cart_to_pure() except +

        long get_shell_type0()
        long get_shell_type1()
        double* get_work()

    cdef cppclass GB2OverlapIntegral:
        GB2OverlapIntegral(long max_shell_type)

    cdef cppclass GB2KineticIntegral:
        GB2KineticIntegral(long max_shell_type)

    void nuclear_attraction_helper(double* work_g, long n0, long n1, double pa, double pb, double cp, double gamma_inv)

    cdef cppclass GB2NuclearAttractionIntegral:
        GB2NuclearAttractionIntegral(long max_shell_type, double* charges, double* centers, long ncharge)

    cdef cppclass GB4Integral:
        long get_nwork()
        long get_max_shell_type()
        long get_max_nbasis()
        void reset(long shell_type0, long shell_type1, long shell_type2, long shell_type3, double* r0, double* r1, double* r2, double* r3)
        void add(double coeff, double alpha0, double alpha1, double alpha2, double alpha3, double* scales0, double* scales1, double* scales2, double* scales3)
        void cart_to_pure() except +

        long get_shell_type0()
        long get_shell_type1()
        long get_shell_type2()
        long get_shell_type3()
        double* get_work()

    cdef cppclass GB4ElectronReuplsionIntegralLibInt:
        GB4ElectronReuplsionIntegralLibInt(long max_shell_type)


cdef extern from "libint2.h":
    void libint2_static_init()
    void libint2_static_cleanup()
