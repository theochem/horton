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


cdef extern from "horton/cell.h":
    cdef cppclass Cell:
        Cell(double* _rvecs, int _nvec) except +
        void mic(double* delta)
        void to_frac(double* cart, double* frac)
        void to_cart(double* frac, double* cart)
        void add_rvec(double* delta, long* r)
        void g_lincomb(double* coeffs, double* gvec)
        void dot_rvecs(double* cart, double* dot_rvecs)

        int get_nvec()
        double get_volume()
        double get_rlength(int i) except +
        double get_glength(int i) except +
        double get_rspacing(int i) except +
        double get_gspacing(int i) except +

        void copy_rvecs(double* _rvecs)
        void copy_gvecs(double* _gvecs)
        void copy_rlengths(double* _rspacings)
        void copy_glengths(double* _gspacings)
        void copy_rspacings(double* _rspacings)
        void copy_gspacings(double* _gspacings)

        void set_ranges_rcut(double* delta, double rcut, long* ranges_begin,
            long* ranges_end)

    long smart_wrap(long i, long shape, long pbc)
