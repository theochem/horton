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


cdef extern from "cell.h":
    cdef cppclass Cell:
        Cell()
        void update(double* _rvecs, double* _gvecs, int _nvec)
        void mic(double* delta)
        void to_center(double* cart, long* center)
        void to_frac(double* cart, double* frac)
        void to_cart(double* frac, double* cart)
        void add_vec(double* delta, long* r)

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

        void set_ranges_rcut(double* delta, double rcut, long mode,
            long* ranges_begin, long* ranges_end)
        long select_inside(double* origin, double* center, double rcut,
            long* ranges_begin, long* ranges_end, long* shape,
            long* pbc_active, long* indexes) except +

    long smart_wrap(long i, long shape, long pbc_active)
