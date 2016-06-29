# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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


cdef extern from "horton/grid/rtransform.h":
    cdef cppclass RTransform:
        RTransform(int npoint) except +
        double radius(double t)
        double deriv(double t)
        double deriv2(double t)
        double deriv3(double t)
        double inv(double t)

        void radius_array(double* t, double* r, int n)
        void deriv_array(double* t, double* d, int n)
        void deriv2_array(double* t, double* d, int n)
        void deriv3_array(double* t, double* d, int n)
        void inv_array(double* r, double* t, int n)
        int get_npoint()


    cdef cppclass IdentityRTransform:
        IdentityRTransform(int npoint) except +


    cdef cppclass LinearRTransform:
        LinearRTransform(double rmin, double rmax, int npoint) except +

        double get_rmin()
        double get_rmax()
        double get_alpha()


    cdef cppclass ExpRTransform:
        ExpRTransform(double rmin, double rmax, int npoint) except +

        double get_rmin()
        double get_rmax()
        double get_alpha()


    cdef cppclass PowerRTransform:
        PowerRTransform(double rmin, double rmax, int npoint) except +

        double get_rmin()
        double get_rmax()
        double get_power()


    cdef cppclass HyperbolicRTransform:
        HyperbolicRTransform(double a, double b, int npoint) except +

        double get_a()
        double get_b()
