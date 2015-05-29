# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2015 The Horton Development Team
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


cimport libcpp

cdef extern from "horton/grid/ode2.h":
    double hermite_overlap2(long xmax, long i0, libcpp.bool deriv0, long i1, libcpp.bool deriv1)
    double hermite_overlap3(long xmax, long i0, libcpp.bool deriv0, long i1, libcpp.bool deriv1, long i2, libcpp.bool deriv2)
    double hermite_node(long x, long center, libcpp.bool kind, libcpp.bool deriv)
    double hermite_product2(long x, long i0, libcpp.bool deriv0, long i1, libcpp.bool deriv1)
    void build_ode2(double* b, double* a, double *f, double** bcs, double* coeffs,
                    double* rhs, long npoint)
