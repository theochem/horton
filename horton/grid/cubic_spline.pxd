# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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


cimport rtransform

cdef extern from "cubic_spline.h":
    void tridiagsym_solve(double* diag_mid, double* diag_up, double* right, double* solution, int n)
    void solve_cubic_spline_system(double* y, double *d, int npoint)
    void compute_cubic_spline_int_weights(double* weights, int npoint)

    cdef cppclass Extrapolation:
        pass


    cdef cppclass CubicSpline:
        CubicSpline(double* y, double* dt, Extrapolation* ep, rtransform.RTransform* rtf, int n) except +
        void eval(double* new_x, double* new_y, int new_n)
        void eval_deriv(double* new_x, double* new_dx, int new_n)

    cdef cppclass ZeroExtrapolation:
        pass

    cdef cppclass CuspExtrapolation:
        pass
