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


from libc.stdint cimport int64_t

cimport rtransform

cdef extern from "horton/grid/cubic_spline.h":
    void tridiagsym_solve(double* diag_mid, double* diag_up, double* right, double* solution, int n)
    void solve_cubic_spline_system(double* y, double *d, int npoint)
    void compute_cubic_spline_int_weights(double* weights, int npoint)

    cdef cppclass Extrapolation:
        double eval_left(double x)
        double eval_right(double x)
        double deriv_left(double x)
        double deriv_right(double x)

    cdef cppclass CubicSpline:
        CubicSpline(double* y, double* dt, Extrapolation* extrapolation, rtransform.RTransform* rtf, int n)
        void eval(double* new_x, double* new_y, int new_n)
        void eval_deriv(double* new_x, double* new_dx, int new_n)

    cdef cppclass ZeroExtrapolation:
        pass

    cdef cppclass CuspExtrapolation:
        pass

    cdef cppclass PowerExtrapolation:
        PowerExtrapolation(double power)
        double get_power()

    cdef cppclass PotentialExtrapolation:
        PotentialExtrapolation(int64_t l) except +
        int64_t get_l()
        double get_amp_left()
        double get_amp_right()
