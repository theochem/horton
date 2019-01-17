// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2019 The HORTON Development Team
//
// This file is part of HORTON.
//
// HORTON is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// HORTON is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

// UPDATELIBDOCTITLE: Auxiliary functions

#ifndef HORTON_GBASIS_COMMON_H
#define HORTON_GBASIS_COMMON_H

#define MAX_SHELL_TYPE 7
#define MAX_NCART_CUMUL ((MAX_SHELL_TYPE+1)*(MAX_SHELL_TYPE+2)*(MAX_SHELL_TYPE+3))/6
#define MAX_NCART_CUMUL_D ((MAX_SHELL_TYPE+2)*(MAX_SHELL_TYPE+3)*(MAX_SHELL_TYPE+4))/6
#define MAX_NCART_CUMUL_DD ((MAX_SHELL_TYPE+3)*(MAX_SHELL_TYPE+4)*(MAX_SHELL_TYPE+5))/6

// Simple math stuff
long fac(long n);
long fac2(long n);
long binom(long n, long m);
long get_shell_nbasis(long shell_type);
long get_max_shell_type();
const double dist_sq(const double* r0, const double* r1);

// Auxiliary functions for Gaussian integrals
void compute_gpt_center(double alpha0, const double* r0, double alpha1, const double* r1, double gamma_inv, double* gpt_center);
double gpt_coeff(long k, long n0, long n1, double pa, double pb);
double gb_overlap_int1d(long n0, long n1, double pa, double pb, double gamma_inv);
void nuclear_attraction_helper(double* work_g, long n0, long n1, double pa, double pb, double pc, double gamma_inv);

// Auxiliary functions for r^alpha integrals

//! cit =t^m (2^(2i)/(2i +1)!)
double cit(int i, double t, int m);

//! j(j-1)(j-2)(j-3)...(j-n)
long jfac(int j, int n);

/** @brief Evaluas the taylor series for r^alpha integrals

    Sum_i t^i * Gamma(i+(alpha+3)/2)*(2^(2i)/(2i +1)!

    @param n
        Angular moment (the order of the derivative of the basic integral Gn in Alhrichs
        Phys. Chem. Chem. Phys., 8, 3072 (2006)). The maximum value implemented is n=10.

    @param alpha
        The power of r in the potential.

    @param t
        rho|p-q|^2

    @param prefac
        e^(-t)/rho^3/2 - This term helps the Taylor series to converge when t is a large
        number, the factor 1/2*sqrt(rho^alpha) was "replaced" and multiplied outside, at
        the end, in the laplace_of_potential function.
*/
double dtaylor(int n, double alpha, double t, double prefac);

#endif
