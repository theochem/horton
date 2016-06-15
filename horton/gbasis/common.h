// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2016 The HORTON Development Team
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

#endif
