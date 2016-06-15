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

// UPDATELIBDOCTITLE: Efficient evaluation of various of polynomials

#ifndef HORTON_MOMENTS_H
#define HORTON_MOMENTS_H

long fill_cartesian_polynomials(double* output, long lmax);
long fill_pure_polynomials(double* output, long lmax);
long fill_pure_polynomials_array(double* output, long lmax, long nrep, long stride);
void fill_radial_polynomials(double* output, long lmax);

#endif
