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

// UPDATELIBDOCTITLE: Efficient evaluation of various of polynomials

#ifndef HORTON_MOMENTS_H
#define HORTON_MOMENTS_H

/** @brief Compute Cartesian mononomials efficiently.

    Computes all Cartesian mononomials (x, y, z) up to a given order in alphabetical
    order. The first three elements of output must contain x, y and z.

    @param output
        The output array, which must be sufficiently large.

    @param lmax
        The highest order of the mononomials.

    @returns The position of the first mononomial of the highest order.
  */
long fill_cartesian_polynomials(double* output, long lmax);

long fill_pure_polynomials(double* output, long lmax);
long fill_pure_polynomials_array(double* output, long lmax, long nrep, long stride);
void fill_radial_polynomials(double* output, long lmax);

#endif
