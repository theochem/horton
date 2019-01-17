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

// UPDATELIBDOCTITLE: The Boys function

#ifndef HORTON_GBASIS_BOYS_H_
#define HORTON_GBASIS_BOYS_H_

#include "horton/gbasis/common.h"
#define BOYS_MAX_M 4*MAX_SHELL_TYPE

/** @brief
        Compute the boys function.

    @param m
        The order parameters.

    @param t
        The rescaled distance between the two centers.
 */
double boys_function(long m, double t);

/** @brief
        Compute the boys function for a range of orders in one go.

    @param mmax
        The highest value of the order, for which the Boys function is to be computed.
        All orders for zero up ot this value (inclusive) are considered.
        The size of the output array has to be at least mmax+1.

    @param t
        The rescaled distance between the two centers.

    @param output
        The output array.
 */
void boys_function_array(long mmax, double t, double *output);

#endif  // HORTON_GBASIS_BOYS_H_
