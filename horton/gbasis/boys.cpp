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


#include <cmath>
#include <stdexcept>
#include "horton/gbasis/boys.h"


// The static arrays in boys_inc.cpp are pre-computed to full precision with
// sympy and mpmath. (See tools/boys.py.)

#include "horton/gbasis/boys_inc.cpp"
#if BOYS_MAX_M + 6!= BOYS_MAX_DATA
#error The file boys_inc.cpp is not consistent with the limits in horton/gbasis
#endif
#define SQRT_PI_D2      8.86226925452758013649e-1 // sqrt(pi)/2

static double boys_function_tail(long m, double t) {
    double result = 1.0;
    for (long i=1; i<=m; i++)
        result *= 0.5*(2*i-1)/t;
    return result*SQRT_PI_D2/sqrt(t);
}


double boys_function(long m, double t) {
    // This implementation of the boys functions is based on the obara saika
    // paper. It is reasonably fast an very accurate. It just consumes a lot of
    // memory with the static arrays.
    if (m < 0 || m > BOYS_MAX_M || t < 0) {
        throw std::domain_error("Arguments to Boys function are outside the valid domain.");
    } else if (round(t*BOYS_RESOLUTION)>=(boys_sizes[m]-1)) {
        return boys_function_tail(m, t);
    } else {
        int i = (int)(round(t*BOYS_RESOLUTION));
        double t_delta = (i - t*BOYS_RESOLUTION)/BOYS_RESOLUTION;
        double result = 0.0;
        double tmp = 1;   result += boys_fn_data[m][i];
        tmp *= t_delta;   result += boys_fn_data[m+1][i]*tmp;
        tmp *= t_delta/2; result += boys_fn_data[m+2][i]*tmp;
        tmp *= t_delta/3; result += boys_fn_data[m+3][i]*tmp;
        tmp *= t_delta/4; result += boys_fn_data[m+4][i]*tmp;
        tmp *= t_delta/5; result += boys_fn_data[m+5][i]*tmp;
        tmp *= t_delta/6; result += boys_fn_data[m+6][i]*tmp;
        return result;
    }
}
