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


#include <cmath>
#include <stdexcept>
#include "horton/gbasis/boys.h"


// The static arrays in boys_inc.cpp are pre-computed to full precision with
// sympy and mpmath. (See tools/boys.py.)

#include "horton/gbasis/boys_inc.cpp"
#if BOYS_MAX_M + 6!= BOYS_MAX_DATA
#error The file boys_inc.cpp is not consistent with the limits in horton/gbasis
#endif
#define SQRT_PI_D2      8.86226925452758013649e-1  // sqrt(pi)/2

static double boys_function_tail(long m, double t) {
    double result = SQRT_PI_D2/sqrt(t);
    for (long i=1; i <= m; i++)
        result *= 0.5*(2*i-1)/t;
    return result;
}


double boys_function(long m, double t) {
    // This implementation of the boys functions is based on the obara saika
    // paper. It is reasonably fast an very accurate. It just consumes a lot of
    // memory with the static arrays.
    if (m < 0 || m > BOYS_MAX_M || t < 0) {
        throw std::domain_error("Arguments to Boys function are outside the valid domain.");
    } else if (round(t*BOYS_RESOLUTION) >= (boys_sizes[m]-1)) {
        return boys_function_tail(m, t);
    } else {
        const int i = static_cast<int>((round(t*BOYS_RESOLUTION)));
        double t_delta = (i - t*BOYS_RESOLUTION)/BOYS_RESOLUTION;
        double result = boys_fn_data[m][i];
        double tmp = t_delta;
        result += boys_fn_data[m+1][i]*tmp;
        tmp *= t_delta/2.0;
        result += boys_fn_data[m+2][i]*tmp;
        tmp *= t_delta/3.0;
        result += boys_fn_data[m+3][i]*tmp;
        tmp *= t_delta/4.0;
        result += boys_fn_data[m+4][i]*tmp;
        tmp *= t_delta/5.0;
        result += boys_fn_data[m+5][i]*tmp;
        tmp *= t_delta/6.0;
        result += boys_fn_data[m+6][i]*tmp;
        return result;
    }
}

void boys_function_array(long mmax, double t, double *output) {
    // for (long m=0; m <= mmax; m++) {
    //   output[m] = boys_function(m, t);
    // }
    // return;

    if (mmax < 0 || mmax > BOYS_MAX_M || t < 0) {
        throw std::domain_error("Arguments to Boys function are outside the valid domain.");
    }
    // Precompute terms in taylor series without coefficients, if needed.
    const int i = static_cast<int>((round(t*BOYS_RESOLUTION)));
    double xs[6];
    if (i < (boys_sizes[mmax]-1)) {
        // The if clause relies on the fact that the size of the pre-computed
        // Boys function arrays increases with m.
        const double t_delta = (i - t*BOYS_RESOLUTION)/BOYS_RESOLUTION;
        xs[0] = t_delta;
        xs[1] = xs[0]*(t_delta/2.0);
        xs[2] = xs[1]*(t_delta/3.0);
        xs[3] = xs[2]*(t_delta/4.0);
        xs[4] = xs[3]*(t_delta/5.0);
        xs[5] = xs[4]*(t_delta/6.0);
    }
    // Start the computation of the approximation of the Boys function and its
    // derivative for large values of t, if needed.
    double tail = NAN;
    if (i >= (boys_sizes[0]-1)) {
        tail = SQRT_PI_D2/sqrt(t);
    }
    // Compute the Boys function and all the requested derivatives.
    for (long m=0; m <= mmax; m++) {
        if (i >= (boys_sizes[m]-1)) {
            if (m > 0) tail *= 0.5*(2*m-1)/t;
            output[m] = tail;
        } else {
            output[m] =
                boys_fn_data[m][i] +
                boys_fn_data[m+1][i]*xs[0] +
                boys_fn_data[m+2][i]*xs[1] +
                boys_fn_data[m+3][i]*xs[2] +
                boys_fn_data[m+4][i]*xs[3] +
                boys_fn_data[m+5][i]*xs[4] +
                boys_fn_data[m+6][i]*xs[5];
        }
    }
}
