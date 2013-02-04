// Horton is a Density Functional Theory program.
// Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
//
// This file is part of Horton.
//
// Horton is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// Horton is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include "evaluate.h"


long index_wrap(long i, long high) {
    // Get around compiler weirdness
    long result = i%high;
    if (result<0) result += high;
    return result;
}


void eval_spline_cube(CubicSpline* spline, double* center, double* output,
                      double* origin, Cell* grid_cell, long* shape,
                      long* pbc_active) {
    // Find the ranges for the triple loop
    double rcut = spline->get_last_x();
    long ranges_low[3], ranges_high[3];
    double oc[3], oc_frac[3];

    for (int i=2; i>=0; i--) {
        oc[i] = center[i] - origin[i];
    }
    grid_cell->to_frac(oc, oc_frac);
    for (int i=2; i>=0; i--) {
        double steps = rcut/grid_cell->get_rspacing(i);
        ranges_low[i] = floor(oc_frac[i] - steps);
        ranges_high[i] = ceil(oc_frac[i] + steps);
    }

#ifdef DEBUG
    printf("shape [%li,%li,%li]\n", shape[0], shape[1], shape[2]);
#endif

#ifdef DEBUG
    printf("ranges [%li,%li,%li] [%li,%li,%li]\n", ranges_low[0], ranges_low[1], ranges_low[2], ranges_high[0], ranges_high[1], ranges_high[2]);
#endif

    // Truncate ranges in case of non-periodic boundary conditions
    for (int i=2; i>=0; i--) {
        if (!pbc_active[i]) {
            if (ranges_low[i] < 0) {
                ranges_low[i] = 0;
            }
            if (ranges_high[i] >= shape[i]) {
                ranges_high[i] = shape[i]-1;
            }
        }
    }

#ifdef DEBUG
    printf("ranges [%li,%li,%li] [%li,%li,%li]\n", ranges_low[0], ranges_low[1], ranges_low[2], ranges_high[0], ranges_high[1], ranges_high[2]);
#endif

    // Run triple loop
    for (long i0=ranges_low[0]; i0 <= ranges_high[0]; i0++) {
        long i0_wrap = index_wrap(i0, shape[0]);
        for (long i1=ranges_low[1]; i1 <= ranges_high[1]; i1++) {
            long i1_wrap = index_wrap(i1, shape[1]);
            for (long i2=ranges_low[2]; i2 <= ranges_high[2]; i2++) {
                long i2_wrap = index_wrap(i2, shape[2]);

                // Compute the distance to the origin
                double frac[3], cart[3];
                frac[0] = i0;
                frac[1] = i1;
                frac[2] = i2;
                grid_cell->to_cart(frac, cart);
                double x = cart[0] - oc[0];
                double y = cart[1] - oc[1];
                double z = cart[2] - oc[2];
                double d = sqrt(x*x+y*y+z*z);

                // Evaluate spline if needed
                if (d < rcut) {
                    double s;
                    spline->eval(&d, &s, 1);

                    // Add to result
                    *(output + (i0_wrap*shape[1] + i1_wrap)*shape[2] + i2_wrap) += s;
                }
            }
        }
    }
}
