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
#include <stdexcept>
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
    // the grid_cell must be 3D
    if (grid_cell->get_nvec() != 3)
        throw std::domain_error("eval_spline_cube only works for 3D grid cells.");

    // Find the ranges for the triple loop
    double rcut = spline->get_last_x();
    double delta[3];
    delta[0] = center[0] - origin[0];
    delta[1] = center[1] - origin[1];
    delta[2] = center[2] - origin[2];
    long ranges_begin[3], ranges_end[3];
    grid_cell->set_ranges_rcut(delta, rcut, 1, ranges_begin, ranges_end);

#ifdef DEBUG
    printf("shape [%li,%li,%li]\n", shape[0], shape[1], shape[2]);
#endif

#ifdef DEBUG
    printf("ranges [%li,%li,%li] [%li,%li,%li]\n", ranges_low[0], ranges_low[1], ranges_low[2], ranges_high[0], ranges_high[1], ranges_high[2]);
#endif

    // Truncate ranges in case of non-periodic boundary conditions
    for (int i=2; i>=0; i--) {
        if (!pbc_active[i]) {
            if (ranges_begin[i] < 0) {
                ranges_begin[i] = 0;
            }
            if (ranges_end[i] > shape[i]) {
                ranges_end[i] = shape[i];
            }
        }
    }

#ifdef DEBUG
    printf("ranges [%li,%li,%li] [%li,%li,%li]\n", ranges_begin[0], ranges_begin[1], ranges_begin[2], ranges_end[0], ranges_end[1], ranges_end[2]);
#endif

    // Run triple loop
    for (long i0 = ranges_begin[0]; i0 < ranges_end[0]; i0++) {
        long i0_wrap = index_wrap(i0, shape[0]);
        for (long i1 = ranges_begin[1]; i1 < ranges_end[1]; i1++) {
            long i1_wrap = index_wrap(i1, shape[1]);
            for (long i2 = ranges_begin[2]; i2 < ranges_end[2]; i2++) {
                long i2_wrap = index_wrap(i2, shape[2]);

                // Compute the distance to the origin
                double frac[3], cart[3];
                frac[0] = i0;
                frac[1] = i1;
                frac[2] = i2;
                grid_cell->to_cart(frac, cart);
                double x = cart[0] - delta[0];
                double y = cart[1] - delta[1];
                double z = cart[2] - delta[2];
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

void eval_spline_grid(CubicSpline* spline, double* center, double* output,
                      double* points, Cell* cell, long npoint) {
    double rcut = spline->get_last_x();

    while (npoint > 0) {
        // Find the ranges for the triple loop
        double delta[3];
        delta[0] = center[0] - points[0];
        delta[1] = center[1] - points[1];
        delta[2] = center[2] - points[2];
        long ranges_begin[3], ranges_end[3];
        cell->set_ranges_rcut(delta, rcut, 1, ranges_begin, ranges_end);

        for (int i=cell->get_nvec(); i < 3; i++) {
            ranges_begin[i] = 0;
            ranges_end[i] = 1;
        }

        // Run the triple loop
        for (long i0 = ranges_begin[0]; i0 < ranges_end[0]; i0++) {
            for (long i1 = ranges_begin[1]; i1 < ranges_end[1]; i1++) {
                for (long i2 = ranges_begin[2]; i2 < ranges_end[2]; i2++) {
                    // Compute the distance between the point and the image of the center
                    double frac[3], cart[3];
                    frac[0] = i0;
                    frac[1] = i1;
                    frac[2] = i2;
                    cell->to_cart(frac, cart);
                    double x = cart[0] - delta[0];
                    double y = cart[1] - delta[1];
                    double z = cart[2] - delta[2];
                    double d = sqrt(x*x+y*y+z*z);

                    // Evaluate spline if needed
                    if (d < rcut) {
                        double s;
                        spline->eval(&d, &s, 1);
                        *output += s;
                    }

                }
            }
        }

        // move on
        points += 3;
        output++;
        npoint--;
    }

#ifdef DEBUG
    printf("\n");
#endif

}
