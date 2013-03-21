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


void eval_spline_cube(CubicSpline* spline, double* center, double* output,
                      UniformIntGrid* ui_grid) {

    // Find the ranges for the triple loop
    double rcut = spline->get_last_x();
    long ranges_begin[3], ranges_end[3];
    ui_grid->set_ranges_rcut(center, rcut, ranges_begin, ranges_end);

    // Run triple loop
    const long* shape = ui_grid->get_shape();
    Range3Iterator r3i = Range3Iterator(ranges_begin, ranges_end, shape);
    long j[3], jwrap[3];
    for (long ipoint=r3i.get_npoint()-1; ipoint >=0; ipoint--) {
        r3i.set_point(ipoint, j, jwrap);
        double d = ui_grid->dist_grid_point(center, j);

        // Evaluate spline if needed
        if (d < rcut) {
            double s;
            spline->eval(&d, &s, 1);

            // Add to result
            *(output + (jwrap[0]*shape[1] + jwrap[1])*shape[2] + jwrap[2]) += s;
        }
    }
}

void eval_spline_grid(CubicSpline* spline, double* center, double* output,
                      double* points, Cell* cell, long npoint) {
    double rcut = spline->get_last_x();

    while (npoint > 0) {
        // Find the ranges for the triple loop
        double delta[3];
        delta[0] = points[0] - center[0];
        delta[1] = points[1] - center[1];
        delta[2] = points[2] - center[2];
        long ranges_begin[3], ranges_end[3];
        cell->set_ranges_rcut(delta, rcut, ranges_begin, ranges_end);

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
                    double x = cart[0] + delta[0];
                    double y = cart[1] + delta[1];
                    double z = cart[2] + delta[2];
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
