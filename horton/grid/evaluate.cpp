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

//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include <stdexcept>
#include "horton/moments.h"
#include "horton/grid/evaluate.h"


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
                    if ((d < rcut) || spline->get_extrapolation()->has_tail()) {
                        double s;
                        spline->eval(&d, &s, 1);
#ifdef DEBUG
                        printf("i=[%li,%li,%li] d=%f s=%f \n", i0, i1, i2, d, s);
#endif
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
}


#define INV_SQRT_4_PI 0.28209479177387814347 // 1/sqrt(4*pi)

void eval_decomposition_grid(CubicSpline** splines, double* center,
                             double* output, double* points, Cell* cell,
                             long nspline, long npoint) {
    long lmax = sqrt(nspline)-1;
    if ((lmax+1)*(lmax+1) != nspline) {
        throw std::domain_error("The number of splines does not match a well-defined lmax.");
    }

    double work[nspline-1];
    double rcut = splines[0]->get_last_x();

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

                    // Evaluate splines if needed
                    if ((d < rcut) || splines[0]->get_extrapolation()->has_tail()) {
                        // l == 0
                        double s;
                        splines[0]->eval(&d, &s, 1);
                        *output += s*INV_SQRT_4_PI;

                        if ((lmax > 0) && (d > 0)) {
                            // l > 0
                            work[0] = z;
                            work[1] = x;
                            work[2] = y;
                            if (lmax > 1) fill_pure_polynomials(work, lmax);

                            long counter = 0;
                            double dpowl = 1.0;
                            for (long l=1; l <= lmax; l++) {
                                dpowl /= d;
                                double factor = sqrt(2*l+1);
                                for (long m=-l; m<=l; m++) {
                                    splines[counter+1]->eval(&d, &s, 1);
                                    *output += s*factor*INV_SQRT_4_PI*dpowl*work[counter];
                                    counter++;
                                }
                            }
                        }
                    }
                }
            }
        }

        // move on
        points += 3;
        output++;
        npoint--;
    }
}
