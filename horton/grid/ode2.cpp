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

//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif
#include <cstddef>
#include <cstdlib>
#include "horton/grid/ode2.h"


double get_value(const double* store, long index) {
    if (index == 0) return 0.0;
    if (index > 0) return store[index-1];
    return -store[-index-1];
}


const double overlap2_store[] = {
    9.0/70, 1.0/2, 13.0/420, 1.0/10, 6.0/5, 1.0/140, 1.0/60, 1.0/30, 13.0/35,
    11.0/210, 1.0/105, 2.0/15
};

const int overlap2_indexes[] = {
    1, 0, -2, 0, 3, 0, -4, 0, 2, 0, -5, 0, 4, 0, -4, 0, -3, 0, 4, 0, -6, 0, 7,
    0, -4, 0, 4, 0, -7, 0, -8, 0, 9, 9, 2, -2, -10, 10, 4, 4, 2, -2, 5, 5, -4,
    -4, -4, 4, -10, 10, -4, -4, 11, 11, 0, 0, 4, 4, -4, 4, 0, 0, 12, 12, 0, 1,
    0, 2, 0, -3, 0, -4, 0, -2, 0, -5, 0, 4, 0, 4, 0, 3, 0, 4, 0, -6, 0, -7, 0,
    -4, 0, -4, 0, 7, 0, -8
};

double hermite_overlap2(long xmax, long i0, bool deriv0, long i1, bool deriv1) {
    // process arguments
    long center0 = i0/2;
    bool kind0 = i0%2;
    long center1 = i1/2;
    bool kind1 = i1%2;

    // separation of the two basis functions
    long delta = center1 - center0;
    if (abs(delta) >= 2) return 0.0;

    // offset in the indexes matrix
    long offset = (((((delta+1)*2 + kind0)*2 + deriv0)*2 + kind1)*2 + deriv1)*2;

    // determine overlap
    double result = 0.0;
    if ((center0 >= 1) & (center0 <= xmax)) {
        // first segment
        result += get_value(overlap2_store, overlap2_indexes[offset]);
#ifdef DEBUG
        printf("first %i %f\n", overlap2_indexes[offset], get_value(overlap2_store, overlap2_indexes[offset]));
#endif
    }
    if ((center0 >= 0) and (center0 <= xmax-1)) {
        // second segment
        result += get_value(overlap2_store, overlap2_indexes[offset+1]);
#ifdef DEBUG
        printf("second %i %f\n", overlap2_indexes[offset+1], get_value(overlap2_store, overlap2_indexes[offset+1]));
#endif
    }
    return result;
}


const double overlap3_store[] = {
    9.0/140, 1.0/6, 1.0/72, 2.0/105, 3.0/5, 17.0/420, 4.0/35, 1.0/315, 1.0/168,
    13.0/420, 1.0/3, 5.0/84, 1.0/70, 54.0/35, 9.0/70, 6.0/35, 1.0/84, 1.0/140,
    3.0/35, 43.0/2520, 1.0/280, 1.0/105, 3.0/140, 1.0/1260, 1.0/840, 17.0/210,
    11.0/840, 1.0/60, 1.0/420, 43.0/140, 97.0/2520, 5.0/42, 2.0/315, 43.0/420,
    1.0/120, 2.0/35
};

const int overlap3_indexes[] = {
    1, 0, -2, 0, 3, 0, -4, 0, -2, 0, 5, 0, -6, 0, 7, 0, 3, 0, -6, 0, 8, 0, -9,
    0, -4, 0, 7, 0, -9, 0, 10, 0, 11, 0, -5, 0, 12, 0, 13, 0, -5, 0, 14, 0, -15,
    0, 16, 0, 12, 0, -15, 0, 17, 0, -18, 0, 13, 0, 16, 0, -18, 0, 19, 0, -20, 0,
    6, 0, -21, 0, 21, 0, 6, 0, -15, 0, 22, 0, -23, 0, -21, 0, 22, 0, -24, 0, 25,
    0, 21, 0, -23, 0, 25, 0, -9, 0, -26, 0, 7, 0, -27, 0, -28, 0, 7, 0, -16, 0,
    23, 0, 13, 0, -27, 0, 23, 0, -29, 0, -25, 0, -28, 0, 13, 0, -25, 0, -22, 0,
    1, 0, 2, 0, -3, 0, -4, 0, -11, 0, -5, 0, 12, 0, -13, 0, 20, 0, 6, 0, -21, 0,
    -21, 0, -26, 0, -7, 0, 27, 0, -28, 0, 2, 0, 5, 0, -6, 0, -7, 0, -5, 0, -14,
    0, 15, 0, 16, 0, 6, 0, 15, 0, -22, 0, -23, 0, -7, 0, -16, 0, 23, 0, -13, 0,
    -3, 0, -6, 0, 8, 0, 9, 0, 12, 0, 15, 0, -17, 0, -18, 0, -21, 0, -22, 0, 24,
    0, 25, 0, 27, 0, 23, 0, -29, 0, 25, 0, -4, 0, -7, 0, 9, 0, 10, 0, -13, 0,
    16, 0, -18, 0, -19, 0, -21, 0, -23, 0, 25, 0, 9, 0, -28, 0, -13, 0, 25, 0,
    -22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, -11, 0, 20, 0, -26, 0, 2, 0, -5, 0, 6, 0, -7, 0, -3, 0,
    12, 0, -21, 0, 27, 0, -4, 0, -13, 0, -21, 0, -28, 0, 2, 0, -5, 0, 6, 0, -7,
    0, 5, 0, -14, 0, 15, 0, -16, 0, -6, 0, 15, 0, -22, 0, 23, 0, -7, 0, 16, 0,
    -23, 0, -13, 0, -3, 0, 12, 0, -21, 0, 27, 0, -6, 0, 15, 0, -22, 0, 23, 0, 8,
    0, -17, 0, 24, 0, -29, 0, 9, 0, -18, 0, 25, 0, 25, 0, -4, 0, -13, 0, -21, 0,
    -28, 0, -7, 0, 16, 0, -23, 0, -13, 0, 9, 0, -18, 0, 25, 0, 25, 0, 10, 0,
    -19, 0, 9, 0, -22, 0, 30, 30, 11, -11, -31, 31, 32, 32, 11, -11, 5, 5, -12,
    -12, 13, -13, -31, 31, -12, -12, 33, 33, -9, 9, 32, 32, 13, -13, -9, 9, 34,
    34, 11, -11, 5, 5, -12, -12, 13, -13, 5, 5, 14, -14, -15, 15, -16, -16, -12,
    -12, -15, 15, 17, -17, 18, 18, 13, -13, -16, -16, 18, 18, 19, -19, -31, 31,
    -12, -12, 33, 33, -9, 9, -12, -12, -15, 15, 17, -17, 18, 18, 33, 33, 17,
    -17, -25, 25, 0, 0, -9, 9, 18, 18, 0, 0, -35, 35, 32, 32, 13, -13, -9, 9,
    34, 34, 13, -13, -16, -16, 18, 18, 19, -19, -9, 9, 18, 18, 0, 0, -35, 35,
    34, 34, 19, -19, -35, 35, 36, 36, 0, 1, 0, 11, 0, -20, 0, -26, 0, -2, 0, -5,
    0, 6, 0, 7, 0, 3, 0, 12, 0, -21, 0, -27, 0, -4, 0, 13, 0, 21, 0, -28, 0, -2,
    0, -5, 0, 6, 0, 7, 0, 5, 0, 14, 0, -15, 0, -16, 0, -6, 0, -15, 0, 22, 0, 23,
    0, 7, 0, 16, 0, -23, 0, 13, 0, 3, 0, 12, 0, -21, 0, -27, 0, -6, 0, -15, 0,
    22, 0, 23, 0, 8, 0, 17, 0, -24, 0, -29, 0, -9, 0, -18, 0, 25, 0, -25, 0, -4,
    0, 13, 0, 21, 0, -28, 0, 7, 0, 16, 0, -23, 0, 13, 0, -9, 0, -18, 0, 25, 0,
    -25, 0, 10, 0, 19, 0, -9, 0, -22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -2, 0, 3, 0, -4, 0, 11,
    0, -5, 0, 12, 0, 13, 0, -20, 0, 6, 0, -21, 0, 21, 0, -26, 0, 7, 0, -27, 0,
    -28, 0, -2, 0, 5, 0, -6, 0, 7, 0, -5, 0, 14, 0, -15, 0, 16, 0, 6, 0, -15, 0,
    22, 0, -23, 0, 7, 0, -16, 0, 23, 0, 13, 0, 3, 0, -6, 0, 8, 0, -9, 0, 12, 0,
    -15, 0, 17, 0, -18, 0, -21, 0, 22, 0, -24, 0, 25, 0, -27, 0, 23, 0, -29, 0,
    -25, 0, -4, 0, 7, 0, -9, 0, 10, 0, 13, 0, 16, 0, -18, 0, 19, 0, 21, 0, -23,
    0, 25, 0, -9, 0, -28, 0, 13, 0, -25, 0, -22, 0, 1, 0, 2, 0, -3, 0, -4, 0, 2,
    0, 5, 0, -6, 0, -7, 0, -3, 0, -6, 0, 8, 0, 9, 0, -4, 0, -7, 0, 9, 0, 10, 0,
    -11, 0, -5, 0, 12, 0, -13, 0, -5, 0, -14, 0, 15, 0, 16, 0, 12, 0, 15, 0,
    -17, 0, -18, 0, -13, 0, 16, 0, -18, 0, -19, 0, 20, 0, 6, 0, -21, 0, -21, 0,
    6, 0, 15, 0, -22, 0, -23, 0, -21, 0, -22, 0, 24, 0, 25, 0, -21, 0, -23, 0,
    25, 0, 9, 0, -26, 0, -7, 0, 27, 0, -28, 0, -7, 0, -16, 0, 23, 0, -13, 0, 27,
    0, 23, 0, -29, 0, 25, 0, -28, 0, -13, 0, 25, 0, -22
};


double hermite_overlap3(long xmax, long i0, bool deriv0, long i1, bool deriv1, long i2, bool deriv2) {
    // process arguments
    long center0 = i0/2;
    bool kind0 = i0%2;
    long center1 = i1/2;
    bool kind1 = i1%2;
    long center2 = i2/2;
    bool kind2 = i2%2;

    // separation of the basis functions
    long delta01 = center1 - center0;
    if (abs(delta01) >= 2) return 0.0;
    long delta02 = center2 - center0;
    if (abs(delta02) >= 2) return 0.0;

    // offset in the indexes matrix
    long offset = ((((((((delta01+1)*3 + (delta02+1))*2 + kind0)*2 + deriv0)*2 + kind1)*2 + deriv1)*2 + kind2)*2 + deriv2)*2;

    // determine overlap
    double result = 0.0;
    if ((center0 >= 1) & (center0 <= xmax)) {
        // first segment
        result += get_value(overlap3_store, overlap3_indexes[offset]);
    }
    if ((center0 >= 0) and (center0 <= xmax-1)) {
        // second segment
        result += get_value(overlap3_store, overlap3_indexes[offset+1]);
    }
    return result;
}


double hermite_node(long x, long center, bool kind, bool deriv) {
#ifdef DEBUG
    printf("node %i %i %i\n", kind, deriv, kind==deriv);
#endif
    if (kind==deriv) return (x==center);
    return 0.0;
}

double hermite_product2(long x, long i0, bool deriv0, long i1, bool deriv1) {
    return hermite_node(x, i0/2, i0%2, deriv0)*
           hermite_node(x, i1/2, i1%2, deriv1);
}


void build_ode2(double* b, double* a, double *f, double** bcs, double* coeffs,
                double* rhs, long npoint) {
    long nfn = 2*npoint;
    long skiprow[2];

    {
        int counter = 0;
        // Fill in the boundary conditions and determine which rows to skip
        if (bcs[0] != NULL) {
            coeffs[0] = 1;
            rhs[0] = *bcs[0];
            skiprow[counter] = 0;
            counter++;
        }
        if (bcs[1] != NULL) {
            coeffs[nfn+1] = 1;
            rhs[1] = *bcs[1];
            skiprow[counter] = 1;
            counter++;
        }
        if (bcs[2] != NULL) {
            coeffs[nfn*(nfn-1)-2] = 1;
            rhs[nfn-2] = *bcs[2];
            skiprow[counter] = nfn-2;
            counter++;
        }
        if (bcs[3] != NULL) {
            coeffs[nfn*nfn-1] = 1;
            rhs[nfn-1] = *bcs[3];
            skiprow[counter] = nfn-1;
            counter++;
        }
    }

    // Loop over all equations, except for boundary contitions
    for (long irow=0; irow<nfn; irow++) {
        // skip rows associated with boundary conditions
        if (irow==skiprow[0]) continue;
        if (irow==skiprow[1]) continue;
        // range of columns to consider. all others are zero due to short-ranged
        // overlap function.
        long begin, end;
        begin = 2*(irow/2-1);
        if (begin < 0) begin = 0;
        end = 2*(irow/2+2);
        if (end > nfn) end = nfn;

        // Loop over all coefficients
        for (long icol=begin; icol<end; icol++) {
            // the output coefficient for this iteration
            double* coeff = coeffs + (irow*nfn + icol);

            // add contribution from second order term
            *coeff -= hermite_overlap2(npoint-1, irow, true, icol, true);

            // Loop over all segments of input function b or a
            for (long ifoo=begin; ifoo<end; ifoo++) {
                // add contribution from first order term
                if (b[ifoo] != 0) {
                    *coeff += b[ifoo]*hermite_overlap3(npoint-1, irow, false, icol, true, ifoo, false);
                }

                // add contribution from zeroth order term
                if (a[ifoo] != 0) {
                    *coeff += a[ifoo]*hermite_overlap3(npoint-1, irow, false, icol, false, ifoo, false);
                }

            }

            // add contribution due to integration by parts
            *coeff += hermite_product2(npoint-1, irow, false, icol, true);
            *coeff -= hermite_product2(0, irow, false, icol, true);
        }

        // Loop over all segments of input function f
        double* output = rhs + irow;
        for (long ifoo=begin; ifoo<end; ifoo++) {
            // add contribution to rhs
            if (f[ifoo] != 0) {
                *output += f[ifoo]*hermite_overlap2(npoint-1, irow, false, ifoo, false);
            }
        }
    }
}
