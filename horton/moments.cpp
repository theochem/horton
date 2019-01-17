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


#include <stdexcept>
#include <cmath>
#include "horton/moments.h"


long fill_cartesian_polynomials(double* output, long lmax) {
    // Shell l=0
    if (lmax <= 0) return -1;

    // Shell l=1
    if (lmax <= 1) return 0;

    // Shell l>1
    int old_offset = 0; // first array index of the moments of the previous shell
    int old_ncart = 3;  // number of moments in previous shell
    for (int l = 2; l <= lmax; l++) {
        // numbers for current iteration
        int new_ncart = old_ncart + l+1;
        int new_offset = old_offset + old_ncart;
        // start by multiplying all the old ones by x
        for (int i=0; i < old_ncart; i++) {
            output[new_offset+i] = output[0]*output[old_offset+i];
        }
        // then multiply the old ones without an x, by y.
        for (int i=0; i < l; i++) {
            output[new_offset+old_ncart+i] = output[1]*output[new_offset-l+i];
        }
        // finally multiply the last old one by z
        output[new_offset+new_ncart-1] = output[2]*output[new_offset-1];
        // translate new to old numbers
        old_ncart = new_ncart;
        old_offset = new_offset;
    }
    return old_offset;
}


long fill_pure_polynomials(double* output, long lmax) {
    // Shell l=0
    if (lmax <= 0) return -1;

    // Shell l=1
    if (lmax <= 1) return 0;

    // Shell l>1

    // auxiliary variables
    double z = output[0];
    double x = output[1];
    double y = output[2];
    double r2 = x*x + y*y + z*z;

    // work arrays in which the PI(z,r) polynomials are stored.
    double pi_old[lmax+1];
    double pi_new[lmax+1];
    double a[lmax+1];
    double b[lmax+1];

    // initialize work arrays
    pi_old[0] = 1;
    pi_new[0] = z;
    pi_new[1] = 1;
    a[1] = x;
    b[1] = y;

    int old_offset = 0; // first array index of the moments of the previous shell
    int old_npure = 3;  // number of moments in previous shell
    for (int l = 2; l <= lmax; l++) {
        // numbers for current iteration
        int new_npure = old_npure + 2;
        int new_offset = old_offset + old_npure;

        // construct polynomials PI(z,r) for current l
        double factor = (2*l-1);
        for (int m=0; m<=l-2; m++) {
            double tmp = pi_old[m];
            pi_old[m] = pi_new[m];
            pi_new[m] = (z*factor*pi_old[m] - r2*(l+m-1)*tmp)/(l-m);
        }
        pi_old[l-1] = pi_new[l-1];
        pi_new[l] = factor*pi_old[l-1];
        pi_new[l-1] = z*pi_new[l];

        // construct new polynomials A(x,y) and B(x,y)
        a[l] = x*a[l-1] - y*b[l-1];
        b[l] = x*b[l-1] + y*a[l-1];

        // construct solid harmonics
        output[new_offset] = pi_new[0];
        factor = M_SQRT2;
        for (int m=1; m<=l; m++) {
            factor /= sqrt((l+m)*(l-m+1));
            output[new_offset+2*m-1] = factor*a[m]*pi_new[m];
            output[new_offset+2*m] = factor*b[m]*pi_new[m];
        }

        // translate new to old numbers
        old_npure = new_npure;
        old_offset = new_offset;
    }
    return old_offset;
}


long fill_pure_polynomials_array(double* output, long lmax, long nrep, long stride) {
    long result = 0;
    for (long irep=0; irep<nrep; irep++) {
        result = fill_pure_polynomials(output, lmax);
        output += stride;
    }
    return result;
}


void fill_radial_polynomials(double* output, long lmax) {
    // Shell l=0
    if (lmax <= 0) return;

    // Shell l=1
    if (lmax <= 1) return;

    // Shell l>1
    int counter = 0;
    while (lmax > 1) {
        output[counter+1] = output[0]*output[counter];
        lmax -= 1;
        counter += 1;
    }
}
