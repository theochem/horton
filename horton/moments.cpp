// Horton is a development platform for electronic structure methods.
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


#include <stdexcept>


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
