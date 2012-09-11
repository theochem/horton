// Horton is a Density Functional Theory program.
// Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


#include "common.h"


long fac(long n) {
    long result = 1;
    while (n > 1) {
        result *= n;
        n--;
    }
    return result;
}


long fac2(long n) {
    long result = 1;
    while (n > 1) {
        result *= n;
        n -= 2;
    }
    return result;
}


long binom(long n, long m) {
    long numer = 1, denom = 1;
    while (n > m) {
        numer *= n;
        denom *= (n-m);
        n--;
    }
    return numer/denom;
}


long get_shell_nbasis(long shell_type) {
    if (shell_type > 0) {
        // Cartesian
        return ((shell_type+1)*(shell_type+2))/2;
    } else if (shell_type == -1) {
        // should not happen.
        return -1;
    } else {
        // Pure
        return -2*shell_type+1;
    }
}


long get_max_shell_type() {
    return MAX_SHELL_TYPE;
}


const double dist_sq(const double* r0, const double* r1) {
    double result, tmp;
    tmp = r0[0] - r1[0];
    result = tmp*tmp;
    tmp = r0[1] - r1[1];
    result += tmp*tmp;
    tmp = r0[2] - r1[2];
    result += tmp*tmp;
    return result;
}
