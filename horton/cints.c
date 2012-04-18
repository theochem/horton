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


#include <math.h>
#include "cints.h"



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


void compute_gpt_center(double exp0, double* r0, double exp1, double* r1, double gamma, double* gpt_center) {
    gpt_center[0] = (exp0*r0[0] + exp1*r1[0])/gamma;
    gpt_center[1] = (exp0*r0[1] + exp1*r1[1])/gamma;
    gpt_center[2] = (exp0*r0[2] + exp1*r1[2])/gamma;
}


double dist_sq(double* r0, double* r1) {
    double result, tmp;
    tmp = r0[0] - r1[0];
    result = tmp*tmp;
    tmp = r0[1] - r1[1];
    result += tmp*tmp;
    tmp = r0[2] - r1[2];
    result += tmp*tmp;
    return result;
}


double gpt_coeff(long k, long n0, long n1, double pa, double pb) {
    long i0, i1;
    double result;

    result = 0;
    i0 = k-n1;
    if (i0<0) i0 = 0;
    i1 = k - i0;
    do {
        result += binom(n0, i0)*binom(n1, i1)*pow(pa, n0-i0)*pow(pb, n1-i1);
        i0++;
        i1--;
    } while ((i1 >= 0)&&(i0 <= n0));
    return result;
}


double gob_overlap_int1d(long n0, long n1, double pa, double pb, double gamma) {
    long k, kmax;
    double result;

    kmax = (n0+n1)/2;
    result = 0.0;
    for (k=0; k<=kmax; k++) {
        result += fac2(2*k-1)*gpt_coeff(2*k, n0, n1, pa, pb)*pow(2*gamma,-k);
    }
    return sqrt(M_PI/gamma)*result;
}


/* TODO
       Generalize this such that any other two-body operator can be computed
       with the same number of arguments.
*/

double gob_overlap(double exp0, long nx0, long ny0, long nz0, double *r0,
                   double exp1, long nx1, long ny1, long nz1, double *r1) {
    double gpt_center[3];
    double gamma;

    gamma = exp0 + exp1;
    compute_gpt_center(exp0, r0, exp1, r1, gamma, gpt_center);
    return (
        exp(-exp0*exp1/gamma*dist_sq(r0, r1))*
        gob_overlap_int1d(nx0, nx1, gpt_center[0] - r0[0], gpt_center[0] - r1[0], gamma)*
        gob_overlap_int1d(ny0, ny1, gpt_center[1] - r0[1], gpt_center[1] - r1[1], gamma)*
        gob_overlap_int1d(nz0, nz1, gpt_center[2] - r0[2], gpt_center[2] - r1[2], gamma)
    );
}
