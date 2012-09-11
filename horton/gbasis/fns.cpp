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

//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif
#include <cmath>
#include <cstring>
#include "cartpure.h"
#include "common.h"
#include "fns.h"
using namespace std;


GB1GridFn::GB1GridFn(long max_shell_type): GBCalculator(max_shell_type) {
    nwork = max_nbasis;
    work_cart = new double[nwork];
    work_pure = new double[nwork];
}

void GB1GridFn::reset(long _shell_type0, const double* _r0, const double* _point) {
    // TODO: add checks (to all reset functions) that _shell_type* is withing bounds
    shell_type0 = _shell_type0;
    r0 = _r0;
    point = _point;
    // We make use of the fact that a floating point zero consists of
    // consecutive zero bytes.
    memset(work_cart, 0, nwork*sizeof(double));
    memset(work_pure, 0, nwork*sizeof(double));
}

void GB1GridFn::add(double coeff, double alpha0, const double* scales0) {
    double pre, poly;
    pre = coeff*exp(-alpha0*dist_sq(r0, point));
    i1p.reset(abs(shell_type0));
    do {
#ifdef DEBUG
        printf("n=[%i,%i,%i]\n", i1p.n0[0], i1p.n0[1], i1p.n0[2]);
#endif
        // For now, simple and inefficient evaluation of polynomial.
        // TODO: make more efficient by moving evaluation of poly to reset
        poly = 1.0;
        for (long j=0; j<3; j++) {
            double tmp = point[j] - r0[j];
            for (long i=0; i<i1p.n0[j]; i++) poly *= tmp;
        }
        work_cart[i1p.ibasis0] += pre*scales0[i1p.ibasis0]*poly;
#ifdef DEBUG
        printf("poly=%f ibasis0=%i scale=%f work[]=%f\n", poly, i1p.ibasis0, scales0[i1p.ibasis0], work_cart[i1p.ibasis0]);
#endif
    } while (i1p.inc());
#ifdef DEBUG
    printf("\n");
#endif
}

void GB1GridFn::cart_to_pure() {
    /*
       The initial results are always stored in work_cart. The projection
       routine always outputs its result in work_pure. Once that is done,
       the pointers to both blocks are swapped such that the final result is
       always back in work_cart.
    */

    if (shell_type0 < -1) {
        cart_to_pure_low(work_cart, work_pure, -shell_type0,
            1, // anterior
            1 // posterior
        );
        swap_work();
    }
}
