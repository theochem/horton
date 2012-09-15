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


//#include <cstdlib>
#include <cmath>
#include <cstring>
#include "cubic_spline.h"

// Note: All cubic spline routines assume h=1.0 and x0=0.0

/*
   Solvers for symmetric and general tridiagonal systems of equations.
*/

void tridiag_solve(double* diag_low, double* diag_mid, double* diag_up,
                   double* right, double* solution, int n) {
    // Simple solver for generic tridiagonal system. The right hand side and the
    // upper diagonal get screwed during the operation.
    double denominator;
    int i;

    diag_up[0] /= diag_mid[0];
    right[0] /= diag_mid[0];
    for(i = 1; i < n; i++){
        denominator = diag_mid[i] - diag_up[i - 1]*diag_low[i-1];
        if (i < n-1) diag_up[i] /= denominator;
        right[i] = (right[i] - right[i - 1]*diag_low[i-1])/denominator;
    }

    solution[n - 1] = right[n - 1];
    for(i = n - 2; i >= 0; i--) {
        solution[i] = right[i] - diag_up[i]*solution[i + 1];
    }
}

void tridiagsym_solve(double* diag_mid, double* diag_up, double* right,
                      double* solution, int n) {
    // Simple solver for symmetric tridiagonal system. The right hand side and the
    // upper diagonal get screwed during the operation.
    double denominator, tmp1, tmp2;
    int i;

    tmp1 = diag_up[0];
    diag_up[0] /= diag_mid[0];
    right[0] /= diag_mid[0];
    for(i = 1; i < n; i++){
        denominator = (diag_mid[i] - diag_up[i - 1]*tmp1);
        if (i < n-1) {
            tmp2 = diag_up[i];
            diag_up[i] /= denominator;
        }
        right[i] = (right[i] - right[i - 1]*tmp1)/denominator;
        tmp1 = tmp2;
    }

    solution[n - 1] = right[n - 1];
    for(i = n - 2; i >= 0; i--) {
        solution[i] = right[i] - diag_up[i]*solution[i + 1];
    }
}

/*
   CubicSpline class.
*/

CubicSpline::CubicSpline(double* _y, double* _d, int _n): y(NULL), d(NULL), n(_n) {
    // Constructs the derivatives (*d) based on the *y values
    // n is the number of nodes, so there are n-1 splines between the nodes.

    // allocate arrays for our own copy of the data
    y = new double[n];
    d = new double[n];

    // always copy the y values
    memcpy((void*)y, (void*)_y, n*sizeof(double));

    if (_d==NULL) {
        // Make our own _d's
        double* diag_mid = new double[n]; // TODO: figure out if this is totally safe with respect to oom.
        double* diag_up = new double[n-1];
        double* right =  new double[n];

        // Setup the tri-diagonal system
        diag_mid[0] = 2.0;
        diag_up[0] = 1.0;
        right[0] = 3.0*(y[1]-y[0]);
        for (int i=1; i<n-1; i++) {
            diag_mid[i] = 4.0;
            diag_up[i] = 1.0;
            right[i] = 3.0*(y[i+1]-y[i-1]);
        }
        diag_mid[n-1] = 2.0;
        right[n-1] = 3.0*(y[n-1]-y[n-2]);

        tridiagsym_solve(diag_mid, diag_up, right, d, n);

        delete[] diag_mid;
        delete[] diag_up;
        delete[] right;
    } else {
        // copy the given d's
        memcpy(d, _d, n*sizeof(double));
    }
}

CubicSpline::~CubicSpline() {
    delete[] y;
    delete[] d;
}


void CubicSpline::eval(double* new_x, double* new_y, int new_n) {
    for (int i=0; i<new_n; i++) {
        if (*new_x < 0.0) {
            *new_y = y[0];
        } else if (*new_x <= n-1) {
            int j;
            double u, z;
            // find the index of the interval in which new_x[i] lies.
            j = (int)floor(*new_x);
            if (j==n - 1) j = n - 2;
            // do the interpolation
            u = *new_x - j;
            z = y[j+1] - y[j];
            *new_y = y[j] + u*(d[j] + u*(3*z - 2*d[j] - d[j+1] + u*(-2*z + d[j] + d[j+1])));
        }
        new_x++;
        new_y++;
    }
}

void CubicSpline::eval_deriv(double* new_x, double* new_d, int new_n) {
    for (int i=0; i<new_n; i++) {
        if (*new_x < 0.0) {
            *new_d = 0.0;
        } else if (*new_x <= n-1) {
            int j;
            double u, z;
            // find the index of the interval in which new_x[i] lies.
            j = (int)floor(*new_x);
            if (j==n - 1) j = n - 2;
            // do the interpolation
            u = *new_x - j;
            z = y[j+1] - y[j];
            *new_d = d[j] + u*(6*z - 4*d[j] - 2*d[j+1] + u*(-6*z + 3*d[j] + 3*d[j+1]));
        }
        new_x++;
        new_d++;
    }
}

void CubicSpline::eval_deriv2(double* new_x, double* new_d2, int new_n) {
    for (int i=0; i<new_n; i++) {
        if (*new_x < 0.0) {
            *new_d2 = 0.0;
        } else if (*new_x <= n-1) {
            int j;
            double u, z;
            // find the index of the interval in which new_x[i] lies.
            j = (int)floor(*new_x);
            if (j==n - 1) j = n - 2;
            // do the interpolation
            u = *new_x - j;
            z = y[j+1] - y[j];
            *new_d2 = 6*z - 4*d[j] - 2*d[j+1] + u*(-12*z + 6*d[j] + 6*d[j+1]);
        }
        new_x++;
        new_d2++;
    }
}

double CubicSpline::integrate() {
    double result;
    double* worky = y;
    double* workd = d;
    result = 0.5*(*worky) + (*workd)/12.0;
    worky++;
    workd++;
    for (int i=1; i<n-1; i++) {
        result += *worky;
        worky++;
        workd++;
    }
    result += 0.5*(*worky) - (*workd)/12.0;
    return result;
}
