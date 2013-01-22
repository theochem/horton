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

/*

   The cubic splines in Horton are a bit special, hence this dedicated
   implementation. Unlike, the general cubic splines, Horton only uses splines
   with a 'standard' equidistant spacing of the grid points, with spacing 1. In
   order to handle data on other grids, one must provide a transformation object
   that defines the relation between that other grid and the standard grid.

   Why, oh why, not just do conventional cubic splines? Answer: more accuracy
   for the same cost. Consider, for example, an exponential decaying function.
   After transformation to a log grid, this function becomes linear and the
   spline becomes exact. This can not be acchieved without the transformation
   trick. Furthermore, cubic splines on equidistant grids can be implemented
   more efficiently than the general ones.

*/


//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include <cstring>
#include <stdexcept>
#include "cubic_spline.h"


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


void solve_cubic_spline_system(double* y, double *d, int npoint) {
    double* diag_mid = new double[npoint];
    double* diag_up = new double[npoint-1];
    double* right =  new double[npoint];

    // Setup the tri-diagonal system
    diag_mid[0] = 2.0;
    diag_up[0] = 1.0;
    right[0] = 3.0*(y[1]-y[0]);
    for (int i=1; i<npoint-1; i++) {
        diag_mid[i] = 4.0;
        diag_up[i] = 1.0;
        right[i] = 3.0*(y[i+1]-y[i-1]);
    }
    diag_mid[npoint-1] = 2.0;
    right[npoint-1] = 3.0*(y[npoint-1]-y[npoint-2]);

    tridiagsym_solve(diag_mid, diag_up, right, d, npoint);

    delete[] diag_mid;
    delete[] diag_up;
    delete[] right;
}


void compute_cubic_spline_int_weights(double* weights, int npoint) {
    // Construct a cubic spline for series 0,0,...,0,1,0...,0. One for every
    // grid point. The integral of such a function is the weight corresponding
    // to the grid point.

    // allocate arrays for the temporary splines
    double* y = new double[npoint];
    double* d = new double[npoint];

    // Set y array to zero. d is initialized in solve_cubic_spline_system.
    memset(y, 0, sizeof(double)*npoint);

    for (int ipoint=0; ipoint<npoint; ipoint++) {
        // setup the input spline
        if (ipoint > 0) y[ipoint-1] = 0.0;
        y[ipoint] = 1.0;
        // solve for the derivatives
        solve_cubic_spline_system(y, d, npoint);
        // compute the integral over the spline, which is heavily simplified...
        double weight = y[ipoint];
        if ((ipoint==0) || (ipoint==npoint-1)) {
            weight *= 0.5;
        }
        weight += (d[0] - d[npoint-1])/12.0;
        weights[ipoint] = weight;
    }

    delete[] y;
    delete[] d;
}

/*
   CubicSpline class.
*/

CubicSpline::CubicSpline(double* _y, double* _d, Extrapolation* _ep, RTransform* _rtf, int _n):
    ep(_ep), own_ep(false), rtf(_rtf), own_rtf(false), first_x(0.0),
    last_x(0.0), y(NULL), d(NULL), n(_n)
{
    // Constructs the derivatives (*d) based on the *y values
    // n is the number of nodes, so there are n-1 splines between the nodes.

    // allocate arrays for our own copy of the data
    y = new double[n];
    d = new double[n];

    // always copy the y values
    memcpy((void*)y, (void*)_y, n*sizeof(double));

    if (rtf == NULL) {
        rtf = new IdentityRTransform(n);
        own_rtf = true;
    } else if (rtf->get_npoint() != _n) {
        throw std::invalid_argument("The length of the cubic spline array data does not match the given rtransform.");
    }
    first_x = rtf->radius(0);
    last_x = rtf->radius(n-1);

    if (_d==NULL) {
        // Make our own d's
        solve_cubic_spline_system(y, d, n);
    } else {
        // Transform the given d's. We get dy/dr, we want dy/dt
        for (int i=0; i<n; i++) {
            d[i] = _d[i]*rtf->deriv(i);
//#ifdef DEBUG
//            printf("TF GIVEN DERIVS d[%i]=%f _d[%i]=%f\n", i, d[i], i, _d[i]);
//#endif
        }
    }

    if (ep == NULL) {
        ep = new ZeroExtrapolation();
        own_ep = true;
    }
    ep->prepare(this);
}

CubicSpline::~CubicSpline() {
    delete[] y;
    delete[] d;
    if (own_ep) delete ep;
    if (own_rtf) delete rtf;
}


void CubicSpline::eval(double* new_x, double* new_y, int new_n) {
    for (int i=0; i<new_n; i++) {
        if (*new_x < first_x) {
            // Left extrapolation
            *new_y = ep->eval_left(*new_x);
        } else if (*new_x <= last_x) {
            // Cubic Spline interpolation
            // 1) transform *new_x to t
            double t = rtf->inv(*new_x);
//#ifdef DEBUG
//            printf("X_TO_T x=%f t=%f\n", *new_x, t);
//#endif
            // 2) find the index of the interval in which t lies.
            int j = (int)floor(t);
            if (j==n - 1) j = n - 2;
            // 3) do the interpolation
            double u = t - j;
            double z = y[j+1] - y[j];
            *new_y = y[j] + u*(d[j] + u*(3*z - 2*d[j] - d[j+1] + u*(-2*z + d[j] + d[j+1])));
        } else {
            // Right extrapolation
            *new_y = ep->eval_right(*new_x);
        }
        new_x++;
        new_y++;
    }
}

void CubicSpline::eval_deriv(double* new_x, double* new_d, int new_n) {
    for (int i=0; i<new_n; i++) {
        if (*new_x < first_x) {
            // Left extrapolation
            *new_d = ep->eval_deriv_left(*new_x);
        } else if (*new_x <= last_x) {
            // Cubic Spline interpolation
            // 1) transform *new_x to t
            double t = rtf->inv(*new_x);
            // 2) find the index of the interval in which t lies.
            int j = (int)floor(t);
            if (j==n - 1) j = n - 2;
            // 3) do the interpolation
            double u = t - j;
            double z = y[j+1] - y[j];
            *new_d = d[j] + u*(6*z - 4*d[j] - 2*d[j+1] + u*(-6*z + 3*d[j] + 3*d[j+1]));
            // 4) transform the derivative, from dy/dt to dy/dr
            *new_d /= rtf->deriv(t);
        } else {
            // Right extrapolation
            *new_d = ep->eval_deriv_right(*new_x);
        }
        new_x++;
        new_d++;
    }
}


/*
   ZeroExtrapolation class
*/

void ZeroExtrapolation::prepare(CubicSpline* cs) {}
double ZeroExtrapolation::eval_left(double x) {return 0.0;}
double ZeroExtrapolation::eval_right(double x) {return 0.0;}
double ZeroExtrapolation::eval_deriv_left(double x) {return 0.0;}
double ZeroExtrapolation::eval_deriv_right(double x) {return 0.0;}


/*
   ExponentialExtrapolation class

   Only extrapolates for values smaller than lowest grid point. Extrapolation
   (of atomic densities) at large values is unreliable and waste of time.
*/

void ExponentialExtrapolation::prepare(CubicSpline* cs) {
    if (cs->d[0] == 0.0) {
        throw std::domain_error("The exponential extrapolation makes no sense when the derivative at the first point is zero.");
    }
    RTransform* rtf = cs->get_rtransform();
    a0 = cs->y[0];
    b0 = cs->d[0]/cs->y[0]/rtf->deriv(0);
    x0 = cs->get_first_x();
#ifdef DEBUG
    printf("PARS EXP EXTRAPOL a0=%f b0=%f x0=%f\n", a0, b0, x0);
#endif
}

double ExponentialExtrapolation::eval_left(double x) {
    return a0*exp(b0*(x-x0));
}

double ExponentialExtrapolation::eval_right(double x) {
    return 0.0;
}

double ExponentialExtrapolation::eval_deriv_left(double x) {
    return a0*b0*exp(b0*(x-x0));
}

double ExponentialExtrapolation::eval_deriv_right(double x) {
    return 0.0;
}
