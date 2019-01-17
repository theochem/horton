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

/*

   The cubic splines in HORTON are a bit special, hence this dedicated
   implementation. Unlike, the general cubic splines, HORTON only uses splines
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


// #define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include <cstring>
#include <stdexcept>
#include "horton/grid/cubic_spline.h"


/*
   Solver for symmetric tridiagonal systems of equations.
*/

void tridiagsym_solve(double* diag_mid, double* diag_up, double* right,
                      double* solution, int n) {
    // Simple solver for symmetric tridiagonal system. The right hand side and the
    // upper diagonal get screwed during the operation.
    double tmp1, tmp2;
    int i;

    tmp1 = diag_up[0];
    diag_up[0] /= diag_mid[0];
    right[0] /= diag_mid[0];
    for (i = 1; i < n; i++) {
        double denominator = (diag_mid[i] - diag_up[i - 1]*tmp1);
        if (i < n-1) {
            tmp2 = diag_up[i];
            diag_up[i] /= denominator;
        }
        right[i] = (right[i] - right[i - 1]*tmp1)/denominator;
        tmp1 = tmp2;
    }

    solution[n - 1] = right[n - 1];
    for (i = n - 2; i >= 0; i--) {
        solution[i] = right[i] - diag_up[i]*solution[i + 1];
    }
}


void solve_cubic_spline_system(double* y, double *dt, int npoint) {
    double* diag_mid = new double[npoint];
    double* diag_up = new double[npoint-1];
    double* right =  new double[npoint];

    // Setup the tri-diagonal system
    diag_mid[0] = 2.0;
    diag_up[0] = 1.0;
    right[0] = 3.0*(y[1]-y[0]);
    for (int i=1; i < npoint-1; i++) {
        diag_mid[i] = 4.0;
        diag_up[i] = 1.0;
        right[i] = 3.0*(y[i+1]-y[i-1]);
    }
    diag_mid[npoint-1] = 2.0;
    right[npoint-1] = 3.0*(y[npoint-1]-y[npoint-2]);

    tridiagsym_solve(diag_mid, diag_up, right, dt, npoint);

    delete[] diag_mid;
    delete[] diag_up;
    delete[] right;
}


void compute_cubic_spline_int_weights(double* weights, int npoint) {
    // Construct a cubic spline for series 0,0,...,0,1,0...,0. One for every
    // grid point. The integral of such a function is the weight corresponding
    // to the grid point.

    // allocate arrays for the temporary splines
    double y[npoint];
    double dt[npoint];

    // Set y array to zero. d is initialized in solve_cubic_spline_system.
    memset(y, 0, sizeof(double)*npoint);

    for (int ipoint=0; ipoint < npoint; ipoint++) {
        // setup the input spline
        if (ipoint > 0) y[ipoint-1] = 0.0;
        y[ipoint] = 1.0;
        // solve for the derivatives
        solve_cubic_spline_system(y, dt, npoint);
        // compute the integral over the spline, which is heavily simplified...
        double weight = y[ipoint];
        if ((ipoint == 0) || (ipoint == npoint-1)) {
            weight *= 0.5;
        }
        weight += (dt[0] - dt[npoint-1])/12.0;
        weights[ipoint] = weight;
    }
}

/*
   CubicSpline class.
*/

CubicSpline::CubicSpline(double* y, double* dt, Extrapolation* extrapolation,
                         RTransform* rtf, int n)
    : extrapolation(extrapolation), rtf(rtf), first_x(0.0), last_x(0.0),
      y(y), dt(dt), n(n) {
  first_x = rtf->radius(0);
  last_x = rtf->radius(n-1);
  extrapolation->prepare(this);
}


void CubicSpline::eval(const double* new_x, double* new_y, int new_n) {
    for (int i=0; i < new_n; i++) {
        if (*new_x < first_x) {
            // Left extrapolation
            *new_y = extrapolation->eval_left(*new_x);
        } else if (*new_x <= last_x) {
            // Cubic Spline interpolation
            // 1) transform *new_x to t
            double t = rtf->inv(*new_x);
            // 2) find the index of the interval in which t lies.
            int j = static_cast<int>(floor(t));
            if (j == n - 1) j = n - 2;
            // 3) do the interpolation
            double u = t - j;
            double z = y[j+1] - y[j];
            *new_y = y[j] + u*(dt[j] + u*(3*z - 2*dt[j] - dt[j+1] + u*(-2*z + dt[j] + dt[j+1])));
        } else {
            // Right extrapolation
            *new_y = extrapolation->eval_right(*new_x);
        }
        new_x++;
        new_y++;
    }
}

void CubicSpline::eval_deriv(const double* new_x, double* new_dx, int new_n) {
    for (int i=0; i < new_n; i++) {
        if (*new_x < first_x) {
            // Left extrapolation
            *new_dx = extrapolation->deriv_left(*new_x);
        } else if (*new_x <= last_x) {
            // Cubic Spline interpolation
            // 1) transform *new_x to t
            double t = rtf->inv(*new_x);
            // 2) find the index of the interval in which t lies.
            int j = static_cast<int>(floor(t));
            if (j == n - 1) j = n - 2;
            // 3) do the interpolation
            double u = t - j;
            double z = y[j+1] - y[j];
            *new_dx = dt[j] + u*(6*z - 4*dt[j] - 2*dt[j+1] + u*(-6*z + 3*dt[j] + 3*dt[j+1]));
            // 4) transform the derivative, from dy/dt to dy/dr
            *new_dx /= rtf->deriv(t);
        } else {
            // Right extrapolation
            *new_dx = extrapolation->deriv_right(*new_x);
        }
        new_x++;
        new_dx++;
    }
}


/*
   ZeroExtrapolation class
*/

void ZeroExtrapolation::prepare(CubicSpline* cs) {}
double ZeroExtrapolation::eval_left(double x) {return 0.0;}
double ZeroExtrapolation::eval_right(double x) {return 0.0;}
double ZeroExtrapolation::deriv_left(double x) {return 0.0;}
double ZeroExtrapolation::deriv_right(double x) {return 0.0;}


/*
   CuspExtrapolation class

   Only extrapolates for values smaller than lowest grid point. Extrapolation
   (of atomic densities) at large values is unreliable and waste of time.
*/

void CuspExtrapolation::prepare(CubicSpline* cs) {
    x0 = cs->get_first_x();
    if (fabs(cs->y[0]) < 0.00001) {
        // If there is no real cusp, don't care
        a0 = 0;
        b0 = 0;
    } else {
        RTransform* rtf = cs->get_rtransform();
        a0 = cs->y[0];
        b0 = cs->dt[0]/cs->y[0]/rtf->deriv(0);
    }
#ifdef DEBUG
    printf("PARS EXP EXTRAPOL a0=%f b0=%f x0=%f\n", a0, b0, x0);
#endif
}

double CuspExtrapolation::eval_left(double x) {
    return a0*exp(b0*(x-x0));
}

double CuspExtrapolation::eval_right(double x) {
    return 0.0;
}

double CuspExtrapolation::deriv_left(double x) {
    return a0*b0*exp(b0*(x-x0));
}

double CuspExtrapolation::deriv_right(double x) {
    return 0.0;
}

/*
   PowerExtrapolation class

   This is used for potentials that obey some power law at large distances.
*/

void PowerExtrapolation::prepare(CubicSpline* cs) {
    double x = cs->get_last_x();
    amp = cs->y[cs->n-1]*pow(x, -power);
}

double PowerExtrapolation::eval_left(double x) {
    return 0.0;
}

double PowerExtrapolation::eval_right(double x) {
    return amp*pow(x, power);
}

double PowerExtrapolation::deriv_left(double x) {
    return 0.0;
}

double PowerExtrapolation::deriv_right(double x) {
    return amp*power*pow(x, power-1);
}


/*
   PotentialExtrapolation class

   This is used for potentials that obey specific trends at short and long distances,
   depending at the angular momentum (l) for which they are computed.
*/

PotentialExtrapolation::PotentialExtrapolation(int64_t l)
    : l(l), amp_left(0.0), amp_right(0.0) {
    if (l < 0) {
        throw std::domain_error("The argument l cannot be negative.");
    }
}


void PotentialExtrapolation::prepare(CubicSpline* cs) {
    amp_left = cs->y[0]/pow(cs->get_first_x(), l);
    amp_right = cs->y[cs->n-1]*pow(cs->get_last_x(), l+1);
}

double PotentialExtrapolation::eval_left(double x) {
    return amp_left*pow(x, l);
}

double PotentialExtrapolation::eval_right(double x) {
    return amp_right/pow(x, l+1);
}

double PotentialExtrapolation::deriv_left(double x) {
    if (l == 0) {
        return 0.0;
    } else {
        return l*amp_left*pow(x, l-1);
    }
}

double PotentialExtrapolation::deriv_right(double x) {
    return -(l+1)*amp_right/pow(x, l+2);
}
