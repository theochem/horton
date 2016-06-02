// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2015 The HORTON Development Team
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

// UPDATELIBDOCTITLE: One-dimensional cubic splines (on uniform grids)

#ifndef HORTON_GRID_CUBIC_SPLINE_H_
#define HORTON_GRID_CUBIC_SPLINE_H_

#include <cstdint>

#include "horton/grid/rtransform.h"


void tridiagsym_solve(double* diag_mid, double* diag_up, double* right, double* solution, int n);
void solve_cubic_spline_system(double* y, double *d, int npoint);
void compute_cubic_spline_int_weights(double* weights, int npoint);


class Extrapolation;


class CubicSpline {
 private:
    Extrapolation* extrapolation;
    RTransform* rtf;
    double first_x, last_x;
 public:
    double* y;
    double* dt;
    int n;

    CubicSpline(double* y, double* dt, Extrapolation* extrapolation, RTransform* rtf, int n);

    void eval(double* new_x, double* new_y, int new_n);
    void eval_deriv(double* new_x, double* new_dx, int new_n);

    RTransform* get_rtransform() {return rtf;}
    double get_first_x() {return first_x;}   // position of first (transformed) grid point
    double get_last_x() {return last_x;}     // position of first (transformed) last point
    Extrapolation* get_extrapolation() {return extrapolation;}
};


class Extrapolation {
 public:
    Extrapolation() {}
    virtual ~Extrapolation() {}
    virtual void prepare(CubicSpline* cs) = 0;
    virtual double eval_left(double x) = 0;
    virtual double eval_right(double x) = 0;
    virtual double deriv_left(double x) = 0;
    virtual double deriv_right(double x) = 0;
    virtual bool has_tail() = 0;
};


class ZeroExtrapolation : public Extrapolation {
 public:
    virtual void prepare(CubicSpline* cs);
    virtual double eval_left(double x);
    virtual double eval_right(double x);
    virtual double deriv_left(double x);
    virtual double deriv_right(double x);
    virtual bool has_tail() {return false;}
};


class CuspExtrapolation : public Extrapolation {
 private:
    double a0, b0, x0;

 public:
    CuspExtrapolation() : a0(0), b0(0), x0(0) {}
    virtual void prepare(CubicSpline* cs);
    virtual double eval_left(double x);
    virtual double eval_right(double x);
    virtual double deriv_left(double x);
    virtual double deriv_right(double x);
    virtual bool has_tail() {return false;}
};


class PowerExtrapolation : public Extrapolation {
 private:
    double amp, power;

 public:
    explicit PowerExtrapolation(double power): amp(0), power(power) {}
    virtual void prepare(CubicSpline* cs);
    virtual double eval_left(double x);
    virtual double eval_right(double x);
    virtual double deriv_left(double x);
    virtual double deriv_right(double x);
    virtual bool has_tail() {return true;}
    double get_power() {return power;}
};


class PotentialExtrapolation : public Extrapolation {
 private:
    int64_t l;
    double amp_left, amp_right;

 public:
    explicit PotentialExtrapolation(int64_t l);
    virtual void prepare(CubicSpline* cs);
    virtual double eval_left(double x);
    virtual double eval_right(double x);
    virtual double deriv_left(double x);
    virtual double deriv_right(double x);
    virtual bool has_tail() {return true;}
    int64_t get_l() {return l;}
    double get_amp_left() {return amp_left;}
    double get_amp_right() {return amp_right;}
};


#endif // HORTON_GRID_CUBIC_SPLINE_H_
