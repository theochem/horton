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

#ifndef HORTON_GRID_CUBIC_SPLINE_H
#define HORTON_GRID_CUBIC_SPLINE_H


#include "rtransform.h"


void tridiag_solve(double* diag_low, double* diag_mid, double* diag_up, double* right, double* solution, int n);
void tridiagsym_solve(double* diag_mid, double* diag_up, double* right, double* solution, int n);
void compute_cubic_spline_int_weights(double* weights, int npoint);


class Extrapolation;


class CubicSpline {
    private:
        Extrapolation* ep;
        bool own_ep;

        RTransform* rtf;
        bool own_rtf;

        double first_x, last_x;
    public:
        double* y;
        double* d;
        int n;

        CubicSpline(double* y, double* d, Extrapolation* ep, RTransform* rtf, int n);
        ~CubicSpline();

        void eval(double* new_x, double* new_y, int new_n);
        void eval_deriv(double* new_x, double* new_d, int new_n);

        RTransform* get_rtransform() {return rtf;}
        double get_first_x() {return first_x;};
        double get_last_x() {return last_x;};
    };


class Extrapolation {
    public:
        Extrapolation() {};
        virtual ~Extrapolation() {};
        virtual void prepare(CubicSpline* cs) = 0;
        virtual double eval_left(double x) = 0;
        virtual double eval_right(double x) = 0;
        virtual double eval_deriv_left(double x) = 0;
        virtual double eval_deriv_right(double x) = 0;
    };


class ZeroExtrapolation : public Extrapolation {
    public:
        virtual void prepare(CubicSpline* cs);
        virtual double eval_left(double x);
        virtual double eval_right(double x);
        virtual double eval_deriv_left(double x);
        virtual double eval_deriv_right(double x);
    };


class ExponentialExtrapolation : public Extrapolation {
    private:
        double a0, b0, x0;

    public:
        virtual void prepare(CubicSpline* cs);
        virtual double eval_left(double x);
        virtual double eval_right(double x);
        virtual double eval_deriv_left(double x);
        virtual double eval_deriv_right(double x);
    };


#endif
