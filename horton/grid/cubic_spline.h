// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2016 The HORTON Development Team
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


/** @brief
        Construct and use cubic splines, with optional transformation of x-axis and
        extrapolation rules beyond spline range.

 */
class CubicSpline {
 private:
    Extrapolation* extrapolation;  //!< Describes how function goes outside spline range.
    RTransform* rtf;               //!< Transformation of the x-axis.
    double first_x;                //!< Transformed first grid point.
    double last_x;                 //!< Transformed last grid point.

 public:
    double* y;                     //!< The y-values of the spline at each grid point.
    double* dt;                    //!< The derivative towards the t-axis (uniform x-axis).
    int n;                         //!< The number of grid points.

    /** @brief
            Create a CubicSpline object.

        @param y
            An array with y-values on the grid.

        @param dt
            An array with derivatives of y toward the uniform x-axis (thus t-axis).

        @param extrapolation
            An instance of an Extrapolation subclass that describes the function outside
            the spline range.

        @param rtf
            The transformation from the t-axis to the x-axis.

        @param n
            The number of grid points.
     */
    CubicSpline(double* y, double* dt, Extrapolation* extrapolation, RTransform* rtf, int n);

    /** @brief
            Evaluate the cubic spline in a set of new x values.

        @param new_x
            An array with new x values.

        @param new_y
            An output array with the corresponding y values.

        @param new_n
            The number of values in new_x.
      */
    void eval(const double* new_x, double* new_y, int new_n);

    /** @brief
            Evaluate the derivative of the cubic spline (dy/dx) in a set of new x values.

        @param new_x
            An array with new x values.

        @param new_dx
            An output array with the corresponding dy/dx values.

        @param new_n
            The number of values in new_x.
      */
    void eval_deriv(const double* new_x, double* new_dx, int new_n);

    //! Return the transformation of the x-axis.
    RTransform* get_rtransform() {return rtf;}

    //! Return the first grid point on the x-axis (transformed coordinates).
    double get_first_x() {return first_x;}

    //! Return the last grid point on the x-axis (transformed coordinates).
    double get_last_x() {return last_x;}

    //! Return the extrapolation function used outside the spline range.
    Extrapolation* get_extrapolation() {return extrapolation;}
};


/** @brief
        Base class for extrapolation of cubic splines outside their range.
  */
class Extrapolation {
 public:
    Extrapolation() {}                          //!< Constructor.
    virtual ~Extrapolation() {}                 //!< Destructor.
    virtual void prepare(CubicSpline* cs) = 0;  //!< Derive parameters from spline.
    virtual double eval_left(double x) = 0;     //!< Compute extrapolation for low x.
    virtual double eval_right(double x) = 0;    //!< Derivative of extrapolation for low x.
    virtual double deriv_left(double x) = 0;    //!< Compute extrapolation for high x.
    virtual double deriv_right(double x) = 0;   //!< Derivative of extrapolation for high x.
    virtual bool has_tail() = 0;                //!< Returns true if right extrapolation is nonzero.
};


/** @brief
        No extrapolation at all. Returns zero outside the spline range.
  */
class ZeroExtrapolation : public Extrapolation {
 public:
    virtual void prepare(CubicSpline* cs);   //!< Derive parameters from spline.
    virtual double eval_left(double x);      //!< Compute extrapolation for low x.
    virtual double eval_right(double x);     //!< Derivative of extrapolation for low x.
    virtual double deriv_left(double x);     //!< Compute extrapolation for high x.
    virtual double deriv_right(double x);    //!< Derivative of extrapolation for high x.
    virtual bool has_tail() {return false;}  //!< Returns false because tail is zero.
};


/** @brief
        An extrapolation suitable for electronic densities. It uses a 1s Slater-type
        density function close to the nucleus.
  */
class CuspExtrapolation : public Extrapolation {
 private:
    double a0;  //!< The value of the spline at the first grid point.
    double b0;  //!< The exponent of the 1s Slater-type function for the cusp.
    double x0;  //!< The position of the first grid point.

 public:
    CuspExtrapolation() : a0(0), b0(0), x0(0) {}  //!< Constructor.
    virtual void prepare(CubicSpline* cs);   //!< Derive parameters from spline.
    virtual double eval_left(double x);      //!< Compute extrapolation for low x.
    virtual double eval_right(double x);     //!< Derivative of extrapolation for low x.
    virtual double deriv_left(double x);     //!< Compute extrapolation for high x.
    virtual double deriv_right(double x);    //!< Derivative of extrapolation for high x.
    virtual bool has_tail() {return false;}  //!< Returns false because tail is zero.
};


/** @brief
        An extrapolation with a polynomial of arbitrary order for high x values.

    The prefactor of the polynomial is such that the function remains continuous at the
    last grid point of the spline.
  */
class PowerExtrapolation : public Extrapolation {
 private:
    double amp;    //!< The prefactor of the polynomial.
    double power;  //!< The power of the polynomial.

 public:
    /** @brief
            Construct a PowerExtrapolation.

        @param power
            The power of the polynomial tail, i.e. extrapolation on the right.
      */
    explicit PowerExtrapolation(double power): amp(0), power(power) {}
    virtual void prepare(CubicSpline* cs);  //!< Derive parameters from spline
    virtual double eval_left(double x);     //!< Compute extrapolation for low x
    virtual double eval_right(double x);    //!< Derivative of extrapolation for low x
    virtual double deriv_left(double x);    //!< Compute extrapolation for high x
    virtual double deriv_right(double x);   //!< Derivative of extrapolation for high x
    virtual bool has_tail() {return true;}  //!< Returns true because if R**power tail.

    //! The power of the polynomial tail.
    double get_power() {return power;}
};


/** @brief
        An extrapolation suitable for solutions of the Poisson equation.

    The prefactor of the left and right polynomial are such that the function remains
    continuous at the last grid point of the spline. The power of R on the left and right
    side is consistent with the boundary conditions of the solutions of the Poisson
    equation.
  */
class PotentialExtrapolation : public Extrapolation {
 private:
    int64_t l;         //!< The angular momentum for which the potential is computed.
    double amp_left;   //!< The prefactor for the polynomial for low x.
    double amp_right;  //!< The prefactor for the polynomial for high x.

 public:
    /** @brief
            Construct a PotentialExtrapolation.

        @param l
            The angular momentum for which the Coulomb potential is generated.
      */
    explicit PotentialExtrapolation(int64_t l);
    virtual void prepare(CubicSpline* cs);  //!< Derive parameters from spline.
    virtual double eval_left(double x);     //!< Compute extrapolation for low x.
    virtual double eval_right(double x);    //!< Derivative of extrapolation for low x.
    virtual double deriv_left(double x);    //!< Compute extrapolation for high x.
    virtual double deriv_right(double x);   //!< Derivative of extrapolation for high x.
    virtual bool has_tail() {return true;}  //!< Returns true because if 1/R**(l+1) tail.

    //! The angular momentum of the Coulomb potential spline.
    int64_t get_l() {return l;}

    //! The prefactor for the polynomial for low x.
    double get_amp_left() {return amp_left;}
    //! The prefactor for the polynomial for high x.
    double get_amp_right() {return amp_right;}
};


#endif  // HORTON_GRID_CUBIC_SPLINE_H_
