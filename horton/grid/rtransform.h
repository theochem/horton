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

// UPDATELIBDOCTITLE: Transformation from uniform 1D to non-uniform 1D grids

#ifndef HORTON_GRID_RTRANSFORM_H_
#define HORTON_GRID_RTRANSFORM_H_


/** @brief
      Transformation of radial grid used for integration in spherical coordinates.

    Some conventions
    - `radius` or `r` are the transformed values, a function of t.
    - `t` is used for the original grid
 */
class RTransform {
 public:
  /** @brief
        Create an RTransform object.

      @param npoint
        The number of points on the radial grid.
   */
  explicit RTransform(int npoint);

  //! Destructor
  virtual ~RTransform() {}

  /** @brief
        Compute the radius for a given t. (Forward transformation.)

      @param t
        The value on the untransformed grid.
   */
  virtual double radius(double t) = 0;

  /** @brief
        Compute the derivative of the radius w.r.t. t.

      @param t
        The value on the untransformed grid.
   */
  virtual double deriv(double t) = 0;

  /** @brief
        Compute the second derivative of the radius w.r.t. t.

      @param t
        The value on the untransformed grid.
   */
  virtual double deriv2(double t) = 0;

  /** @brief
        Compute the third derivative of the radius w.r.t. t.

      @param t
        The value on the untransformed grid.
   */
  virtual double deriv3(double t) = 0;

  /** @brief
        Compute t for a given radius. (Reverse transform.)

      @param r
        The value on the transformed grid.
   */
  virtual double inv(double r) = 0;

  /** @brief
        Compute the radii for a given array of t values. (Forward transformation.)

      @param t
        The array with values on the untransformed grid.

      @param r
        The output array with radii.

      @param n
        The number of elements in the array.
   */
  void radius_array(double* t, double* r, int n);

  /** @brief
        Compute derivatives of the radius for a given array of t values.

      @param t
        The array with values on the untransformed grid.

      @param d
        The output array with (d r)/(d t).

      @param n
        The number of elements in the array.
   */
  void deriv_array(double* t, double* d, int n);

  /** @brief
        Compute second derivatives of the radius for a given array of t values.

      @param t
        The array with values on the untransformed grid.

      @param d
        The output array with (d^2 r)/(d t^2).

      @param n
        The number of elements in the array.
   */
  void deriv2_array(double* t, double* d, int n);

  /** @brief
        Compute third derivatives of the radius for a given array of t values.

      @param t
        The array with values on the untransformed grid.

      @param d
        The output array with (d^3 r)/(d t^3).

      @param n
        The number of elements in the array.
   */
  void deriv3_array(double* t, double* d, int n);

  /** @brief
        Compute the t values for a given array radii. (Reverse transformation.)

      @param r
        The array with radii.

      @param t
        The output array with values on the untransformed grid.

      @param n
        The number of elements in the array.
   */
  void inv_array(double* r, double* t, int n);

  int get_npoint() {return npoint;}  //!< The number of grid points .

 protected:
  int npoint;  //!< The number of grid points.
};


/** @brief
      Identity transformation, mostly used for debugging. r(t) = t
 */
class IdentityRTransform : public RTransform {
 public:
  /** @brief
        Create an IdentityRTransform object.

      @param npoint
        The number of points on the radial grid.
   */
  explicit IdentityRTransform(int npoint): RTransform(npoint) {}

  virtual double radius(double t);  //!< Forward transformation
  virtual double deriv(double t);   //!< (d r)/(d t)
  virtual double deriv2(double t);  //!< (d^2 r)/(d t^2)
  virtual double deriv3(double t);  //!< (d^3 r)/(d t^3)
  virtual double inv(double r);     //!< Reverse transformation
};


/** @brief
      Linear transformation. r(t) = rmin + t*(rmax - rmin)/(npoint - 1)
 */
class LinearRTransform : public RTransform {
 public:
  /** @brief
        Create an LinearRTransform object.

      @param rmin
        First radius

      @param rmax
        Last radius

      @param npoint
        The number of points on the radial grid.
   */
  LinearRTransform(double rmin, double rmax, int npoint);

  virtual double radius(double t);    //!< Forward transformation
  virtual double deriv(double t);     //!< (d r)/(d t)
  virtual double deriv2(double t);    //!< (d^2 r)/(d t^2)
  virtual double deriv3(double t);    //!< (d^3 r)/(d t^3)
  virtual double inv(double r);       //!< Reverse transformation

  double get_rmin() {return rmin;}    //!< First radius
  double get_rmax() {return rmax;}    //!< Last radius
  double get_alpha() {return alpha;}  //!< Slope

 private:
  double rmin;   //!< First radius
  double rmax;   //!< Last radius
  double alpha;  //!< Slope
};


/** @brief
      Exponential transformation. r(t) = rmin*exp(t*log(rmax/rmin)/(npoint - 1))
 */
class ExpRTransform : public RTransform {
 public:
  /** @brief
        Create an ExpRTransform object.

      @param rmin
        First radius

      @param rmax
        Last radius

      @param npoint
        The number of points on the radial grid.
   */
  ExpRTransform(double rmin, double rmax, int npoint);

  virtual double radius(double t);   //!< Forward transformation
  virtual double deriv(double t);    //!< (d r)/(d t)
  virtual double deriv2(double t);   //!< (d^2 r)/(d t^2)
  virtual double deriv3(double t);   //!< (d^3 r)/(d t^3)
  virtual double inv(double r);      //!< Reverse transformation

  double get_rmin() {return rmin;}    //!< First radius
  double get_rmax() {return rmax;}    //!< Last radius
  double get_alpha() {return alpha;}  //!< Exponent

 private:
  double rmin;   //!< First radius
  double rmax;   //!< Last radius
  double alpha;  //!< Exponent
};


/** @brief
      Power transformation. r(t) = rmin*t^power  with  power = log(rmax/rmin)/log(npoint)
 */
class PowerRTransform : public RTransform {
 public:
  /** @brief
        Create a PowerRTransform object.

      @param rmin
        First radius

      @param rmax
        Last radius

      @param npoint
        The number of points on the radial grid.
   */
  PowerRTransform(double rmin, double rmax, int npoint);

  virtual double radius(double t);  //!< Forward transformation
  virtual double deriv(double t);   //!< (d r)/(d t)
  virtual double deriv2(double t);  //!< (d^2 r)/(d t^2)
  virtual double deriv3(double t);  //!< (d^3 r)/(d t^3)
  virtual double inv(double r);     //!< Reverse transformation

  double get_rmin() {return rmin;}    //!< First radius
  double get_rmax() {return rmax;}    //!< Last radius
  double get_power() {return power;}  //!< The power of t

 private:
  double rmin;   //!< First radius
  double rmax;   //!< Last radius
  double power;  //!< The power of t
};


/** @brief
      Hyperbolic transformation. r(t) = a*t/(1 - b*t)
 */
class HyperbolicRTransform : public RTransform {
 public:
  /** @brief
        Create a HyperbolicRTransform object.

      @param a
        Prefactor in r(t).

      @param b
        Parameter in denominator of r(t).

      @param npoint
        The number of points on the radial grid.
   */
  HyperbolicRTransform(double a, double b, int npoint);

  virtual double radius(double t);  //!< Forward transformation
  virtual double deriv(double t);   //!< (d r)/(d t)
  virtual double deriv2(double t);  //!< (d^2 r)/(d t^2)
  virtual double deriv3(double t);  //!< (d^3 r)/(d t^3)
  virtual double inv(double r);     //!< Reverse transformation

  double get_a() {return a;}  //!< Prefactor in r(t)
  double get_b() {return b;}  //!< Parameter in denominator of r(t)

 private:
  double a;  //!< Prefactor in r(t)
  double b;  //!< Parameter in denominator of r(t)
};

#endif  // HORTON_GRID_RTRANSFORM_H_
