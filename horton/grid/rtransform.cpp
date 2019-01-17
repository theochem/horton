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


//  #define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include <stdexcept>
#include "horton/grid/cubic_spline.h"
#include "horton/grid/rtransform.h"


/*
   RTransform
*/

RTransform::RTransform(int npoint) : npoint(npoint) {
  if (npoint < 2)
    throw std::domain_error("A radial grid consists of at least two points.");
}

void RTransform::radius_array(double* t, double* r, int n) {
  while (n > 0) {
    *r = radius(*t);
    n--;
    t++;
    r++;
  }
}

void RTransform::deriv_array(double* t, double* d, int n) {
  while (n > 0) {
    *d = deriv(*t);
    n--;
    t++;
    d++;
  }
}

void RTransform::deriv2_array(double* t, double* d, int n) {
  while (n > 0) {
    *d = deriv2(*t);
    n--;
    t++;
    d++;
  }
}

void RTransform::deriv3_array(double* t, double* d, int n) {
  while (n > 0) {
    *d = deriv3(*t);
    n--;
    t++;
    d++;
  }
}

void RTransform::inv_array(double* r, double* t, int n) {
  while (n > 0) {
    *t = inv(*r);
    n--;
    r++;
    t++;
  }
}


/*
   IdentityRTransform
*/

double IdentityRTransform::radius(double t) {
  return t;
}

double IdentityRTransform::deriv(double t) {
  return 1.0;
}

double IdentityRTransform::deriv2(double t) {
  return 0.0;
}

double IdentityRTransform::deriv3(double t) {
  return 0.0;
}

double IdentityRTransform::inv(double r) {
  return r;
}


/*
   LinearRTransform
*/


LinearRTransform::LinearRTransform(double rmin, double rmax, int npoint)
    : RTransform(npoint), rmin(rmin), rmax(rmax) {
  if (rmin >= rmax)
    throw std::domain_error("rmin must be below rmax.");
  alpha = (rmax-rmin)/(npoint-1);
}

double LinearRTransform::radius(double t) {
  return alpha*t + rmin;
}

double LinearRTransform::deriv(double t) {
  return alpha;
}

double LinearRTransform::deriv2(double t) {
  return 0.0;
}

double LinearRTransform::deriv3(double t) {
  return 0.0;
}

double LinearRTransform::inv(double r) {
  return (r-rmin)/alpha;
}


/*
   ExpRTransform
*/

ExpRTransform::ExpRTransform(double rmin, double rmax, int npoint)
    : RTransform(npoint), rmin(rmin), rmax(rmax) {
  if (rmin >= rmax)
    throw std::domain_error("rmin must be below rmax.");
  if ((rmin <= 0.0) || (rmax <= 0.0))
    throw std::domain_error("The minimum and maximum radii must be positive.");
  alpha = log(rmax/rmin)/(npoint-1);
}

double ExpRTransform::radius(double t) {
  return rmin*exp(t*alpha);
}

double ExpRTransform::deriv(double t) {
  return rmin*alpha*exp(t*alpha);
}

double ExpRTransform::deriv2(double t) {
  return rmin*alpha*alpha*exp(t*alpha);
}

double ExpRTransform::deriv3(double t) {
  return rmin*alpha*alpha*alpha*exp(t*alpha);
}

double ExpRTransform::inv(double r) {
  return log(r/rmin)/alpha;
}


/*
   PowerRTransform
*/

PowerRTransform::PowerRTransform(double rmin, double rmax, int npoint)
    : RTransform(npoint), rmin(rmin), rmax(rmax) {
  if (rmin >= rmax)
    throw std::domain_error("rmin must be below rmax.");
  if ((rmin <= 0.0) || (rmax <= 0.0))
    throw std::domain_error("The minimum and maximum radii must be positive.");
  power = (log(rmax) - log(rmin))/log(npoint);
  if (power < 2.0)
    throw std::domain_error("Power must be at least two for a decent intgration.");
}

double PowerRTransform::radius(double t) {
  return rmin*pow(t+1, power);
}

double PowerRTransform::deriv(double t) {
  return power*rmin*pow(t+1, power-1);
}

double PowerRTransform::deriv2(double t) {
  return power*(power-1)*rmin*pow(t+1, power-2);
}

double PowerRTransform::deriv3(double t) {
  return power*(power-1)*(power-2)*rmin*pow(t+1, power-3);
}

double PowerRTransform::inv(double r) {
  return pow(r/rmin, 1.0/power)-1;
}


/*
   HyperbolicRTransform
*/

HyperbolicRTransform::HyperbolicRTransform(double a, double b, int npoint)
    : RTransform(npoint), a(a), b(b) {
  if (a <= 0.0)
    throw std::domain_error("a must be strictly positive.");
  if (b <= 0.0)
    throw std::domain_error("b must be strictly positive.");
  if (b*(npoint-1) >= 1.0)
    throw std::domain_error("b*(npoint-1) must be smaller than one.");
}

double HyperbolicRTransform::radius(double t) {
  return a*t/(1-b*t);
}

double HyperbolicRTransform::deriv(double t) {
  double x = 1.0/(1-b*t);
  return a*x*x;
}

double HyperbolicRTransform::deriv2(double t) {
  double x = 1.0/(1-b*t);
  return 2.0*a*b*x*x*x;
}

double HyperbolicRTransform::deriv3(double t) {
  double x = 1.0/(1-b*t);
  return 6.0*a*b*b*x*x*x*x;
}

double HyperbolicRTransform::inv(double r) {
  return r/(a+b*r);
}
