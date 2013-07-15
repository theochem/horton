// Horton is a development platform for electronic structure methods.
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


//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include <stdexcept>
#include "cubic_spline.h"
#include "rtransform.h"


/*
   RTransform
*/

RTransform::RTransform(int npoint): npoint(npoint) {
    if (npoint < 2)
        throw std::domain_error("A radial grid consists of at least two points.");
}

void RTransform::radius_array(double* t, double* r, int n) {
    while (n>0) {
        *r = radius(*t);
        n--;
        t++;
        r++;
    }
}

void RTransform::deriv_array(double* t, double* d, int n) {
    while (n>0) {
        *d = deriv(*t);
        n--;
        t++;
        d++;
    }
}

void RTransform::inv_array(double* r, double* t, int n) {
    while (n>0) {
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

double IdentityRTransform::inv(double r) {
    return r;
}


/*
   LinearRTransform
*/


LinearRTransform::LinearRTransform(double rmin, double rmax, int npoint):
    RTransform(npoint), rmin(rmin), rmax(rmax)
{
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

double LinearRTransform::inv(double r) {
    return (r-rmin)/alpha;
}


/*
   ExpRTransform
*/

ExpRTransform::ExpRTransform(double rmin, double rmax, int npoint):
    RTransform(npoint), rmin(rmin), rmax(rmax)
{
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

double ExpRTransform::inv(double r) {
    return log(r/rmin)/alpha;
}


/*
   ShiftedExpRTransform
*/

ShiftedExpRTransform::ShiftedExpRTransform(double rmin, double rshift, double rmax, int npoint):
    RTransform(npoint), rmin(rmin), rshift(rshift), rmax(rmax)
{
    if (rmin >= rmax)
        throw std::domain_error("rmin must be below rmax.");
    if ((rmin <= 0.0) || (rmax <= 0.0))
        throw std::domain_error("The minimum and maximum radii must be positive.");
    r0 = rmin + rshift;
    if (r0 <= 0.0)
        throw std::domain_error("The parameter r0 must be positive.");
    alpha = log((rmax+rshift)/r0)/(npoint-1);
}

double ShiftedExpRTransform::radius(double t) {
    return r0*exp(t*alpha) - rshift;
}

double ShiftedExpRTransform::deriv(double t) {
    return r0*alpha*exp(t*alpha);
}

double ShiftedExpRTransform::inv(double r) {
    return log((r + rshift)/r0)/alpha;
}


/*
   PowerExpRTransform
*/

PowerExpRTransform::PowerExpRTransform(double alpha, double rmax, double power, int npoint):
    RTransform(npoint), alpha(alpha), rmax(rmax), power(power)
{
    if (rmax <= 0.0)
        throw std::domain_error("The maximum radius must be positive.");
    if (power < 2.0)
        throw std::domain_error("Power must be at least two for a decent intgration.");
    if (alpha != 0) {
        amp = rmax/pow(fabs(expm1(alpha*npoint)), power);
    } else {
        amp = rmax/pow(npoint, power);
    }
}

double PowerExpRTransform::radius(double t) {
    if (fabs(alpha) < 1e-15) {
        return amp*pow(t+1, power);
    } else {
        return amp*pow(fabs(expm1(alpha*(t+1))), power);
    }
}

#define sign(x) ((x>0)-(x<0))

double PowerExpRTransform::deriv(double t) {
    if (fabs(alpha) < 1e-15) {
        return power*amp*pow(t+1, power-1);
    } else {
        double tmp = expm1(alpha*(t+1));
        return power*amp*pow(fabs(tmp), power-1)*sign(tmp)*exp(alpha*(t+1))*alpha;
    }
}

double PowerExpRTransform::inv(double r) {
    if (fabs(alpha) < 1e-15) {
        return pow(r/amp, 1.0/power)-1;
    } else {
        return log1p(pow(r/amp, 1.0/power))/alpha-1;
    }
}


/*
   BakerRTransform
*/

BakerRTransform::BakerRTransform(double rmax, int npoint):
    RTransform(npoint), rmax(rmax)
{
    if (rmax <= 0.0)
        throw std::domain_error("The maximum radius must be positive.");
    scale = npoint/(npoint+1.0);
    scale = rmax/log(1-scale*scale);
}

double BakerRTransform::radius(double t) {
    double tmp = (t+1)/(npoint+1);
    return scale*log(1-tmp*tmp);
}

double BakerRTransform::deriv(double t) {
    double tmp = (t+1)/(npoint+1);
    return -scale*2.0*tmp/(1-tmp*tmp)/(npoint+1);
}

double BakerRTransform::inv(double r) {
    return (npoint+1)*sqrt(1-exp(r/scale))-1;
}
