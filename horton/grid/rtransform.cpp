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


#include <cmath>
#include <stdexcept>
#include "cubic_spline.h"
#include "rtransform.h"


// TODO: get rid of Base prefix.

/*
   BaseRTransform
*/

BaseRTransform::BaseRTransform(int npoint): npoint(npoint) {
    if (npoint < 2)
        throw std::domain_error("A radial grid consists of at least two points.");
}

void BaseRTransform::radius_array(double* t, double* r, int n) {
    while (n>0) {
        *r = radius(*t);
        n--;
        t++;
        r++;
    }
}

void BaseRTransform::deriv_array(double* t, double* d, int n) {
    while (n>0) {
        *d = deriv(*t);
        n--;
        t++;
        d++;
    }
}

void BaseRTransform::inv_array(double* r, double* t, int n) {
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
    BaseRTransform(npoint), rmin(rmin), rmax(rmax)
{
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
   LogRTransform
*/

LogRTransform::LogRTransform(double rmin, double rmax, int npoint):
    BaseRTransform(npoint), rmin(rmin), rmax(rmax)
{
    if ((rmin <= 0.0) || (rmax <= 0.0))
        throw std::domain_error("The minimum and maximum radii of a log grid must be positive.");
    alpha = log(rmax/rmin)/(npoint-1);
}

double LogRTransform::radius(double t) {
    return rmin*exp(t*alpha);
}

double LogRTransform::deriv(double t) {
    return rmin*alpha*exp(t*alpha);
}

double LogRTransform::inv(double r) {
    return log(r/rmin)/alpha;
}
