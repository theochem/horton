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
    if ((rmin <= 0.0) || (rmax <= 0.0))
        throw std::domain_error("The minimum and maximum radii of a log grid must be positive.");
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
   BakerRTransform
*/

BakerRTransform::BakerRTransform(double rmax, int npoint):
    RTransform(npoint), rmax(rmax)
{
    if (rmax <= 0.0)
        throw std::domain_error("The minimum and maximum radii of a log grid must be positive.");
    scale = (npoint-1.0)/npoint;
    scale = rmax/log(1-scale*scale);
}

double BakerRTransform::radius(double t) {
    double tmp = t/npoint;
    return scale*log(1-tmp*tmp);
}

double BakerRTransform::deriv(double t) {
    double tmp = t/npoint;
    return -scale*2.0*tmp/(1-tmp*tmp)/npoint;
}

double BakerRTransform::inv(double r) {
    return npoint*sqrt(1-exp(r/scale));
}
