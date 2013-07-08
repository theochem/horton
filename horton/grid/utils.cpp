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
#include "utils.h"
#include "moments.h"


double dot_multi(long npoint, long nvector, double** data) {
    double result = 0.0;

    #pragma omp parallel for reduction(+:result)
    for (long ipoint=npoint-1; ipoint>=0; ipoint--) {
        double tmp = data[nvector-1][ipoint];
        for (long ivector=nvector-2; ivector>=0; ivector--) {
           tmp *= data[ivector][ipoint];
        }
#ifdef DEBUG
        printf("i=%i  tmp=%f\n", ipoint, tmp);
#endif
        result += tmp;
    }
    return result;
}


double intexp(double base, long exp) {
    double result = 1.0;
    while (exp > 0) {
        if (exp%2 == 1) result *= base;
        base *= base;
        exp /= 2;
    }
    return result;
}


// TODO: eliminate duplicate code

void dot_multi_moments_cube(long nvector, double** data, UniformGrid* ugrid, double* center, long lmax, long mtype, double* output, long nmoment) {
    if (ugrid->get_cell()->get_nvec() != 0) {
        throw std::domain_error("dot_multi_moments_cube only works for non-periodic grids.");
    }
    if ((lmax<0) || (lmax>7)) { // TODO: refer to some constant for maximum lmax
        throw std::domain_error("lmax should be in the range [0,7].");
    }
    if ((mtype != 1) && (mtype != 3)) {
        throw std::domain_error("mtype should be 1 or 3.");
    }

    // reset the output to zero

    Cube3Iterator c3i = Cube3Iterator(NULL, ugrid->get_shape());
    for (long ipoint=c3i.get_npoint()-1; ipoint >= 0; ipoint--) {
        // do the usual product of integranda
        double term = data[nvector-1][ipoint];
        for (long ivector=nvector-2; ivector>=0; ivector--)
           term *= data[ivector][ipoint];
        output[0] += term;

        if (lmax > 0) {
            // construct relative vector
            long j[3];
            c3i.set_point(ipoint, j);
            double delta[3];
            delta[0] = center[0];
            delta[1] = center[1];
            delta[2] = center[2];
            ugrid->delta_grid_point(delta, j);

            // evaluate polynomials in work array
            double work[nmoment-1];
            if (mtype==1) {
                work[0] = delta[0];
                work[1] = delta[1];
                work[2] = delta[2];
                fill_cartesian_polynomials(work, lmax);
            } else if (mtype==3) {
                work[0] = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
                fill_radial_polynomials(work, lmax);
            }

            // add product of polynomial and integrand to output
            for (long imoment=1; imoment < nmoment; imoment++) {
                output[imoment] += term*work[imoment-1];
            }
        }
    }
}

void dot_multi_moments(long npoint, long nvector, double** data, double* points, double* center, long lmax, long mtype, double* output, long nmoment) {
    if ((lmax<0) || (lmax>7)) { // TODO: refer to some constant for maximum lmax
        throw std::domain_error("lmax should be in the range [0,7].");
    }
    if ((mtype != 1) && (mtype != 3)) {
        throw std::domain_error("mtype should be 1 or 3.");
    }

    for (long ipoint=npoint-1; ipoint >= 0; ipoint--) {
        // do the usual product of integranda
        double term = data[nvector-1][ipoint];
        for (long ivector=nvector-2; ivector>=0; ivector--)
           term *= data[ivector][ipoint];

        output[0] += term;

        if (lmax > 0) {
            // construct relative vector
            double delta[3];
            delta[0] = points[ipoint*3  ] - center[0];
            delta[1] = points[ipoint*3+1] - center[1];
            delta[2] = points[ipoint*3+2] - center[2];

            // evaluate polynomials in work array
            double work[nmoment-1];
            if (mtype==1) {
                work[0] = delta[0];
                work[1] = delta[1];
                work[2] = delta[2];
                fill_cartesian_polynomials(work, lmax);
            } else if (mtype==3) {
                work[0] = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
                fill_radial_polynomials(work, lmax);
            }

            // add product of polynomial and integrand to output
            for (long imoment=1; imoment < nmoment; imoment++) {
                output[imoment] += term*work[imoment-1];
            }
        }
    }
}

void dot_multi_parts(long npoint, long nvector, long noutput, double** data, long* sizes, double* output) {
    long begin=0;
    for (long ioutput=0; ioutput<noutput; ioutput++) {
        double total = 0.0;
        long end = begin + sizes[ioutput];

        for (long ipoint=begin; ipoint<end; ipoint++) {
            double tmp = data[nvector-1][ipoint];
            for (long ivector=nvector-2; ivector>=0; ivector--) {
                tmp *= data[ivector][ipoint];
            }
            total += tmp;
        }
        output[ioutput] = total;

        begin = end;
    }
}

void grid_distances(double *points, double *center, double *distances, long n) {
  double d, tmp;
  for (long i=0; i<n; i++) {
    // x
    d = *points - center[0];
    tmp = d*d;
    points++;
    // y
    d = *points - center[1];
    tmp += d*d;
    points++;
    // z
    d = *points - center[2];
    tmp += d*d;
    *distances = sqrt(tmp);
    points++;
    distances++;
  }
}
