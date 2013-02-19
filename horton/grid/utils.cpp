// Horton is a Density Functional Theory program.
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


double dot_multi(long npoint, long nvector, double** data) {
    double result = 0.0;
    for (long ipoint=npoint; ipoint>0; ipoint--) {
        double tmp = *(data[nvector-1]);
        data[nvector-1]++;
        for (long ivector=nvector-2; ivector>=0; ivector--) {
           tmp *= *(data[ivector]);
           data[ivector]++;
        }
        result += tmp;
    }
    return result;
}

double dot_multi_poly_cube(long npoint, long nvector, double** data,
    double* origin, Cell* grid_cell, long* shape, long* pbc_active,
    Cell* cell, double* center, long mask, double powx, double powy,
    double powz, double powr) {

    // relative vector from origin to center
    double oc[3];
    oc[0] = center[0] - origin[0];
    oc[1] = center[1] - origin[1];
    oc[2] = center[2] - origin[2];
    // dot product of relative vector with grid_cell vectors
    double dot_oc[3];
    grid_cell->dot_cart(oc, dot_oc);

    // compute the cutoff radii for the mask function
    double rcut1 = -1;
    for (long i=grid_cell->get_nvec()-1; i >= 0; i--) {
        if (pbc_active[i]) {
            double tmp = 0.5*grid_cell->get_rspacing(i)*shape[i];
            if ((rcut1 == -1) || (tmp < rcut1)) rcut1 = tmp;
        } else {
            double tmp, l;
            l = grid_cell->get_rlength(i);
            tmp = dot_oc[i]/l;
            if (tmp < 0) throw std::domain_error("The center must lie inside the cell (1)");
            if ((rcut1 == -1) || (tmp < rcut1)) rcut1 = tmp;
#ifdef DEBUG
            printf("l=%f dot_oc[i]=%f\n", l, dot_oc[i]);
            printf("tmp1=%f\n", tmp);
#endif
            tmp = l*(shape[i]-1) - tmp;
            if (tmp < 0) throw std::domain_error("The center must lie inside the cell (2)");
#ifdef DEBUG
            printf("tmp2=%f\n", tmp);
#endif
            if ((rcut1 == -1) || (tmp < rcut1)) rcut1 = tmp;
        }
    }
    double rcut0 = 0.9*rcut1;
#ifdef DEBUG
    printf("rcut1=%f\n", rcut1);
    printf("rcut0=%f\n", rcut0);
#endif

    double result = 0.0;
    long i0=0;
    long i1=0;
    long i2=0;
    for (long ipoint=npoint; ipoint>0; ipoint--) {
        // Compute the products of the integranda
        double work = *(data[nvector-1]);
        data[nvector-1]++;
        for (long ivector=nvector-2; ivector>=0; ivector--) {
           work *= *(data[ivector]);
           data[ivector]++;
        }

        if (work != 0.0) {
            // compute the relative vector in the mic
            double frac[3];
            frac[0] = i0;
            frac[1] = i1;
            frac[2] = i2;
            double delta[3];
            grid_cell->to_cart(frac, delta);
            delta[0] -= oc[0];
            delta[1] -= oc[1];
            delta[2] -= oc[2];
            cell->mic(delta);

            // compute the distance
            double dist = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);

#ifdef DEBUG
            printf("i0=%i i1=%i i2=%i  |  %10.5f %10.5f %10.5f  |  %10.5f\n", i0, i1, i2, delta[0], delta[1], delta[2], dist);
#endif
            // multiply with a masking function and/or a user defined polynomial
            if (mask!=0) {
                if (dist > rcut1) {
                    work = 0.0;
                } else if (dist > rcut0) {
                    double x = (rcut1-dist)/(rcut1-rcut0);
                    work *= (3 - 2*x)*x*x;
                }
            }

            if (work != 0.0) {
                // multiply with four types of polynomials
                if (powx != 0) work *= pow(delta[0], powx);
                if (powy != 0) work *= pow(delta[1], powy);
                if (powz != 0) work *= pow(delta[2], powz);
                if (powr != 0) work *= pow(dist, powr);

                // Add to the multi-dot product
                result += work;
            }
        }

        // increment the indexes i0, i1 and i2
        if (i2<shape[2]-1) {
            i2++;
        } else {
            i2 = 0;
            if (i1<shape[1]-1) {
                i1++;
            } else {
                i1 = 0;
                if (i0<shape[0]-1) {
                    i0++;
                } else {
                    i0=0;
                }
            }
        }
    }
    return result;
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
