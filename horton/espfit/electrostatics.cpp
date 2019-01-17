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


#include <cstdio>
#include <cmath>
#include <cstddef>
#include "horton/espfit/electrostatics.h"


double pair_electrostatics(double* delta, const Cell* cell, double rcut,
    double alpha, double gcut) {

    double result = 0.0;
    int nvec = cell->get_nvec();

    if (nvec == 0) {
        // 0D pbc. Do not use cutoff. All pairs are included.
        result = 1.0/sqrt(delta[0]*delta[0] + delta[1]*delta[1] +
                          delta[2]*delta[2]);
    }
    else if (nvec == 3) {
        // 3D pbc
        result = pair_ewald3d(delta, cell, rcut, alpha, gcut);
    }
    return result;
}


double pair_ewald3d(double* delta, const Cell* cell, double rcut, double alpha,
    double gcut) {

    double result = 0;

    // Determine the ranges of the real sum
    long rbegin[3], rend[3];
    cell->set_ranges_rcut(delta, rcut, rbegin, rend);

    // the real-space terms
    double j[3];
    for (int j0=rbegin[0]; j0 < rend[0]; j0++) {
        j[0] = j0;
        for (int j1=rbegin[1]; j1 < rend[1]; j1++) {
            j[1] = j1;
            for (int j2=rbegin[2]; j2 < rend[2]; j2++) {
                j[2] = j2;
                double tmp[3];
                cell->to_cart(j, tmp);
                tmp[0] += delta[0];
                tmp[1] += delta[1];
                tmp[2] += delta[2];
                double d = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
                result += erfc(alpha*d)/d;
            }
        }
    }

    // Precomput some factors
    double fac1 = 4.0*M_PI/cell->get_volume();
    double fac2 = 0.25/alpha/alpha;

    // Determine the ranges of the reciprocal sum
    double kmax[3];
    kmax[0] = ceil(gcut/cell->get_gspacing(0));
    kmax[1] = ceil(gcut/cell->get_gspacing(1));
    kmax[2] = ceil(gcut/cell->get_gspacing(2));

    // add the reciprocal-space terms
    for (int j0=-kmax[0]; j0 <= kmax[0]; j0++) {
        j[0] = 2*M_PI*j0;
        for (int j1=-kmax[1]; j1 <= kmax[1]; j1++) {
            j[1] = 2*M_PI*j1;
            for (int j2=0; j2 <= kmax[2]; j2++) {
                if ((j0==0) && (j1==0) && (j2==0)) continue;
                j[2] = 2*M_PI*j2;
                double k[3];
                cell->g_lincomb(j, k);
                double ksq = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
                result += (1+(j2>0))*fac1*exp(-ksq*fac2)/ksq*cos(
                    k[0]*delta[0] + k[1]*delta[1] + k[2]*delta[2]
                );
            }
        }
    }

    // background correction (needed to make result indepdent of alpha
    result -= M_PI/cell->get_volume()/alpha/alpha;

    return result;
}


void setup_esp_cost_cube(UniformGrid* ugrid, double* vref,
    double* weights, double* centers, double* A, double* B, double* C,
    long ncenter, double rcut, double alpha, double gcut) {

    Cell* cell = ugrid->get_cell();
    Cell* grid_cell = ugrid->get_grid_cell();
    double gvol = grid_cell->get_volume();
    bool is3d = (cell->get_nvec() == 3);
    long neq = ncenter + is3d;
    double* work = new double[neq];
    double grid_cart[3];

    Cube3Iterator c3i = Cube3Iterator(NULL, ugrid->shape);
    long i[3];
    long npoint = c3i.get_npoint();

    for (long ipoint=0; ipoint<npoint; ipoint++) {
        c3i.set_point(ipoint, i);
        grid_cart[0] = 0;
        grid_cart[1] = 0;
        grid_cart[2] = 0;
        ugrid->delta_grid_point(grid_cart, i);

        if (*weights > 0) {
            double sqrtw = sqrt((*weights)*gvol);

            // Do some electrostatics
            for (long icenter=0; icenter<ncenter; icenter++) {
                double delta[3];
                delta[0] = centers[3*icenter]   - grid_cart[0];
                delta[1] = centers[3*icenter+1] - grid_cart[1];
                delta[2] = centers[3*icenter+2] - grid_cart[2];

                work[icenter] = sqrtw*pair_electrostatics(delta, cell, rcut, alpha, gcut);
            }
            if (is3d) work[ncenter] = sqrtw;

            // Add to the quadratic cost function
            double vrefw = (*vref)*sqrtw;
            for (long ic0=0; ic0<neq; ic0++) {
                for (long ic1=0; ic1<neq; ic1++) {
                    A[ic0+neq*ic1] += work[ic0]*work[ic1];
                }
                B[ic0] += vrefw*work[ic0];
            }
            *C += vrefw*vrefw;
        }

        // move on
        vref++;
        weights++;
    }

    delete[] work;
    delete cell;
    delete grid_cell;
}


void compute_esp_cube(UniformGrid* ugrid, double* esp,
    double* centers, double* charges, long ncenter, double rcut, double alpha,
    double gcut) {

    double grid_cart[3];
    Cube3Iterator c3i = Cube3Iterator(NULL, ugrid->shape);
    long i[3];
    long npoint = c3i.get_npoint();
    Cell* cell = ugrid->get_cell();

    for (long ipoint=0; ipoint<npoint; ipoint++) {
        // Compute the position of the grid point
        c3i.set_point(ipoint, i);
        grid_cart[0] = 0;
        grid_cart[1] = 0;
        grid_cart[2] = 0;
        ugrid->delta_grid_point(grid_cart, i);

        // Do some ewald stuff
        double tmp = 0;
        for (long icenter=0; icenter<ncenter; icenter++) {
            double delta[3];
            delta[0] = centers[3*icenter]   - grid_cart[0];
            delta[1] = centers[3*icenter+1] - grid_cart[1];
            delta[2] = centers[3*icenter+2] - grid_cart[2];
            tmp += charges[icenter]*pair_electrostatics(delta, cell, rcut,
                                                        alpha, gcut);
        }
        (*esp) = tmp;

        // move on
        esp++;
    }

    delete cell;
}
