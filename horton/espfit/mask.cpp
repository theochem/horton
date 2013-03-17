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


#include <cmath>
#include "mask.h"


double switch_fn(double x, double a) {
    x /= a;
    return 0.5*(x*(3-x*x));
}


void multiply_dens_mask(double* rho, double rho0, double alpha, double* weights, long npoint) {
    while (npoint > 0) {
        *weights *= switch_fn(log10(*rho) - log10(rho0), alpha);

        // move on
        npoint --;
        rho++;
        weights++;
    }
}

/*
void multiply_near_mask(double* center, double* origin, Cell* grid_cell,
    long* shape, Cell* cell, double r0, double gamma, double* weights) {

    double delta[3]
    long ranges_begin[3], ranges_end[3]

    grid_cell->set_ranges_rcut(
}*/
