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

#include <stdexcept>
#include <cmath>
#include "uniform.h"

UniformIntGrid::UniformIntGrid(double* _origin, Cell* _grid_cell, long* _shape, long* _pbc)
    : grid_cell(_grid_cell) {

    if (_grid_cell->get_nvec() != 3) {
        throw std::domain_error("The grid_cell of a UniformIntGrid must be 3D.");
    }

    for (int i=0; i<3; i++) {
        origin[i] = _origin[i];
        shape[i] = _shape[i];
        pbc[i] = _pbc[i];
    }
}

void UniformIntGrid::copy_origin(double* output) {
    output[0] = origin[0];
    output[1] = origin[1];
    output[2] = origin[2];
}

void UniformIntGrid::copy_shape(long* output) {
    output[0] = shape[0];
    output[1] = shape[1];
    output[2] = shape[2];
}

void UniformIntGrid::copy_pbc(long* output) {
    output[0] = pbc[0];
    output[1] = pbc[1];
    output[2] = pbc[2];
}



void UniformIntGrid::set_ranges_rcut(double* center, double rcut, long* ranges_begin, long* ranges_end) {
    double delta[3];
    delta[0] = origin[0] - center[0];
    delta[1] = origin[1] - center[1];
    delta[2] = origin[2] - center[2];
    grid_cell->set_ranges_rcut(delta, rcut, -1, ranges_begin, ranges_end);

    // Truncate ranges in case of non-periodic boundary conditions
    for (int i=2; i>=0; i--) {
        if (!pbc[i]) {
            if (ranges_begin[i] < 0) {
                ranges_begin[i] = 0;
            }
            if (ranges_end[i] > shape[i]) {
                ranges_end[i] = shape[i];
            }
        }
    }

#ifdef DEBUG
    printf("ranges [%li,%li,%li] [%li,%li,%li]\n", ranges_begin[0], ranges_begin[1], ranges_begin[2], ranges_end[0], ranges_end[1], ranges_end[2]);
#endif
}

double UniformIntGrid::dist_grid_point(double* center, long* i) {
    double delta[3];
    delta[0] = origin[0] - center[0];
    delta[1] = origin[1] - center[1];
    delta[2] = origin[2] - center[2];
    grid_cell->add_rvec(delta, i);
    return sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
}
