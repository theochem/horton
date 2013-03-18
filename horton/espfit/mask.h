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


#ifndef HORTON_ESPFIT_MASK
#define HORTON_ESPFIT_MASK

#include "cell.h"
#include "uniform.h"

void multiply_dens_mask(double* rho, double rho0, double alpha, double* weights, long npoint);

void multiply_near_mask(double* center, UniformIntGrid* ui_grid, double r0,
    double gamma, double* weights);

void multiply_far_mask(double* centers, long ncenter, UniformIntGrid* ui_grid,
    double r0, double gamma, double* weights);

#endif
