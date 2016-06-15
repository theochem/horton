// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2016 The HORTON Development Team
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

// UPDATELIBDOCTITLE: Computation of weights for ESP cost functions

#ifndef HORTON_ESPFIT_MASK
#define HORTON_ESPFIT_MASK

#include "horton/cell.h"
#include "horton/grid/uniform.h"

void multiply_dens_mask(double* rho, double lnrho0, double sigma, double* weights, long npoint);

void multiply_near_mask(double* center, UniformGrid* ugrid, double r0,
    double gamma, double* weights);

void multiply_far_mask(double* centers, long ncenter, UniformGrid* ugrid,
    double r0, double gamma, double* weights);

#endif
