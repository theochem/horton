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

// UPDATELIBDOCTITLE: Computation of ESP cost functions

#ifndef HORTON_ESPFIT_ELECTROSTATICS_H
#define HORTON_ESPFIT_ELECTROSTATICS_H

#include "horton/cell.h"
#include "horton/grid/uniform.h"

double pair_electrostatics(double* delta, const Cell* cell, double rcut, double alpha,
    double gcut);

double pair_ewald3d(double* delta, const Cell* cell, double rcut, double alpha,
    double gcut);

void setup_esp_cost_cube(UniformGrid* ugrid, double* vref,
    double* weights, double* centers, double* A, double* B, double* C,
    long ncenter, double rcut, double alpha, double gcut);

void compute_esp_cube(UniformGrid* ugrid, double* esp,
    double* centers, double* charges, long ncenter, double rcut, double alpha,
    double gcut);

#endif
