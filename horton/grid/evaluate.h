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



#ifndef HORTON_GRID_EVALUATE_H
#define HORTON_GRID_EVALUATE_H


#include "cell.h"
#include "cubic_spline.h"
#include "uniform.h"


void eval_spline_cube(CubicSpline* spline, double* center, double* output,
                      UniformIntGrid* ui_grid);

void eval_spline_grid(CubicSpline* spline, double* center, double* output,
                      double* points, Cell* cell, long npoint);

#endif
