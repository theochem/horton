# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

cimport cell
cimport uniform

cdef extern from "ewald3d.h":
    double pair_ewald3d(double* delta, cell.Cell* cell, double rcut,
        double alpha, double gcut)

    void setup_esp_cost_cube_ewald3d(uniform.UniformIntGrid* ui_grid,
        double* vref, double* weights, double* centers, double* A, double* B,
        double* C, long ncenter, double rcut, double alpha, double gcut) except +
