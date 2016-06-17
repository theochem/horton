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

// UPDATELIBDOCTITLE: Uniform 3D grids

#ifndef HORTON_GRID_UNIFORM_H
#define HORTON_GRID_UNIFORM_H

#include "horton/cell.h"
#include "horton/grid/cubic_spline.h"


class UniformGrid {
    public:
        double origin[3];
        double grid_rvecs[9];
        long shape[3];
        long pbc[3];

        UniformGrid(double* _origin, double* _grid_rvecs, long* _shape, long* _pbc);

        Cell* get_cell();
        Cell* get_grid_cell();

        void set_ranges_rcut(double* center, double rcut, long* ranges_begin, long* ranges_end);
        double dist_grid_point(double* center, long* i);
        void delta_grid_point(double* center, long* i);
        double* get_pointer(double* array, long* i);
};


long index_wrap(long i, long high);


class Range3Iterator {
    private:
        const long* ranges_begin;
        const long* ranges_end;
        const long* shape;
        long loop_shape[3];
        long npoint;
    public:
        Range3Iterator(const long* ranges_begin, const long* ranges_end, const long* shape);

        long get_npoint() const { return npoint; };
        void set_point(long ipoint, long* i, long* iwrap);
};


class Cube3Iterator {
    private:
        const long* begin;
        const long* end;
        long shape[3];
        long npoint;
    public:
        Cube3Iterator(const long* begin, const long* end);

        long get_npoint() const { return npoint; };
        void set_point(long ipoint, long* j);
};


#endif
