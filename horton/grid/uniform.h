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


#ifndef HORTON_GRID_UNIFORM_H
#define HORTON_GRID_UNIFORM_H

#include "cell.h"
#include "cubic_spline.h"


class UniformIntGrid {
    private:
        double origin[3];
        Cell* grid_cell;
        Cell* cell;
        long shape[3], pbc[3];
    public:
        UniformIntGrid(double* _origin, Cell* _grid_cell, long* _shape, long* _pbc, Cell* _cell);

        void copy_origin(double* output);
        void copy_shape(long* output);
        void copy_pbc(long* output);

        const long* get_shape() { return shape; };
        const Cell* get_cell() { return cell; };

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

        // for progress info
        const bool progress;
        long imarker, idot, icall;
    public:
        Range3Iterator(const long* ranges_begin, const long* ranges_end, const long* shape, const bool progress);

        long get_npoint() const { return npoint; };
        void set_point(long ipoint, long* i, long* iwrap);
};


#endif
