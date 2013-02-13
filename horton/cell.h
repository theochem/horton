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


#ifndef HORTON_CELL_H
#define HORTON_CELL_H


class Cell {
    private:
        double rvecs[9], gvecs[9];
        double rspacings[3], gspacings[3];
        double volume;
        int nvec;
    public:
        Cell(): nvec(0) {};

        void update(double* _rvecs, double* _gvecs, int _nvec);
        void mic(double* delta);
        void to_center(double* car, long* center);
        void to_frac(double* cart, double* frac);
        void to_cart(double* frac, double* cart);
        void add_vec(double* delta, long* r);

        int get_nvec() {return nvec;};
        double get_volume() {return volume;};
        double get_rspacing(int i);
        double get_gspacing(int i);

        void copy_rvecs(double* _rvecs);
        void copy_gvecs(double* _gvecs);
        void copy_rspacings(double* _rspacings);
        void copy_gspacings(double* _gspacings);

        void set_ranges_rcut(double* origin, double* center, double rcut,
                             long* ranges_begin, long* ranges_end);
};


#endif
