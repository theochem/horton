// Horton is a development platform for electronic structure methods.
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
        double rlengths[3], glengths[3]; // TODO: test these
        double rspacings[3], gspacings[3];
        double volume;
        int nvec;
    public:
        Cell(): nvec(0) {};

        void update(double* _rvecs, double* _gvecs, int _nvec);
        void mic(double* delta) const;
        void to_frac(double* cart, double* frac) const;
        void to_cart(double* frac, double* cart) const;
        void g_lincomb(double* coeffs, double* gvec) const;
        void dot_rvecs(double* cart, double* dot_rvecs) const;
        void add_rvec(double* delta, long* r) const;

        int get_nvec() const {return nvec;};
        double get_volume() const {return volume;};
        double get_rspacing(int i) const;
        double get_gspacing(int i) const;
        double get_rlength(int i) const;
        double get_glength(int i) const;

        void copy_rvecs(double* _rvecs) const;
        void copy_gvecs(double* _gvecs) const;
        void copy_rlengths(double* _rlengths) const;
        void copy_glengths(double* _glengths) const;
        void copy_rspacings(double* _rspacings) const;
        void copy_gspacings(double* _gspacings) const;

        void set_ranges_rcut(double* delta, double rcut, long* ranges_begin,
            long* ranges_end) const;

        long select_inside(double* origin, double* center, double rcut,
            long* ranges_begin, long* ranges_end, long* shape,
            long* pbc, long* indexes) const;
};

long smart_wrap(long i, long shape, long pbc);

#endif
