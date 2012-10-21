// Horton is a Density Functional Theory program.
// Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


#ifndef HORTON_GBASIS_FNS_H
#define HORTON_GBASIS_FNS_H

#include "calc.h"
#include "iter_pow.h"


class GB1GridFn : public GBCalculator  {
    protected:
        long shell_type0;
        const double *r0;
        const double *point;
        IterPow1 i1p;
    public:
        GB1GridFn(long max_shell_type);

        void reset(long shell_type0, const double* r0, const double* point);
        void cart_to_pure();
        const long get_shell_type0() const {return shell_type0;};

        virtual long get_dim_work() = 0;
        virtual long get_dim_output() = 0;
        virtual void add(double coeff, double alpha0, const double* scales0) = 0;
        virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output) = 0;
        virtual void compute_fock_from_dm(double factor, double* work_basis, long nbasis, double* output) = 0;
    };


class GB1GridDensityFn : public GB1GridFn  {
    public:
        GB1GridDensityFn(long max_shell_type): GB1GridFn(max_shell_type) {};

        long get_dim_work() {return 1;};
        long get_dim_output() {return 1;};
        void add(double coeff, double alpha0, const double* scales0);
        void compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output);
        void compute_fock_from_dm(double factor, double* work_basis, long nbasis, double* output);
    };


#endif
