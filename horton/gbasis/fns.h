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

// UPDATELIBDOCTITLE: Evaluation of functions expanded in a Gaussian basis

#ifndef HORTON_GBASIS_FNS_H
#define HORTON_GBASIS_FNS_H

#include "horton/gbasis/calc.h"
#include "horton/gbasis/common.h"
#include "horton/gbasis/iter_pow.h"


class GB1GridFn : public GBCalculator  {
    protected:
        long shell_type0;
        const long dim_work, dim_output;
        const double *r0;    // The center of the basis function
        const double *point; // The grid point at which the fn is evaluated
        IterPow1 i1p;
    public:
        GB1GridFn(long max_shell_type, long dim_work, long dim_output);

        virtual void reset(long shell_type0, const double* r0, const double* point);
        void cart_to_pure();
        const long get_shell_type0() const {return shell_type0;};

        long get_dim_work() {return dim_work;};
        long get_dim_output() {return dim_output;};
        virtual void add(double coeff, double alpha0, const double* scales0) = 0;
    };




class GB1ExpGridFn : public GB1GridFn  {
    protected:
        long nfn;
    public:
        GB1ExpGridFn(long max_shell_type, long nfn, long dim_work, long dim_output) : GB1GridFn(max_shell_type, dim_work, dim_output), nfn(nfn) {};
        virtual void compute_point_from_exp(double* work_basis, double* coeffs, long nbasis, double* output) = 0;
    };


class GB1ExpGridOrbitalFn : public GB1ExpGridFn  {
    protected:
        long* iorbs;
        long norb;
    public:
        GB1ExpGridOrbitalFn(long max_shell_type, long nfn, long* iorbs, long norb) : GB1ExpGridFn(max_shell_type, nfn, 1, norb), iorbs(iorbs), norb(norb) {};
        virtual void add(double coeff, double alpha0, const double* scales0);
        virtual void compute_point_from_exp(double* work_basis, double* coeffs, long nbasis, double* output);
    };




class GB1DMGridFn : public GB1GridFn  {
    public:
        GB1DMGridFn(long max_shell_type, long dim_work, long dim_output) : GB1GridFn(max_shell_type, dim_work, dim_output) {};
        virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output, double epsilon, double* dmmaxrow) = 0;
        virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis, double* output) = 0;
    };


class GB1DMGridDensityFn : public GB1DMGridFn  {
    private:
        double poly_work[MAX_NCART_CUMUL-1];
        long offset;
    public:
        GB1DMGridDensityFn(long max_shell_type): GB1DMGridFn(max_shell_type, 1, 1) {};

        virtual void reset(long _shell_type0, const double* _r0, const double* _point);
        virtual void add(double coeff, double alpha0, const double* scales0);
        virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output, double epsilon, double* dmmaxrow);
        virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis, double* output);
    };


class GB1DMGridGradientFn : public GB1DMGridFn  {
    public:
        GB1DMGridGradientFn(long max_shell_type): GB1DMGridFn(max_shell_type, 4, 3) {};

        virtual void add(double coeff, double alpha0, const double* scales0);
        virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output, double epsilon, double* dmmaxrow);
        virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis, double* output);
    };


class GB1DMGridKineticFn : public GB1DMGridFn  {
    public:
        GB1DMGridKineticFn(long max_shell_type): GB1DMGridFn(max_shell_type, 3, 1) {};

        virtual void add(double coeff, double alpha0, const double* scales0);
        virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output, double epsilon, double* dmmaxrow);
        virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis, double* output);
    };




class GB2DMGridFn : public GBCalculator  {
    protected:
        long shell_type0, shell_type1;
        const double *r0, *r1;  // The centers of the basis functions
        const double *point;    // The grid point at which the fn is evaluated
        IterPow2 i2p;
    public:
        GB2DMGridFn(long max_shell_type);

        void reset(long shell_type0, long shell_type1, const double* r0, const double* r1, const double* point);
        void cart_to_pure();
        const long get_shell_type0() const {return shell_type0;};
        const long get_shell_type1() const {return shell_type1;};

        virtual void add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1) = 0;
    };


class GB2DMGridHartreeFn : public GB2DMGridFn  {
    private:
        double* work_g0;
        double* work_g1;
        double* work_g2;
        double* work_boys;
    public:
        GB2DMGridHartreeFn(long max_shell_type);
        ~GB2DMGridHartreeFn();

        virtual void add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1);
    };




#endif
