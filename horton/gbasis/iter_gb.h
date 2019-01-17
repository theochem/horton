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

// UPDATELIBDOCTITLE: Iterators over Gaussian basis functions

#ifndef HORTON_GBASIS_ITER_GB_H
#define HORTON_GBASIS_ITER_GB_H

#include "horton/gbasis/gbasis.h"

class IterGB1 {
    private:
        // input fields
        const GBasis* gbasis;
        const long* basis_offsets;
    public:
        IterGB1(GBasis* gbasis);

        int inc_shell();
        void update_shell();
        int inc_prim();
        void update_prim();
        void store(const double* work, double* output, long dim);

        // 'public' iterator fields
        long shell_type0;
        double con_coeff, alpha0;
        const double* r0; // current centers
        const double* scales0; // normalization constants
        long ibasis0; // basis function counters (for output storage)

        // 'private' iterator fields
        long ishell0; // shell counters
        long nprim0, iprim0, oprim0; // primitive counters
    };


class IterGB2 {
    private:
        // input fields
        const GBasis* gbasis;
        const long* basis_offsets;
    public:
        IterGB2(GBasis* gbasis);

        int inc_shell();
        void update_shell();
        int inc_prim();
        void update_prim();
        void store(const double* work, double* output);
        double dot(const double* work, const double* dm);

        // 'public' iterator fields
        long shell_type0, shell_type1;
        double con_coeff, alpha0, alpha1;
        const double* r0;// current centers
        const double* r1; // current centers
        const double* scales0; // normalization constants
        const double* scales1; // normalization constants
        long ibasis0, ibasis1; // basis function counters (for output storage)

        // 'private' iterator fields
        long ishell0, ishell1; // shell counters
        long nprim0, nprim1, iprim0, iprim1, oprim0, oprim1; // primitive counters
    };


class IterGB4 {
    private:
        // input fields
        const GBasis* gbasis;
        const long* basis_offsets;
    public:
        IterGB4(GBasis* gbasis);

        int inc_shell();
        void update_shell();
        int inc_prim();
        void update_prim();
        void store(const double* work, double* output);

        // 'public' iterator fields
        long shell_type0, shell_type1, shell_type2, shell_type3;
        double con_coeff, alpha0, alpha1, alpha2, alpha3;
        const double* r0; // current centers
        const double* r1; // current centers
        const double* r2; // current centers
        const double* r3; // current centers
        const double* scales0; // normalization constants
        const double* scales1; // normalization constants
        const double* scales2; // normalization constants
        const double* scales3; // normalization constants
        long ibasis0, ibasis1, ibasis2, ibasis3; // basis function counters (for output storage)

        // 'private' iterator fields
        long ishell0, ishell1, ishell2, ishell3; // shell counters
        long ishell3_max;
        long nprim0, nprim1, nprim2, nprim3; // primitive counters
        long iprim0, iprim1, iprim2, iprim3;
        long oprim0, oprim1, oprim2, oprim3;
    };

#endif
