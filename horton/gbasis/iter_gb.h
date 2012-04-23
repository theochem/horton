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


#ifndef HORTON_GBASIS_ITERGB_H
#define HORTON_GBASIS_ITERGB_H

#include "gbasis.h"

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
        void store(const double *work, double *output);

        // 'public' iterator fields
        long shell_type0, shell_type1;
        double con_coeff, alpha0, alpha1;
        const double *r0, *r1; // current centers
        const double *scales0, *scales1; // normalization constants
        long ibasis0, ibasis1; // basis function counters (for output storage)

        // 'private' iterator fields
        long ishell0, ishell1; // shell counters
        long nprim0, nprim1, iprim0, iprim1, oprim0, oprim1; // primitive counters
    };

#endif
