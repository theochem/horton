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


#ifndef HORTON_GBASIS_GBASIS_H
#define HORTON_GBASIS_GBASIS_H

#include "ints.h"

class GBasis {
    private:
        // Auxiliary arrays that contain convenient derived information.
        long* basis_offsets;
        long* scales_offsets;
        double* scales; // pre-computed normalization constants.
        long nbasis, nscales;
        long max_shell_nbasis;

    public:
        // Arrays that fully describe the basis set.
        const double* centers;
        const long* shell_map;
        const long* nprims;
        const long* shell_types;
        const double* alphas;
        const double* con_coeffs;
        const long ncenter, nshell, nprim_total;

        GBasis(const double* centers, const long* shell_map, const long* nprims,
               const long* shell_types, const double* alphas, const double* con_coeffs,
               const long ncenter, const long nshell, const long nprim_total);
        ~GBasis();
        virtual const double normalization(const double alpha, const long* n) const =0;
        void init_scales();
        void compute_one_body(double* output, GB2Integral* integral);

        const long get_nbasis() const {return nbasis;};
        const long get_nscales() const {return nscales;};
        const long get_max_shell_nbasis() const {return max_shell_nbasis;};
        const long* get_basis_offsets() const {return basis_offsets;};
        const double* get_scales(long iprim) const {return scales + scales_offsets[iprim];};
    };


class GOBasis : public GBasis {
    public:
        //double compute_normalization(double alpha, long nx, long ny, long nz);
        GOBasis(const double* centers, const long* shell_map, const long* nprims,
                const long* shell_types, const double* alphas, const double* con_coeffs,
                const long ncenter, const long nshell, const long nprim_total);
        const double normalization(const double alpha, const long* n) const;

        void compute_overlap(double* output);
        void compute_kinetic(double* output);
        void compute_nuclear_attraction(double* charges, double* centers, long ncharge, double* output);
    };

#endif
