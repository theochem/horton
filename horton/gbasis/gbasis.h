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

// UPDATELIBDOCTITLE: Gaussian basis set classes

#ifndef HORTON_GBASIS_GBASIS_H
#define HORTON_GBASIS_GBASIS_H

#include "horton/gbasis/ints.h"
#include "horton/gbasis/fns.h"


const double gob_cart_normalization(const double alpha, const long* n);
const double gob_pure_normalization(const double alpha, const long l);


class GBasis {
    private:
        // Auxiliary arrays that contain convenient derived information.
        long* basis_offsets;
        long* prim_offsets;
        long* scales_offsets;
        long* shell_lookup;
        double* scales;  // pre-computed normalization constants.
        long nbasis, nscales;
        long max_shell_type;

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
        virtual ~GBasis();
        virtual const double normalization(const double alpha, const long* n) const = 0;
        void init_scales();
        void compute_two_index(double* output, GB2Integral* integral);
        void compute_four_index(double* output, GB4Integral* integral);
        void compute_grid_point1(double* output, double* point, GB1GridFn* grid_fn);
        double compute_grid_point2(double* dm, double* point, GB2DMGridFn* grid_fn);

        const long get_nbasis() const {return nbasis;}
        const long get_nscales() const {return nscales;}
        const long get_max_shell_type() const {return max_shell_type;}
        const long* get_basis_offsets() const {return basis_offsets;}
        const long* get_prim_offsets() const {return prim_offsets;}
        const long* get_shell_lookup() const {return shell_lookup;}
        const double* get_scales(long iprim) const {return scales + scales_offsets[iprim];}
};


class GOBasis : public GBasis {
    public:
         // double compute_normalization(double alpha, long nx, long ny, long nz);
        GOBasis(const double* centers, const long* shell_map, const long* nprims,
                const long* shell_types, const double* alphas, const double* con_coeffs,
                const long ncenter, const long nshell, const long nprim_total);
        const double normalization(const double alpha, const long* n) const;

        /** @brief
                Computes the overlap integrals.

            @param output
                The output array with the integrals.
         */
        void compute_overlap(double* output);

        /** @brief
                Computes the kinetic integrals.

            @param output
                The output array with the integrals.
         */
        void compute_kinetic(double* output);

        /** @brief
                Computes the nuclear attraction integrals.

            @param charges
                The array with values on the nuclear charges.

            @param centers
                The array with location of the nuclear charges.

            @param ncharge
                The number of nuclear charges.

            @param output
                The output array with the integrals.
         */
        void compute_nuclear_attraction(double* charges, double* centers, long ncharge,
                                        double* output);

        /** @brief
                Computes the nuclear attraction integrals.

            @param charges
                The array with values on the nuclear charges.

            @param centers
                The array with location of the nuclear charges.

            @param ncharge
                The number of nuclear charges.

            @param output
                The output array with the integrals.

            @param mu
                The range-separation parameter.
         */
        void compute_erf_attraction(double* charges, double* centers, long ncharge,
                                    double* output, double mu);

        /** @brief
                Computes the nuclear attraction integrals.

            @param charges
                The array with values on the nuclear charges.

            @param centers
                The array with location of the nuclear charges.

            @param ncharge
                The number of nuclear charges.

            @param output
                The output array with the integrals.
            @param c
                Coefficient of the gaussian.

            @param alpha
                Exponential parameter of the gaussian.
         */
        void compute_gauss_attraction(double* charges, double* centers, long ncharge,
                                      double* output, double c, double alpha);

        /** @brief
                Computes the electron repulsion integrals.

            @param output
                The output array with the integrals.
         */
        void compute_electron_repulsion(double* output);

        /** @brief
                Computes the ERF electron repulsion integrals.

            @param output
                The output array with the integrals.

            @param mu
                The range-separation parameter.
         */
        void compute_erf_repulsion(double* output, double mu);

        /** @brief
                Computes the Gaussian electron repulsion integrals.

            @param output
                The output array with the integrals.

            @param c
                Coefficient of the gaussian.

            @param alpha
                Exponential parameter of the gaussian.
         */
        void compute_gauss_repulsion(double* output, double c, double alpha);

        /** @brief
                Computes the r^alpha electron repulsion integrals.

            @param output
                The output array with the integrals.

            @param alpha
                The power of r in the potential.
         */
        void compute_ralpha_repulsion(double* output, double alpha);

        /** @brief
                Computes the (multipole) moment integrals.

            @param xyz
                The powers of xyz in the integrals.

            @param center
                The location around which the moment integrals are computed.

            @param output
                The output array with the integrals.
         */
        void compute_multipole_moment(long* xyz, double* center, double* output);

        void compute_grid1_exp(long nfn, double* coeffs, long npoint, double* points,
                               long norb, long* iorbs, double* output);

        /** @brief
                Computes the gradient of the molecular orbital on a grid.

            @param nfn
                The number of functions.

            @param coeffs
                The coefficients for the basisfunction expanion.

            @param npoint
                The number of grid points to be calculated.

            @param points
                The coordinates of grid points to be calculated.

            @param norb
                The number of orbitals to be calculated.

            @param iorbs
                The orbitals to be calculated.

            @param output
                The output array with the integrals.
     */
        void compute_grid1_grad_exp(long nfn, double* coeffs, long npoint,
                                    double* points, long norb, long* iorbs, double* output);
        void compute_grid1_dm(double* dm, long npoint, double* points,
                              GB1DMGridFn* grid_fn, double* output,
                              double epsilon, double* dmmaxrow);
        void compute_grid2_dm(double* dm, long npoint, double* points, double* output);
        void compute_grid1_fock(long npoint, double* points, double* weights,
                                long pot_stride, double* pots,
                                GB1DMGridFn* grid_fn, double* output);
};

#endif  // HORTON_GBASIS_GBASIS_H_
