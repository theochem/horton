# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


cimport fns

cdef extern from "horton/gbasis/gbasis.h":
    double gob_cart_normalization(double alpha, long* n)
    double gob_pure_normalization(double alpha, long l)

    cdef cppclass GBasis:
        # Arrays that fully describe the basis set.
        double* centers
        long* shell_map
        long* nprims
        long* shell_types
        double* alphas
        double* con_coeffs
        long ncenter, nshell, nprim_total

        # Auxiliary bits
        long get_nbasis()
        long get_nscales()
        long get_max_shell_type()
        double* get_scales(long iprim)
        long* get_shell_lookup()
        long* get_basis_offsets()

        # low-level compute routines
        void compute_grid_point1(double* output, double* point, fns.GB1DMGridFn* grid_fn)
        double compute_grid_point2(double* dm, double* point, fns.GB2DMGridFn* grid_fn)

    cdef cppclass GOBasis:
        GOBasis(double* centers, long* shell_map, long* nprims,
                long* shell_types, double* alphas, double* con_coeffs,
                long ncenter, long nshell, long nprim_total) except +

        void compute_overlap(double* output)
        void compute_kinetic(double* output)
        void compute_nuclear_attraction(double* charges, double* centers, long ncharge, double* output)
        void compute_electron_repulsion(double* output)
        void compute_grid1_exp(long nfn, double* coeffs, long npoint, double* points, long norb, long* iorbs, double* output)
        void compute_grid1_dm(double* dm, long npoint, double* points, fns.GB1DMGridFn* grid_fn, double* output, double epsilon, double* dmmaxrow)
        void compute_grid2_dm(double* dm, long npoint, double* points, double* output)
        void compute_grid1_fock(long npoint, double* points, double* weights, long pot_stride, double* pots, fns.GB1DMGridFn* grid_fn, double* output)
