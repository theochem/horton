# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


cdef extern from "contraction.h":
    long get_shell_nbasis(long shell_type)

    int compute_gobasis_overlap(double* centers, long* shell_map, long* nprims,
        long* shell_types, double* alphas, double* con_coeffs, long nshell,
        long ncenter, long nprim_total, double* output, long nbasis)

    double gob_normalization(double exponent, int nx, int ny, int nz)

    ctypedef struct i2gob_type:
        long nshell, nbasis

        long max_nbasis, shell_type0, shell_type1
        double con_coeff, alpha0, alpha1
        double r0[3], r1[3]
        long ibasis0, ibasis1

        long ishell0, ishell1
        long nprim0, nprim1, oprim0, oprim1, iprim0, iprim1


    i2gob_type* i2gob_new()
    void i2gob_free(i2gob_type* i2)
    int i2gob_init(i2gob_type* i2, double* centers, long* shell_map, long* nprims,
        long* shell_types, double* alphas, double* con_coeffs, long nshell,
        long ncenter, long nprim_total, long nbasis)
    bint i2gob_inc_shell(i2gob_type* i2)
    void i2gob_update_shell(i2gob_type* i2)
    bint i2gob_inc_prim(i2gob_type* i2)
    void i2gob_update_prim(i2gob_type* i2)
    void i2gob_store(i2gob_type* i2, double *work_pure, double *output)

    ctypedef struct i2pow_type:
        int nx0, ny0, nz0, nx1, ny1, nz1, offset

    i2pow_type* i2pow_new()
    void i2pow_free(i2pow_type* i2p)
    void i2pow_init(i2pow_type* i2p, long shell_type0, long shell_type1, long max_nbasis)
    bint i2pow_inc(i2pow_type* i2p)

    bint i1pow_inc(int* nx, int* ny, int* nz)
