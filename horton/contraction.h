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


#ifndef HORTON_CONTRACTION_H
#define HORTON_CONTRACTION_H


long get_shell_nbasis(long shell_type);

int compute_gobasis_overlap(double* centers, long* shell_map, long* nprims,
    long* shell_types, double* alphas, double* con_coeffs, long nshell,
    long ncenter, long nprim_total, double* output, long nbasis);

double gob_normalization(double exponent, int nx, int ny, int nz);


typedef struct {
    // input data
    double* centers;
    long* shell_map;
    long* nprims;
    long* shell_types;
    long* basis_offsets;
    double* alphas;
    double* con_coeffs;
    long nshell, nbasis;

    // public iterator fields
    long max_nbasis, shell_type0, shell_type1;
    double con_coeff, alpha0, alpha1;
    double *r0, *r1; // current centers
    long ibasis0, ibasis1; // basis function counters (for output storage)

    // private iterator fields
    long ishell0, ishell1; // shell counters
    long nprim0, nprim1, iprim0, iprim1, oprim0, oprim1; // primitive counters
} i2gob_type;


i2gob_type* i2gob_new(void);
void i2gob_free(i2gob_type* i2);
int i2gob_init(i2gob_type* i2, double* centers, long* shell_map, long* nprims,
    long* shell_types, double* alphas, double* con_coeffs, long nshell,
    long ncenter, long nprim_total, long nbasis);
int i2gob_inc_shell(i2gob_type* i2);
void i2gob_update_shell(i2gob_type* i2);
int i2gob_inc_prim(i2gob_type* i2);
void i2gob_update_prim(i2gob_type* i2);
void i2gob_store(i2gob_type* i2, double *work_pure, double *output);


typedef struct {
    // input data
    long shell_type0, shell_type1, skip;

    // public iterator fields
    int nx0, ny0, nz0, nx1, ny1, nz1, offset;
} i2pow_type;

i2pow_type* i2pow_new(void);
void i2pow_free(i2pow_type* i2p);
void i2pow_init(i2pow_type* i2p, long shell_type0, long shell_type1, long max_nbasis);
int i2pow_inc(i2pow_type* i2p);

int i1pow_inc(int* nx, int* ny, int* nz);


#endif
