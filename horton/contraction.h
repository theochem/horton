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


long get_con_nbasis(long con_type);

int compute_gobasis_overlap(double* centers, long* shell_map,
    long* ncons, long* nexps, long* con_types, double* exponents, double*
    con_coeffs, long nshell, long ncenter, long ncon_total, long nexp_total,
    long ncon_coeff_total, double* output, long nbasis);



typedef struct {
    // input data
    double* centers;
    long* shell_map;
    long* ncons;
    long* nexps;
    long* con_types;
    double* exponents;
    double* con_coeffs;
    long nshell;

    // public iterator fields
    long max_nbasis, con_type0, con_type1;
    double con_coeff, exp0, exp1, x0, y0, z0, x1, y1, z1;

    // private iterator fields
    long ishell0, ishell1, ncon0, ncon1, icon0, icon1, ocon0, ocon1;
} i2gob_type;


i2gob_type* i2gob_new(void);
void i2gob_free(i2gob_type* i2);
int i2gob_init(i2gob_type* i2, double* centers, long* shell_map,
    long* ncons, long* nexps, long* con_types, double* exponents,
    double* con_coeffs, long nshell, long ncenter, long ncon_total,
    long nexp_total, long ncon_coeff_total, long nbasis);
int i2gob_check(i2gob_type* i2, long nshell, long ncenter, long ncon_total,
    long nexp_total, long ncon_coeff_total, long nbasis);
int i2gob_inc_shell(i2gob_type* i2);
void i2gob_update_shell(i2gob_type* i2);
int i2gob_inc_con(i2gob_type* i2);
void i2gob_update_con(i2gob_type* i2);
int i2gob_inc_exp(i2gob_type* i2);
void i2gob_store(i2gob_type* i2, double *work_pure, double *output);


typedef struct {
    // input data
    long con_type0, con_type1, skip;

    // public iterator fields
    int l0, m0, n0, l1, m1, n1, offset;
} i2pow_type;

i2pow_type* i2pow_new(void);
void i2pow_free(i2pow_type* i2p);
void i2pow_init(i2pow_type* i2p, long con_type0, long con_type1, long max_nbasis);
int i2pow_inc(i2pow_type* i2p);

int i1pow_inc(int* l, int* m, int* n);

void project_cartesian_to_pure(double *work_cart, double* work_pure, long
    con_type, long stride, long spacing, long count);


#endif
