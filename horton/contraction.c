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


#include "contraction.h"
#include "cints.h"

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define safe_exit(x) {result=x; goto exit;}


long get_con_nbasis(long con_type) {
    if (con_type > 0) {
        // Cartesian
        return (con_type+1)*(con_type+2)/2;
    } else if (con_type == -1) {
        // should not happen.
        return -1;
    } else {
        // Pure
        return -2*con_type+1;
    }
}

void swap_ptr(double **x, double **y) {
    double *t = *x;
    *x = *y;
    *y = t;
}


int compute_gobasis_overlap(double* centers, long* shell_map,
    long* ncons, long* nexps, long* con_types, double* exponents, double*
    con_coeffs, long nshell, long ncenter, long ncon_total, long nexp_total,
    long ncon_coeff_total, double* output, long nbasis) {

    /*

    This routine has the following structure:

    ----------------------------------------------------------------------------

    (A) Allocated working memory, based on largest angular momentum in the basis
        set.

    (B) Double loop over the shells:

        (C) Double loop over the contractions:

            (D) Clear working memory for storage of Cartesian Gaussian integrals.

            (E) Double loop over the exponents:

                (F) Double loop over basis functions associated with the angular
                    momenta of each contraction:

                    (G) Add contribution from a pair of primitives to the
                        working memory. This part contains a call to a function
                        from cints.c.

            (H) If needed, perform projection from Cartesian to pure basis
                functions. (A second piece of working memory is needed for this
                transformation.)

            (I) Store the result in the output array.

    (J) Free memory and return.

    ----------------------------------------------------------------------------

    The structure above does not make use of generalized contractions. For that
    we would also need an alternative Gaussian integral code, which should be
    able compute integrals for different angular momenta at the same time,
    making use of shared intermediate results.

    TODO: The following steps are needed and can be tested independently.

        - Write iterators for the double loops.
        - Write transformation routines from Cartesian to pure basis functions.
        - Routine to copy results in working memory to output array.

    */

    int nwork, result;
    double* work_cart;
    double* work_pure;
    i2gob_type* i2;
    i2pow_type* i2p;

    // initialize some variables

    result = 0;
    work_cart = NULL;
    work_pure = NULL;
    i2 = NULL;
    i2p = NULL;

    // (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A)

    // Allocate iterators
    i2 = i2gob_new();
    if (i2==NULL) safe_exit(-1);
    i2p = i2pow_new();
    if (i2p==NULL) safe_exit(-1);

    // i2gob stands for double iterator for gaussian orbital basis sets. For
    // now there is just one c-ish implementation of this iterator. We may
    // need to make a c++-ish iterator in future if we to generalize this
    // routine for different storage types. Note that the iterator is designed
    // for the nesting of the loops below.
    result = i2gob_init(i2, centers, shell_map, ncons, nexps, con_types,
                        exponents, con_coeffs, nshell, ncenter, ncon_total,
                        nexp_total, ncon_coeff_total, nbasis);
    if (result != 0) goto exit;

    // Allocate work arrays.
    nwork = (*i2).max_nbasis*(*i2).max_nbasis;
    work_cart = malloc(nwork*sizeof(double));
    if (work_cart==NULL) safe_exit(-1);
    work_pure = malloc(nwork*sizeof(double));
    if (work_pure==NULL) safe_exit(-1);


    // (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B)
    i2gob_update_shell(i2);
    do {

        // (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C)
        i2gob_update_con(i2);
        do {
            // (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D)
            // We make use of the fact that a floating point zero consists of
            // consecutive zero bytes.
            memset(work_cart, 0, nwork*sizeof(double));
            memset(work_pure, 0, nwork*sizeof(double));

            // (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E)
            i2gob_update_exp(i2);
            do {

                // (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F)
                i2pow_init(i2p, abs((*i2).con_type0), abs((*i2).con_type1), (*i2).max_nbasis);
                do {
                    work_cart[(*i2p).offset] += (*i2).con_coeff*overlap(
                        (*i2).exp0, (*i2p).l0, (*i2p).m0, (*i2p).n0,
                        (*i2).x0, (*i2).y0, (*i2).z0,
                        (*i2).exp1, (*i2p).l1, (*i2p).m1, (*i2p).n1,
                        (*i2).x1, (*i2).y1, (*i2).z1
                    );
                } while (i2pow_inc(i2p));
            } while (i2gob_inc_exp(i2));

            // (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H)
            if ((*i2).con_type0 < -1) {
                project_cartesian_to_pure(work_cart, work_pure, (*i2).con_type0,
                    1, // stride
                    (*i2).max_nbasis, // spacing
                    get_con_nbasis(abs((*i2).con_type1)) // count
                );
                swap_ptr(&work_cart, &work_pure);
            }

            if ((*i2).con_type1 < -1) {
                project_cartesian_to_pure(work_cart, work_pure, (*i2).con_type1,
                    (*i2).max_nbasis, // stride
                    1, // spacing
                    get_con_nbasis((*i2).con_type0) // count
                );
                swap_ptr(&work_cart, &work_pure);
            }

            // make sure the final result after projections is in work_pure.
            swap_ptr(&work_cart, &work_pure);

            // (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I)
            i2gob_store(i2, work_pure, output);
        } while (i2gob_inc_con(i2));
    } while (i2gob_inc_shell(i2));

exit:
    // (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J)

    free(work_cart);
    free(work_pure);
    i2gob_free(i2);
    i2pow_free(i2p);
    return result;
}



i2gob_type* i2gob_new(void) {
    return malloc(sizeof(i2gob_type));
}


void i2gob_free(i2gob_type* i2) {
    free(i2);
}


int i2gob_init(i2gob_type* i2, double* centers, long* shell_map,
    long* ncons, long* nexps, long* con_types, double* exponents,
    double* con_coeffs, long nshell, long ncenter, long ncon_total,
    long nexp_total, long ncon_coeff_total, long nbasis) {

    int result;

    // assign input variables
    (*i2).centers = centers;
    (*i2).shell_map = shell_map;
    (*i2).ncons = ncons;
    (*i2).nexps = nexps;
    (*i2).con_types = con_types;
    (*i2).exponents = exponents;
    (*i2).con_coeffs = con_coeffs;
    (*i2).nshell = nshell;

    result = i2gob_check(i2, nshell, ncenter, ncon_total, nexp_total,
                         ncon_coeff_total, nbasis);
    if (result < 0) return result;

    // reset public fields
    (*i2).max_nbasis = result;
    (*i2).con_type0 = 0;
    (*i2).con_type1 = 0;
    (*i2).con_coeff = 0.0;
    (*i2).exp0 = 0.0;
    (*i2).exp1 = 0.0;
    (*i2).x0 = 0.0;
    (*i2).y0 = 0.0;
    (*i2).z0 = 0.0;
    (*i2).x1 = 0.0;
    (*i2).y1 = 0.0;
    (*i2).z1 = 0.0;

    // reset internal fields
    (*i2).ishell0 = 0;
    (*i2).ishell1 = 0;

    (*i2).ncon0 = 0; // number of contractions in a shell
    (*i2).ncon1 = 0;
    (*i2).icon0 = 0; // index of the current con_type within the current shell
    (*i2).icon1 = 0;
    (*i2).ocon0 = 0; // offset of the first con_type within the current shell
    (*i2).ocon1 = 0;

    (*i2).nexp0 = 0; // number of exponents in a shell
    (*i2).nexp1 = 0;
    (*i2).iexp0 = 0; // index of the current exponent within the current shell
    (*i2).iexp1 = 0;
    (*i2).oexp0 = 0; // offset of the first exponent within the current shell
    (*i2).oexp1 = 0;

    (*i2).occ0 = 0; // offset of the first contraction coefficient within the current shell
    (*i2).occ1 = 0;

    return 0;
}


int i2gob_check(i2gob_type* i2, long nshell, long ncenter, long ncon_total,
    long nexp_total, long ncon_coeff_total, long nbasis) {

    /*

    This routine computes the maximum number of basis functions associated with
    a contraction over the entire basis set.

    Besides this principal purpose, this routine also checks for all possible
    inconsistencies in the basis set description:

        -2: the center indexes in the shell_map are out of range.
        -3: an element of ncons is less than one.
        -4: ncon_total does not match the sum of all values in ncons.
        -5: an element of nexps is less than one.
        -6: nexp_total does not match the sum of all values in nexps.
        -7: con_type -1 occurs.
        -8: ncon_coeff_total does not match the total number of expected
            con_coeffs.
        -9: the total number of basis functions is incorrect.

    */

    long i, j, max_angmom, i_con, i_exp, i_con_coeff, i_basis, con_type;

    max_angmom = 0;
    i_con = 0;
    i_exp = 0;
    i_con_coeff = 0;
    i_basis = 0;
    for (i=0; i<nshell; i++) {
        if (((*i2).shell_map[i] < 0) || ((*i2).shell_map[i] >= ncenter)) return -2;
        if ((*i2).ncons[i] <= 0) return -3;
        for (j=0; j < (*i2).ncons[i]; j++) {
            con_type = (*i2).con_types[i_con];
            if (con_type == -1) return -7;
            if (abs(con_type) > max_angmom) max_angmom = abs(con_type);
            i_con++;
            if (i_con > ncon_total) return -4;
            i_basis += get_con_nbasis(con_type);
            if (i_basis > nbasis) return -9;
        }
        if ((*i2).nexps[i] <= 0) return -5;
        i_exp += (*i2).nexps[i];
        if (i_exp > nexp_total) return -6;
        i_con_coeff += (*i2).nexps[i]*(*i2).ncons[i];
        if (i_con_coeff > ncon_coeff_total) return -8;
    }
    if (i_con != ncon_total) return -4;
    if (i_exp != nexp_total) return -6;
    if (i_con_coeff != ncon_coeff_total) return -8;
    if (i_basis != nbasis) return -9;
    return ((max_angmom+1)*(max_angmom+2))/2;
}


int i2gob_inc_shell(i2gob_type* i2) {
    // Increment shell and related counters.
    if ((*i2).ishell0 < (*i2).ishell1) {
        (*i2).ishell0++;
        (*i2).ocon0 += (*i2).ncon0;
        (*i2).oexp0 += (*i2).nexp0;
        (*i2).occ0 += (*i2).ncon0 * (*i2).nexp0;
        i2gob_update_shell(i2);
        return 1;
    } else if ((*i2).ishell1 < (*i2).nshell-1) {
        (*i2).ishell0 = 0;
        (*i2).ocon0 = 0;
        (*i2).oexp0 = 0;
        (*i2).occ0 = 0;
        (*i2).ocon1 += (*i2).ncon1;
        (*i2).oexp1 += (*i2).nexp1;
        (*i2).occ1 += (*i2).ncon1 * (*i2).nexp1;
        (*i2).ishell1++;
        i2gob_update_shell(i2);
        return 1;
    } else {
        (*i2).ishell0 = 0;
        (*i2).ishell1 = 0;
        (*i2).ocon0 = 0;
        (*i2).ocon1 = 0;
        (*i2).oexp0 = 0;
        (*i2).oexp1 = 0;
        (*i2).occ0 = 0;
        (*i2).occ1 = 0;
        i2gob_update_shell(i2);
        return 0;
    }
}


void i2gob_update_shell(i2gob_type* i2) {
    double* tmp;
    // Update fields that depend on shell and related counters.
    (*i2).ncon0 = (*i2).ncons[(*i2).ishell0];
    (*i2).ncon1 = (*i2).ncons[(*i2).ishell1];
    (*i2).nexp0 = (*i2).nexps[(*i2).ishell0];
    (*i2).nexp1 = (*i2).nexps[(*i2).ishell1];
    // update center 0
    tmp = (*i2).centers + 3*(*i2).shell_map[(*i2).ishell0];
    (*i2).x0 = *tmp;
    tmp++;
    (*i2).y0 = *tmp;
    tmp++;
    (*i2).z0 = *tmp;
    // update center 1
    tmp = (*i2).centers + 3*(*i2).shell_map[(*i2).ishell1];
    (*i2).x1 = *tmp;
    tmp++;
    (*i2).y1 = *tmp;
    tmp++;
    (*i2).z1 = *tmp;
    // reset contraction counters
    (*i2).icon0 = 0;
    (*i2).icon1 = 0;
}


int i2gob_inc_con(i2gob_type* i2) {
    // Increment contraction counters.
    if ((*i2).icon0 < (*i2).ncon0-1) {
        (*i2).icon0++;
        i2gob_update_con(i2);
        return 1;
    } else if ((*i2).icon1 < (*i2).ncon1-1) {
        (*i2).icon1++;
        (*i2).icon0 = 0;
        i2gob_update_con(i2);
        return 1;
    } else {
        (*i2).icon0 = 0;
        (*i2).icon1 = 0;
        i2gob_update_con(i2);
        return 0;
    }
}


void i2gob_update_con(i2gob_type* i2) {
    // Update fields that depend on contraction counters.
    (*i2).con_type0 = (*i2).con_types[(*i2).ocon0 + (*i2).icon0];
    (*i2).con_type1 = (*i2).con_types[(*i2).ocon1 + (*i2).icon1];
    // Reset exponent counters.
    (*i2).iexp0 = 0;
    (*i2).iexp1 = 0;
}


int i2gob_inc_exp(i2gob_type* i2) {
    // Increment exponent counters.
    if ((*i2).iexp0 < (*i2).nexp0-1) {
        (*i2).iexp0++;
        i2gob_update_exp(i2);
        return 1;
    } else if ((*i2).iexp1 < (*i2).nexp1-1) {
        (*i2).iexp1++;
        (*i2).iexp0 = 0;
        i2gob_update_exp(i2);
        return 1;
    } else {
        (*i2).iexp0 = 0;
        (*i2).iexp1 = 0;
        i2gob_update_exp(i2);
        return 0;
    }
}


void i2gob_update_exp(i2gob_type* i2) {
    // Update fields that depend on exponent counters.
    (*i2).exp0 = (*i2).exponents[(*i2).oexp0 + (*i2).iexp0];
    (*i2).exp1 = (*i2).exponents[(*i2).oexp1 + (*i2).iexp1];
    (*i2).con_coeff = (*i2).con_coeffs[(*i2).occ0 + (*i2).ncon0*(*i2).iexp0 + (*i2).icon0]*
                      (*i2).con_coeffs[(*i2).occ1 + (*i2).ncon1*(*i2).iexp1 + (*i2).icon1];
}

void i2gob_store(i2gob_type* i2, double *work_pure, double *output) {
}





i2pow_type* i2pow_new(void) {
    return malloc(sizeof(i2pow_type));
}


void i2pow_free(i2pow_type* i2p) {
    free(i2p);
}


void i2pow_init(i2pow_type* i2p, long con_type0, long con_type1, long max_nbasis) {
    assert(con_type0 >= 0);
    assert(con_type1 >= 0);
    (*i2p).con_type0 = con_type0;
    (*i2p).con_type1 = con_type1;
    (*i2p).skip = max_nbasis - get_con_nbasis(con_type0) + 1;
    assert((*i2p).skip >= 0);
    assert(max_nbasis >= get_con_nbasis(con_type1));
    (*i2p).l0 = con_type0;
    (*i2p).m0 = 0;
    (*i2p).n0 = 0;
    (*i2p).l1 = con_type1;
    (*i2p).m1 = 0;
    (*i2p).n1 = 0;
    (*i2p).offset = 0;
}


int i2pow_inc(i2pow_type* i2p) {
    // Increment indexes of shell 0
    int result;
    result = i1pow_inc(&(*i2p).l0, &(*i2p).m0, &(*i2p).n0);
    if (result) {
        (*i2p).offset++;
    } else {
        result = i1pow_inc(&(*i2p).l1, &(*i2p).m1, &(*i2p).n1);
        if (result) {
            (*i2p).offset += (*i2p).skip;
        } else {
            (*i2p).offset = 0;
        }
    }
    return result;
}


int i1pow_inc(int* l, int* m, int* n) {
    // Modify the indexes in place as to move to the next combination of powers
    // withing one angular momentum.
    // Note: con_type = (*l) + (*m) + (*n);
    if (*m == 0) {
      if (*l == 0) {
        *l = *n;
        *n = 0;
        return 0;
      } else {
        *m = *n + 1;
        *n = 0;
        (*l)--;
      }
    } else {
      (*m)--;
      (*n)++;
    }
    return 1;
}


void project_cartesian_to_pure(double *work_cart, double* work_pure, long
    con_type, long stride, long spacing, long count) {
}
