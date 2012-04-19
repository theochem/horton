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
#include "cartpure.h"
#include "cints.h"

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define safe_exit(x) {result=x; goto exit;}


long get_shell_nbasis(long shell_type) {
    if (shell_type > 0) {
        // Cartesian
        return (shell_type+1)*(shell_type+2)/2;
    } else if (shell_type == -1) {
        // should not happen.
        return -1;
    } else {
        // Pure
        return -2*shell_type+1;
    }
}

void swap_ptr(double **x, double **y) {
    double *t = *x;
    *x = *y;
    *y = t;
}


int compute_gobasis_overlap(double* centers, long* shell_map, long* nprims,
    long* shell_types, double* alphas, double* con_coeffs, long nshell,
    long ncenter, long nprim_total, double* output, long nbasis) {

    /* TODO
             This routine will eventually become a general contraction routine
             for most one-body operators. For now, it just does the overlap
             matrix. In future, it will be renamed and an argument will be added
             (datatype or class) that contains the operator specific
             information.
    */

    /*

    This routine has the following structure:

    ----------------------------------------------------------------------------

    (A) Allocated working memory, based on largest angular momentum in the basis
        set.

    (B) Double loop over the shells:

        (C) Clear working memory for storage of Cartesian Gaussian integrals.

        (D) Double loop over the alphas:

            (E) Double loop over basis functions associated with the angular
                momenta of each contraction:

                (F) Add contribution from a pair of primitives to the
                    working memory. This part contains a call to a function
                    from cints.c.

        (G) If needed, perform projection from Cartesian to pure basis
            functions. (A second piece of working memory is needed for this
            transformation.)

        (H) Store the result in the output array.

    (I) Free memory and return.

    ----------------------------------------------------------------------------

    The structure above does not make use of generalized contractions. For that
    we would also need an alternative Gaussian integral code, which should be
    able compute integrals for different angular momenta at the same time,
    making use of shared intermediate results.

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
    result = i2gob_init(i2, centers, shell_map, nprims, shell_types, alphas,
                        con_coeffs, nshell, ncenter, nprim_total, nbasis);
    if (result != 0) goto exit;

    // Allocate work arrays.
    nwork = i2->max_nbasis*i2->max_nbasis;
    work_cart = malloc(nwork*sizeof(double));
    if (work_cart==NULL) safe_exit(-1);
    work_pure = malloc(nwork*sizeof(double));
    if (work_pure==NULL) safe_exit(-1);


    // (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B) (B)
    i2gob_update_shell(i2);
    do {
        // (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C) (C)
        // We make use of the fact that a floating point zero consists of
        // consecutive zero bytes.
        memset(work_cart, 0, nwork*sizeof(double));
        memset(work_pure, 0, nwork*sizeof(double));

        // (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D) (D)
        i2gob_update_prim(i2);
        do {

            // (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E) (E)
            /* TODO:
               This loop should go into the integral evaluation routine
               to avoid trivially redundant floating point operations.
               It would be yet more useful to have a data structure
               for each type of Gaussian integral, which can be initialized
               prior to this loop.
            */
            i2pow_init(i2p, abs(i2->shell_type0), abs(i2->shell_type1), i2->max_nbasis);
            do {
                // (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F) (F)
                work_cart[i2p->offset] += i2->con_coeff*gob_overlap(
                    i2->alpha0, i2p->nx0, i2p->ny0, i2p->nz0, i2->r0,
                    i2->alpha1, i2p->nx1, i2p->ny1, i2p->nz1, i2->r1
                )*gob_normalization(i2->alpha0, i2p->nx0, i2p->ny0, i2p->nz0)
                *gob_normalization(i2->alpha1, i2p->nx1, i2p->ny1, i2p->nz1);
            } while (i2pow_inc(i2p));
        } while (i2gob_inc_prim(i2));

        // (G) (G) (G) (G) (G) (G) (G) (G) (G) (G) (G) (G) (G) (G) (G) (G)  (G)
        if (i2->shell_type0 < -1) {
            project_cartesian_to_pure(work_cart, work_pure, -i2->shell_type0,
                1, // stride
                i2->max_nbasis, // spacing
                get_shell_nbasis(abs(i2->shell_type1)) // count
            );
            swap_ptr(&work_cart, &work_pure);
        }

        if (i2->shell_type1 < -1) {
            project_cartesian_to_pure(work_cart, work_pure, -i2->shell_type1,
                i2->max_nbasis, // stride
                1, // spacing
                get_shell_nbasis(i2->shell_type0) // count
            );
            swap_ptr(&work_cart, &work_pure);
        }

        // make sure the final result after projections is in work_pure.
        swap_ptr(&work_cart, &work_pure);

        // (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H) (H)
        i2gob_store(i2, work_pure, output);
    } while (i2gob_inc_shell(i2));

exit:
    // (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I) (I)

    free(work_cart);
    free(work_pure);
    i2gob_free(i2);
    i2pow_free(i2p);
    return result;
}


double gob_normalization(double exponent, int nx, int ny, int nz) {
    return sqrt(pow(4.0*exponent, nx+ny+nz)*pow(2.0*exponent/M_PI, 1.5)
           /(fac2(2*nx-1)*fac2(2*ny-1)*fac2(2*nz-1)));
}



i2gob_type* i2gob_new(void) {
    return malloc(sizeof(i2gob_type));
}


void i2gob_free(i2gob_type* i2) {
    free(i2->basis_offsets);
    free(i2);
}


int i2gob_init(i2gob_type* i2, double* centers, long* shell_map, long* nprims,
    long* shell_types, double* alphas, double* con_coeffs, long nshell,
    long ncenter, long nprim_total, long nbasis) {

    long ishell, shell_nbasis;

    // assign input variables
    i2->centers = centers;
    i2->shell_map = shell_map;
    i2->nprims = nprims;
    i2->shell_types = shell_types;
    i2->basis_offsets = NULL;
    i2->alphas = alphas;
    i2->con_coeffs = con_coeffs;
    i2->nshell = nshell;
    i2->nbasis = nbasis;

    // allocate and set up auxilliary array
    i2->basis_offsets = malloc(nshell*sizeof(long));
    if (i2->basis_offsets == NULL) return -1;
    i2->basis_offsets[0] = 0;
    i2->max_nbasis = 0;
    for (ishell=1; ishell<nshell; ishell++) {
        shell_nbasis = get_shell_nbasis(i2->shell_types[ishell-1]);
        i2->basis_offsets[ishell] = i2->basis_offsets[ishell-1] + shell_nbasis;
        shell_nbasis = get_shell_nbasis(abs(i2->shell_types[ishell-1]));
        if (i2->max_nbasis < shell_nbasis) {
            i2->max_nbasis = shell_nbasis;
        }
    }

    // reset public fields
    i2->shell_type0 = 0;
    i2->shell_type1 = 0;
    i2->con_coeff = 0.0;
    i2->alpha0 = 0.0;
    i2->alpha1 = 0.0;
    i2->r0 = NULL;
    i2->r1 = NULL;
    i2->ibasis0 = 0;
    i2->ibasis1 = 0;

    // reset internal fields
    i2->ishell0 = 0;
    i2->ishell1 = 0;

    i2->nprim0 = 0; // number of alphas in a shell
    i2->nprim1 = 0;
    i2->iprim0 = 0; // index of the current exponent within the current shell
    i2->iprim1 = 0;
    i2->oprim0 = 0; // offset of the first exponent within the current shell
    i2->oprim1 = 0;

    return 0;
}


int i2gob_inc_shell(i2gob_type* i2) {
    // Increment shell and related counters.
    if (i2->ishell1 < i2->ishell0) {
        i2->ishell1++;
        i2->oprim1 += i2->nprim1;
        i2gob_update_shell(i2);
        return 1;
    } else if (i2->ishell0 < i2->nshell-1) {
        i2->ishell1 = 0;
        i2->oprim1 = 0;
        i2->oprim0 += i2->nprim0;
        i2->ishell0++;
        i2gob_update_shell(i2);
        return 1;
    } else {
        i2->ishell1 = 0;
        i2->ishell0 = 0;
        i2->oprim0 = 0;
        i2->oprim1 = 0;
        i2gob_update_shell(i2);
        return 0;
    }
}


void i2gob_update_shell(i2gob_type* i2) {
    // Update fields that depend on shell and related counters.
    i2->nprim0 = i2->nprims[i2->ishell0];
    i2->nprim1 = i2->nprims[i2->ishell1];
    // Update indexes in output array
    i2->ibasis0 = i2->basis_offsets[i2->ishell0];
    i2->ibasis1 = i2->basis_offsets[i2->ishell1];
    // update centers
    i2->r0 = i2->centers + 3*i2->shell_map[i2->ishell0];
    i2->r1 = i2->centers + 3*i2->shell_map[i2->ishell1];
    // update shell types
    i2->shell_type0 = i2->shell_types[i2->ishell0];
    i2->shell_type1 = i2->shell_types[i2->ishell1];
    // reset contraction counters
    i2->iprim0 = 0;
    i2->iprim1 = 0;
}


int i2gob_inc_prim(i2gob_type* i2) {
    // Increment exponent counters.
    if (i2->iprim1 < i2->nprim1-1) {
        i2->iprim1++;
        i2gob_update_prim(i2);
        return 1;
    } else if (i2->iprim0 < i2->nprim0-1) {
        i2->iprim0++;
        i2->iprim1 = 0;
        i2gob_update_prim(i2);
        return 1;
    } else {
        i2->iprim0 = 0;
        i2->iprim1 = 0;
        i2gob_update_prim(i2);
        return 0;
    }
}


void i2gob_update_prim(i2gob_type* i2) {
    // Update fields that depend on exponent counters.
    i2->alpha0 = i2->alphas[i2->oprim0 + i2->iprim0];
    i2->alpha1 = i2->alphas[i2->oprim1 + i2->iprim1];
    i2->con_coeff = i2->con_coeffs[i2->oprim0 + i2->iprim0]*
                      i2->con_coeffs[i2->oprim1 + i2->iprim1];
}


void i2gob_store(i2gob_type* i2, double *work_pure, double *output) {
    // TODO: this routine is hardwired to work only for the dense storage
    long i0, i1, n0, n1;
    n0 = get_shell_nbasis(i2->shell_type0);
    n1 = get_shell_nbasis(i2->shell_type1);
    for (i0=0; i0<n0; i0++) {
        for (i1=0; i1<n1; i1++) {
            output[(i0+i2->ibasis0)*i2->nbasis+i1+i2->ibasis1] = work_pure[i0*i2->max_nbasis+i1];
        }
    }
    if (i2->ibasis0 != i2->ibasis1) {
        for (i0=0; i0<n0; i0++) {
            for (i1=0; i1<n1; i1++) {
                output[(i1+i2->ibasis1)*i2->nbasis+i0+i2->ibasis0] = work_pure[i0*i2->max_nbasis+i1];
            }
        }
    }
}



i2pow_type* i2pow_new(void) {
    return malloc(sizeof(i2pow_type));
}


void i2pow_free(i2pow_type* i2p) {
    free(i2p);
}


void i2pow_init(i2pow_type* i2p, long shell_type0, long shell_type1, long max_nbasis) {
    assert(shell_type0 >= 0);
    assert(shell_type1 >= 0);
    i2p->shell_type0 = shell_type0;
    i2p->shell_type1 = shell_type1;
    i2p->skip = max_nbasis - get_shell_nbasis(shell_type1) + 1;
    assert(i2p->skip >= 0);
    assert(max_nbasis >= get_shell_nbasis(shell_type1));
    i2p->nx0 = shell_type0;
    i2p->ny0 = 0;
    i2p->nz0 = 0;
    i2p->nx1 = shell_type1;
    i2p->ny1 = 0;
    i2p->nz1 = 0;
    i2p->offset = 0;
}


int i2pow_inc(i2pow_type* i2p) {
    // Increment indexes of shell 0
    int result;
    result = i1pow_inc(&i2p->nx1, &i2p->ny1, &i2p->nz1);
    if (result) {
        i2p->offset++;
    } else {
        result = i1pow_inc(&i2p->nx0, &i2p->ny0, &i2p->nz0);
        if (result) {
            i2p->offset += i2p->skip;
        } else {
            i2p->offset = 0;
        }
    }
    return result;
}


int i1pow_inc(int* nx, int* ny, int* nz) {
    // Modify the indexes in place as to move to the next combination of powers
    // withing one angular momentum.
    // Note: shell_type = (*nx) + (*ny) + (*nz);
    if (*ny == 0) {
      if (*nx == 0) {
        *nx = *nz;
        *nz = 0;
        return 0;
      } else {
        *ny = *nz + 1;
        *nz = 0;
        (*nx)--;
      }
    } else {
      (*ny)--;
      (*nz)++;
    }
    return 1;
}
