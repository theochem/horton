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
#include <math.h>
#include <stdlib.h>


#define safe_exit(x) {result=x; goto exit;}


int compute_gobasis_overlap(double* centers, long* shell_map,
    long* num_exponents, long* num_contractions, long* con_types,
    double* exponents, double* con_coeffs, long num_shells,
    long tot_num_con, double* output) {

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
                    momenta of each contraction

                    (G) Add contribution from a pair of primitives to the
                        working memory. This part contains a call to a function
                        from cints.c

            (H) If needed, perform projection from Cartesian to pure basis
                functions. (A second piece of working memory is needed for this
                transformation.)

            (I) Store the result in the output array.

    (J) Free memory and return

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

    int max_nbasis, result;
    double* work_cart;
    double* work_pure;

    result = 0;
    work_cart = NULL;
    work_pure = NULL;

    // (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A) (A)

    // Get the maximum number of basis functions associated with a contraction.
    max_nbasis = get_max_nbasis(num_contractions, con_types, num_shells,
                                tot_num_con);
    if (max_nbasis < 0) safe_exit(max_nbasis);

    // Allocate work arrays.
    work_cart = (double *)malloc(max_nbasis*max_nbasis*sizeof(double));
    if (work_cart==NULL) safe_exit(-1);
    work_pure = (double *)malloc(max_nbasis*max_nbasis*sizeof(double));
    if (work_pure==NULL) safe_exit(-1);

exit:
    // (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J) (J)

    free(work_cart);
    free(work_pure);
    return result;
}


int get_max_nbasis(long* num_contractions, long* con_types, long num_shells,
    long tot_num_con) {

    /*

    This routine computes the maximum number of basis functions associated with
    a contraction over the entire basis set.

    Besides this principal purpose, this routine also asserts two expected
    properties of the input data:

        -2: each element of num_contractions is at least 1
        -3: con_type -1 does not occur
        -4: tot_num_con matches the sum of all values in num_contractions

    */

    long i, j, max_angmom, counter;

    max_angmom = 0;
    counter = 0;
    for (i=0; i<num_shells; i++) {
        if ((*num_contractions) <= 0) return -2;
        for (j=0; j < (*num_contractions); j++) {
            if ((*con_types) == -1) return -3;
            if (abs(*con_types) > max_angmom) max_angmom = abs(*con_types);
            con_types++;
            counter++;
            if (counter > tot_num_con) return -4;
        }
        num_contractions++;
    }
    return ((max_angmom+1)*(max_angmom+2))/2;
}
