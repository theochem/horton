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


int compute_gobasis_overlap(double* centers, long* shell_map,
    long* num_exponents, long* num_contractions, long* con_types,
    double* exponents, double* con_coeffs, double* output) {

    /*

    This routine has the following structure:

    ----------------------------------------------------------------------------

    (A) Double loop over the shells:

        (B) Double loop over the contractions:

            (C) Clear working memory for storage of Cartesian Gaussian integrals.

            (D) Double loop over the exponents:

                (E) Double loop over basis functions associated with the angular
                    momenta of each contraction

                    (F) Add contribution from a pair of primitives to the
                        working memory. This part contains a call to a function
                        from cints.c

            (G) If needed, perform projection from Cartesian to pure basis
                functions. (A second piece of working memory is needed for this
                transformation.)

            (H) Store the result in the output array.

    ----------------------------------------------------------------------------

    The structure above does not make use of generalized contractions. For that
    we would also need an alternative Gaussian integral code, which should be
    able compute integrals for different angular momenta at the same time,
    making use of shared intermediate results.

    TODO: everything. :) The following steps are needed and can be tested
    independently.

        - Define a maximum angular momentum such that the size of the working
          memory can be fixed.
        - Allocate working memory in advance. Use goto statement to an
          exit clause that frees all memory in case of some error.
        - Write iterators for the double loops.
        - Write transformation routines from Cartesian to pure basis functions.
        - Routine to copy results in working memory to output array.

    */

    return 0;
}
