// Horton is a development platform for electronic structure methods.
// Copyright (C) 2011-2015 The Horton Development Team
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

// UPDATELIBDOCTITLE: Cholesky decomposition of (any) four-center integrals

#ifndef CHOLESKY_H
#define CHOLESKY_H

#include <vector>
#include "horton/gbasis/gbw.h"

// Include the CBLAS headers
#ifdef BLAS_MKL
#include <mkl.h>
#else
extern "C"
{
#include <cblas.h>
}
#endif

/**
    @brief
        Computes Cholesky vectors for a four-index object

    Only the 4-center integrals relevant for the decomposition are actually
    computed. This implementation computes slices of the four-index object for
    a pair of shells at a time. (This is because most implementations of a
    four-center work like that.)

    @param gbw4
        A wrapper around a definition of the 4-center integral. See gbw.h

    @param uninit_result
        An output pointer. The Cholesky vectors will be allocated as part of
        this routine and the pointer to the Cholesky vectors is assigned to this
        output argument.

    @param threshold
        A threshold for the error on the (double) diagonal of the four-center
        object. The Cholesky decomposition stops when sufficient vectors are
        generated such that the error on the diagonal falls below this
        threshold.
*/
long cholesky(GB4IntegralWrapper* gbw4, double** uninit_result,
    double threshold);

#endif
