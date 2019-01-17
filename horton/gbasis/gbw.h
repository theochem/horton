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

// UPDATELIBDOCTITLE: A four-center integral wrapper for the Cholesky code

#ifndef GBW_H
#define GBW_H

#include "horton/gbasis/gbasis.h"
#include "horton/gbasis/ints.h"

/**
    @brief
        A wrapper around a four-center integral implementation that is suitable
        for a Cholesky algorithm.
*/
class GB4IntegralWrapper {
    private:
        GOBasis* gobasis;
        GB4Integral* gb4int;
        long max_shell_size;
        long slice_size;
        double* integrals;

        long ishell0;
        long ishell2;
        long begin0; // beginning of the basis indexes for shell0
        long begin2; // beginning of the basis indexes for shell2

        /**
            @brief
                Compute four-center integrals for a quadruplet of shells
                (defined by ishell0, ishell1, ishell2 and ishell3).
        */
        void compute_shell(long ishell1, long ishell3);
    public:
        /**
            @brief
                Construct the wrapper.

            @param gobasis
                A Gaussian basis set for which the four-center integrals are
                to be computed.

            @param gb4int
                A definition/implementation of a four-center integral.
        */
        GB4IntegralWrapper(GOBasis* gobasis, GB4Integral* gb4int);
        ~GB4IntegralWrapper();

        /**
            @brief
                The number of basis functions.
        */
        long get_nbasis() {return gobasis->get_nbasis();}

        /**
            @brief
                Select a pair of shells to which a pair of basis indexes belong.

            The output of this function is also stored internally for the next
            call to the compute method.

            @param index0
                The first basis index.

            @param index2
                The third basis index. (Physicist notation!)

            @param pbegin0
                Output argument. The begin of the range for the first shell.

            @param pend0
                Output argument. The end (non-inclusive) of the range for the
                first shell.

            @param pbegin2
                Output argument. The begin of the range for the third shell.
                (Physicist notation!)

            @param pend2
                Output argument. The end (non-inclusive) of the range for the
                third shell. (Physicist notation!)
        */
        void select_2index(long index0, long index2,
                           long* pbegin0, long* pend0,
                           long* pbegin2, long* pend2);

        /**
            @brief
                Compute four-center integrals for the slices selected with
                the select_2index method.
        */
        void compute();

        /**
            @brief
                Compute the (double) diagonal of the four-index object. This
                is usually needed in the initialization of a Cholesky algorithm.

            @param diagonal
                The output array to which the result is written
                (size=nbasis*basis).
        */
        void compute_diagonal(double* diagonal);

        /**
            @brief
                Get a slice from the four-center integral computed with the last
                call to the compute method.

            @param index0
                The first index of the slice. (Absolute indexes, not relative
                within the currently selected shell.)

            @param index2
                The third index of the slice. (Physicist notation!) (Absolute
                indexes, not relative within the currently selected shell.)
        */
        double* get_2index_slice(long index0, long index2);
};

#endif
