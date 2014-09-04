// Horton is a development platform for electronic structure methods.
// Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

#ifndef GBW4_H
#define GBW4_H

#include "gbasis.h"
#include "ints.h"

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
    public:
        GB4IntegralWrapper(GOBasis* gobasis, GB4Integral* gb4int);
        ~GB4IntegralWrapper();
        long get_nbasis() {return gobasis->get_nbasis();}
        void select_2index(long index0, long index2,
                                    long* pbegin0, long* pend0,
                                    long* pbegin2, long* pend2);
        //void compute(bool* mask);
        void compute();
        void compute_diagonal(double* diagonal);
        double* get_2index_slice(long index0, long index2);
};

#endif
