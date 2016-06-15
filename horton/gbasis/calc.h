// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2016 The HORTON Development Team
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

// UPDATELIBDOCTITLE: Base class for any integral/evaluation of Gaussian functions

#ifndef HORTON_GBASIS_CALC_H
#define HORTON_GBASIS_CALC_H


/*

    Base class for anything that computes stuff by iterating over Gaussian
    basis functions

*/


class GBCalculator {
    protected:
        long nwork, max_shell_type, max_nbasis;
        double *work_pure, *work_cart; // contiguous work arrays sufficiently large for max_shell_type
        void swap_work();
    public:
        GBCalculator(long max_shell_type);
        virtual ~GBCalculator();
        const long get_nwork() const {return nwork;};
        const long get_max_shell_type() const {return max_shell_type;};
        const long get_max_nbasis() const {return max_nbasis;};
        const double* get_work() const {return work_cart;};
    };


#endif
