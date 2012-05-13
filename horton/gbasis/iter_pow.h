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


#ifndef HORTON_GBASIS_ITERPOW_H
#define HORTON_GBASIS_ITERPOW_H

int iter_pow1_inc(long* n);

class IterPow2 {
    private:
        // TODO: remove skip and offset
        long shell_type0, shell_type1, skip;
    public:
        void reset(long shell_type0, long shell_type1, long max_nbasis);
        int inc();
        long n0[3], n1[3], offset, ibasis0, ibasis1;
    };

#endif
