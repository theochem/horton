// Horton is a Density Functional Theory program.
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


#include <cmath>
#include "nucpot.h"

void compute_grid_nucpot(long* numbers, double* coordinates, long natom,
                         double* points, double* output, long npoint) {
    for (long ipoint=npoint-1; ipoint >= 0; ipoint--) {
        double tmp = 0.0;
        for (long iatom=natom-1; iatom >= 0; iatom--) {
            double delta = points[0]-coordinates[3*iatom];
            double dsq = delta*delta;
            delta = points[1]-coordinates[3*iatom+1];
            dsq += delta*delta;
            delta = points[2]-coordinates[3*iatom+2];
            dsq += delta*delta;
            tmp += numbers[iatom]/sqrt(dsq);
        }
        *output += tmp;

        points += 3;
        output++;
    }
}
