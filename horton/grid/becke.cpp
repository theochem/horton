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


//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif
#include <cmath>
#include <stdexcept>

#include "becke.h"


/* dist

   compute the Euclidian distance between two points.
*/
static double dist(double* p0, double* p1) {
    double tmp, result;
    tmp = (p0[0] - p1[0]);
    result = tmp*tmp;
    tmp = (p0[1] - p1[1]);
    result += tmp*tmp;
    tmp = (p0[2] - p1[2]);
    result += tmp*tmp;
    return sqrt(result);
}


/* becke_helper_atom

   Computes the Becke weighting function for every point in the grid

   npoint
        The number of grid points

   points
        The grid points. row-major storate of (npoint,3) array.

   weights
        The output, i.e. the Becke weights for the grid points.

   natom
        The number of atoms

   radii
        The radii of the atoms

   centers
        The positions of the atoms.

   select
        The selected atom for which the Becke weights are to be computed.

   order
        The order of the switching function in the Becke scheme.

   return value
*/
void becke_helper_atom(int npoint, double* points, double* weights, int natom,
                       double* radii, double* centers, int select, int order)
{
    double nom, denom; // The nominator and the denominator in the weight definition
    double alpha, p, s; // Used to build up the value of the switching function

    for (int ipoint = npoint-1; ipoint>=0; ipoint--) {
        nom = 0;
        denom = 0;
        for (int iatom0 = 0; iatom0 < natom; iatom0++) {
            p = 1;
            for (int iatom1 = 0; iatom1 < natom; iatom1++) {
                if (iatom0 == iatom1) continue;

                // TODO: move the following block out of the loops
                alpha = (radii[iatom0] - radii[iatom1])/(radii[iatom0] + radii[iatom1]);
                alpha = alpha/(alpha*alpha-1);
                if (alpha > 0.45) {
                    alpha = 0.45;
                } else if (alpha < -0.45) {
                    alpha = -0.45;
                }

                // TODO: move the constant parts, independent of grid point, out of the loops
                s = (dist(points, &centers[3*iatom0]) - dist(points, &centers[3*iatom1]))/dist(&centers[3*iatom0], &centers[3*iatom1]);
                s = s + alpha*(1-s*s);

                for (int k=1; k <= order; k++) {
                    s = 0.5*s*(3-s*s);
                }
                s = 0.5*(1-s);

                p *= s;
#ifdef DEBUG
                printf("iatom0=%i  iatom1=%i s=%f p=%f\n", iatom0, iatom1, s, p);
#endif
            }

            if (iatom0 == select) nom = p;
            denom += p;
        }
#ifdef DEBUG
        printf("nom=%f  denom=%f\n", nom, denom);
#endif

        *weights = nom/denom;

        // go to next point
        points += 3;
        weights++;
    }
}
