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


//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif
#include <cmath>
#include <stdexcept>

#include "horton/grid/becke.h"


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
        The output, i.e. the Becke weights for the grid points. Note that the
        becke weight is **multiplied** with the original contents of the array!

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

   See Becke's paper for the details:
   A. D. Becke, The Journal of Chemical Physics 88, 2547 (1988)
   URL http://dx.doi.org/10.1063/1.454033.
*/
void becke_helper_atom(int npoint, double* points, double* weights, int natom,
                       double* radii, double* centers, int select, int order)
{
    double nom, denom; // The nominator and the denominator in the weight definition
    double p, s; // Used to build up the value of the switching function

    // precompute the the alpha parameters for each atom pair
    double alphas[(natom*(natom+1))/2];
    long offset = 0;
    for (int iatom0 = 0; iatom0 < natom; iatom0++) {
        for (int iatom1 = 0; iatom1 <= iatom0; iatom1++) {
            // Heteronuclear assignment of the boundary. (Appendix in Becke's paper.)
            double alpha = (radii[iatom0] - radii[iatom1])/(radii[iatom0] + radii[iatom1]); // Eq. (A6)
            alpha = alpha/(alpha*alpha-1); // Eq. (A5)
            // Eq. (A3), except that we use some safe margin (0.45 instead of 0.5)
            // to stay away from a ridiculous imbalance.
            if (alpha > 0.45) {
                alpha = 0.45;
            } else if (alpha < -0.45) {
                alpha = -0.45;
            }
            alphas[offset] = alpha;
            offset += 1;
        }
    }

    // precompute interatomic distances
    double atomic_dists[(natom*(natom+1))/2];
    offset = 0;
    for (int iatom0 = 0; iatom0 < natom; iatom0++) {
        for (int iatom1 = 0; iatom1 <= iatom0; iatom1++) {
            atomic_dists[offset] = dist(&centers[3*iatom0], &centers[3*iatom1]);
            offset += 1;
        }
    }

    // actual computations of Becke weights
    for (int ipoint = npoint-1; ipoint>=0; ipoint--) {
        nom = 0;
        denom = 0;
        for (int iatom0 = 0; iatom0 < natom; iatom0++) {
            p = 1;
            for (int iatom1 = 0; iatom1 < natom; iatom1++) {
                if (iatom0 == iatom1) continue;

                // compute offset for alpha and interatomic distance
                if (iatom0 < iatom1) {
                    offset = (iatom1*(iatom1+1))/2+iatom0;
                } else {
                    offset = (iatom0*(iatom0+1))/2+iatom1;
                }

                // Diatomic switching function
                s = (dist(points, &centers[3*iatom0])
                     -dist(points, &centers[3*iatom1]))
                    /atomic_dists[offset]; // Eq. (11)
                s = s + alphas[offset]*(1 - 2*(iatom0<iatom1))*(1-s*s); // Eq. (A2)

                for (int k=1; k <= order; k++) { // Eq. (19) and (20)
                    s = 0.5*s*(3-s*s);
                }
                s = 0.5*(1-s); // Eq. (18)

                p *= s; // Eq. (13)
#ifdef DEBUG
                printf("iatom0=%i  iatom1=%i s=%f p=%f\n", iatom0, iatom1, s, p);
#endif
            }

            if (iatom0 == select) nom = p;
            denom += p; // Eq. (22)
        }
#ifdef DEBUG
        printf("nom=%f  denom=%f\n", nom, denom);
#endif

        // Weight function at this grid point:
        *weights *= nom/denom; // Eq. (22)

        // go to next point
        points += 3;
        weights++;
    }
}
