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
#include "horton/cell.h"


Cell::Cell(double* _rvecs, int _nvec) {
    // check if nvec is sensible
    if ((_nvec < 0) || (_nvec > 3)) {
        throw std::domain_error("The number of cell vectors must be 0, 1, 2 or 3.");
    }

    // copy the given _rvecs and _nvec:
    nvec = _nvec;
    for (int i=nvec*3-1; i>=0; i--) {
        rvecs[i] = _rvecs[i];
    }

    // compute the volume
    switch(nvec) {
        case 0:
            volume = 0.0;
            break;
        case 1:
            volume = sqrt(
                rvecs[0]*rvecs[0]+rvecs[1]*rvecs[1]+rvecs[2]*rvecs[2]
            );
            break;
        case 2:
            double tmp;
            tmp = rvecs[0]*rvecs[3]+rvecs[1]*rvecs[4]+rvecs[2]*rvecs[5];
            tmp = (rvecs[0]*rvecs[0]+rvecs[1]*rvecs[1]+rvecs[2]*rvecs[2])*
                  (rvecs[3]*rvecs[3]+rvecs[4]*rvecs[4]+rvecs[5]*rvecs[5]) - tmp*tmp;
            if (tmp > 0) {
                volume = sqrt(tmp);
            } else {
                volume = 0.0;
            }
            break;
        case 3:
            volume = fabs(
                rvecs[0]*(rvecs[4]*rvecs[8]-rvecs[5]*rvecs[7])+
                rvecs[1]*(rvecs[5]*rvecs[6]-rvecs[3]*rvecs[8])+
                rvecs[2]*(rvecs[3]*rvecs[7]-rvecs[4]*rvecs[6])
            );
            break;
    }

    // If the volume is zero and nvec > 0, raise an error. In this case, the
    // reciprocal cell vectors can not be computed.
    if ((volume == 0.0) && (nvec > 0)) {
        throw std::domain_error("The cell vectors are degenerate");
    }

    // complete the list of rvecs in case nvec < 3
    switch(nvec) {
        case 0:
            // Just put in the identity matrix.
            rvecs[0] = 1.0;
            rvecs[1] = 0.0;
            rvecs[2] = 0.0;
            rvecs[3] = 0.0;
            rvecs[4] = 1.0;
            rvecs[5] = 0.0;
            rvecs[6] = 0.0;
            rvecs[7] = 0.0;
            rvecs[8] = 1.0;
            break;
        case 1: {
            // Add two rvecs that are orthogonal to the given rvec, orthogonal
            // to each other and normalized. The three vectors will be
            // right-handed.
            // 1) find the component of the given vector with the smallest
            // absolute value
            int ismall = 0;
            if (fabs(rvecs[1]) < fabs(rvecs[0])) {
                ismall = 1;
                if (fabs(rvecs[2]) < fabs(rvecs[1])) {
                    ismall = 2;
                }
            } else if (fabs(rvecs[2]) < fabs(rvecs[0])) {
                ismall = 2;
            }
            // 2) store a temporary vector in position 3
            rvecs[6] = 0.0;
            rvecs[7] = 0.0;
            rvecs[8] = 0.0;
            rvecs[ismall+6] = 1.0;
            // 3) compute the cross product of vector 1 and 3
            rvecs[3] = rvecs[1]*rvecs[8] - rvecs[2]*rvecs[7];
            rvecs[4] = rvecs[2]*rvecs[6] - rvecs[0]*rvecs[8];
            rvecs[5] = rvecs[0]*rvecs[7] - rvecs[1]*rvecs[6];
            // 4) normalize
            double norm = sqrt(rvecs[3]*rvecs[3] + rvecs[4]*rvecs[4] + rvecs[5]*rvecs[5]);
            rvecs[3] /= norm;
            rvecs[4] /= norm;
            rvecs[5] /= norm;
            // the rest is done in case 2, so no break here!
        }
        case 2:
            // Add one rvec that is normalized and orthogonal to the two given
            // rvecs. The three vectors will be right-handed.
            // 1) compute the cross product of vector 1 and 2
            rvecs[6] = rvecs[1]*rvecs[5] - rvecs[2]*rvecs[4];
            rvecs[7] = rvecs[2]*rvecs[3] - rvecs[0]*rvecs[5];
            rvecs[8] = rvecs[0]*rvecs[4] - rvecs[1]*rvecs[3];
            // 2) normalize
            double norm = sqrt(rvecs[6]*rvecs[6] + rvecs[7]*rvecs[7] + rvecs[8]*rvecs[8]);
            rvecs[6] /= norm;
            rvecs[7] /= norm;
            rvecs[8] /= norm;
    }

    // Now we assume that rvecs contains a set of three well-behaved
    // non-degenerate vectors. Cramer's rule is used to compute the reciprocal
    // space vectors. This is fairly ugly in terms of numerical stability but
    // it keeps things simple.
    gvecs[0] = rvecs[4]*rvecs[8] - rvecs[5]*rvecs[7];
    gvecs[1] = rvecs[5]*rvecs[6] - rvecs[3]*rvecs[8];
    gvecs[2] = rvecs[3]*rvecs[7] - rvecs[4]*rvecs[6];
    gvecs[3] = rvecs[7]*rvecs[2] - rvecs[8]*rvecs[1];
    gvecs[4] = rvecs[8]*rvecs[0] - rvecs[6]*rvecs[2];
    gvecs[5] = rvecs[6]*rvecs[1] - rvecs[7]*rvecs[0];
    gvecs[6] = rvecs[1]*rvecs[5] - rvecs[2]*rvecs[4];
    gvecs[7] = rvecs[2]*rvecs[3] - rvecs[0]*rvecs[5];
    gvecs[8] = rvecs[0]*rvecs[4] - rvecs[1]*rvecs[3];
    // determinant
    double det = gvecs[0]*rvecs[0] + gvecs[1]*rvecs[1] + gvecs[2]*rvecs[2];
    // inverse
    gvecs[0] /= det;
    gvecs[1] /= det;
    gvecs[2] /= det;
    gvecs[3] /= det;
    gvecs[4] /= det;
    gvecs[5] /= det;
    gvecs[6] /= det;
    gvecs[7] /= det;
    gvecs[8] /= det;

    // compute the spacings and the lengths of the cell vectors
    for (int i=2; i>=0; i--) {
        rlengths[i] = sqrt(rvecs[3*i]*rvecs[3*i] + rvecs[3*i+1]*rvecs[3*i+1] + rvecs[3*i+2]*rvecs[3*i+2]);
        glengths[i] = sqrt(gvecs[3*i]*gvecs[3*i] + gvecs[3*i+1]*gvecs[3*i+1] + gvecs[3*i+2]*gvecs[3*i+2]);
        rspacings[i] = 1.0/glengths[i];
        gspacings[i] = 1.0/rlengths[i];
    }
}


void Cell::mic(double* delta) const {
    // Applies the Minimum Image Convention. Well, sort of. It does not always work like this.
    // This function contains an unrolled loop for speed.
    double x;
    if (nvec == 0) return;
    x = ceil(gvecs[0]*delta[0] + gvecs[1]*delta[1] + gvecs[2]*delta[2] - 0.5);
    delta[0] -= x*rvecs[0];
    delta[1] -= x*rvecs[1];
    delta[2] -= x*rvecs[2];
    if (nvec == 1) return;
    x = ceil(gvecs[3]*delta[0] + gvecs[4]*delta[1] + gvecs[5]*delta[2] - 0.5);
    delta[0] -= x*rvecs[3];
    delta[1] -= x*rvecs[4];
    delta[2] -= x*rvecs[5];
    if (nvec == 2) return;
    x = ceil(gvecs[6]*delta[0] + gvecs[7]*delta[1] + gvecs[8]*delta[2] - 0.5);
    delta[0] -= x*rvecs[6];
    delta[1] -= x*rvecs[7];
    delta[2] -= x*rvecs[8];
}


void Cell::to_frac(double* cart, double* frac) const {
    // Transfroms to fractional coordinates
    frac[0] = gvecs[0]*cart[0] + gvecs[1]*cart[1] + gvecs[2]*cart[2];
    frac[1] = gvecs[3]*cart[0] + gvecs[4]*cart[1] + gvecs[5]*cart[2];
    frac[2] = gvecs[6]*cart[0] + gvecs[7]*cart[1] + gvecs[8]*cart[2];
}

void Cell::to_cart(double* frac, double* cart) const {
    // Transfroms to Cartesian coordinates
    cart[0] = rvecs[0]*frac[0] + rvecs[3]*frac[1] + rvecs[6]*frac[2];
    cart[1] = rvecs[1]*frac[0] + rvecs[4]*frac[1] + rvecs[7]*frac[2];
    cart[2] = rvecs[2]*frac[0] + rvecs[5]*frac[1] + rvecs[8]*frac[2];
}


void Cell::g_lincomb(double* coeffs, double* gvec) const {
    // Make a linear combination of reciprocal cell vectors
    gvec[0] = gvecs[0]*coeffs[0] + gvecs[3]*coeffs[1] + gvecs[6]*coeffs[2];
    gvec[1] = gvecs[1]*coeffs[0] + gvecs[4]*coeffs[1] + gvecs[7]*coeffs[2];
    gvec[2] = gvecs[2]*coeffs[0] + gvecs[5]*coeffs[1] + gvecs[8]*coeffs[2];
}


void Cell::dot_rvecs(double* cart, double* dots) const {
    // Take dot product with real cell vectors
    dots[0] = rvecs[0]*cart[0] + rvecs[1]*cart[1] + rvecs[2]*cart[2];
    dots[1] = rvecs[3]*cart[0] + rvecs[4]*cart[1] + rvecs[5]*cart[2];
    dots[2] = rvecs[6]*cart[0] + rvecs[7]*cart[1] + rvecs[8]*cart[2];
}


void Cell::add_rvec(double* delta, long* coeffs) const {
    // Simply adds an linear combination of real cell vectors to delta.
    // This function contains an unrolled loop for speed.
    if (nvec == 0) return;
    delta[0] += coeffs[0]*rvecs[0];
    delta[1] += coeffs[0]*rvecs[1];
    delta[2] += coeffs[0]*rvecs[2];
    if (nvec == 1) return;
    delta[0] += coeffs[1]*rvecs[3];
    delta[1] += coeffs[1]*rvecs[4];
    delta[2] += coeffs[1]*rvecs[5];
    if (nvec == 2) return;
    delta[0] += coeffs[2]*rvecs[6];
    delta[1] += coeffs[2]*rvecs[7];
    delta[2] += coeffs[2]*rvecs[8];
}


double Cell::get_rspacing(int i) const {
    if ((i < 0) || (i > 3)) {
        throw std::domain_error("Index must be 0, 1 or 2.");
    }
    return rspacings[i];
}


double Cell::get_gspacing(int i) const {
    if ((i < 0) || (i > 3)) {
        throw std::domain_error("Index must be 0, 1 or 2.");
    }
    return gspacings[i];
}

double Cell::get_rlength(int i) const {
    if ((i < 0) || (i > 3)) {
        throw std::domain_error("Index must be 0, 1 or 2.");
    }
    return rlengths[i];
}

double Cell::get_glength(int i) const {
    if ((i < 0) || (i > 3)) {
        throw std::domain_error("Index must be 0, 1 or 2.");
    }
    return glengths[i];
}


void Cell::copy_rvecs(double* _rvecs) const {
    for (int i=nvec*3-1; i>=0; i--) _rvecs[i] = rvecs[i];
}


void Cell::copy_gvecs(double* _gvecs) const {
    for (int i=nvec*3-1; i>=0; i--) _gvecs[i] = gvecs[i];
}


void Cell::copy_rlengths(double* _rlengths) const {
    for (int i=nvec-1; i>=0; i--) _rlengths[i] = rlengths[i];
}


void Cell::copy_glengths(double* _glengths) const {
    for (int i=nvec-1; i>=0; i--) _glengths[i] = glengths[i];
}


void Cell::copy_rspacings(double* _rspacings) const {
    for (int i=nvec-1; i>=0; i--) _rspacings[i] = rspacings[i];
}


void Cell::copy_gspacings(double* _gspacings) const {
    for (int i=nvec-1; i>=0; i--) _gspacings[i] = gspacings[i];
}


void Cell::set_ranges_rcut(double* center, double rcut,  long* ranges_begin,
    long* ranges_end) const {
    double frac[3];
    to_frac(center, frac);
    for (int i=nvec-1; i>=0; i--) {
        double step = rcut/rspacings[i];
        ranges_begin[i] = ceil(-frac[i]-step);
        ranges_end[i] = ceil(-frac[i]+step);
    }
}

long smart_wrap(long i, long shape, long pbc) {
    if ((i < 0) || (i >= shape)) {
        if (pbc) {
            long j = i%shape;
            if (j < 0) j += shape; // just to make sure that this works on all compilers.
            return j;
        } else {
            return -1;
        }
    }
    return i;
}
