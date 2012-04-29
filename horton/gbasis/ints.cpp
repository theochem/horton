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


#include <cmath>
#include <cstring>
#include "boys.h"
#include "cartpure.h"
#include "common.h"
#include "ints.h"
using namespace std;

/*

    The implementation of the one-body operators in this module are based on
    the paper "Gaussian-Expansion Methods for Molecular Integrals", H. Taketa,
    S. Huzinga, K. O-ohata, Journal of the Physical Society of Japan, vol. 21,
    p. 2313, y. 1996. Be aware that there are some misprints in the paper:

    - eq. 2.18: CP_x should be p_x, i.e. change of sign.
    - list of expressions at the bottom of p. 2319:
        case l_1+l_2=2, G_1 = -f_1*p - f_2/(2*gamma)

*/

/*

   Auxiliary functions

*/

void compute_gpt_center(double alpha0, const double* r0, double alpha1, const double* r1, double gamma_inv, double* gpt_center) {
    gpt_center[0] = (alpha0*r0[0] + alpha1*r1[0])*gamma_inv;
    gpt_center[1] = (alpha0*r0[1] + alpha1*r1[1])*gamma_inv;
    gpt_center[2] = (alpha0*r0[2] + alpha1*r1[2])*gamma_inv;
}


const double dist_sq(const double* r0, const double* r1) {
    double result, tmp;
    tmp = r0[0] - r1[0];
    result = tmp*tmp;
    tmp = r0[1] - r1[1];
    result += tmp*tmp;
    tmp = r0[2] - r1[2];
    result += tmp*tmp;
    return result;
}


double gpt_coeff(long k, long n0, long n1, double pa, double pb) {
    long i0, i1;
    double result;

    result = 0;
    i0 = k-n1;
    if (i0<0) i0 = 0;
    i1 = k - i0;
    do {
        result += binom(n0, i0)*binom(n1, i1)*pow(pa, n0-i0)*pow(pb, n1-i1);
        i0++;
        i1--;
    } while ((i1 >= 0)&&(i0 <= n0));
    return result;
}


double gb_overlap_int1d(long n0, long n1, double pa, double pb, double gamma_inv) {
    long k, kmax;
    double result;

    kmax = (n0+n1)/2;
    result = 0.0;
    for (k=0; k<=kmax; k++) {
        result += fac2(2*k-1)*gpt_coeff(2*k, n0, n1, pa, pb)*pow(0.5*gamma_inv,k);
    }
    return sqrt(M_PI*gamma_inv)*result;
}


const double gob_normalization(const double alpha, const long* n) {
    return sqrt(pow(4.0*alpha, n[0]+n[1]+n[2])*pow(2.0*alpha/M_PI, 1.5)
           /(fac2(2*n[0]-1)*fac2(2*n[1]-1)*fac2(2*n[2]-1)));
}



/*

   GB2Integral

*/


GB2Integral::GB2Integral(long max_nbasis): max_nbasis(max_nbasis) {
    nwork = max_nbasis*max_nbasis;
    work_cart = new double[nwork*sizeof(double)];
    work_pure = new double[nwork*sizeof(double)];
}

GB2Integral::~GB2Integral() {
    delete[] work_cart;
    delete[] work_pure;
}

long GB2Integral::get_max_nbasis() {
    return max_nbasis;
}

void GB2Integral::reset(long _shell_type0, long _shell_type1, const double* _r0, const double* _r1) {
    shell_type0 = _shell_type0;
    shell_type1 = _shell_type1;
    r0 = _r0;
    r1 = _r1;
    // We make use of the fact that a floating point zero consists of
    // consecutive zero bytes.
    memset(work_cart, 0, nwork*sizeof(double));
    memset(work_pure, 0, nwork*sizeof(double));
}

void GB2Integral::cart_to_pure() {
    /*
       The initial results are always stored in work_cart. The projection
       routine always outputs its result in work_pure. Once that is done,
       the pointers to both blocks are swapped such that the final result is
       always back in work_cart.
    */

    // Project along index 0 (rows)
    if (shell_type0 < -1) {
        cart_to_pure_low(work_cart, work_pure, -shell_type0,
            1, // stride
            max_nbasis, // spacing
            get_shell_nbasis(abs(shell_type1)) // count
        );
        swap_work();
    }

    // Project along index 1 (rows)
    if (shell_type1 < -1) {
        cart_to_pure_low(work_cart, work_pure, -shell_type1,
            max_nbasis, // stride
            1, // spacing
            get_shell_nbasis(shell_type0) // count
        );
        swap_work();
    }
}

void GB2Integral::swap_work() {
    double* tmp;
    tmp = work_cart;
    work_cart = work_pure;
    work_pure = tmp;
}


/*

   GB2OverlapIntegral

*/


void GB2OverlapIntegral::add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1) {
    double pre, gamma_inv;
    double gpt_center[3];

    gamma_inv = 1.0/(alpha0 + alpha1);
    pre = coeff*exp(-alpha0*alpha1*gamma_inv*dist_sq(r0, r1));
    compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
    i2p.reset(abs(shell_type0), abs(shell_type1), max_nbasis);
    do {
        work_cart[i2p.offset] += pre*(
            gb_overlap_int1d(i2p.n0[0], i2p.n1[0], gpt_center[0] - r0[0], gpt_center[0] - r1[0], gamma_inv)*
            gb_overlap_int1d(i2p.n0[1], i2p.n1[1], gpt_center[1] - r0[1], gpt_center[1] - r1[1], gamma_inv)*
            gb_overlap_int1d(i2p.n0[2], i2p.n1[2], gpt_center[2] - r0[2], gpt_center[2] - r1[2], gamma_inv)*
            scales0[i2p.ibasis0]*scales1[i2p.ibasis1]
        );
    } while (i2p.inc());
}


/*

   GB2KineticIntegral

*/


double kinetic_helper(double alpha0, double alpha1, long n0, long n1, double pa, double pb, double gamma_inv) {
    double poly = 0;

    if (n0 > 0) {
        if (n1 > 0) {
            // <-1|-1>
            poly += 0.5*n0*n1*gb_overlap_int1d(n0-1, n1-1, pa, pb, gamma_inv);
        }
        // <-1|+1>
        poly -= alpha1*n0*gb_overlap_int1d(n0-1, n1+1, pa, pb, gamma_inv);
    }
    if (n1 > 0) {
        // <+1|-1>
        poly -= alpha0*n1*gb_overlap_int1d(n0+1, n1-1, pa, pb, gamma_inv);
    }
    // <+1|+1>
    poly += 2.0*alpha0*alpha1*gb_overlap_int1d(n0+1, n1+1, pa, pb, gamma_inv);

    return poly;
}

void GB2KineticIntegral::add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1) {
    double pre, gamma_inv, poly, fx0, fy0, fz0;
    double gpt_center[3], pa[3], pb[3];

    gamma_inv = 1.0/(alpha0 + alpha1);
    pre = coeff*exp(-alpha0*alpha1*gamma_inv*dist_sq(r0, r1));
    compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
    pa[0] = gpt_center[0] - r0[0];
    pa[1] = gpt_center[1] - r0[1];
    pa[2] = gpt_center[2] - r0[2];
    pb[0] = gpt_center[0] - r1[0];
    pb[1] = gpt_center[1] - r1[1];
    pb[2] = gpt_center[2] - r1[2];

    i2p.reset(abs(shell_type0), abs(shell_type1), max_nbasis);
    do {
        fx0 = gb_overlap_int1d(i2p.n0[0], i2p.n1[0], pa[0], pb[0], gamma_inv);
        fy0 = gb_overlap_int1d(i2p.n0[1], i2p.n1[1], pa[1], pb[1], gamma_inv);
        fz0 = gb_overlap_int1d(i2p.n0[2], i2p.n1[2], pa[2], pb[2], gamma_inv);
        poly = fy0*fz0*kinetic_helper(alpha0, alpha1, i2p.n0[0], i2p.n1[0], pa[0], pb[0], gamma_inv) +
               fz0*fx0*kinetic_helper(alpha0, alpha1, i2p.n0[1], i2p.n1[1], pa[1], pb[1], gamma_inv) +
               fx0*fy0*kinetic_helper(alpha0, alpha1, i2p.n0[2], i2p.n1[2], pa[2], pb[2], gamma_inv);
        work_cart[i2p.offset] += pre*scales0[i2p.ibasis0]*scales1[i2p.ibasis1]*poly;
    } while (i2p.inc());
}


/*

   GB2NuclearAttractionIntegral

*/


GB2NuclearAttractionIntegral::GB2NuclearAttractionIntegral(long max_nbasis, double* charges, double* centers, long ncharge) :
            GB2Integral(max_nbasis), charges(charges), centers(centers), ncharge(ncharge) {
    // inverse of nbasis = ((max_shell_type+1)*(max_shell_type+2))/2
    // TODO. Replace max_nbasis by max_shell_type in all Integral classes.
    max_shell_type = (long)((sqrt(1+8*max_nbasis) - 3)/2);
    work_g0 = new double[2*max_shell_type];
    work_g1 = new double[2*max_shell_type];
    work_g2 = new double[2*max_shell_type];
    work_boys = new double[2*max_shell_type];
}


GB2NuclearAttractionIntegral::~GB2NuclearAttractionIntegral() {
    delete[] work_g0;
    delete[] work_g1;
    delete[] work_g2;
    delete[] work_boys;
}


void nuclear_attraction_helper(double* work_g, long n0, long n1, double pa, double pb, double pc, double gamma_inv) {
    for (long index=n0+n1; index>=0; index--) {
        double tmp=0;
        for (long i=n0+n1; i>=index; i--) {
            long rmin = (i+1)/2-index;
            if (rmin<=0) rmin = 0;
            for (long r=(i-index)/2; r>=rmin; r--) {
                long u=i-2*r-index;
                tmp += (
                    (1-2*(i%2))* // (-1)**i
                    gpt_coeff(i, n0, n1, pa, pb)*
                    (1-2*(u%2))* // (-1)**u
                    fac(i)*
                    pow(pc, i-2*r-2*u)*
                    pow(0.25*gamma_inv, r+u)/
                    (fac(r)*fac(u)*fac(i-2*r-2*u))
                );
            }
        }
        work_g[index] = tmp;
    }
}


void GB2NuclearAttractionIntegral::add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1) {
    double pre, gamma, gamma_inv, arg;
    double gpt_center[3], pa[3], pb[3], pc[3];

    gamma = alpha0 + alpha1;
    gamma_inv = 1.0/gamma;
    pre = 2*M_PI*gamma_inv*coeff*exp(-alpha0*alpha1*gamma_inv*dist_sq(r0, r1));
    compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
    pa[0] = gpt_center[0] - r0[0];
    pa[1] = gpt_center[1] - r0[1];
    pa[2] = gpt_center[2] - r0[2];
    pb[0] = gpt_center[0] - r1[0];
    pb[1] = gpt_center[1] - r1[1];
    pb[2] = gpt_center[2] - r1[2];

    for (long icharge=0; icharge < ncharge; icharge++) {
        // thrid center for the current charge
        pc[0] = gpt_center[0] - centers[icharge*3  ];
        pc[1] = gpt_center[1] - centers[icharge*3+1];
        pc[2] = gpt_center[2] - centers[icharge*3+2];

        // Fill the work array with the Boys function values
        arg = gamma*(pc[0]*pc[0] + pc[1]*pc[1] + pc[2]*pc[2]);
        for (long nu=abs(shell_type0)+abs(shell_type1); nu>=0; nu--) {
            work_boys[nu] = boys_function(nu, arg);
        }

        // Iterate over all combinations of Cartesian exponents
        i2p.reset(abs(shell_type0), abs(shell_type1), max_nbasis);
        do {
            // Fill the work arrays with the polynomials
            nuclear_attraction_helper(work_g0, i2p.n0[0], i2p.n1[0], pa[0], pb[0], pc[0], gamma_inv);
            nuclear_attraction_helper(work_g1, i2p.n0[1], i2p.n1[1], pa[1], pb[1], pc[1], gamma_inv);
            nuclear_attraction_helper(work_g2, i2p.n0[2], i2p.n1[2], pa[2], pb[2], pc[2], gamma_inv);

            // Take the product
            arg = 0;
            for (long i0=i2p.n0[0]+i2p.n1[0]; i0>=0; i0--)
                for (long i1=i2p.n0[1]+i2p.n1[1]; i1>=0; i1--)
                    for (long i2=i2p.n0[2]+i2p.n1[2]; i2>=0; i2--)
                        arg += work_g0[i0]*work_g1[i1]*work_g2[i2]*work_boys[i0+i1+i2];

            // Finally add to the work array
            work_cart[i2p.offset] += pre*scales0[i2p.ibasis0]*scales1[i2p.ibasis1]*arg*charges[icharge];
        } while (i2p.inc());
    }
}
