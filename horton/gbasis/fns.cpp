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

//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include "horton/moments.h"
#include "horton/gbasis/boys.h"
#include "horton/gbasis/cartpure.h"
#include "horton/gbasis/fns.h"
using namespace std;


/*
    GB1GridFn
*/

GB1GridFn::GB1GridFn(long max_shell_type, long dim_work, long dim_output):
    GBCalculator(max_shell_type), dim_work(dim_work), dim_output(dim_output)
{
    nwork = max_nbasis*dim_work;
    work_cart = new double[nwork];
    work_pure = new double[nwork];
}

void GB1GridFn::reset(long _shell_type0, const double* _r0, const double* _point) {
    if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type)) {
      throw domain_error("shell_type0 out of range.");
    }
    shell_type0 = _shell_type0;
    r0 = _r0;
    point = _point;
    // We make use of the fact that a floating point zero consists of
    // consecutive zero bytes.
    memset(work_cart, 0, nwork*sizeof(double));
    memset(work_pure, 0, nwork*sizeof(double));
}

void GB1GridFn::cart_to_pure() {
    /*
       The initial results are always stored in work_cart. The projection
       routine always outputs its result in work_pure. Once that is done,
       the pointers to both blocks are swapped such that the final result is
       always back in work_cart.
    */

    if (shell_type0 < -1) {
        cart_to_pure_low(work_cart, work_pure, -shell_type0,
            1, // anterior
            dim_work // posterior
        );
        swap_work();
    }
}


/*
    GB1ExpGridOrbitalFn
*/


void GB1ExpGridOrbitalFn::add(double coeff, double alpha0, const double* scales0) {
    double pre, poly;
    pre = coeff*exp(-alpha0*dist_sq(r0, point));
    i1p.reset(abs(shell_type0));
    do {
#ifdef DEBUG
        printf("n=[%i,%i,%i]\n", i1p.n0[0], i1p.n0[1], i1p.n0[2]);
#endif
        // For now, simple and inefficient evaluation of polynomial.
        poly = 1.0;
        for (long j=0; j<3; j++) {
            double tmp = point[j] - r0[j];
            for (long i=0; i<i1p.n0[j]; i++) poly *= tmp;
        }
        work_cart[i1p.ibasis0] += pre*scales0[i1p.ibasis0]*poly;
#ifdef DEBUG
        printf("poly=%f ibasis0=%i scale=%f work[]=%f\n", poly, i1p.ibasis0, scales0[i1p.ibasis0], work_cart[i1p.ibasis0]);
#endif
    } while (i1p.inc());
#ifdef DEBUG
    printf("\n");
#endif
}

void GB1ExpGridOrbitalFn::compute_point_from_exp(double* work_basis, double* coeffs, long nbasis, double* output) {
    for (long i=0; i < norb; i++) {
        long iorb = iorbs[i];
        for (long ibasis=0; ibasis < nbasis; ibasis++) {
            output[i] += coeffs[ibasis*nfn + iorb]*work_basis[ibasis];
        }
    }
}


/*
    GB1DMGridDensityFn
*/

void GB1DMGridDensityFn::reset(long _shell_type0, const double* _r0, const double* _point) {
    GB1GridFn::reset(_shell_type0, _r0, _point);
    if (shell_type0!=0) {
        poly_work[0] = point[0] - r0[0];
        poly_work[1] = point[1] - r0[1];
        poly_work[2] = point[2] - r0[2];
        offset = fill_cartesian_polynomials(poly_work, abs(shell_type0));
    }
}

void GB1DMGridDensityFn::add(double coeff, double alpha0, const double* scales0) {
    double pre;
    pre = coeff*exp(-alpha0*dist_sq(r0, point));
    if (shell_type0==0) {
        work_cart[0] += pre*scales0[0];
    } else {
        for (long ibasis0=get_shell_nbasis(abs(shell_type0))-1; ibasis0>=0; ibasis0--) {
            work_cart[ibasis0] += pre*scales0[ibasis0]*poly_work[ibasis0+offset];
        }
    }
}

void GB1DMGridDensityFn::compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output, double epsilon, double* dmmaxrow) {
    double absmax_basis = 0.0;
    if (epsilon > 0) {
        // compute the maximum basis function
        for (long ibasis=0; ibasis<nbasis; ibasis++) {
            double tmp = fabs(work_basis[ibasis]);
            if (tmp > absmax_basis) absmax_basis = tmp;
        }
        // upper estimate of the density
        double rho_upper = 0.0;
        for (long ibasis=0; ibasis<nbasis; ibasis++) {
            rho_upper += fabs(work_basis[ibasis])*dmmaxrow[ibasis];
        }
        rho_upper *= nbasis*absmax_basis;

        // if the upper bound is too low, do not compute density.
        if (rho_upper < epsilon) return;

        // modify epsilon to avoid recomputation
        epsilon /= absmax_basis*nbasis*nbasis;
    }

    // Loop over all basis functions and add significant contributions
    double rho = 0.0;
    for (long ibasis0=0; ibasis0<nbasis; ibasis0++) {
        // if the contribution of this loop is smaller than epsilon/nbasis, skipt it.
        if (epsilon>0) {
            if (fabs(work_basis[ibasis0])*dmmaxrow[ibasis0] < epsilon)
                continue;
        }
        double tmp = 0;
        for (long ibasis1=ibasis0-1; ibasis1>=0; ibasis1--) {
            tmp += work_basis[ibasis1]*dm[ibasis0*nbasis+ibasis1];
        }
        rho += (2*tmp+dm[ibasis0*(nbasis+1)]*work_basis[ibasis0])*work_basis[ibasis0];
    }
    *output += rho;
}

void GB1DMGridDensityFn::compute_fock_from_pot(double* pot, double* work_basis, long nbasis, double* output) {
    for (long ibasis0=0; ibasis0<nbasis; ibasis0++) {
        double tmp1 = (*pot)*work_basis[ibasis0];
        for (long ibasis1=0; ibasis1<=ibasis0; ibasis1++) {
            double tmp2 = tmp1*work_basis[ibasis1];
            output[ibasis0*nbasis+ibasis1] += tmp2;
            if (ibasis0 != ibasis1) output[ibasis1*nbasis+ibasis0] += tmp2;
        }
    }
}



/*
    GB1DMGridGradientFn
*/

static void poly_helper(double x, long n, double* poly, double* poly1) {
    for (long i=n; i>0; i--) {
        *poly *= x;
        if (i==2) {
            *poly1 = *poly;
        }
    }
}


void GB1DMGridGradientFn::add(double coeff, double alpha0, const double* scales0) {
    double x = point[0] - r0[0];
    double y = point[1] - r0[1];
    double z = point[2] - r0[2];
    double pre = coeff*exp(-alpha0*(x*x+y*y+z*z));
    i1p.reset(abs(shell_type0));
    do {
        double pre0 = pre*scales0[i1p.ibasis0];

        // For now, simple and inefficient evaluation of polynomial.
        double poly_x = 1.0;
        double poly_1x = 1.0;
        poly_helper(x, i1p.n0[0], &poly_x, &poly_1x);
        double poly_y = 1.0;
        double poly_1y = 1.0;
        poly_helper(y, i1p.n0[1], &poly_y, &poly_1y);
        double poly_z = 1.0;
        double poly_1z = 1.0;
        poly_helper(z, i1p.n0[2], &poly_z, &poly_1z);

        double tmp0 = pre0*poly_x*poly_y*poly_z;
        double tmp1;
        // Basis function value
        work_cart[i1p.ibasis0*4] += tmp0;
        // Basis function derivative towards x
        tmp0 *= -2.0*alpha0;
        tmp1 = x*tmp0;
        if (i1p.n0[0] > 0) tmp1 += i1p.n0[0]*pre0*poly_1x*poly_y*poly_z;
        work_cart[i1p.ibasis0*4+1] += tmp1;
        // Basis function derivative towards y
        tmp1 = y*tmp0;
        if (i1p.n0[1] > 0) tmp1 += i1p.n0[1]*pre0*poly_x*poly_1y*poly_z;
        work_cart[i1p.ibasis0*4+2] += tmp1;
        // Basis function derivative towards z
        tmp1 = z*tmp0;
        if (i1p.n0[2] > 0) tmp1 += i1p.n0[2]*pre0*poly_x*poly_y*poly_1z;
        work_cart[i1p.ibasis0*4+3] += tmp1;
    } while (i1p.inc());
}

void GB1DMGridGradientFn::compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output, double epsilon, double* dmmaxrow) {
    double rho_x = 0, rho_y = 0, rho_z = 0;
    for (long ibasis0=0; ibasis0<nbasis; ibasis0++) {
        double row = 0;
        for (long ibasis1=0; ibasis1<nbasis; ibasis1++) {
            row += work_basis[ibasis1*4]*dm[ibasis0*nbasis+ibasis1];
        }
        rho_x += row*work_basis[ibasis0*4+1];
        rho_y += row*work_basis[ibasis0*4+2];
        rho_z += row*work_basis[ibasis0*4+3];
    }
    output[0] += 2*rho_x;
    output[1] += 2*rho_y;
    output[2] += 2*rho_z;
}

void GB1DMGridGradientFn::compute_fock_from_pot(double* pot, double* work_basis, long nbasis, double* output) {
    for (long ibasis0=0; ibasis0<nbasis; ibasis0++) {
        double tmp0 = work_basis[ibasis0*4];
        double tmp1 = pot[0]*work_basis[ibasis0*4+1] +
                      pot[1]*work_basis[ibasis0*4+2] +
                      pot[2]*work_basis[ibasis0*4+3];
        for (long ibasis1=0; ibasis1<=ibasis0; ibasis1++) {
            double result = tmp0*(pot[0]*work_basis[ibasis1*4+1] +
                                  pot[1]*work_basis[ibasis1*4+2] +
                                  pot[2]*work_basis[ibasis1*4+3]) +
                            tmp1*work_basis[ibasis1*4];
            output[ibasis1*nbasis+ibasis0] += result;
            if (ibasis1!=ibasis0) {
                // Enforce symmetry
                output[ibasis0*nbasis+ibasis1] += result;
            }
        }
    }
}


/*
    GB1DMGridKineticFn
*/


void GB1DMGridKineticFn::add(double coeff, double alpha0, const double* scales0) {
    double x = point[0] - r0[0];
    double y = point[1] - r0[1];
    double z = point[2] - r0[2];
    double pre = coeff*exp(-alpha0*(x*x+y*y+z*z));
    i1p.reset(abs(shell_type0));
    do {
        double pre0 = pre*scales0[i1p.ibasis0];

        // For now, simple and inefficient evaluation of polynomial.
        double poly_x = 1.0;
        double poly_1x = 1.0;
        poly_helper(x, i1p.n0[0], &poly_x, &poly_1x);
        double poly_y = 1.0;
        double poly_1y = 1.0;
        poly_helper(y, i1p.n0[1], &poly_y, &poly_1y);
        double poly_z = 1.0;
        double poly_1z = 1.0;
        poly_helper(z, i1p.n0[2], &poly_z, &poly_1z);

        double tmp0 = pre0*poly_x*poly_y*poly_z;
        double tmp1;
        // Basis function derivative towards x
        tmp0 *= -2.0*alpha0;
        tmp1 = x*tmp0;
        if (i1p.n0[0] > 0) tmp1 += i1p.n0[0]*pre0*poly_1x*poly_y*poly_z;
        work_cart[i1p.ibasis0*3  ] += tmp1;
        // Basis function derivative towards y
        tmp1 = y*tmp0;
        if (i1p.n0[1] > 0) tmp1 += i1p.n0[1]*pre0*poly_x*poly_1y*poly_z;
        work_cart[i1p.ibasis0*3+1] += tmp1;
        // Basis function derivative towards z
        tmp1 = z*tmp0;
        if (i1p.n0[2] > 0) tmp1 += i1p.n0[2]*pre0*poly_x*poly_y*poly_1z;
        work_cart[i1p.ibasis0*3+2] += tmp1;
    } while (i1p.inc());
}

void GB1DMGridKineticFn::compute_point_from_dm(double* work_basis, double* dm, long nbasis, double* output, double epsilon, double* dmmaxrow) {
    double tau = 0.0;
    for (long ibasis0=0; ibasis0<nbasis; ibasis0++) {
        double tmp_x = 0;
        double tmp_y = 0;
        double tmp_z = 0;
        for (long ibasis1=0; ibasis1<ibasis0; ibasis1++) {
            tmp_x += work_basis[ibasis1*3  ]*dm[ibasis0*nbasis+ibasis1];
            tmp_y += work_basis[ibasis1*3+1]*dm[ibasis0*nbasis+ibasis1];
            tmp_z += work_basis[ibasis1*3+2]*dm[ibasis0*nbasis+ibasis1];
        }
        tmp_x += 0.5*work_basis[ibasis0*3  ]*dm[ibasis0*nbasis+ibasis0];
        tmp_y += 0.5*work_basis[ibasis0*3+1]*dm[ibasis0*nbasis+ibasis0];
        tmp_z += 0.5*work_basis[ibasis0*3+2]*dm[ibasis0*nbasis+ibasis0];
        tau += tmp_x*work_basis[ibasis0*3  ] +
               tmp_y*work_basis[ibasis0*3+1] +
               tmp_z*work_basis[ibasis0*3+2];
    }
    *output += tau;
}

void GB1DMGridKineticFn::compute_fock_from_pot(double* pot, double* work_basis, long nbasis, double* output) {
    for (long ibasis0=0; ibasis0<nbasis; ibasis0++) {
        double tmp_x = 0.5*(*pot)*work_basis[ibasis0*3  ];
        double tmp_y = 0.5*(*pot)*work_basis[ibasis0*3+1];
        double tmp_z = 0.5*(*pot)*work_basis[ibasis0*3+2];
        for (long ibasis1=0; ibasis1<=ibasis0; ibasis1++) {
            double result = tmp_x*work_basis[ibasis1*3  ] +
                            tmp_y*work_basis[ibasis1*3+1] +
                            tmp_z*work_basis[ibasis1*3+2];
            output[ibasis1*nbasis+ibasis0] += result;
            if (ibasis0 != ibasis1) {
                output[ibasis0*nbasis+ibasis1] += result;
            }
        }
    }
}


/*
    GB2DMGridFn
*/

GB2DMGridFn::GB2DMGridFn(long max_shell_type):
    GBCalculator(max_shell_type)
{
    nwork = max_nbasis*max_nbasis;
    work_cart = new double[nwork];
    work_pure = new double[nwork];
}

void GB2DMGridFn::reset(long _shell_type0, long _shell_type1, const double* _r0, const double* _r1, const double* _point) {
    if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type)) {
      throw domain_error("shell_type0 out of range.");
    }
    if ((_shell_type1 < -max_shell_type) || (_shell_type1 > max_shell_type)) {
      throw domain_error("shell_type0 out of range.");
    }
    shell_type0 = _shell_type0;
    shell_type1 = _shell_type1;
    r0 = _r0;
    r1 = _r1;
    point = _point;
    // We make use of the fact that a floating point zero consists of
    // consecutive zero bytes.
    memset(work_cart, 0, nwork*sizeof(double));
    memset(work_pure, 0, nwork*sizeof(double));
}

void GB2DMGridFn::cart_to_pure() {
    /*
       The initial results are always stored in work_cart. The projection
       routine always outputs its result in work_pure. Once that is done,
       the pointers to both blocks are swapped such that the final result is
       always back in work_cart.
    */

    // For now, this is just a copy from GB2Integral. It must be changed when electrical fields are implemented.

    // Project along index 0 (rows)
    if (shell_type0 < -1) {
        cart_to_pure_low(work_cart, work_pure, -shell_type0,
            1, // anterior
            get_shell_nbasis(abs(shell_type1)) // posterior
        );
        swap_work();
    }

    // Project along index 1 (cols)
    if (shell_type1 < -1) {
        cart_to_pure_low(work_cart, work_pure, -shell_type1,
            get_shell_nbasis(shell_type0), // anterior
            1 // posterior
        );
        swap_work();
    }
}


/*
    GB2DMGridHartreeFn
*/

GB2DMGridHartreeFn::GB2DMGridHartreeFn(long max_shell_type): GB2DMGridFn(max_shell_type) {
    work_g0 = new double[2*max_shell_type+1];
    work_g1 = new double[2*max_shell_type+1];
    work_g2 = new double[2*max_shell_type+1];
    work_boys = new double[2*max_shell_type+1];
}


GB2DMGridHartreeFn::~GB2DMGridHartreeFn() {
    delete[] work_g0;
    delete[] work_g1;
    delete[] work_g2;
    delete[] work_boys;
}


void GB2DMGridHartreeFn::add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1) {
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

    // thrid center for the current charge
    pc[0] = gpt_center[0] - point[0];
    pc[1] = gpt_center[1] - point[1];
    pc[2] = gpt_center[2] - point[2];

    // Fill the work array with the Boys function values
    arg = gamma*(pc[0]*pc[0] + pc[1]*pc[1] + pc[2]*pc[2]);
    for (long nu=abs(shell_type0)+abs(shell_type1); nu>=0; nu--) {
        work_boys[nu] = boys_function(nu, arg);
    }

    // Iterate over all combinations of Cartesian exponents
    i2p.reset(abs(shell_type0), abs(shell_type1));
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
        work_cart[i2p.offset] += pre*scales0[i2p.ibasis0]*scales1[i2p.ibasis1]*arg;
    } while (i2p.inc());
}
