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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <cstdlib>
#include "horton/gbasis/boys.h"
#include "horton/gbasis/cartpure.h"
#include "horton/gbasis/common.h"
#include "horton/gbasis/ints.h"
using std::abs;

/*

    The implementation of the two-index operators in this module are based on
    the paper "Gaussian-Expansion Methods for Molecular Integrals", H. Taketa,
    S. Huzinga, K. O-ohata, Journal of the Physical Society of Japan, vol. 21,
    p. 2313, y. 1966. Be aware that there are some misprints in the paper:

    - eq. 2.18: CP_x should be p_x, i.e. change of sign.
    - list of expressions at the bottom of p. 2319:
        case l_1+l_2=2, G_1 = -f_1*p - f_2/(2*gamma)

    The implementation of the (multipole) moment integrals in this module are
    based on the paper "Efficient recursive computation of molecular integrals
    over Cartesian Gaussian functions", S. Obara, A. Saika, Journal of Chemical
    Physics, vol. 84, p. 3963, y. 1986.

*/


/*

   GB2Integral

*/


GB2Integral::GB2Integral(long max_shell_type): GBCalculator(max_shell_type) {
    nwork = max_nbasis*max_nbasis;
    work_cart = new double[nwork];
    work_pure = new double[nwork];
}

void GB2Integral::reset(long _shell_type0, long _shell_type1, const double* _r0, const double* _r1) {
    if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type)) {
      throw std::domain_error("shell_type0 out of range.");
    }
    if ((_shell_type1 < -max_shell_type) || (_shell_type1 > max_shell_type)) {
      throw std::domain_error("shell_type1 out of range.");
    }
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

   GB2OverlapIntegral

*/


void GB2OverlapIntegral::add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1) {
    double pre, gamma_inv;
    double gpt_center[3];

    gamma_inv = 1.0/(alpha0 + alpha1);
    pre = coeff*exp(-alpha0*alpha1*gamma_inv*dist_sq(r0, r1));
    compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
    i2p.reset(abs(shell_type0), abs(shell_type1));
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

    i2p.reset(abs(shell_type0), abs(shell_type1));
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

   GB2NuclearAttractionIntegrals

*/


GB2AttractionIntegral::GB2AttractionIntegral(long max_shell_type, double* charges,
                                             double* centers, long ncharge) :
            GB2Integral(max_shell_type), charges(charges), centers(centers), ncharge(ncharge) {
    work_g0 = new double[2*max_shell_type+1];
    work_g1 = new double[2*max_shell_type+1];
    work_g2 = new double[2*max_shell_type+1];
    work_boys = new double[2*max_shell_type+1];
}


GB2AttractionIntegral::~GB2AttractionIntegral() {
    delete[] work_g0;
    delete[] work_g1;
    delete[] work_g2;
    delete[] work_boys;
}


void GB2AttractionIntegral::add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1) {
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
        // for (long nu=abs(shell_type0)+abs(shell_type1); nu>=0; nu--) {
        //    work_boys[nu] = boys_function(nu, arg);
        //}
        /*
           Laplace transform of the potential
        */

        int mmax = abs(shell_type0) + abs(shell_type1);
        laplace_of_potential(gamma, arg, mmax, work_boys);

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

            // Finally add to the work array, accounting for opposite charge of electron and nucleus
            work_cart[i2p.offset] -= pre*scales0[i2p.ibasis0]*scales1[i2p.ibasis1]*arg*charges[icharge];
        } while (i2p.inc());
    }
}


void GB2NuclearAttractionIntegral::laplace_of_potential(double gamma, double arg, long mmax,
                                                        double* output) {
  boys_function_array(mmax, arg, output);
}


void GB2ErfAttractionIntegral::laplace_of_potential(double gamma, double arg, long mmax,
                                                    double* output) {
  double efac = mu*mu/(mu*mu + gamma);
  boys_function_array(mmax, arg*efac, output);
  double prefac = sqrt(efac);
  for (long m=0; m <= mmax; m++) {
    output[m] *= prefac;
    prefac *= efac;
  }
}


void GB2GaussAttractionIntegral::laplace_of_potential(double gamma, double arg, long mmax,
                                                      double* output) {
  double afac = alpha/(gamma+alpha);
  double prefac = (M_PI/(gamma+alpha))*sqrt(M_PI/(gamma+alpha))*c*exp(-arg*afac)*(gamma/(2*M_PI));
  for (long m=0; m <= mmax; m++) {
    output[m] = prefac;
    prefac *= afac;
  }
}


/*

 GB2MomentIntegral

 */


double moment_helper(long n0, long n1, long n2, double pa, double pb, double pc, double gamma_inv) {
    /*
     The Obara Saika Scheme, equations A7 and A8 in "Efficient recursive computation of
     molecular integrals over Cartesian Gaussian functions", S. Obara, A. Saika, Journal
     of Chemical Physics, vol. 84, p. 3963, y. 1986.

     Note that we use a shift in the indices:
     (0_A|R(0)|0_B) in the paper is stored in work_mm[1][1][1].
     In this allows to set all terms in with the angular momentum index becomes -1 to 0,
     avoiding a lot of if-statements.
    */

    double result;
    // if n2 == 0, we just need the overlap.
    if (n2 == 0) {
        result = gb_overlap_int1d(n0, n1, pa, pb, (2.0e0 * gamma_inv));
    } else {
        long m, l, n;
        double work_mm[n0+2][n1+2][n2+2];

        for (l=0; l < (n0+2); l++) {
            for (m=0; m < (n1+2); m++) {
                for (n=0; n < (n2+2); n++) {
                    work_mm[l][m][n] = 0.0e0;
                }
            }
        }

        // The auxiliary overlap (0_A|R(0)|0_B) is stored in work_mm[1][1][1]:
        work_mm[1][1][1] = gb_overlap_int1d(0, 0, pa, pb, (2.0e0 * gamma_inv));

        // Equation A8 in the Obara-Saika paper for (0_A|R(mu + 1)|0_B):
        for (n=1; n < (n2+1); n++) {
            work_mm[1][1][n+1] = pc * work_mm[1][1][n] + (gamma_inv * (n-1) * work_mm[1][1][n-1]);
        }

        // Equation A7 in the Obara-Saika paper for (0_A|R(mu)|b + 1):
        for (m=1; m < (n1+1); m++) {
            for (n=1; n <= (n2+1); n++) {
                work_mm[1][m+1][n] = pb * work_mm[1][m][n]
                                                + (gamma_inv * (m-1) * work_mm[1][m-1][n])
                                                + (gamma_inv * (n-1) * work_mm[1][m][n-1]);
            }
        }

        // Equation A7 in the Obara-Saika paper for (a + 1|R(mu)|b):
        for (l=1; l < (n0+1); l++) {
            for (m=1; m <= (n1+1); m++) {
                for (n=1; n <= (n2+1); n++) {
                    work_mm[l+1][m][n] = pa * work_mm[l][m][n]
                                                + (gamma_inv * (l-1) * work_mm[l-1][m][n])
                                                + (gamma_inv * (m-1) * work_mm[l][m-1][n])
                                                + (gamma_inv * (n-1) * work_mm[l][m][n-1]);
                }
            }
        }

        result = work_mm[n0+1][n1+1][n2+1];
    }

    return result;
}


GB2MomentIntegral::GB2MomentIntegral(long max_shell_type, long* xyz, double* center)
: GB2Integral(max_shell_type), xyz(xyz), center(center) {
    if (xyz[0] + xyz[1] + xyz[2] < 0)
    throw std::domain_error(" sum < 0");
    if ((xyz[0] < 0) || (xyz[1] < 0) || (xyz[2] < 0))
    throw std::domain_error(" all elements of xyz must be >= 0");
}


void GB2MomentIntegral::add(double coeff, double alpha0, double alpha1,
                            const double* scales0, const double* scales1) {
    double pre, gamma_inv, twogamma_inv;
    double gpt_center[3], pa[3], pb[3], pc[3];

    gamma_inv = 1.0/(alpha0 + alpha1);
    twogamma_inv = 0.50/(alpha0 + alpha1);
    pre = coeff*exp(-alpha0*alpha1*gamma_inv*dist_sq(r0, r1));
    compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
    pa[0] = gpt_center[0] - r0[0];
    pa[1] = gpt_center[1] - r0[1];
    pa[2] = gpt_center[2] - r0[2];
    pb[0] = gpt_center[0] - r1[0];
    pb[1] = gpt_center[1] - r1[1];
    pb[2] = gpt_center[2] - r1[2];
    pc[0] = gpt_center[0] - center[0];
    pc[1] = gpt_center[1] - center[1];
    pc[2] = gpt_center[2] - center[2];
    i2p.reset(abs(shell_type0), abs(shell_type1));

    do {
        work_cart[i2p.offset] += pre*(
            moment_helper(i2p.n0[0], i2p.n1[0], xyz[0], pa[0], pb[0], pc[0], twogamma_inv)*
            moment_helper(i2p.n0[1], i2p.n1[1], xyz[1], pa[1], pb[1], pc[1], twogamma_inv)*
            moment_helper(i2p.n0[2], i2p.n1[2], xyz[2], pa[2], pb[2], pc[2], twogamma_inv)*
            scales0[i2p.ibasis0]*scales1[i2p.ibasis1]);
    } while (i2p.inc());
}

/*

   GB4Integral

*/


GB4Integral::GB4Integral(long max_shell_type)
    : GBCalculator(max_shell_type), shell_type0(0), shell_type1(0), shell_type2(0),
      shell_type3(0), r0(NULL), r1(NULL), r2(NULL), r3(NULL) {
  nwork = max_nbasis*max_nbasis*max_nbasis*max_nbasis;
  work_cart = new double[nwork];
  work_pure = new double[nwork];
}

void GB4Integral::reset(
    long _shell_type0, long _shell_type1, long _shell_type2, long _shell_type3,
    const double* _r0, const double* _r1, const double* _r2, const double* _r3) {
  if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type))
    throw std::domain_error("shell_type0 out of range.");
  if ((_shell_type1 < -max_shell_type) || (_shell_type1 > max_shell_type))
    throw std::domain_error("shell_type1 out of range.");
  if ((_shell_type2 < -max_shell_type) || (_shell_type2 > max_shell_type))
    throw std::domain_error("shell_type2 out of range.");
  if ((_shell_type3 < -max_shell_type) || (_shell_type3 > max_shell_type))
    throw std::domain_error("shell_type3 out of range.");
  shell_type0 = _shell_type0;
  shell_type1 = _shell_type1;
  shell_type2 = _shell_type2;
  shell_type3 = _shell_type3;
  r0 = _r0;
  r1 = _r1;
  r2 = _r2;
  r3 = _r3;
  // We make use of the fact that a floating point zero consists of
  // consecutive zero bytes.
  memset(work_cart, 0, nwork*sizeof(double));
  memset(work_pure, 0, nwork*sizeof(double));
}

void GB4Integral::cart_to_pure() {
  /* The initial results are always stored in work_cart. The projection routine always
     outputs its result in work_pure. Once that is done, the pointers to both blocks are
     swapped such that the final result is always back in work_cart.
  */

  // Transform along index 0
  if (shell_type0 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type0,
                     1,                                    // anterior
                     get_shell_nbasis(abs(shell_type1))*
                     get_shell_nbasis(abs(shell_type2))*
                     get_shell_nbasis(abs(shell_type3)));  // posterior
    swap_work();
  }

  // Project along index 1
  if (shell_type1 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type1,
                     get_shell_nbasis(shell_type0),        // anterior
                     get_shell_nbasis(abs(shell_type2))*
                     get_shell_nbasis(abs(shell_type3)));  // posterior
    swap_work();
  }

  // Project along index 2
  if (shell_type2 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type2,
                     get_shell_nbasis(shell_type0)*
                     get_shell_nbasis(shell_type1),        // anterior
                     get_shell_nbasis(abs(shell_type3)));  // posterior
    swap_work();
  }

  // Project along index 3
  if (shell_type3 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type3,
                     get_shell_nbasis(shell_type0)*
                     get_shell_nbasis(shell_type1)*
                     get_shell_nbasis(shell_type2),  // anterior
                     1);                             // posterior
    swap_work();
  }
}


/*

   GB4IntegralLibInt

*/

#if LIBINT2_REALTYPE != double
#error LibInt must be compiled with REALTYPE=double.
#endif

#if LIBINT2_SUPPORT_ERI == 0
#error LibInt must be compiled with support for electron repulsion integrals.
#endif

#if LIBINT2_MAX_AM_ERI < MAX_SHELL_TYPE
#error LibInt must be compiled with an angular momentum limit of at least MAX_SHELL_TYPE.
#endif

GB4IntegralLibInt::GB4IntegralLibInt(long max_shell_type)
    : GB4Integral(max_shell_type),
      libint_args{{0, NULL, 0.0}, {0, NULL, 0.0}, {0, NULL, 0.0}, {0, NULL, 0.0}},
      order{0, 0, 0, 0}, ab{0.0, 0.0, 0.0}, cd{0.0, 0.0, 0.0}, ab2(0.0), cd2(0.0) {
  libint2_init_eri(&erieval, max_shell_type, 0);
  erieval.contrdepth = 1;
}

GB4IntegralLibInt::~GB4IntegralLibInt() {
  libint2_cleanup_eri(&erieval);
}

void GB4IntegralLibInt::reset(
    long _shell_type0, long _shell_type1, long _shell_type2, long _shell_type3,
    const double* _r0, const double* _r1, const double* _r2, const double* _r3) {
  GB4Integral::reset(_shell_type0, _shell_type1, _shell_type2, _shell_type3, _r0, _r1, _r2, _r3);

  // Store the arguments for libint such that they can be reordered conveniently.

  libint_args[0].am = abs(shell_type0);
  libint_args[1].am = abs(shell_type1);
  libint_args[2].am = abs(shell_type2);
  libint_args[3].am = abs(shell_type3);

  libint_args[0].r = r0;
  libint_args[1].r = r1;
  libint_args[2].r = r2;
  libint_args[3].r = r3;

  /* Figure out the ordering of the shell quartet that is compatible with libint. Note
     that libint uses the chemist's notation, while HORTON uses the physicists notation
     for the indexes of a four-index operator.

     The arguments must be reordered such that the following conditions are met:
        abs(shell_type0) >= abs(shell_type2)
        abs(shell_type1) >= abs(shell_type3)
        abs(shell_type1+abs(shell_type3) >= abs(shell_type0)+abs(shell_type2)
     using the following permutations respectively:
        (2,1,0,3)
        (0,3,2,1)
        (1,0,3,2)
     The first two permutations are set directly below. The last one is applied after the
     first two are set.
  */

  if (libint_args[0].am >= libint_args[2].am) {
    order[0] = 0;
    order[2] = 2;
  } else {
    order[0] = 2;
    order[2] = 0;
  }

  if (libint_args[1].am >= libint_args[3].am) {
    order[1] = 1;
    order[3] = 3;
  } else {
    order[1] = 3;
    order[3] = 1;
  }

  if (libint_args[1].am + libint_args[3].am < libint_args[0].am + libint_args[2].am) {
    long tmp;
    tmp = order[0];
    order[0] = order[1];
    order[1] = tmp;
    tmp = order[2];
    order[2] = order[3];
    order[3] = tmp;
  }

  /* Compute the distances squared of AB and CD.
     AB corresponds to libint_args[order[0]].r - libint_args[order[2]].r
     CD corresponds to libint_args[order[1]].r - libint_args[order[3]].r
  */

  ab[0] = libint_args[order[0]].r[0] - libint_args[order[2]].r[0];
  ab[1] = libint_args[order[0]].r[1] - libint_args[order[2]].r[1];
  ab[2] = libint_args[order[0]].r[2] - libint_args[order[2]].r[2];
  ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
#if LIBINT2_DEFINED(eri, AB_x)
  erieval.AB_x[0] = ab[0];
#endif
#if LIBINT2_DEFINED(eri, AB_y)
  erieval.AB_y[0] = ab[1];
#endif
#if LIBINT2_DEFINED(eri, AB_z)
  erieval.AB_z[0] = ab[2];
#endif

  cd[0] = libint_args[order[1]].r[0] - libint_args[order[3]].r[0];
  cd[1] = libint_args[order[1]].r[1] - libint_args[order[3]].r[1];
  cd[2] = libint_args[order[1]].r[2] - libint_args[order[3]].r[2];
  cd2 = cd[0] * cd[0] + cd[1] * cd[1] + cd[2] * cd[2];
#if LIBINT2_DEFINED(eri, CD_x)
  erieval.CD_x[0] = cd[0];
#endif
#if LIBINT2_DEFINED(eri, CD_y)
  erieval.CD_y[0] = cd[1];
#endif
#if LIBINT2_DEFINED(eri, CD_z)
  erieval.CD_z[0] = cd[2];
#endif
}


void GB4IntegralLibInt::add(
    double coeff, double alpha0, double alpha1, double alpha2, double alpha3,
    const double* scales0, const double* scales1, const double* scales2, const double* scales3) {
  /*
      Store the arguments for libint such that they can be reordered
      conveniently.
  */

  libint_args[0].alpha = alpha0;
  libint_args[1].alpha = alpha1;
  libint_args[2].alpha = alpha2;
  libint_args[3].alpha = alpha3;

  /*
      Precompute some variables for libint. The approach here is not
      super-efficient, but why bother...
  */

  const double gammap = libint_args[order[0]].alpha + libint_args[order[2]].alpha;
  const double gammap_inv = 1.0/gammap;
  double p[3];
  compute_gpt_center(libint_args[order[0]].alpha, libint_args[order[0]].r,
                     libint_args[order[2]].alpha, libint_args[order[2]].r,
                     gammap_inv, p);
  const double pa[3] = {
      p[0] - libint_args[order[0]].r[0],
      p[1] - libint_args[order[0]].r[1],
      p[2] - libint_args[order[0]].r[2]
  };
  const double pb[3] = {
      p[0] - libint_args[order[2]].r[0],
      p[1] - libint_args[order[2]].r[1],
      p[2] - libint_args[order[2]].r[2]
  };

#if LIBINT2_DEFINED(eri, PA_x)
  erieval.PA_x[0] = pa[0];
#endif
#if LIBINT2_DEFINED(eri, PA_y)
  erieval.PA_y[0] = pa[1];
#endif
#if LIBINT2_DEFINED(eri, PA_z)
  erieval.PA_z[0] = pa[2];
#endif
#if LIBINT2_DEFINED(eri, PB_x)
  erieval.PB_x[0] = pb[0];
#endif
#if LIBINT2_DEFINED(eri, PB_y)
  erieval.PB_y[0] = pb[1];
#endif
#if LIBINT2_DEFINED(eri, PB_z)
  erieval.PB_z[0] = pb[2];
#endif
#if LIBINT2_DEFINED(eri, oo2z)
  erieval.oo2z[0] = 0.5*gammap_inv;
#endif

  const double gammaq = libint_args[order[1]].alpha + libint_args[order[3]].alpha;
  const double gammaq_inv = 1.0/gammaq;
  double q[3];
  compute_gpt_center(libint_args[order[1]].alpha, libint_args[order[1]].r,
                     libint_args[order[3]].alpha, libint_args[order[3]].r,
                     gammaq_inv, q);
  const double qc[3] = {
    q[0] - libint_args[order[1]].r[0],
    q[1] - libint_args[order[1]].r[1],
    q[2] - libint_args[order[1]].r[2]
  };
  const double qd[3] = {
    q[0] - libint_args[order[3]].r[0],
    q[1] - libint_args[order[3]].r[1],
    q[2] - libint_args[order[3]].r[2]
  };

#if LIBINT2_DEFINED(eri, QC_x)
  erieval.QC_x[0] = qc[0];
#endif
#if LIBINT2_DEFINED(eri, QC_y)
  erieval.QC_y[0] = qc[1];
#endif
#if LIBINT2_DEFINED(eri, QC_z)
  erieval.QC_z[0] = qc[2];
#endif
#if LIBINT2_DEFINED(eri, QD_x)
  erieval.QD_x[0] = qd[0];
#endif
#if LIBINT2_DEFINED(eri, QD_y)
  erieval.QD_y[0] = qd[1];
#endif
#if LIBINT2_DEFINED(eri, QD_z)
  erieval.QD_z[0] = qd[2];
#endif
#if LIBINT2_DEFINED(eri, oo2e)
  erieval.oo2e[0] = 0.5*gammaq_inv;
#endif

  // Prefactors for inter-electron transfer relation
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_x)
  erieval.TwoPRepITR_pfac0_0_x[0] = -(libint_args[order[2]].alpha*ab[0] +
                                      libint_args[order[3]].alpha*cd[0])/gammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_y)
  erieval.TwoPRepITR_pfac0_0_y[0] = -(libint_args[order[2]].alpha*ab[1] +
                                      libint_args[order[3]].alpha*cd[1])/gammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_z)
  erieval.TwoPRepITR_pfac0_0_z[0] = -(libint_args[order[2]].alpha*ab[2] +
                                      libint_args[order[3]].alpha*cd[2])/gammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_x)
  erieval.TwoPRepITR_pfac0_1_x[0] = -(libint_args[order[2]].alpha*ab[0] +
                                      libint_args[order[3]].alpha*cd[0])/gammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_y)
  erieval.TwoPRepITR_pfac0_1_y[0] = -(libint_args[order[2]].alpha*ab[1] +
                                      libint_args[order[3]].alpha*cd[1])/gammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_z)
  erieval.TwoPRepITR_pfac0_1_z[0] = -(libint_args[order[2]].alpha*ab[2] +
                                      libint_args[order[3]].alpha*cd[2])/gammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac1_0)
  erieval.TwoPRepITR_pfac1_0[0] = -gammaq*gammap_inv;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac1_1)
  erieval.TwoPRepITR_pfac1_1[0] = -gammap*gammaq_inv;
#endif

  const double pq[3] = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
  const double pq2 = pq[0]*pq[0] + pq[1]*pq[1] + pq[2]*pq[2];
  const double eta = gammap + gammaq;
  const double eta_inv = 1.0/eta;
  double w[3];
  compute_gpt_center(gammap, p, gammaq, q, eta_inv, w);

#if LIBINT2_DEFINED(eri, WP_x)
  erieval.WP_x[0] = w[0] - p[0];
#endif
#if LIBINT2_DEFINED(eri, WP_y)
  erieval.WP_y[0] = w[1] - p[1];
#endif
#if LIBINT2_DEFINED(eri, WP_z)
  erieval.WP_z[0] = w[2] - p[2];
#endif
#if LIBINT2_DEFINED(eri, WQ_x)
  erieval.WQ_x[0] = w[0] - q[0];
#endif
#if LIBINT2_DEFINED(eri, WQ_y)
  erieval.WQ_y[0] = w[1] - q[1];
#endif
#if LIBINT2_DEFINED(eri, WQ_z)
  erieval.WQ_z[0] = w[2] - q[2];
#endif
#if LIBINT2_DEFINED(eri, oo2ze)
  erieval.oo2ze[0] = 0.5*eta_inv;
#endif
#if LIBINT2_DEFINED(eri, roz)
  erieval.roz[0] = gammaq*eta_inv;
#endif
#if LIBINT2_DEFINED(eri, roe)
  erieval.roe[0] = gammap*eta_inv;
#endif

  /*
      Arguments for the kernel (using Boy's function or something else)
  */

  const double k1 = exp(-libint_args[order[0]].alpha * libint_args[order[2]].alpha *
                        ab2 * gammap_inv);
  const double k2 = exp(-libint_args[order[1]].alpha * libint_args[order[3]].alpha *
                        cd2 * gammaq_inv);
#define PI_POW_3_2 5.5683279968317078
  const double pfac = PI_POW_3_2*k1*k2*eta_inv*sqrt(eta_inv)*coeff;
  const double rho = 1.0/(gammaq_inv + gammap_inv);
  const double t = pq2*rho;

  /*
      Laplace transform of the potential
  */

  int mmax = libint_args[0].am + libint_args[1].am + libint_args[2].am + libint_args[3].am;
  double kernel[4*MAX_SHELL_TYPE+1];
  laplace_of_potential(pfac, rho, t, mmax, kernel);

#define TEST_END_BOYS mmax--; if (mmax < 0) goto end_boys;
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(0))
  erieval.LIBINT_T_SS_EREP_SS(0)[0] = kernel[0]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(1))
  erieval.LIBINT_T_SS_EREP_SS(1)[0] = kernel[1]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(2))
  erieval.LIBINT_T_SS_EREP_SS(2)[0] = kernel[2]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(3))
  erieval.LIBINT_T_SS_EREP_SS(3)[0] = kernel[3]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(4))
  erieval.LIBINT_T_SS_EREP_SS(4)[0] = kernel[4]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(5))
  erieval.LIBINT_T_SS_EREP_SS(5)[0] = kernel[5]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(6))
  erieval.LIBINT_T_SS_EREP_SS(6)[0] = kernel[6]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(7))
  erieval.LIBINT_T_SS_EREP_SS(7)[0] = kernel[7]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(8))
  erieval.LIBINT_T_SS_EREP_SS(8)[0] = kernel[8]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(9))
  erieval.LIBINT_T_SS_EREP_SS(9)[0] = kernel[9]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(10))
  erieval.LIBINT_T_SS_EREP_SS(10)[0] = kernel[10]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(11))
  erieval.LIBINT_T_SS_EREP_SS(11)[0] = kernel[11]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(12))
  erieval.LIBINT_T_SS_EREP_SS(12)[0] = kernel[12]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(13))
  erieval.LIBINT_T_SS_EREP_SS(13)[0] = kernel[13]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(14))
  erieval.LIBINT_T_SS_EREP_SS(14)[0] = kernel[14]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(15))
  erieval.LIBINT_T_SS_EREP_SS(15)[0] = kernel[15]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(16))
  erieval.LIBINT_T_SS_EREP_SS(16)[0] = kernel[16]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(17))
  erieval.LIBINT_T_SS_EREP_SS(17)[0] = kernel[17]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(18))
  erieval.LIBINT_T_SS_EREP_SS(18)[0] = kernel[18]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(19))
  erieval.LIBINT_T_SS_EREP_SS(19)[0] = kernel[19]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(20))
  erieval.LIBINT_T_SS_EREP_SS(20)[0] = kernel[20]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(21))
  erieval.LIBINT_T_SS_EREP_SS(21)[0] = kernel[21]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(22))
  erieval.LIBINT_T_SS_EREP_SS(22)[0] = kernel[22]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(23))
  erieval.LIBINT_T_SS_EREP_SS(23)[0] = kernel[23]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(24))
  erieval.LIBINT_T_SS_EREP_SS(24)[0] = kernel[24]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(25))
  erieval.LIBINT_T_SS_EREP_SS(25)[0] = kernel[25]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(26))
  erieval.LIBINT_T_SS_EREP_SS(26)[0] = kernel[26]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(27))
  erieval.LIBINT_T_SS_EREP_SS(27)[0] = kernel[27]; TEST_END_BOYS;
#endif
#if LIBINT2_DEFINED(eri, LIBINT_T_SS_EREP_SS(28))
  erieval.LIBINT_T_SS_EREP_SS(28)[0] = kernel[28]; TEST_END_BOYS;
#endif
end_boys:

  /*
      Actual computation of all the integrals in this shell-set by libint
  */

  if ((libint_args[0].am == 0) && (libint_args[1].am == 0) &&
      (libint_args[2].am == 0) && (libint_args[3].am == 0)) {
    work_cart[0] += erieval.LIBINT_T_SS_EREP_SS(0)[0] *
                    scales0[0] * scales1[0] * scales2[0] * scales3[0];
  } else {
    libint2_build_eri
        [libint_args[order[0]].am]
        [libint_args[order[2]].am]
        [libint_args[order[1]].am]
        [libint_args[order[3]].am]
        (&erieval);

    /*
        Extract the integrals from the erieval target array, taking into account
        the reordering of the indexes.
    */

    double strides[4];
    strides[3] = 1;
    strides[2] = get_shell_nbasis(libint_args[3].am);
    strides[1] = strides[2]*get_shell_nbasis(libint_args[2].am);
    strides[0] = strides[1]*get_shell_nbasis(libint_args[1].am);
    const int nbasis0 = get_shell_nbasis(libint_args[order[0]].am);
    const int nbasis1 = get_shell_nbasis(libint_args[order[1]].am);
    const int nbasis2 = get_shell_nbasis(libint_args[order[2]].am);
    const int nbasis3 = get_shell_nbasis(libint_args[order[3]].am);
    const double* scales[4];
    scales[0] = scales0;
    scales[1] = scales1;
    scales[2] = scales2;
    scales[3] = scales3;
    long counter = 0;
    for (int i0=0; i0 < nbasis0; i0++) {
      for (int i2=0; i2 < nbasis2; i2++) {
        for (int i1=0; i1 < nbasis1; i1++) {
          for (int i3=0; i3 < nbasis3; i3++) {
            const long offset = i0*strides[order[0]] + i1*strides[order[1]] +
                                i2*strides[order[2]] + i3*strides[order[3]];
            work_cart[offset] += erieval.targets[0][counter] *
                                 scales[order[0]][i0] * scales[order[1]][i1] *
                                 scales[order[2]][i2] * scales[order[3]][i3];
            counter++;
          }
        }
      }
    }
  }
}


void GB4ElectronRepulsionIntegralLibInt::laplace_of_potential(
    double prefac, double rho, double t, long mmax, double* output) {
  boys_function_array(mmax, t, output);
  prefac *= 2.0*M_PI/rho;
  for (long m=0; m <= mmax; m++) {
    output[m] *= prefac;
  }
}


void GB4ErfIntegralLibInt::laplace_of_potential(
    double prefac, double rho, double t, long mmax, double* output) {
  double efac = mu*mu/(mu*mu + rho);
  boys_function_array(mmax, t*efac, output);
  prefac *= 2.0*M_PI/rho*sqrt(efac);
  for (long m=0; m <= mmax; m++) {
    output[m] *= prefac;
    prefac *= efac;
  }
}


void GB4GaussIntegralLibInt::laplace_of_potential(double prefac, double rho, double t,
                                                  long mmax, double* output) {
  double afac = alpha/(rho+alpha);
  prefac *= (sqrt(M_PI*M_PI*M_PI)/(rho+alpha))*sqrt(1.0/(rho+alpha))*c*exp(-t*afac);
  for (long m=0; m <= mmax; m++) {
    output[m] = prefac;
    prefac *= afac;
  }
}


void GB4RAlphaIntegralLibInt::laplace_of_potential(double prefac, double rho, double t,
                                                   long mmax, double* output) {
  if (mmax > 10)
    throw std::domain_error("mmax > 10, highest cartesian angular value implemented is 10");
  prefac *= exp(-t)/(rho*sqrt(rho));
  double tfactor = ((4.0*M_PI)/(2.0*sqrt(pow(rho, alpha))));
  for (long m=0; m <= mmax; m++) {
    output[m] = tfactor*dtaylor(m, alpha, t, prefac);
  }
}
