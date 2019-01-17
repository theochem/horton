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


#include "horton/gbasis/common.h"

#include <cmath>
#include <iostream>


long fac(long n) {
    long result = 1;
    while (n > 1) {
        result *= n;
        n--;
    }
    return result;
}


long fac2(long n) {
    long result = 1;
    while (n > 1) {
        result *= n;
        n -= 2;
    }
    return result;
}


long binom(long n, long m) {
    long numer = 1, denom = 1;
    while (n > m) {
        numer *= n;
        denom *= (n-m);
        n--;
    }
    return numer/denom;
}


long get_shell_nbasis(long shell_type) {
    if (shell_type > 0) {
        // Cartesian
        return ((shell_type+1)*(shell_type+2))/2;
    } else if (shell_type == -1) {
        // should not happen.
        return -1;
    } else {
        // Pure
        return -2*shell_type+1;
    }
}


long get_max_shell_type() {
    return MAX_SHELL_TYPE;
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


/*

   Auxiliary functions for Gaussian integrals

*/

void compute_gpt_center(double alpha0, const double* r0, double alpha1, const double* r1, double gamma_inv, double* gpt_center) {
    gpt_center[0] = (alpha0*r0[0] + alpha1*r1[0])*gamma_inv;
    gpt_center[1] = (alpha0*r0[1] + alpha1*r1[1])*gamma_inv;
    gpt_center[2] = (alpha0*r0[2] + alpha1*r1[2])*gamma_inv;
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


/*

    Auxiliary functions for r^alpha integrals

*/


double cit(int i, double t, int m) {
  double result = pow(t, m);
  for (int j=1; j <= i; j++)
    result *= 4.0/((2.0*j + 1.0)*(2.0*j));
  return result;
}


long jfac(int j, int n) {
  double result = j;
  for (int i=1; i < n; i++)
    result *= j-i;
  return result;
}


double dtaylor(int n, double alpha, double t, double prefac) {
  double taylor0, taylor1;
  int j = 0;
  double gj = (alpha+3.0)/2.0;
  if (n == 0) {
    // s type orbitals
    taylor0 = prefac*tgamma(j+gj)*cit(j, t, j);
    j += 1;
    taylor1 = taylor0 + prefac*tgamma(j+gj)*cit(j, t, j);
    while (fabs(taylor1-taylor0) > 1e-9) {
      taylor0 = taylor1;
      j += 1;
      taylor1 = taylor0 + prefac*tgamma(j+gj)*cit(j, t, j);
    }
  } else {
    // higher angular moment
    // matrix with the coefficients that come up from the derivatives (up to 10)
    double matrix[10][11] =
      {{1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 3.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 4.0, 6.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 5.0, 10.0, 10.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 7.0, 21.0, 35.0, 35.0, 21.0, 7.0, 1.0, 0.0, 0.0, 0.0},
       {1.0, 8.0, 28.0, 56.0, 70.0, 56.0, 28.0, 8.0, 1.0, 0.0, 0.0},
       {1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0, 0.0},
       {1.0, 10.0, 45.0, 120.0, 210.0, 252.0, 210.0, 120.0, 45.0, 10.0, 1.0}};
    // first cycle
    taylor0 = prefac*tgamma(j+gj)*cit(j, t, j);
    j += 1;
    taylor1 = taylor0;
    if (j-n >= 0)
      taylor1 += pow(-1, n)*(jfac(j, n))*prefac*tgamma(j+gj)*cit(j, t, j-n);
    taylor1 += prefac*tgamma(j+gj)*cit(j, t, j);
    for (int k = n-1; k >= 1; k--) {
      if (j-k >= 0)
        taylor1 += pow(-1, k)*matrix[n-1][k]*(jfac(j, k))*prefac*tgamma(j+gj)*cit(j, t, j-k);
    }
    // other cycles
    while (fabs(taylor1-taylor0) > 1e-10) {
      taylor0 = taylor1;
      j += 1;
      if (j-n >= 0)
        taylor1 += pow(-1, n)*(jfac(j, n))*prefac*tgamma(j+gj)*cit(j, t, j-n);
      taylor1 += prefac*tgamma(j+gj)*cit(j, t, j);
      for (int k = n-1; k >= 1; k--) {
        if (j-k >= 0)
          taylor1 += pow(-1, k)*matrix[n-1][k]*(jfac(j, k))*prefac*tgamma(j+gj)*cit(j, t, j-k);
      }
    }
  }
  return taylor1;
}
