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

// #define DEBUG

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


/*
    GB1GridFn
*/

GB1GridFn::GB1GridFn(long max_shell_type, long dim_work, long dim_output)
    : GBCalculator(max_shell_type), dim_work(dim_work), dim_output(dim_output),
      shell_type0(0), r0(NULL), point(NULL), i1p() {
  nwork = max_nbasis*dim_work;
  work_cart = new double[nwork];
  work_pure = new double[nwork];
}

void GB1GridFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type)) {
    throw std::domain_error("shell_type0 out of range.");
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
      1,          // anterior
      dim_work);  // posterior
    swap_work();
  }
}


/*
    GB1ExpGridOrbitalFn
*/

void GB1ExpGridOrbitalFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  if (shell_type0 != 0) {
    poly_work[1] = point[0] - r0[0];
    poly_work[2] = point[1] - r0[1];
    poly_work[3] = point[2] - r0[2];
    offset = fill_cartesian_polynomials(poly_work+1, abs(shell_type0))+1;
  } else {
    offset = 0;
  }
}

void GB1ExpGridOrbitalFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The (Cartesian) basis function evaluated in `point` (see reset method) is added to
  // work_cart, where ibasis0 is an index for the primitive in the current Cartesian
  // shell.
  for (long ibasis0=get_shell_nbasis(abs(shell_type0))-1; ibasis0 >= 0; ibasis0--) {
    work_cart[ibasis0] += pre*scales0[ibasis0]*poly_work[ibasis0+offset];
  }
}

void GB1ExpGridOrbitalFn::compute_point_from_exp(double* work_basis, double* coeffs,
                                                 long nbasis, double* output) {
  for (long i=0; i < norb; i++) {
    long iorb = iorbs[i];
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
      // Just evaluate the contribution of each basis function to an orbital.
      // The values of the basis functions at `point` (see reset method) are available in
      // work_basis.
      output[i] += coeffs[ibasis*nfn + iorb]*work_basis[ibasis];
    }
  }
}


/*
 GB1ExpGridOrbGradientFn
 */

void GB1ExpGridOrbGradientFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher and lower polynomials are required because of first derivative.
  offset_h1 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+1)+1;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
}

void GB1ExpGridOrbGradientFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // For every primitive, work_cart contains four values computed at `point` (see reset
  // method): the basis function and its derivatives toward x, y and z.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function derived toward x
    work_cart[3*i1p.ibasis0+0] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[3*i1p.ibasis0+0] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
      work_cart[3*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[3*i1p.ibasis0+1] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[3*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[3*i1p.ibasis0+2] += i1p.n0[2]*pre0*poly_work[i1p.ibasis0-nnotx-1+offset_l1];
    } while (i1p.inc());
}

void GB1ExpGridOrbGradientFn::compute_point_from_exp(double* work_basis, double* coeffs,
                                                     long nbasis, double* output) {
  for (long i=0; i < norb; i++) {
    double g_x = 0, g_y = 0, g_z = 0;
    long iorb = iorbs[i];
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
      g_x += coeffs[ibasis*nfn + iorb]*work_basis[ibasis*3+0];
      g_y += coeffs[ibasis*nfn + iorb]*work_basis[ibasis*3+1];
      g_z += coeffs[ibasis*nfn + iorb]*work_basis[ibasis*3+2];
    }
    output[i*3+0] += g_x;
    output[i*3+1] += g_y;
    output[i*3+2] += g_z;
  }
}


/*
    GB1DMGridDensityFn
*/

void GB1DMGridDensityFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  if (shell_type0 != 0) {
    poly_work[1] = point[0] - r0[0];
    poly_work[2] = point[1] - r0[1];
    poly_work[3] = point[2] - r0[2];
    offset = fill_cartesian_polynomials(poly_work+1, abs(shell_type0))+1;
  } else {
    offset = 0;
  }
}

void GB1DMGridDensityFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The (Cartesian) basis function evaluated in `point` (see reset method) is added to
  // work_cart, where ibasis0 is an index for the primitive in the current Cartesian
  // shell.
  for (long ibasis0=get_shell_nbasis(abs(shell_type0))-1; ibasis0 >= 0; ibasis0--) {
    work_cart[ibasis0] += pre*scales0[ibasis0]*poly_work[ibasis0+offset];
  }
}

void GB1DMGridDensityFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The values of the basis functions at `point` (see reset method) are available in
  // work_basis.

  // The epsilon argument allows one to skip part of the calculation where the density is
  // low.
  if (epsilon > 0) {
    double absmax_basis = 0.0;
    // compute the maximum basis function
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
      double tmp = fabs(work_basis[ibasis]);
      if (tmp > absmax_basis) absmax_basis = tmp;
    }
    // upper estimate of the density
    double rho_upper = 0.0;
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
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
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    // if the contribution of this loop is smaller than epsilon/nbasis, skip it.
    if (epsilon > 0) {
      if (fabs(work_basis[ibasis0])*dmmaxrow[ibasis0] < epsilon)
          continue;
    }
    double tmp = 0;
    // Loop for off-diagonal contributions of density matrix.
    for (long ibasis1=ibasis0-1; ibasis1 >= 0; ibasis1--) {
      tmp += work_basis[ibasis1]*dm[ibasis0*nbasis+ibasis1];
    }
    // Finally, also include diagonal contribution
    rho += (2*tmp +  // off-diagonal
            dm[ibasis0*(nbasis+1)]*work_basis[ibasis0])  // diagonal
           *work_basis[ibasis0];
  }
  *output += rho;
}

void GB1DMGridDensityFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The potential in `point` (see reset method) is given in `*pot`. It is the functional
  // derivative of the energy w.r.t. to the density in `point`. This gets transformed to
  // a contribution to the Fock matrix, i.e. the derivative of the energy w.r.t. the 1RDM
  // elements, using the chain rule.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp1 = (*pot)*work_basis[ibasis0];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double tmp2 = tmp1*work_basis[ibasis1];
      fock[ibasis0*nbasis+ibasis1] += tmp2;
      if (ibasis0 != ibasis1) fock[ibasis1*nbasis+ibasis0] += tmp2;
    }
  }
}



/*
    GB1DMGridGradientFn
*/

void GB1DMGridGradientFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher and lower polynomials are required because of first derivative.
  offset_h1 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+1)+1;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
}

void GB1DMGridGradientFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // For every primitive, work_cart contains four values computed at `point` (see reset
  // method): the basis function and its derivatives toward x, y and z.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[4*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[4*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[4*i1p.ibasis0+1] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
    work_cart[4*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[4*i1p.ibasis0+2] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[4*i1p.ibasis0+3] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[4*i1p.ibasis0+3] += i1p.n0[2]*pre0*poly_work[i1p.ibasis0-nnotx-1+offset_l1];
  } while (i1p.inc());
}

void GB1DMGridGradientFn::compute_point_from_dm(double* work_basis, double* dm,
                                                long nbasis, double* output,
                                                double epsilon, double* dmmaxrow) {
  // The value of the basis function and its derivatives toward x, y and z (at `point`)
  // are available in work_basis. These are used, together with the 1RDM coefficients,
  // to compute the density gradient at `point`.
  double rho_x = 0, rho_y = 0, rho_z = 0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
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

void GB1DMGridGradientFn::compute_fock_from_pot(double* pot, double* work_basis,
                                                long nbasis, double* fock) {
  // The functional derivative of the energy w.r.t. to the density gradient is given in
  // the argument *pot (three components: x, y, and z). The chain rule is used to
  // transform these into a contribution to Fock matrix.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp0 = work_basis[ibasis0*4];
    double tmp1 = pot[0]*work_basis[ibasis0*4+1] +
                  pot[1]*work_basis[ibasis0*4+2] +
                  pot[2]*work_basis[ibasis0*4+3];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double result = tmp0*(pot[0]*work_basis[ibasis1*4+1] +
                            pot[1]*work_basis[ibasis1*4+2] +
                            pot[2]*work_basis[ibasis1*4+3]) +
                      tmp1*work_basis[ibasis1*4];
      fock[ibasis1*nbasis+ibasis0] += result;
      if (ibasis1 != ibasis0) {
        // Enforce symmetry
        fock[ibasis0*nbasis+ibasis1] += result;
      }
    }
  }
}


/*
    GB1DMGridGGAFn
*/

void GB1DMGridGGAFn::compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                           double* output, double epsilon,
                                           double* dmmaxrow) {
  // The value of the basis function and its derivatives toward x, y and z (at `point`)
  // are available in work_basis. These are used, together with the 1RDM coefficients,
  // to compute the density and its gradient at `point`.
  double rho = 0, rho_x = 0, rho_y = 0, rho_z = 0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      row += work_basis[ibasis1*4]*dm[ibasis0*nbasis+ibasis1];
    }
    rho += row*work_basis[ibasis0*4];
    rho_x += row*work_basis[ibasis0*4+1];
    rho_y += row*work_basis[ibasis0*4+2];
    rho_z += row*work_basis[ibasis0*4+3];
  }
  output[0] += rho;
  output[1] += 2*rho_x;
  output[2] += 2*rho_y;
  output[3] += 2*rho_z;
}

void GB1DMGridGGAFn::compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                           double* fock) {
  // The functional derivative of the energy w.r.t. to the density and its gradient are
  // given inthe argument *pot (four components). The chain rule is used to transform
  // these into a contribution to Fock matrix.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp0 = pot[0]*work_basis[ibasis0*4] +
                  pot[1]*work_basis[ibasis0*4+1] +
                  pot[2]*work_basis[ibasis0*4+2] +
                  pot[3]*work_basis[ibasis0*4+3];
    double tmp1 = pot[1]*work_basis[ibasis0*4];
    double tmp2 = pot[2]*work_basis[ibasis0*4];
    double tmp3 = pot[3]*work_basis[ibasis0*4];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double result = tmp0*work_basis[ibasis1*4] +
                      tmp1*work_basis[ibasis1*4+1] +
                      tmp2*work_basis[ibasis1*4+2] +
                      tmp3*work_basis[ibasis1*4+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      if (ibasis1 != ibasis0) {
        // Enforce symmetry
        fock[ibasis0*nbasis+ibasis1] += result;
      }
    }
  }
}


/*
    GB1DMGridKineticFn
*/

void GB1DMGridKineticFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h1 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+1)+1;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
}

void GB1DMGridKineticFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // For every primitive, work_cart contains three values computed at `point` (see reset
  // method): the derivatives of the basis function toward x, y and z.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function derived toward x
    work_cart[3*i1p.ibasis0] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[3*i1p.ibasis0] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
    work_cart[3*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[3*i1p.ibasis0+1] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[3*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[3*i1p.ibasis0+2] += i1p.n0[2]*pre0*poly_work[i1p.ibasis0-nnotx-1+offset_l1];
  } while (i1p.inc());
}

void GB1DMGridKineticFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The derivatives of the basis function w.r.t. x, y and z, at `point` (see reset
  // method) are given in work_basis. Together with the 1RDM coefficients, dm, these are
  // used to compute the kinetic energy density in `point`.
  double tau = 0.0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;
    for (long ibasis1=0; ibasis1 < ibasis0; ibasis1++) {
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

void GB1DMGridKineticFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the kinetic energy density in `point`, is given
  // in *pot. This is used to compute a contribution to the Fock matrix from `point`.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp_x = 0.5*(*pot)*work_basis[ibasis0*3  ];
    double tmp_y = 0.5*(*pot)*work_basis[ibasis0*3+1];
    double tmp_z = 0.5*(*pot)*work_basis[ibasis0*3+2];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double result = tmp_x*work_basis[ibasis1*3  ] +
                      tmp_y*work_basis[ibasis1*3+1] +
                      tmp_z*work_basis[ibasis1*3+2];
      fock[ibasis1*nbasis+ibasis0] += result;
      if (ibasis0 != ibasis1)
        fock[ibasis0*nbasis+ibasis1] += result;
    }
  }
}


/*
    GB1DMGridHessianFn
*/

void GB1DMGridHessianFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h2 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+2)+1;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
}

void GB1DMGridHessianFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The value of the basis functions in the Cartesian primitive shell are added to
  // work_basis, together with its first and second derivatives, all evaluated at `point`.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    double pre0_hh = -pre0_h*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[10*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[10*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[10*i1p.ibasis0+1] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
    work_cart[10*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[10*i1p.ibasis0+2] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[10*i1p.ibasis0+3] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[10*i1p.ibasis0+3] += i1p.n0[2]*pre0*poly_work[i1p.ibasis0-nnotx-1+offset_l1];

    // Basis function derived toward xx
    work_cart[10*i1p.ibasis0+4] += pre0_hh*poly_work[i1p.ibasis0+offset_h2];
    work_cart[10*i1p.ibasis0+4] += (2*i1p.n0[0]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[0] > 1)
      work_cart[10*i1p.ibasis0+4] += i1p.n0[0]*(i1p.n0[0]-1)*pre0*
                                      poly_work[i1p.ibasis0+offset_l2];

    // Basis function derived toward xy
    work_cart[10*i1p.ibasis0+5] += pre0_hh*poly_work[i1p.ibasis0+1+nnotx+offset_h2];
    if (i1p.n0[0] > 0)
      work_cart[10*i1p.ibasis0+5] += i1p.n0[0]*pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset];
    if (i1p.n0[1] > 0)
      work_cart[10*i1p.ibasis0+5] += i1p.n0[1]*pre0_h*poly_work[i1p.ibasis0-nnotx+offset];
    if ((i1p.n0[0] > 0) && (i1p.n0[1] > 0))
      work_cart[10*i1p.ibasis0+5] += i1p.n0[0]*i1p.n0[1]*pre0*
                                     poly_work[i1p.ibasis0-nnotx+offset_l2];

    // Basis function derived toward xz
    work_cart[10*i1p.ibasis0+6] += pre0_hh*poly_work[i1p.ibasis0+2+nnotx+offset_h2];
    if (i1p.n0[0] > 0)
      work_cart[10*i1p.ibasis0+6] += i1p.n0[0]*pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset];
    if (i1p.n0[2] > 0)
      work_cart[10*i1p.ibasis0+6] += i1p.n0[2]*pre0_h*poly_work[i1p.ibasis0-1-nnotx+offset];
    if ((i1p.n0[0] > 0) && (i1p.n0[2] > 0))
      work_cart[10*i1p.ibasis0+6] += i1p.n0[0]*i1p.n0[2]*pre0*
                                     poly_work[i1p.ibasis0-1-nnotx+offset_l2];

    // Basis function derived toward yy
    work_cart[10*i1p.ibasis0+7] += pre0_hh*poly_work[i1p.ibasis0+3+2*nnotx+offset_h2];
    work_cart[10*i1p.ibasis0+7] += (2*i1p.n0[1]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[1] > 1)
      work_cart[10*i1p.ibasis0+7] += i1p.n0[1]*(i1p.n0[1]-1)*pre0*
                                     poly_work[i1p.ibasis0+1-2*nnotx+offset_l2];

    // Basis function derived toward yz
    work_cart[10*i1p.ibasis0+8] += pre0_hh*poly_work[i1p.ibasis0+4+2*nnotx+offset_h2];
    if (i1p.n0[1] > 0)
      work_cart[10*i1p.ibasis0+8] += i1p.n0[1]*pre0_h*poly_work[i1p.ibasis0+1+offset];
    if (i1p.n0[2] > 0)
      work_cart[10*i1p.ibasis0+8] += i1p.n0[2]*pre0_h*poly_work[i1p.ibasis0-1+offset];
    if ((i1p.n0[1] > 0) && (i1p.n0[2] > 0))
      work_cart[10*i1p.ibasis0+8] += i1p.n0[1]*i1p.n0[2]*pre0*
                                     poly_work[i1p.ibasis0-2*nnotx+offset_l2];

    // Basis function derived toward zz
    work_cart[10*i1p.ibasis0+9] += pre0_hh*poly_work[i1p.ibasis0+5+2*nnotx+offset_h2];
    work_cart[10*i1p.ibasis0+9] += (2*i1p.n0[2]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[2] > 1)
      work_cart[10*i1p.ibasis0+9] += i1p.n0[2]*(i1p.n0[2]-1)*pre0*
                                     poly_work[i1p.ibasis0-1-2*nnotx+offset_l2];
  } while (i1p.inc());
}

void GB1DMGridHessianFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The density Hessian is computed in `point` using the results in work_basis (basis
  // function values and their first and second derivatives).
  double rho_xx = 0, rho_xy = 0, rho_xz = 0;
  double rho_yy = 0, rho_yz = 0, rho_zz = 0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      row += work_basis[ibasis1*10]*dm[ibasis0*nbasis+ibasis1];
      tmp_x += work_basis[ibasis1*10+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*10+2]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*10+3]*dm[ibasis0*nbasis+ibasis1];
    }
    rho_xx += row*work_basis[ibasis0*10+4] + tmp_x*work_basis[ibasis0*10+1];
    rho_xy += row*work_basis[ibasis0*10+5] + tmp_x*work_basis[ibasis0*10+2];
    rho_xz += row*work_basis[ibasis0*10+6] + tmp_x*work_basis[ibasis0*10+3];
    rho_yy += row*work_basis[ibasis0*10+7] + tmp_y*work_basis[ibasis0*10+2];
    rho_yz += row*work_basis[ibasis0*10+8] + tmp_y*work_basis[ibasis0*10+3];
    rho_zz += row*work_basis[ibasis0*10+9] + tmp_z*work_basis[ibasis0*10+3];
  }
  output[0] += 2*rho_xx;
  output[1] += 2*rho_xy;
  output[2] += 2*rho_xz;
  output[3] += 2*rho_yy;
  output[4] += 2*rho_yz;
  output[5] += 2*rho_zz;
}

void GB1DMGridHessianFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the density Hessian matrix elements, evaluated in
  // `point` is given in `*pot` (six elements). These are transformed with the chain rule
  // to a contribution to the Fock matrix from `point`.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp = pot[0]*work_basis[ibasis0*10+4]
                +pot[1]*work_basis[ibasis0*10+5]
                +pot[2]*work_basis[ibasis0*10+6]
                +pot[3]*work_basis[ibasis0*10+7]
                +pot[4]*work_basis[ibasis0*10+8]
                +pot[5]*work_basis[ibasis0*10+9];
    double tmp_x = pot[0]*work_basis[ibasis0*10+1]
                  +pot[1]*work_basis[ibasis0*10+2]
                  +pot[2]*work_basis[ibasis0*10+3];
    double tmp_y = pot[1]*work_basis[ibasis0*10+1]
                  +pot[3]*work_basis[ibasis0*10+2]
                  +pot[4]*work_basis[ibasis0*10+3];
    double tmp_z = pot[2]*work_basis[ibasis0*10+1]
                  +pot[4]*work_basis[ibasis0*10+2]
                  +pot[5]*work_basis[ibasis0*10+3];
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      double result = tmp*work_basis[ibasis1*10]
                     +tmp_x*work_basis[ibasis1*10+1]
                     +tmp_y*work_basis[ibasis1*10+2]
                     +tmp_z*work_basis[ibasis1*10+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      fock[ibasis0*nbasis+ibasis1] += result;
    }
  }
}


/*
    GB1DMGridMGGAFn
*/

void GB1DMGridMGGAFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h2 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+2)+1;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
}

void GB1DMGridMGGAFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The primitive basis functions in one Cartesian shell are added to work_cart,
  // including the first derivatives and the Laplacian. (Five elements in total per basis
  // function.)
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    double pre0_hh = -pre0_h*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[5*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[5*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[5*i1p.ibasis0+1] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
    work_cart[5*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[5*i1p.ibasis0+2] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[5*i1p.ibasis0+3] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[5*i1p.ibasis0+3] += i1p.n0[2]*pre0*
                                    poly_work[i1p.ibasis0-nnotx-1+offset_l1];

    // Laplacian of the Basis function
    work_cart[5*i1p.ibasis0+4] += pre0_hh*poly_work[i1p.ibasis0+offset_h2];
    work_cart[5*i1p.ibasis0+4] += (2*i1p.n0[0]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[0] > 1)
      work_cart[5*i1p.ibasis0+4] += i1p.n0[0]*(i1p.n0[0]-1)*
                                    pre0*poly_work[i1p.ibasis0+offset_l2];
    work_cart[5*i1p.ibasis0+4] += pre0_hh*poly_work[i1p.ibasis0+3+2*nnotx+offset_h2];
    work_cart[5*i1p.ibasis0+4] += (2*i1p.n0[1]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[1] > 1)
      work_cart[5*i1p.ibasis0+4] += i1p.n0[1]*(i1p.n0[1]-1)*pre0*
                                    poly_work[i1p.ibasis0+1-2*nnotx+offset_l2];
    work_cart[5*i1p.ibasis0+4] += pre0_hh*poly_work[i1p.ibasis0+5+2*nnotx+offset_h2];
    work_cart[5*i1p.ibasis0+4] += (2*i1p.n0[2]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[2] > 1)
      work_cart[5*i1p.ibasis0+4] += i1p.n0[2]*(i1p.n0[2]-1)*pre0*
                                    poly_work[i1p.ibasis0-1-2*nnotx+offset_l2];
  } while (i1p.inc());
}

void GB1DMGridMGGAFn::compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                            double* output, double epsilon,
                                            double* dmmaxrow) {
  // The density, its gradient, the kinetic energy density and the Laplacian, are computed
  // in `point` (see reset method), using the results in work_basis and the 1RDM
  // coefficients in dm.
  double rho = 0, rho_x = 0, rho_y = 0, rho_z = 0;
  double lapl_part = 0, tau2 = 0.0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      row += work_basis[ibasis1*5]*dm[ibasis0*nbasis+ibasis1];
      tmp_x += work_basis[ibasis1*5+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*5+2]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*5+3]*dm[ibasis0*nbasis+ibasis1];
    }
    rho += row*work_basis[ibasis0*5];
    rho_x += row*work_basis[ibasis0*5+1];
    rho_y += row*work_basis[ibasis0*5+2];
    rho_z += row*work_basis[ibasis0*5+3];
    lapl_part += row*work_basis[ibasis0*5+4];
    tau2 += tmp_x*work_basis[ibasis0*5+1] +
            tmp_y*work_basis[ibasis0*5+2] +
            tmp_z*work_basis[ibasis0*5+3];
  }
  output[0] += rho;
  output[1] += 2*rho_x;
  output[2] += 2*rho_y;
  output[3] += 2*rho_z;
  output[4] += 2*lapl_part + 2*tau2;
  output[5] += 0.5*tau2;
}

void GB1DMGridMGGAFn::compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                            double* fock) {
  // The functional derivative of the energy w.r.t. density, its gradient, the kinetic
  // energy density and the Laplacian, are given in *pot. These are transformed to a
  // contribution to the Fock matrix using the chain rule.
  double auxpot = 0.5*pot[5] + 2.0*pot[4];
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp0 = pot[0]*work_basis[ibasis0*5] +
                  pot[1]*work_basis[ibasis0*5+1] +
                  pot[2]*work_basis[ibasis0*5+2] +
                  pot[3]*work_basis[ibasis0*5+3] +
                  pot[4]*work_basis[ibasis0*5+4];
    double tmp1 = pot[1]*work_basis[ibasis0*5] + auxpot*work_basis[ibasis0*5+1];
    double tmp2 = pot[2]*work_basis[ibasis0*5] + auxpot*work_basis[ibasis0*5+2];
    double tmp3 = pot[3]*work_basis[ibasis0*5] + auxpot*work_basis[ibasis0*5+3];
    double tmp4 = pot[4]*work_basis[ibasis0*5];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double result = tmp0*work_basis[ibasis1*5] +
                      tmp1*work_basis[ibasis1*5+1] +
                      tmp2*work_basis[ibasis1*5+2] +
                      tmp3*work_basis[ibasis1*5+3] +
                      tmp4*work_basis[ibasis1*5+4];
      fock[ibasis1*nbasis+ibasis0] += result;
      if (ibasis1 != ibasis0) {
        // Enforce symmetry
        fock[ibasis0*nbasis+ibasis1] += result;
      }
    }
  }
}


/*
    GB2DMGridFn
*/

GB2DMGridFn::GB2DMGridFn(long max_shell_type)
    : GBCalculator(max_shell_type), shell_type0(0), shell_type1(0), r0(NULL), r1(NULL),
      point(NULL), i2p() {
  nwork = max_nbasis*max_nbasis;
  work_cart = new double[nwork];
  work_pure = new double[nwork];
}

void GB2DMGridFn::reset(long _shell_type0, long _shell_type1, const double* _r0,
                        const double* _r1, const double* _point) {
  if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type)) {
    throw std::domain_error("shell_type0 out of range.");
  }
  if ((_shell_type1 < -max_shell_type) || (_shell_type1 > max_shell_type)) {
    throw std::domain_error("shell_type0 out of range.");
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

  // For now, this is just a copy from GB2Integral.
  // (This must be changed when electrical fields are implemented.)

  // Project along index 0 (rows)
  if (shell_type0 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type0,
      1,                                     // anterior
      get_shell_nbasis(abs(shell_type1)) );  // posterior
    swap_work();
  }

  // Project along index 1 (cols)
  if (shell_type1 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type1,
      get_shell_nbasis(shell_type0),  // anterior
      1);                             // posterior
    swap_work();
  }
}


/*
    GB2DMGridHartreeFn
*/

GB2DMGridHartreeFn::GB2DMGridHartreeFn(long max_shell_type) : GB2DMGridFn(max_shell_type) {
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


void GB2DMGridHartreeFn::add(double coeff, double alpha0, double alpha1,
                            const double* scales0, const double* scales1) {
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
  for (long nu=abs(shell_type0)+abs(shell_type1); nu >= 0; nu--) {
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
    for (long i0 = i2p.n0[0]+i2p.n1[0]; i0 >= 0; i0--)
      for (long i1 = i2p.n0[1]+i2p.n1[1]; i1 >= 0; i1--)
        for (long i2 = i2p.n0[2]+i2p.n1[2]; i2 >= 0; i2--)
          arg += work_g0[i0]*work_g1[i1]*work_g2[i2]*work_boys[i0+i1+i2];

    // Finally add to the work array
    work_cart[i2p.offset] += pre*scales0[i2p.ibasis0]*scales1[i2p.ibasis1]*arg;
  } while (i2p.inc());
}
