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
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include "horton/gbasis/gbasis.h"
#include "horton/gbasis/common.h"
#include "horton/gbasis/iter_gb.h"
using std::abs;

/*

  Auxiliary routines

*/

const double gob_cart_normalization(const double alpha, const long* n) {
    return sqrt(pow(4.0*alpha, n[0]+n[1]+n[2])*pow(2.0*alpha/M_PI, 1.5)
           /(fac2(2*n[0]-1)*fac2(2*n[1]-1)*fac2(2*n[2]-1)));
}


const double gob_pure_normalization(const double alpha, const long l) {
    return sqrt(pow(4.0*alpha, l)*pow(2.0*alpha/M_PI, 1.5)
           /fac2(2*l-1));
}

/*
    GBasis

    The constructor assumes that the arrays are allocated externally and will
    not be deallocated prematurely. It is also assumed that the arrays are not
    changed once the constructor is called.
*/

GBasis::GBasis(const double* centers, const long* shell_map, const long* nprims,
               const long* shell_types, const double* alphas, const double* con_coeffs,
               const long ncenter, const long nshell, const long nprim_total) :
    nbasis(0), nscales(0), max_shell_type(0),
    centers(centers), shell_map(shell_map), nprims(nprims),
    shell_types(shell_types), alphas(alphas), con_coeffs(con_coeffs),
    ncenter(ncenter), nshell(nshell), nprim_total(nprim_total)
{
    long shell_nbasis, shell_type;

    // check for maximum shell type
    for (long ishell=0; ishell < nshell; ishell++) {
        shell_type = abs(shell_types[ishell]);
        if (shell_type > MAX_SHELL_TYPE) {
            throw std::domain_error("Exceeded the maximum shell type.");
        }
        if (max_shell_type < shell_type) {
            max_shell_type = shell_type;
        }
    }

    // basis_offsets: the index of the first basis function for every shell of
    // contracted Gaussians.
    basis_offsets = new long[nshell];
    basis_offsets[0] = 0;
    for (long ishell=1; ishell < nshell; ishell++) {
        shell_nbasis = get_shell_nbasis(shell_types[ishell-1]);
        basis_offsets[ishell] = basis_offsets[ishell-1] + shell_nbasis;
    }
    shell_nbasis = get_shell_nbasis(shell_types[nshell-1]);
    nbasis = basis_offsets[nshell-1] + shell_nbasis;

    // prim_offsets: the index of the first primitive for every shell of
    // contracted Gaussians.
    prim_offsets = new long[nshell];
    prim_offsets[0] = 0;
    for (long ishell=1; ishell < nshell; ishell++) {
        prim_offsets[ishell] = prim_offsets[ishell-1] + nprims[ishell-1];
    }

    // shell_lookup: the index of the contracted shell of Gaussians for every
    // basis function.
    shell_lookup = new long[nbasis];
    long ishell = 0;
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
        shell_lookup[ibasis] = ishell;
        if ((ishell < nshell-1) && (ibasis+1 == basis_offsets[ishell+1])) {
            ishell++;
        }
    }

    // nscales
    for (long ishell=0; ishell < nshell; ishell++) {
        shell_nbasis = get_shell_nbasis(abs(shell_types[ishell]));
        nscales += shell_nbasis*nprims[ishell];
    }

    // scales
    scales = new double[nscales];
    scales_offsets = new long[nprim_total];
}

GBasis::~GBasis() {
    delete[] basis_offsets;
    delete[] prim_offsets;
    delete[] shell_lookup;
    delete[] scales;
    delete[] scales_offsets;
}

void GBasis::init_scales() {
    long n[3], counter=0, oprim=0;
    double alpha;
    for (long ishell=0; ishell < nshell; ishell++) {
        for (long iprim=0; iprim < nprims[ishell]; iprim++) {
            scales_offsets[oprim + iprim] = counter;
            alpha = alphas[oprim + iprim];
            n[0] = abs(shell_types[ishell]);
            n[1] = 0;
            n[2] = 0;
            do {
                scales[counter] = normalization(alpha, n);
                counter += 1;
            } while (iter_pow1_inc(n));
        }
        oprim += nprims[ishell];
    }
}

void GBasis::compute_two_index(double* output, GB2Integral* integral) {
    IterGB2 iter = IterGB2(this);
    iter.update_shell();
    do {
        integral->reset(iter.shell_type0, iter.shell_type1, iter.r0, iter.r1);
        iter.update_prim();
        do {
            integral->add(iter.con_coeff, iter.alpha0, iter.alpha1, iter.scales0, iter.scales1);
        } while (iter.inc_prim());
        integral->cart_to_pure();
        iter.store(integral->get_work(), output);
    } while (iter.inc_shell());
}

void GBasis::compute_four_index(double* output, GB4Integral* integral) {
    IterGB4 iter = IterGB4(this);
    iter.update_shell();
    do {
        integral->reset(iter.shell_type0, iter.shell_type1, iter.shell_type2, iter.shell_type3,
                        iter.r0, iter.r1, iter.r2, iter.r3);
        iter.update_prim();
        do {
            integral->add(iter.con_coeff, iter.alpha0, iter.alpha1, iter.alpha2, iter.alpha3,
                          iter.scales0, iter.scales1, iter.scales2, iter.scales3);
        } while (iter.inc_prim());
        integral->cart_to_pure();
        iter.store(integral->get_work(), output);
    } while (iter.inc_shell());
}

void GBasis::compute_grid_point1(double* output, double* point, GB1GridFn* grid_fn) {
    IterGB1 iter = IterGB1(this);
    iter.update_shell();
    do {
        grid_fn->reset(iter.shell_type0, iter.r0, point);
        iter.update_prim();
        do {
            grid_fn->add(iter.con_coeff, iter.alpha0, iter.scales0);
        } while (iter.inc_prim());
        grid_fn->cart_to_pure();
        iter.store(grid_fn->get_work(), output, grid_fn->get_dim_work());
    } while (iter.inc_shell());
}

double GBasis::compute_grid_point2(double* dm, double* point, GB2DMGridFn* grid_fn) {
    double result = 0.0;
    IterGB2 iter = IterGB2(this);
    iter.update_shell();
    do {
        grid_fn->reset(iter.shell_type0, iter.shell_type1, iter.r0, iter.r1, point);
        iter.update_prim();
        do {
            grid_fn->add(iter.con_coeff, iter.alpha0, iter.alpha1, iter.scales0, iter.scales1);
        } while (iter.inc_prim());
        grid_fn->cart_to_pure();
        result += iter.dot(grid_fn->get_work(), dm);
    } while (iter.inc_shell());
    return result;
}

GOBasis::GOBasis(const double* centers, const long* shell_map, const long* nprims,
                 const long* shell_types, const double* alphas, const double* con_coeffs,
                 const long ncenter, const long nshell, const long nprim_total) :
    GBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs,
    ncenter, nshell, nprim_total) {
    init_scales();
}

const double GOBasis::normalization(const double alpha, const long* n) const {
    return gob_cart_normalization(alpha, n);
}

void GOBasis::compute_overlap(double* output) {
    GB2OverlapIntegral integral = GB2OverlapIntegral(get_max_shell_type());
    compute_two_index(output, &integral);
}

void GOBasis::compute_kinetic(double* output) {
    GB2KineticIntegral integral = GB2KineticIntegral(get_max_shell_type());
    compute_two_index(output, &integral);
}

void GOBasis::compute_nuclear_attraction(double* charges, double* centers, long ncharge,
                                         double* output) {
    GB2NuclearAttractionIntegral integral = GB2NuclearAttractionIntegral(get_max_shell_type(),
                                            charges, centers, ncharge);
    compute_two_index(output, &integral);
}

void GOBasis::compute_erf_attraction(double* charges, double* centers, long ncharge, double* output,
                                     double mu) {
    GB2ErfAttractionIntegral integral = GB2ErfAttractionIntegral(get_max_shell_type(),
                                        charges, centers, ncharge, mu);
    compute_two_index(output, &integral);
}

void GOBasis::compute_gauss_attraction(double* charges, double* centers, long ncharge,
                                       double* output, double c, double alpha) {
    GB2GaussAttractionIntegral integral = GB2GaussAttractionIntegral(get_max_shell_type(),
                                          charges, centers, ncharge, c, alpha);
    compute_two_index(output, &integral);
}

void GOBasis::compute_multipole_moment(long* xyz, double* center, double* output) {
    GB2MomentIntegral integral = GB2MomentIntegral(get_max_shell_type(), xyz, center);
    compute_two_index(output, &integral);
}

void GOBasis::compute_electron_repulsion(double* output) {
  GB4ElectronRepulsionIntegralLibInt integral =
    GB4ElectronRepulsionIntegralLibInt(get_max_shell_type());
  compute_four_index(output, &integral);
}

void GOBasis::compute_erf_repulsion(double* output, double mu) {
    GB4ErfIntegralLibInt integral = GB4ErfIntegralLibInt(get_max_shell_type(), mu);
    compute_four_index(output, &integral);
}

void GOBasis::compute_gauss_repulsion(double* output, double c, double alpha) {
    GB4GaussIntegralLibInt integral = GB4GaussIntegralLibInt(get_max_shell_type(), c, alpha);
    compute_four_index(output, &integral);
}

void GOBasis::compute_ralpha_repulsion(double* output, double alpha) {
    GB4RAlphaIntegralLibInt integral = GB4RAlphaIntegralLibInt(get_max_shell_type(), alpha);
    compute_four_index(output, &integral);
}

void GOBasis::compute_grid1_exp(long nfn, double* coeffs, long npoint, double* points,
                                long norb, long* iorbs, double* output) {
    // The work array contains the basis functions evaluated at the grid point,
    // and optionally some of its derivatives.
    GB1ExpGridOrbitalFn grid_fn = GB1ExpGridOrbitalFn(get_max_shell_type(), nfn, iorbs, norb);

    long nwork = get_nbasis()*grid_fn.get_dim_work();
    long dim_output = grid_fn.get_dim_output();
    double* work_basis = new double[nwork];

    for (long ipoint=0; ipoint < npoint; ipoint++) {
        // A) clear the basis functions.
        memset(work_basis, 0, nwork*sizeof(double));

        // B) evaluate the basis functions in the current point.
        compute_grid_point1(work_basis, points, &grid_fn);

        // C) Use the basis function results and the density matrix to evaluate
        // the function at the grid point. The result is added to the output.
        grid_fn.compute_point_from_exp(work_basis, coeffs, get_nbasis(), output);

        // D) Prepare for next iteration
        output += dim_output;
        points += 3;
    }

    delete[] work_basis;
}

void GOBasis::compute_grid1_grad_exp(long nfn, double* coeffs, long npoint,
                                     double* points, long norb, long* iorbs, double* output) {
    // The work array contains the basis functions evaluated at the grid point,
    // and optionally some of its derivatives.
    GB1ExpGridOrbGradientFn grid_fn = GB1ExpGridOrbGradientFn(get_max_shell_type(),
                                                              nfn, iorbs, norb);

    long nwork = get_nbasis()*grid_fn.get_dim_work();
    long dim_output = grid_fn.get_dim_output();
    double* work_basis = new double[nwork];

    for (long ipoint=0; ipoint < npoint; ipoint++) {
        // A) clear the basis functions.
        memset(work_basis, 0, nwork*sizeof(double));

        // B) evaluate the basis functions in the current point.
        compute_grid_point1(work_basis, points, &grid_fn);

        // C) Use the basis function results and the density matrix to evaluate
        // the function at the grid point. The result is added to the output.
        grid_fn.compute_point_from_exp(work_basis, coeffs, get_nbasis(), output);

        // D) Prepare for next iteration
        output += dim_output;
        points += 3;
    }

    delete[] work_basis;
}

void GOBasis::compute_grid1_dm(double* dm, long npoint, double* points,
                               GB1DMGridFn* grid_fn, double* output,
                               double epsilon, double* dmmaxrow) {
    // The work array contains the basis functions evaluated at the grid point,
    // and optionally some of its derivatives.
    long nwork = get_nbasis()*grid_fn->get_dim_work();
    long dim_output = grid_fn->get_dim_output();
    double* work_basis = new double[nwork];

    for (long ipoint=0; ipoint < npoint; ipoint++) {
        // A) clear the basis functions.
        memset(work_basis, 0, nwork*sizeof(double));

        // B) evaluate the basis functions in the current point.
        compute_grid_point1(work_basis, points, grid_fn);
#ifdef DEBUG
        for (int i=0; i<nwork; i++) printf("%f ", work_basis[i]);
        printf("\n");
#endif

        // C) Use the basis function results and the density matrix to evaluate
        // the function at the grid point. The result is added to the output.
        grid_fn->compute_point_from_dm(work_basis, dm, get_nbasis(), output, epsilon, dmmaxrow);

        // D) Prepare for next iteration
        output += dim_output;
        points += 3;
    }

    delete[] work_basis;
}

void GOBasis::compute_grid2_dm(double* dm, long npoint, double* points, double* output) {
    // For the moment, it is only possible to compute the Hartree potential on
    // a grid with this routine. Generalizations with electrical field and
    // other things are for later.
    GB2DMGridHartreeFn grid_fn = GB2DMGridHartreeFn(get_max_shell_type());

    for (long ipoint=0; ipoint < npoint; ipoint++) {
        *output += compute_grid_point2(dm, points, &grid_fn);
        output++;
        points += 3;
    }
}

void GOBasis::compute_grid1_fock(long npoint, double* points, double* weights, long pot_stride, double* pots, GB1DMGridFn* grid_fn, double* output) {
    // The work array contains the basis functions evaluated at the grid point,
    // and optionally some of its derivatives.
    long nwork = get_nbasis()*grid_fn->get_dim_work();
    double* work_basis = new double[nwork];
    long dim_output = grid_fn->get_dim_output();
    double* work_pot = new double[dim_output];

    for (long ipoint=0; ipoint < npoint; ipoint++) {
        // A) clear the work array.
        memset(work_basis, 0, nwork*sizeof(double));

        // B) evaluate the basis functions in the current point.
        compute_grid_point1(work_basis, points, grid_fn);

        // C) Add the contribution from this grid point to the operator
        for (long i=dim_output-1; i >= 0; i--) {
            work_pot[i] = (*weights)*pots[i];
        }
        grid_fn->compute_fock_from_pot(work_pot, work_basis, get_nbasis(), output);

        // D) Prepare for next iteration
        points += 3;
        weights++;
        pots += pot_stride;
    }

    delete[] work_basis;
    delete[] work_pot;
}
