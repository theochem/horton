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
#include <cstdio>
#include <cmath>
#include <stdexcept>
#include "gbasis.h"
#include "common.h"
#include "iter_gb.h"
using namespace std;

/*

  Auxiliary routines

*/

const double gob_normalization(const double alpha, const long* n) {
    return sqrt(pow(4.0*alpha, n[0]+n[1]+n[2])*pow(2.0*alpha/M_PI, 1.5)
           /(fac2(2*n[0]-1)*fac2(2*n[1]-1)*fac2(2*n[2]-1)));
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
    for (long ishell=0; ishell<nshell; ishell++) {
        shell_type = abs(shell_types[ishell]);
        if (shell_type > MAX_SHELL_TYPE) {
            throw std::domain_error("Exceeded the maximum shell type.");
        }
        if (max_shell_type < shell_type) {
            max_shell_type = shell_type;
        }
    }

    // basis_offsets
    basis_offsets = new long[nshell];
    basis_offsets[0] = 0;
    for (long ishell=1; ishell<nshell; ishell++) {
        shell_nbasis = get_shell_nbasis(shell_types[ishell-1]);
        basis_offsets[ishell] = basis_offsets[ishell-1] + shell_nbasis;
    }
    shell_nbasis = get_shell_nbasis(shell_types[nshell-1]);
    nbasis = basis_offsets[nshell-1] + shell_nbasis;

    // nscales
    for (long ishell=0; ishell<nshell; ishell++) {
        shell_nbasis = get_shell_nbasis(abs(shell_types[ishell]));
        nscales += shell_nbasis*nprims[ishell];
    }

    // scales
    scales = new double[nscales];
    scales_offsets = new long[nprim_total];
}

GBasis::~GBasis() {
    delete[] basis_offsets;
    delete[] scales;
    delete[] scales_offsets;
}

void GBasis::init_scales() {
    long n[3], counter=0, oprim=0;
    double alpha;
    for (long ishell=0; ishell<nshell; ishell++) {
        for (long iprim=0; iprim<nprims[ishell]; iprim++) {
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

void GBasis::compute_one_body(double* output, GB2Integral* integral) {
    /*
        TODO
             When multiple different memory storage schemes are implemented for
             the operators, the iterator must also become an argument for this
             function

    */
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

void GBasis::compute_two_body(double* output, GB4Integral* integral) {
    /*
        TODO
             When multiple different memory storage schemes are implemented for
             the operators, the iterator must also become an argument for this
             function

    */
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

void GBasis::compute_grid(double* output, double* point, GB1GridFn* grid_fn) {
    /*
        TODO
             When multiple different memory storage schemes are implemented for
             the operators, the iterator must also become an argument for this
             function

    */
    IterGB1 iter = IterGB1(this);
    iter.update_shell();
    do {
        grid_fn->reset(iter.shell_type0, iter.r0, point);
        iter.update_prim();
        do {
            grid_fn->add(iter.con_coeff, iter.alpha0, iter.scales0);
        } while (iter.inc_prim());
        grid_fn->cart_to_pure();
        iter.store(grid_fn->get_work(), output);
    } while (iter.inc_shell());
}



GOBasis::GOBasis(const double* centers, const long* shell_map, const long* nprims,
                 const long* shell_types, const double* alphas, const double* con_coeffs,
                 const long ncenter, const long nshell, const long nprim_total) :
    GBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs,
    ncenter, nshell, nprim_total) {
    init_scales();
}

const double GOBasis::normalization(const double alpha, const long* n) const {
    return gob_normalization(alpha, n);
}

void GOBasis::compute_overlap(double* output) {
    GB2OverlapIntegral integral = GB2OverlapIntegral(get_max_shell_type());
    compute_one_body(output, &integral);
}

void GOBasis::compute_kinetic(double* output) {
    GB2KineticIntegral integral = GB2KineticIntegral(get_max_shell_type());
    compute_one_body(output, &integral);
}

void GOBasis::compute_nuclear_attraction(double* charges, double* centers, long ncharge, double* output) {
    GB2NuclearAttractionIntegral integral = GB2NuclearAttractionIntegral(get_max_shell_type(), charges, centers, ncharge);
    compute_one_body(output, &integral);
}

void GOBasis::compute_electron_repulsion(double* output) {
    GB4ElectronReuplsionIntegralLibInt integral = GB4ElectronReuplsionIntegralLibInt(get_max_shell_type());
    compute_two_body(output, &integral);
}

void GOBasis::compute_density_grid_dm(double* dm, long npoint, double* points, double* rhos) {
    double basis_fns[get_nbasis()];
    GB1GridFn grid_fn = GB1GridFn(get_max_shell_type());

    for (long ipoint=0; ipoint<npoint; ipoint++) {

        // A) clear the basis functions.
        for (long ibasis=0; ibasis<get_nbasis(); ibasis++) {
            basis_fns[ibasis] = 0.0;
        }

        // B) evaluate the basis functions in the current point.
        compute_grid(basis_fns, points, &grid_fn);
#ifdef DEBUG
        for (long ibasis=0; ibasis<get_nbasis(); ibasis++) {
            printf("basis_fns[%i] = %f\n", ibasis, basis_fns[ibasis]);
        }
        printf("\n");
#endif

        // C) Make dot product of basis functions with density matrix. This is the
        // bottle neck. Note that the result is added to the output array!
        double rho = 0, row;
        for (long ibasis0=0; ibasis0<get_nbasis(); ibasis0++) {
            row = 0;
            for (long ibasis1=0; ibasis1<get_nbasis(); ibasis1++) {
                row += basis_fns[ibasis1]*dm[ibasis0*get_nbasis()+ibasis1];
            }
            rho += row*basis_fns[ibasis0];
        }
        *rhos += rho;

        // D) Prepare for next iteration
        rhos++;
        points += 3;
    }

}

void GOBasis::compute_density_grid_orb(double* orbs, long nocc, long npoint, double* points, double* rhos) {
    // TODO
}
