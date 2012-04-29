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
#include <stdexcept>
#include "gbasis.h"
#include "common.h"
#include "iter_gb.h"
using namespace std;

/*
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
