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


#include <cstdlib>
#include <cstring>
#include "horton/gbasis/common.h"
#include "horton/gbasis/iter_gb.h"
using namespace std;



IterGB1::IterGB1(GBasis* gbasis) :
    gbasis(gbasis),
    basis_offsets(gbasis->get_basis_offsets()),
    // public fields
    shell_type0(0),
    con_coeff(0.0),
    alpha0(0.0),
    r0(NULL),
    ibasis0(0),
    // internal fields
    ishell0(0),
    nprim0(0), iprim0(0), oprim0(0)
{}


int IterGB1::inc_shell() {
    // Increment shell and related counters.
    if (ishell0 < gbasis->nshell-1) {
        oprim0 += nprim0;
        ishell0++;
        update_shell();
        return 1;
    } else {
        ishell0 = 0;
        oprim0 = 0;
        update_shell();
        return 0;
    }
}


void IterGB1::update_shell() {
    // Update fields that depend on shell and related counters.
    nprim0 = gbasis->nprims[ishell0];
    // Update indexes in output array
    ibasis0 = basis_offsets[ishell0];
    // update centers
    r0 = gbasis->centers + 3*gbasis->shell_map[ishell0];
    // update shell types
    shell_type0 = gbasis->shell_types[ishell0];
    // reset contraction counters
    iprim0 = 0;
}


int IterGB1::inc_prim() {
    // Increment primitive counters.
    if (iprim0 < nprim0-1) {
        iprim0++;
        update_prim();
        return 1;
    } else {
        iprim0 = 0;
        update_prim();
        return 0;
    }
}


void IterGB1::update_prim() {
    // Update fields that depend on primitive counters.
    alpha0 = gbasis->alphas[oprim0 + iprim0];
    con_coeff = gbasis->con_coeffs[oprim0 + iprim0];
    scales0 = gbasis->get_scales(oprim0 + iprim0);
}


void IterGB1::store(const double *work, double *output, long dim) {
    // This routine is hardwired to work only for the dense storage
    const long n0 = get_shell_nbasis(shell_type0);
    memcpy(output+ibasis0*dim, work, n0*dim*sizeof(double));
}






IterGB2::IterGB2(GBasis* gbasis) :
    gbasis(gbasis),
    basis_offsets(gbasis->get_basis_offsets()),
    // public fields
    shell_type0(0), shell_type1(0),
    con_coeff(0.0),
    alpha0(0.0), alpha1(0.0),
    r0(NULL), r1(NULL),
    ibasis0(0), ibasis1(0),
    // internal fields
    ishell0(0), ishell1(0),
    nprim0(0), nprim1(0), iprim0(0), iprim1(0), oprim0(0), oprim1(0)
{}


int IterGB2::inc_shell() {
    // Increment shell and related counters. This loop takes into account the
    // symmetry of the two-index operator
    if (ishell1 < ishell0) {
        ishell1++;
        oprim1 += nprim1;
        update_shell();
        return 1;
    } else if (ishell0 < gbasis->nshell-1) {
        ishell1 = 0;
        oprim1 = 0;
        oprim0 += nprim0;
        ishell0++;
        update_shell();
        return 1;
    } else {
        ishell1 = 0;
        ishell0 = 0;
        oprim0 = 0;
        oprim1 = 0;
        update_shell();
        return 0;
    }
}


void IterGB2::update_shell() {
    // Update fields that depend on shell and related counters.
    nprim0 = gbasis->nprims[ishell0];
    nprim1 = gbasis->nprims[ishell1];
    // Update indexes in output array
    ibasis0 = basis_offsets[ishell0];
    ibasis1 = basis_offsets[ishell1];
    // update centers
    r0 = gbasis->centers + 3*gbasis->shell_map[ishell0];
    r1 = gbasis->centers + 3*gbasis->shell_map[ishell1];
    // update shell types
    shell_type0 = gbasis->shell_types[ishell0];
    shell_type1 = gbasis->shell_types[ishell1];
    // reset contraction counters
    iprim0 = 0;
    iprim1 = 0;
}


int IterGB2::inc_prim() {
    // Increment primitive counters.
    if (iprim1 < nprim1-1) {
        iprim1++;
        update_prim();
        return 1;
    } else if (iprim0 < nprim0-1) {
        iprim0++;
        iprim1 = 0;
        update_prim();
        return 1;
    } else {
        iprim0 = 0;
        iprim1 = 0;
        update_prim();
        return 0;
    }
}


void IterGB2::update_prim() {
    // Update fields that depend on primitive counters.
    alpha0 = gbasis->alphas[oprim0 + iprim0];
    alpha1 = gbasis->alphas[oprim1 + iprim1];
    con_coeff = gbasis->con_coeffs[oprim0 + iprim0]*
                gbasis->con_coeffs[oprim1 + iprim1];
    scales0 = gbasis->get_scales(oprim0 + iprim0);
    scales1 = gbasis->get_scales(oprim1 + iprim1);
}


void IterGB2::store(const double *work, double *output) {
    // This routine is hardwired to work only for the dense storage
    long i0, i1;
    const long n0 = get_shell_nbasis(shell_type0);
    const long n1 = get_shell_nbasis(shell_type1);
    const long nbasis = gbasis->get_nbasis();
    const double* tmp;
    tmp = work;
    // diagonal or first off-diagonal block
    for (i0=0; i0<n0; i0++) {
        for (i1=0; i1<n1; i1++) {
            output[(i0+ibasis0)*nbasis+i1+ibasis1] = *tmp;
            tmp++;
        }
    }
    // optional second off-diagonal block
    if (ibasis0 != ibasis1) {
        tmp = work;
        for (i0=0; i0<n0; i0++) {
            for (i1=0; i1<n1; i1++) {
                output[(i1+ibasis1)*nbasis+i0+ibasis0] = *tmp;
                tmp++;
            }
        }
    }
}

double IterGB2::dot(const double *work, const double *dm) {
    // This routine is hardwired to work only for the dense storage
    long i0, i1;
    const long n0 = get_shell_nbasis(shell_type0);
    const long n1 = get_shell_nbasis(shell_type1);
    const long nbasis = gbasis->get_nbasis();
    const double* tmp;
    double result = 0.0;
    tmp = work;
    // diagonal or first off-diagonal block
    for (i0=0; i0<n0; i0++) {
        for (i1=0; i1<n1; i1++) {
            result += dm[(i0+ibasis0)*nbasis+i1+ibasis1] * (*tmp);
            tmp++;
        }
    }
    // optional second off-diagonal block
    if (ibasis0 != ibasis1) {
        tmp = work;
        for (i0=0; i0<n0; i0++) {
            for (i1=0; i1<n1; i1++) {
                result += dm[(i1+ibasis1)*nbasis+i0+ibasis0] * (*tmp);
                tmp++;
            }
        }
    }
    return result;
}





IterGB4::IterGB4(GBasis* gbasis) :
    gbasis(gbasis),
    basis_offsets(gbasis->get_basis_offsets()),
    // public fields
    shell_type0(0), shell_type1(0), shell_type2(0), shell_type3(0),
    con_coeff(0.0),
    alpha0(0.0), alpha1(0.0), alpha2(0.0), alpha3(0.0),
    r0(NULL), r1(NULL), r2(NULL), r3(NULL),
    ibasis0(0), ibasis1(0), ibasis2(0), ibasis3(0),
    // internal fields
    ishell0(0), ishell1(0), ishell2(0), ishell3(0),
    ishell3_max(0),
    nprim0(0), nprim1(0), nprim2(0), nprim3(0),
    iprim0(0), iprim1(0), iprim2(0), iprim3(0),
    oprim0(0), oprim1(0), oprim2(0), oprim3(0)
{}


int IterGB4::inc_shell() {
    // Increment shell and related counters. This loop takes into account the
    // symmetry of the four-index operator:
    //   <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
    //   <il|kj> = <jk|li> = <kj|il> = <li|jk>
    if (ishell0==ishell1) {
        ishell3_max = ishell2;
    } else {
        ishell3_max = ishell1;
    }
    if (ishell3 < ishell3_max) {
        ishell3++;
        oprim3 += nprim3;
        update_shell();
        return 1;
    } else if (ishell2 < ishell0) {
        ishell3 = 0;
        oprim3 = 0;
        ishell2++;
        oprim2 += nprim2;
        update_shell();
        return 1;
    } else if (ishell1 < ishell0) {
        ishell3 = 0;
        oprim3 = 0;
        ishell2 = 0;
        oprim2 = 0;
        ishell1++;
        oprim1 += nprim1;
        update_shell();
        return 1;
    } else if (ishell0 < gbasis->nshell-1) {
        ishell3 = 0;
        oprim3 = 0;
        ishell2 = 0;
        oprim2 = 0;
        ishell1 = 0;
        oprim1 = 0;
        ishell0++;
        oprim0 += nprim0;
        update_shell();
        return 1;
    } else {
        ishell3 = 0;
        ishell2 = 0;
        ishell1 = 0;
        ishell0 = 0;
        oprim3 = 0;
        oprim2 = 0;
        oprim1 = 0;
        oprim0 = 0;
        update_shell();
        return 0;
    }
}


void IterGB4::update_shell() {
    // Update fields that depend on shell and related counters.
    nprim0 = gbasis->nprims[ishell0];
    nprim1 = gbasis->nprims[ishell1];
    nprim2 = gbasis->nprims[ishell2];
    nprim3 = gbasis->nprims[ishell3];
    // Update indexes in output array
    ibasis0 = basis_offsets[ishell0];
    ibasis1 = basis_offsets[ishell1];
    ibasis2 = basis_offsets[ishell2];
    ibasis3 = basis_offsets[ishell3];
    // update centers
    r0 = gbasis->centers + 3*gbasis->shell_map[ishell0];
    r1 = gbasis->centers + 3*gbasis->shell_map[ishell1];
    r2 = gbasis->centers + 3*gbasis->shell_map[ishell2];
    r3 = gbasis->centers + 3*gbasis->shell_map[ishell3];
    // update shell types
    shell_type0 = gbasis->shell_types[ishell0];
    shell_type1 = gbasis->shell_types[ishell1];
    shell_type2 = gbasis->shell_types[ishell2];
    shell_type3 = gbasis->shell_types[ishell3];
    // reset contraction counters
    iprim0 = 0;
    iprim1 = 0;
    iprim2 = 0;
    iprim3 = 0;
}


int IterGB4::inc_prim() {
    // Increment primitive counters.
    if (iprim3 < nprim3-1) {
        iprim3++;
        update_prim();
        return 1;
    } else if (iprim2 < nprim2-1) {
        iprim2++;
        iprim3 = 0;
        update_prim();
        return 1;
    } else if (iprim1 < nprim1-1) {
        iprim1++;
        iprim3 = 0;
        iprim2 = 0;
        update_prim();
        return 1;
    } else if (iprim0 < nprim0-1) {
        iprim0++;
        iprim3 = 0;
        iprim2 = 0;
        iprim1 = 0;
        update_prim();
        return 1;
    } else {
        iprim0 = 0;
        iprim1 = 0;
        iprim2 = 0;
        iprim3 = 0;
        update_prim();
        return 0;
    }
}


void IterGB4::update_prim() {
    // Update fields that depend on primitive counters.
    alpha0 = gbasis->alphas[oprim0 + iprim0];
    alpha1 = gbasis->alphas[oprim1 + iprim1];
    alpha2 = gbasis->alphas[oprim2 + iprim2];
    alpha3 = gbasis->alphas[oprim3 + iprim3];
    con_coeff = gbasis->con_coeffs[oprim0 + iprim0]*
                gbasis->con_coeffs[oprim1 + iprim1]*
                gbasis->con_coeffs[oprim2 + iprim2]*
                gbasis->con_coeffs[oprim3 + iprim3];
    scales0 = gbasis->get_scales(oprim0 + iprim0);
    scales1 = gbasis->get_scales(oprim1 + iprim1);
    scales2 = gbasis->get_scales(oprim2 + iprim2);
    scales3 = gbasis->get_scales(oprim3 + iprim3);
}


#define ITERGB4_ASSIGN(idx0,idx1,idx2,idx3) output[(((i##idx0 + ibasis##idx0)*nbasis + i##idx1 + ibasis##idx1)*nbasis + i##idx2 + ibasis##idx2)*nbasis + i##idx3 + ibasis##idx3] = *tmp;

void IterGB4::store(const double *work, double *output) {
    // This routine is hardwired to work only for the dense storage
    // The code is horribly inefficient, but we don't care.
    long i0, i1, i2, i3;
    const long n0 = get_shell_nbasis(shell_type0);
    const long n1 = get_shell_nbasis(shell_type1);
    const long n2 = get_shell_nbasis(shell_type2);
    const long n3 = get_shell_nbasis(shell_type3);
    const long nbasis = gbasis->get_nbasis();
    const double* tmp = work;
    for (i0=0; i0<n0; i0++) {
        for (i1=0; i1<n1; i1++) {
            for (i2=0; i2<n2; i2++) {
                for (i3=0; i3<n3; i3++) {
                    ITERGB4_ASSIGN(0,1,2,3);
                    ITERGB4_ASSIGN(1,0,3,2);
                    ITERGB4_ASSIGN(2,3,0,1);
                    ITERGB4_ASSIGN(3,2,1,0);
                    ITERGB4_ASSIGN(0,3,2,1);
                    ITERGB4_ASSIGN(1,2,3,0);
                    ITERGB4_ASSIGN(2,1,0,3);
                    ITERGB4_ASSIGN(3,0,1,2);
                    tmp++;
                }
            }
        }
    }
}
