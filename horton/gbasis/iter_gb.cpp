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


#include <cstdlib>
#include "common.h"
#include "iter_gb.h"
using namespace std;


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
    // Increment shell and related counters.
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
    // Update indexes in output array (TODO: get them directly from the basis.)
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
    long i0, i1, n0, n1;
    const long nbasis = gbasis->get_nbasis();
    const long max_shell_nbasis = gbasis->get_max_shell_nbasis();
    n0 = get_shell_nbasis(shell_type0);
    n1 = get_shell_nbasis(shell_type1);
    for (i0=0; i0<n0; i0++) {
        for (i1=0; i1<n1; i1++) {
            output[(i0+ibasis0)*nbasis+i1+ibasis1] = work[i0*max_shell_nbasis+i1];
        }
    }
    if (ibasis0 != ibasis1) {
        for (i0=0; i0<n0; i0++) {
            for (i1=0; i1<n1; i1++) {
                output[(i1+ibasis1)*nbasis+i0+ibasis0] = work[i0*max_shell_nbasis+i1];
            }
        }
    }
}
