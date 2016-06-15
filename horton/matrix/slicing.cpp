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

#include <cstring>
#include "slicing.h"

void slice_to_three_abbc_abc(double* inp, double* inp2, double* out, double factor, bool clear, long nbasis, long nvec) {
    if (clear) {
        memset(out, 0, sizeof(double)*nbasis*nbasis*nbasis);
    }
    for (long k=0; k<nvec; k++){
        for (long a=0; a<nbasis; a++){
            for (long c=0; c<nbasis; c++){
                for (long b=0; b<nbasis; b++){
                    out[(a*nbasis + b)*nbasis + c] +=
                        factor *
                        inp[(k*nbasis + a)*nbasis + b] *
                        inp2[(k*nbasis + b)*nbasis + c];
                }
            }
        }
    }
}


void slice_to_three_abcc_bac(double* inp, double* inp2, double* out, double factor, bool clear, long nbasis, long nvec) {
    if (clear) {
        memset(out, 0, sizeof(double)*nbasis*nbasis*nbasis);
    }
    for (long k=0; k<nvec; k++){
        for (long a=0; a<nbasis; a++){
            for (long b=0; b<nbasis; b++){
                for (long c=0; c<nbasis; c++){
                    out[(b*nbasis + a)*nbasis + c] +=
                        factor *
                        inp[(k*nbasis + a)*nbasis + c] *
                        inp2[(k*nbasis + b)*nbasis + c];
                }
            }
        }
    }
}


void slice_to_three_abcc_abc(double* inp, double* inp2, double* out, double factor, bool clear, long nbasis, long nvec) {
    if (clear) {
        memset(out, 0, sizeof(double)*nbasis*nbasis*nbasis);
    }
    for (long k=0; k<nvec; k++){
        for (long a=0; a<nbasis; a++){
            for (long b=0; b<nbasis; b++){
                for (long c=0; c<nbasis; c++){
                    out[(a*nbasis + b)*nbasis + c] +=
                        factor *
                        inp[(k*nbasis + a)*nbasis + c] *
                        inp2[(k*nbasis + b)*nbasis + c];
                }
            }
        }
    }
}
