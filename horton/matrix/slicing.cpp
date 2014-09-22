// Horton is a development platform for electronic structure methods.
// Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

#include <iostream>
#include "slicing.h"

void get_slice_abcc(double* inp, double* inp2, double* out, long nbasis, long nvec){
    for (long k=0; k<nvec; k++){
        for (long a=0; a<nbasis; a++){
            for (long b=0; b<nbasis; b++){
                for (long c=0; c<nbasis; c++){
                    out[a*nbasis*nbasis + b*nbasis + c] += inp[k*nbasis*nbasis + a*nbasis + c] * inp2[k*nbasis*nbasis + b*nbasis + c];
                }
            }
        }
    }
}
void get_slice_abbc(double* inp, double* inp2, double* out, long nbasis, long nvec){
    for (long k=0; k<nvec; k++){
        for (long a=0; a<nbasis; a++){
            for (long c=0; c<nbasis; c++){
                for (long b=0; b<nbasis; b++){
                    out[a*nbasis*nbasis + b*nbasis + c] += inp[k*nbasis*nbasis + a*nbasis + b] * inp2[k*nbasis*nbasis + b*nbasis + c];
                }
            }
        }
    }

}

void sub_slice_abbc(double* inp, double* inp2, double* out, long nbasis, long nvec){
    for (long k=0; k<nvec; k++){
        for (long a=0; a<nbasis; a++){
            for (long c=0; c<nbasis; c++){
                for (long b=0; b<nbasis; b++){
                    out[a*nbasis*nbasis + b*nbasis + c] -= inp[k*nbasis*nbasis + a*nbasis + b] * inp2[k*nbasis*nbasis + b*nbasis + c];
                }
            }
        }
    }

}
