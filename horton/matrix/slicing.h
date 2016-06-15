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

// UPDATELIBDOCTITLE: Low-level routines for ``CholeskyLinalgFactory``

#ifndef HORTON_MATRIX_SLICING_H
#define HORTON_MATRIX_SLICING_H

void slice_to_three_abbc_abc(double* inp, double* inp2, double* out, double factor, bool clear, long nbasis, long nvec);
void slice_to_three_abcc_bac(double* inp, double* inp2, double* out, double factor, bool clear, long nbasis, long nvec);
void slice_to_three_abcc_abc(double* inp, double* inp2, double* out, double factor, bool clear, long nbasis, long nvec);

#endif
