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


#include "horton/gbasis/gbw.h"
#include "horton/gbasis/common.h"

GB4IntegralWrapper::GB4IntegralWrapper(GOBasis* gobasis, GB4Integral* gb4int) :
    gobasis(gobasis), gb4int(gb4int)
{
  max_shell_size = get_shell_nbasis(gobasis->get_max_shell_type());
  slice_size = gobasis->get_nbasis()*gobasis->get_nbasis();
  /*
     Allocate temporary storage for four-center integrals:
       - The maximum number of 2-index objects computed in one go with libint is
         max_shell_size*max_shell_size.
       - The size of a 2-index object is nbasis*nbasis.
  */
  integrals = new double[max_shell_size*max_shell_size*slice_size];
}

GB4IntegralWrapper::~GB4IntegralWrapper() {
  delete[] integrals;
}

void GB4IntegralWrapper::compute_shell(long ishell1, long ishell3)
{
  // Configure the four-center integral with the right input for this
  // quadruple of shells. (index0 and index2 are fixed.)
  gb4int->reset(gobasis->shell_types[ishell0], gobasis->shell_types[ishell1],
                gobasis->shell_types[ishell2], gobasis->shell_types[ishell3],
                gobasis->centers + gobasis->shell_map[ishell0]*3, gobasis->centers + gobasis->shell_map[ishell1]*3,
                gobasis->centers + gobasis->shell_map[ishell2]*3, gobasis->centers + gobasis->shell_map[ishell3]*3);

  // Quadruple loop over all primitives for these four shells.
  for (long iprim0 = 0; iprim0 < gobasis->nprims[ishell0]; iprim0++) {
    for (long iprim1 = 0; iprim1 < gobasis->nprims[ishell1]; iprim1++) {
      for (long iprim2 = 0; iprim2 < gobasis->nprims[ishell2]; iprim2++) {
        for (long iprim3 = 0; iprim3 < gobasis->nprims[ishell3]; iprim3++) {
          gb4int->add(gobasis->con_coeffs[gobasis->get_prim_offsets()[ishell0] + iprim0]*
                      gobasis->con_coeffs[gobasis->get_prim_offsets()[ishell1] + iprim1]*
                      gobasis->con_coeffs[gobasis->get_prim_offsets()[ishell2] + iprim2]*
                      gobasis->con_coeffs[gobasis->get_prim_offsets()[ishell3] + iprim3],
                      gobasis->alphas[gobasis->get_prim_offsets()[ishell0] + iprim0],
                      gobasis->alphas[gobasis->get_prim_offsets()[ishell1] + iprim1],
                      gobasis->alphas[gobasis->get_prim_offsets()[ishell2] + iprim2],
                      gobasis->alphas[gobasis->get_prim_offsets()[ishell3] + iprim3],
                      gobasis->get_scales(gobasis->get_prim_offsets()[ishell0] + iprim0),
                      gobasis->get_scales(gobasis->get_prim_offsets()[ishell1] + iprim1),
                      gobasis->get_scales(gobasis->get_prim_offsets()[ishell2] + iprim2),
                      gobasis->get_scales(gobasis->get_prim_offsets()[ishell3] + iprim3));
        }
      }
    }
  }

  // Convert to pure functions if needed.
  gb4int->cart_to_pure();
}


void GB4IntegralWrapper::select_2index(long index0, long index2,
                            long* pbegin0, long* pend0,
                            long* pbegin2, long* pend2) {
  ishell0 = gobasis->get_shell_lookup()[index0];
  begin0 = gobasis->get_basis_offsets()[ishell0];
  *pbegin0 = begin0;
  *pend0 = *pbegin0 + get_shell_nbasis(gobasis->shell_types[ishell0]);

  ishell2 = gobasis->get_shell_lookup()[index2];
  begin2 = gobasis->get_basis_offsets()[ishell2];
  *pbegin2 = begin2;
  *pend2 = *pbegin2 + get_shell_nbasis(gobasis->shell_types[ishell2]);
}

void GB4IntegralWrapper::compute() {
  // Double loop over second and fourth shell of the four-index object. The
  // entire range over these two indexes is included in the 2-index slices.
  for (long ishell1 = 0; ishell1 < gobasis->nshell; ishell1++) {
    for (long ishell3 = 0; ishell3 < gobasis->nshell; ishell3++) {
      // Compute integrals for the given combination of shells.
      compute_shell(ishell1, ishell3);

      // Copy data from work array to ``integrals``, the temporary storage of
      // this wrapper.
      const double* tmp = gb4int->get_work();
      const long n0 = get_shell_nbasis(gobasis->shell_types[ishell0]);
      const long n1 = get_shell_nbasis(gobasis->shell_types[ishell1]);
      const long n2 = get_shell_nbasis(gobasis->shell_types[ishell2]);
      const long n3 = get_shell_nbasis(gobasis->shell_types[ishell3]);
      for (long i0=0; i0<n0; i0++) {
        for (long i1=0; i1<n1; i1++) {
          for (long i2=0; i2<n2; i2++) {
            for (long i3=0; i3<n3; i3++) {
              integrals[((i0)*max_shell_size + i2)*slice_size +
                        (i1+gobasis->get_basis_offsets()[ishell1])*gobasis->get_nbasis() +
                        (i3+gobasis->get_basis_offsets()[ishell3])] = *tmp;
              tmp++;
            }
          }
        }
      }
    }
  }
}

void GB4IntegralWrapper::compute_diagonal(double* diagonal) {
  // Double loop over second and fourth shell of the four-index object. The
  // entire range over these two indexes is included in the 2-index slices.
  for (long ishell1 = 0; ishell1 < gobasis->nshell; ishell1++) {
    for (long ishell3 = 0; ishell3 < gobasis->nshell; ishell3++) {
      // Compute integrals for the given combination of shells.
      ishell0 = ishell1;
      ishell2 = ishell3;
      compute_shell(ishell1, ishell3);

      // copy data from work array to the output array.
      const double* tmp = gb4int->get_work();
      const long n1 = get_shell_nbasis(gobasis->shell_types[ishell1]);
      const long n3 = get_shell_nbasis(gobasis->shell_types[ishell3]);
      for (long i1=0; i1<n1; i1++) {
        for (long i3=0; i3<n3; i3++) {
          diagonal[(i1+gobasis->get_basis_offsets()[ishell1])*gobasis->get_nbasis() +
                   (i3+gobasis->get_basis_offsets()[ishell3])] = tmp[(n1*i1+i1)*n3*n3 + n3*i3+i3];
        }
      }
    }
  }
}

double* GB4IntegralWrapper::get_2index_slice(long index0, long index2) {
    return integrals + ((index0-begin0)*max_shell_size + (index2-begin2))*slice_size;
}
