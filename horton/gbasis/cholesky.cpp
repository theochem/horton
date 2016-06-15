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

#include <cstddef>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include "horton/gbasis/cholesky.h"

/**
    Find the maximum diagonal error, used several times in the cholesky routine.
*/
double find_maxdiag(double* diagerr, long nbasis, long begin1, long end1,
    long begin2, long end2, long &index1, long &index2)
{
  double maxdiag = 0;
  index1 = -1; //safety
  index2 = -1; //safety
  for (long i1=begin1; i1<end1; i1++){
    for (long i2=begin2; i2<end2; i2++){
       if (maxdiag < diagerr[i1*nbasis + i2]){
         maxdiag = diagerr[i1*nbasis + i2];
         index1 = i1;
         index2 = i2;
      }
    }
  }
  return maxdiag;
}


long cholesky(GB4IntegralWrapper* gbw4, double** uninit_result,
    double threshold)
{
  if (threshold <= 0) {
    // The algorithm below may go crazy with a non-positive threshold.
    throw std::domain_error("Cholesky threshold must be strictly positive.");
  }

  long nbasis = gbw4->get_nbasis();
  double* diagonal = new double[nbasis*nbasis];  // allocate 2 index object
  double* diagerr = new double[nbasis*nbasis];   //  "
  double* slice = NULL;                          //  "
  // storage for 2-index cholesky vectors, start with allocation of 15% of the
  // full 4-center matrix
  std::vector<double>* vectors = new std::vector<double>;
  vectors->reserve(nbasis*nbasis*nbasis*nbasis*0.15);

  /*
    In some future version, we'll introduce a mask to only compute integrals
    that will be relevant after some pre-screening test. These will be flagged
    with a mask array.
  */
  //bool* mask = NULL;        // 2-index object

  /*
    Initialize the diagonal and set the diagerr equal to the diagonal (because
    we start with zero Cholesky vectors).
  */
  gbw4->compute_diagonal(diagonal);
  memcpy(diagerr, diagonal, sizeof(double)*nbasis*nbasis);

  /*
    This is extra stuff in the wrapper we'll need in future:
    gb4w->compute_centers(centers);
        centers of products 3 two 2-index objects
        (one for x, one for y and one for z)
    gb4w->compute_extents(extents);
        2-index object (should be 2-index for shells)
    gb4w->compute_qs(qs);
        2-index object (should be 2-index for shells) ?????
  */

  // Locate the maximum of diagerr -> 2 indexes. This determines the first
  // Cholesky vector that will be computed.
  long index1;
  long index2;
  double maxdiag = find_maxdiag(diagerr, nbasis, 0, nbasis, 0, nbasis,
                                index1, index2);

  // std::cout << "initial maxdiag " << maxdiag << " " << index1 << " " << index2 << std::endl;
  unsigned long nvec=0;
  do {
    // call wrapper to let it select a pair of shells for the given variables
    // index1 and index2.
    long begin1;
    long begin2;
    long end1;
    long end2;
    // index 2 is least significant
    gbw4->select_2index(index1, index2, &begin1, &end1, &begin2, &end2);
    // All integrals are computed for the selected pair of shells.
    // gbw4->compute(mask);
    gbw4->compute();

    do {
      // Get the the slice of computed four-center integrals that correspond
      // to the pair index1,index2.
      slice = gbw4->get_2index_slice(index1, index2);

      // construct new cholesky vector
      double* pastvector_sum = new double[nbasis*nbasis];
      memset(pastvector_sum, 0, sizeof(double)*nbasis*nbasis);
      maxdiag = 1.0 / sqrt(maxdiag);

      //compute sum of past Ls
      for (unsigned long l=0; l<nvec; l++){
        /*
        for (long i=0; i<nbasis; i++){
          for (long j=0; j<nbasis; j++){
            pastvector_sum[i*nbasis + j] +=
                (*vectors)[(l*nbasis*nbasis) + index1*nbasis + index2]
                *(*vectors)[(l*nbasis*nbasis) + i*nbasis + j];
          }
        }
        */
        cblas_daxpy(nbasis*nbasis,
                    (*vectors)[(l*nbasis*nbasis) + index1*nbasis + index2],
                    &(*vectors)[l*nbasis*nbasis], 1, pastvector_sum, 1);
      }

      //compute current L
      for (long i=0; i<nbasis; i++){
        for (long j=0; j<nbasis; j++){
          vectors->push_back(maxdiag * (slice[i*nbasis + j] -
                             pastvector_sum[i*nbasis + j]));
        }
      }

      // update diagerr
      for (long i=0; i<nbasis; i++){
        for (long j=0; j<nbasis; j++){
          diagerr[i*nbasis + j] -=
              (*vectors)[(nvec*nbasis*nbasis) + i*nbasis + j] *
              (*vectors)[(nvec*nbasis*nbasis) + i*nbasis + j];
        }
      }

      // We've just added one vector.
      nvec++;

      // Decide which 2-index within the shell, i.e. largest error.
      // Here, begin1, end1, begin2 and end2 are used to limit the search for
      // the maximum to the selected shells.
      maxdiag = find_maxdiag(diagerr, nbasis, begin1, end1, begin2, end2,
                             index1, index2);
      // std::cout << "current maxdiag " << maxdiag << " " << index1 << " " << index2 << std::endl;
    } while (maxdiag > threshold*1000);

    // Look for the new maximum error on the diagonal
    maxdiag = find_maxdiag(diagerr, nbasis, 0, nbasis, 0, nbasis,
                           index1, index2);

    // std::cout << "shell maxdiag " << maxdiag << " " << index1 << " " << index2 << std::endl;
  } while (maxdiag > threshold);

  // Concatenate results array. Not efficient! Who cares?
  *uninit_result = &((*vectors)[0]);

  return nvec;
}
