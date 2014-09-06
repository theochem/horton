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

#include "cholesky.h"
#include <cstddef>
#include <cstring>
#include <cmath>


long cholesky(GB4IntegralWrapper* gbw4, double** uninit_result, double threshold) {
  long nbasis = gbw4->get_nbasis();
  double* diagonal = new double[nbasis*nbasis];  // allocate 2 index object
  double* diagerr = new double[nbasis*nbasis];  //  "
//  double* slice = new double[nbasis*nbasis];   //  "
  double* slice = NULL;   //  "
  std::vector<double>* vectors = new std::vector<double>;
  vectors->reserve(nbasis*nbasis*nbasis*nbasis*0.15);
  //bool* mask = NULL;        // 2-index object


  gbw4->compute_diagonal(diagonal);
  memcpy(diagerr, diagonal, sizeof(double)*nbasis*nbasis);

  /*
      gb4w->compute_centers(centers); // centers of products 3 two 2-index objects (one for x, one for y and one for z)
      gb4w->compute_extents(extents); // 2-index object (should be 2-index for shells)
      gb4w->compute_qs(qs);           // 2-index object (should be 2-index for shells) ?????
  */

  // locate the maximum of diagerr -> 2 indexes
  double maxdiag =0;
  long index1 = -1; //safety
  long index2 = -1; //safety
  for (long i=0; i<nbasis; i++){
      for (long j=0; j<nbasis; j++){
          if (maxdiag < diagerr[i*nbasis + j]){
              maxdiag = diagerr[i*nbasis + j];
              index1 = i;
              index2 = j;
          }
      }
  }
  //std::cout << std::endl << std::endl << "initial maxdiag " << maxdiag << std::endl;
  unsigned long nvec=0;
  do {
    // call wrapper to let it select a pair of shells.
    long begin1;
    long begin2;
    long end1;
    long end2;
    gbw4->select_2index(index1, index2, &begin1, &end1, &begin2, &end2); // index 2 is least significant
    // All integrals (not masked) are computed for the selected pair of shells.
    // gbw4->compute(mask);  // will replace folowing line at some point ...
    gbw4->compute();

    do {
        // break if there are no slices left in the selected shell
        slice = gbw4->get_2index_slice(index1, index2);
        // construct new cholesky vector
        double* pastvector_sum = new double[nbasis*nbasis];
        memset(pastvector_sum, 0, sizeof(double)*nbasis*nbasis);
        maxdiag = 1.0 / sqrt(maxdiag);


    //std::cout << "indices " << index1 << " " << index2 << std::endl;

    //compute sum of past Ls
    for (unsigned long l=0; l<nvec; l++){
            for (long i=0; i<nbasis; i++){
                for (long j=0; j<nbasis; j++){
                    pastvector_sum[i*nbasis + j] += (*vectors)[(l*nbasis*nbasis) + index1*nbasis + index2]
                                    * (*vectors)[(l*nbasis*nbasis) + i*nbasis + j];
                }
            }
    }

    //compute current L
        for (long i=0; i<nbasis; i++){
            for (long j=0; j<nbasis; j++){
                vectors->push_back(maxdiag * (slice[i*nbasis + j] - pastvector_sum[i*nbasis + j]));
            }
        }

        // update diagerr
        for (long i=0; i<nbasis; i++){
            for (long j=0; j<nbasis; j++){
                diagerr[i*nbasis + j] -= (*vectors)[(nvec*nbasis*nbasis) + i*nbasis + j] * (*vectors)[(nvec*nbasis*nbasis) + i*nbasis + j];
            }
        }
    nvec++;
        // Decide which 2-index within the shell, i.e. largest error.
        // Here, begin1, end1, begin2 and end2 are used to limit the search for
        // the maximum to the select shells.
        maxdiag=0;
    index1 = -1; //safety
        index2 = -1; //safety
        for (long i=begin1; i<end1; i++){
            for (long j=begin2; j<end2; j++){
                if (maxdiag < diagerr[i*nbasis + j]){
                    maxdiag = diagerr[i*nbasis + j];
                    index1 = i;
                    index2 = j;
                }
            }
        }
    } while (maxdiag > threshold*1000);
    maxdiag=0;
    index1 = -1; //safety
    index2 = -1; //safety
    for (long i=0; i<nbasis; i++){
        for (long j=0; j<nbasis; j++){
            if (maxdiag < diagerr[i*nbasis + j]){
                maxdiag = diagerr[i*nbasis + j];
                index1 = i;
                index2 = j;
            }
        }
    }
  //std::cout << std::endl << std::endl << "maxdiag " << maxdiag << std::endl;
  } while (maxdiag > threshold);

  // Concatenate results array. Not efficient!
  *uninit_result = &((*vectors)[0]);

  return nvec;
}

/*
void compute(bool* mask) {
    mask = NULL; //not implemented
    for (long shell3=0; shell3<nshell; shell3++) {
        for (long shell4=0; shell4<nshell; shell4++) {
            // compute shell

            // copy work array
        }
    }
}
*/
