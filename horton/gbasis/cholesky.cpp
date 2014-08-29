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


std::vector<double*> cholesky(GB4IntegralWrapper* gbw4, double threshold) {
  long nbasis = gbw4->get_nbasis();
  double* diagonal = new double[nbasis*nbasis];  // allocate 2 index object
  double* diagerr = new double[nbasis*nbasis];  //  "
  double* slice = new double[nbasis*nbasis];   //  "
  std::vector<double*> vectors;
  //bool* mask = NULL;        // 2-index object
  long maxdiag = 0;


  gbw4->compute_diagonal(diagonal);
  memcpy(diagerr, diagonal, sizeof(double)*nbasis*nbasis);

  /*
      gb4w->compute_centers(centers); // centers of products 3 two 2-index objects (one for x, one for y and one for z)
      gb4w->compute_extents(extents); // 2-index object (should be 2-index for shells)
      gb4w->compute_qs(qs);           // 2-index object (should be 2-index for shells) ?????
  */

  do {
    // locate the maximum of diagerr -> 2 indexes

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
    // call wrapper to let it select a pair of shells.
    long begin1;
    long begin2;
    long end1;
    long end2;
    gbw4->select_2index(index1, index2, &begin1, &end1, &begin2, &end2); // index 2 is least significant
    /*
        // Decide which parts need to be computed according to some inequality
        // Here, begin1, end1, begin2 and end2 are used.
        for (long index3=0; index3<basis; index3++) {
            for (long index4=0; index4<basis; index4++) {
                bool mask34 = false;
                for (long index1=begin1; index1<end1; index1++) {
                    for (long index2=begin2; index2<end2; index2++) {
                        distance = norm(centers[index1, index2] - centers[index3, index4]);
                        extentssum = extents[index1, index2] + [index3, index4];
                        if (distances < extentssum) {
                            mask34 = true;
                            break;
                        }
                        qqr = qs[index1,index2]*qs[index3,index4]/(distance - extentssum)**3;
                        if (qqr > threshold) {
                            mask34 = true;
                            break;
                        }
                    }
                    if (mask34) break;
                }
                mask[index3, index4] = mask34;
            }
        }
    */
    // All integrals (not masked) are computed for the selected pair of shells.
    // gbw4->compute(mask);  // will replace folowing line at some point ...
    gbw4->compute();

    do {
        // break if there are no slices left in the selected shell
        gbw4->get_2index_slice(index1, index2, slice);
        // construct new cholesky vector
        double* vector = new double[nbasis*nbasis];
        double* pastvector_sum = new double[nbasis*nbasis];
        memset(pastvector_sum, 0, sizeof(double)*nbasis*nbasis);
        maxdiag = 1.0 / sqrt(maxdiag);


        //OPTION A:
        //store old l
        for (long l=0; l<vectors.size(); l++){
            for (long i=0; i<nbasis; i++){
                for (long j=0; j<nbasis; j++){
                    pastvector_sum[i*nbasis + j] += vectors[l][index1*nbasis + index2]
                                                        * vectors[l][i*nbasis + j];
                }
            }
        }

        for (long i=0; i<nbasis; i++){
            for (long j=0; j<nbasis; j++){
                vector[i*nbasis + j] = maxdiag * (slice[index1*nbasis + index2]
                                                    - pastvector_sum[i*nbasis + j]);
            }
        }

        /*
        //OPTION B:
        for (long i=0; i<nbasis; i++){
            for (long j=0; j<nbasis; j++){
                double temp = 0;
                for (long l=0; l<nbasis; l++){
                    temp += vectors[l][i*nbasis + j];
                }
                vector[i*nbasis + j] = maxdiag * (slice[index1*nbasis + index2] - temp);
            }
        }
        */

        // update diagerr
        for (long i=0; i<nbasis; i++){
            for (long j=0; j<nbasis; j++){
                diagerr[i*nbasis + j] -= vector[i*nbasis + j]*vector[i*nbasis + j];
            }
        }
        //store vector
        vectors.push_back(vector);
        // Decide which 2-index within the shell, i.e. largest error.
        // Here, begin1, end1, begin2 and end2 are used to limit the search for
        // the maximum to the select shells.
        maxdiag = 0;
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

  } while (maxdiag > threshold);

  return vectors;
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
