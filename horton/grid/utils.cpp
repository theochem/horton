// Horton is a Density Functional Theory program.
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


//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include <stdexcept>
#include "utils.h"


double dot_multi(long npoint, long nvector, double** data) {
    double result = 0.0;
    for (long ipoint=npoint; ipoint>0; ipoint--) {
        double tmp = *(data[nvector-1]);
        data[nvector-1]++;
        for (long ivector=nvector-2; ivector>=0; ivector--) {
           tmp *= *(data[ivector]);
           data[ivector]++;
        }
        result += tmp;
    }
    return result;
}


void grid_distances(double *points, double *center, double *distances, long n) {
  double d, tmp;
  for (long i=0; i<n; i++) {
    // x
    d = *points - center[0];
    tmp = d*d;
    points++;
    // y
    d = *points - center[1];
    tmp += d*d;
    points++;
    // z
    d = *points - center[2];
    tmp += d*d;
    *distances = sqrt(tmp);
    points++;
    distances++;
  }
}


long index_wrap(long i, long high) {
    // Get around compiler weirdness
    long result = i%high;
    if (result<0) result += high;
    return result;
}


Range3Iterator::Range3Iterator(const long* ranges_begin, const long* ranges_end, const long* shape) :
    first(true), ranges_begin(ranges_begin), ranges_end(ranges_end), shape(shape) {};


bool next_one(long* i, long* iwrap, const long* ranges_begin, const long* ranges_end, const long* shape, const int axis) {
    bool result;
    i[axis] += 1;
    result = (i[axis] != ranges_end[axis]);
    if (not result) i[axis] = ranges_begin[axis];
    if (iwrap!=NULL) iwrap[axis] = index_wrap(i[axis], shape[axis]);
    return result;
}

bool Range3Iterator::next(long* i, long* iwrap) {
    if (first) {
        first = false;
        i[0] = ranges_begin[0];
        i[1] = ranges_begin[1];
        i[2] = ranges_begin[2];
        if (iwrap!=NULL) {
            iwrap[0] = index_wrap(i[0], shape[0]);
            iwrap[1] = index_wrap(i[1], shape[1]);
            iwrap[2] = index_wrap(i[2], shape[2]);
        }
        return true;
    } else {
        if (next_one(i, iwrap, ranges_begin, ranges_end, shape, 0)) return true;
        if (next_one(i, iwrap, ranges_begin, ranges_end, shape, 1)) return true;
        if (next_one(i, iwrap, ranges_begin, ranges_end, shape, 2)) return true;
        return false;
    }
}
