// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2017 The HORTON Development Team
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

// UPDATELIBDOCTITLE: Iterators over Cartesian polynomials in one shell

#ifndef HORTON_GBASIS_ITER_POW_H_
#define HORTON_GBASIS_ITER_POW_H_

int iter_pow1_inc(long *n);

class IterPow1 {
 private:
  long shell_type0;
 public:
  IterPow1() : shell_type0(0), n0({0, 0, 0}), ibasis0(0) {}
  IterPow1(const IterPow1& other) = delete;

  void reset(long shell_type0);

  int inc();

  long n0[3];
  long ibasis0;
};

class IterPow2 {
 private:
  long shell_type0, shell_type1;
 public:
  IterPow2() : shell_type0(0), shell_type1(0), n0({0, 0, 0}), n1({0, 0, 0}),
               offset(0), ibasis0(0), ibasis1(0) {}
  IterPow2(const IterPow2& other) = delete;

  void reset(long shell_type0, long shell_type1);

  int inc();

  long n0[3];
  long n1[3];
  long offset, ibasis0, ibasis1;
};

#endif  // HORTON_GBASIS_ITER_POW_H_
