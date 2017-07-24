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

// #define DEBUG

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include "calc.h"
#include "common.h"

GBCalculator::GBCalculator(long max_shell_type, long dim_work, int basis_work)
    : max_shell_type(max_shell_type), work_pure(NULL), work_cart(NULL) {
  if (max_shell_type < 0) {
    throw std::domain_error("max_shell_type must be positive.");
  }
  max_nbasis = get_shell_nbasis(max_shell_type);
  nwork = dim_work;
  for (int i = 0; i < basis_work; i++) {
    nwork *= max_nbasis;
  }
  work_cart = new double[nwork];
  work_pure = new double[nwork];
}

GBCalculator::~GBCalculator() {
  delete[] work_cart;
  delete[] work_pure;
}

void GBCalculator::swap_work() {
  double *tmp;
  tmp = work_cart;
  work_cart = work_pure;
  work_pure = tmp;
}
