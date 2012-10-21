// Horton is a Density Functional Theory program.
// Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
#include <cstdlib>
#include <stdexcept>
#include "calc.h"
#include "cartpure.h"
#include "common.h"
using namespace std;



GBCalculator::GBCalculator(long max_shell_type): max_shell_type(max_shell_type), work_pure(NULL), work_cart(NULL) {
    if (max_shell_type < 0) {
        throw std::domain_error("max_shell_type must be positive.");
    }
    max_nbasis = get_shell_nbasis(max_shell_type);
}

GBCalculator::~GBCalculator() {
#ifdef DEBUG
    printf("%x %x\n", work_cart, work_pure);
#endif
    delete[] work_cart;
    delete[] work_pure;
}

void GBCalculator::swap_work() {
    double* tmp;
    tmp = work_cart;
    work_cart = work_pure;
    work_pure = tmp;
}
