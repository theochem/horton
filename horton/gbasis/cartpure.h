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


#ifndef HORTON_GBASIS_CARTPURE_H
#define HORTON_GBASIS_CARTPURE_H


#include <exception>
#define MAX_SHELL_TYPE 9

class bad_shell_type_exception: public std::exception
{
  virtual const char* what() const throw()
  {
    return "shell_type not supported by cart_to_pure_low.";
  }
};

void cart_to_pure_low(double *work_cart, double* work_pure, long
    shell_type, long stride, long spacing, long count);


#endif
