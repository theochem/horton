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


long fac2(long n);
long binom(long n, long m);

double gpt_coeff(long k, long n0, long n1, double pa, double pb);
double gob_overlap_int1d(long n0, long n1, double pa, double pb, double gamma);

double gob_overlap(double exp0, long nx0, long ny0, long nz0, double *r0,
                   double exp1, long nx1, long ny1, long nz1, double *r1);
