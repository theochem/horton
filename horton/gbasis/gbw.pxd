# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
cimport gbasis
cimport ints

cdef extern from "horton/gbasis/gbw.h":
    cdef cppclass GB4IntegralWrapper:
        GB4IntegralWrapper(gbasis.GOBasis* gobasis, ints.GB4Integral* gb4int)
        long get_nbasis()
        void select_2index(long index0, long index2,
                            long* pbegin0, long* pend0,
                            long* pbegin2, long* pend2)
        void compute()
        void compute_diagonal(double* diagonal)
        double* get_2index_slice(long index0, long index2)
