# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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


cdef extern from "horton/gbasis/iter_pow.h":
    bint iter_pow1_inc(long* n)

    cdef cppclass IterPow1:
        void reset(long shell_type0)
        bint inc()
        long n0[3]
        long ibasis0

    cdef cppclass IterPow2:
        void reset(long shell_type0, long shell_type1)
        bint inc()
        long n0[3]
        long n1[3]
        long offset, ibasis0, ibasis1
