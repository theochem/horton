# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


cimport gbasis


cdef extern from "iter_gb.h":
    cdef cppclass IterGB1:
        IterGB1(gbasis.GBasis* gbasis)

        bint inc_shell()
        void update_shell()
        bint inc_prim()
        void update_prim()
        void store(double *work, double *output, long dim)

        # 'public' iterator fields
        long shell_type0
        double con_coeff, alpha0
        double *r0
        double *scales0
        long ibasis0

        # 'private' iterator fields
        long ishell0
        long nprim0, iprim0, oprim0

    cdef cppclass IterGB2:
        IterGB2(gbasis.GBasis* gbasis)

        bint inc_shell()
        void update_shell()
        bint inc_prim()
        void update_prim()
        void store(double *work, double *output)

        # 'public' iterator fields
        long shell_type0, shell_type1
        double con_coeff, alpha0, alpha1
        double *r0, *r1
        double *scales0, *scales1
        long ibasis0, ibasis1

        # 'private' iterator fields
        long ishell0, ishell1
        long nprim0, nprim1, iprim0, iprim1, oprim0, oprim1

    cdef cppclass IterGB4:
        IterGB4(gbasis.GBasis* gbasis)

        bint inc_shell()
        void update_shell()
        bint inc_prim()
        void update_prim()
        void store(double *work, double *output)

        # 'public' iterator fields
        long shell_type0, shell_type1, shell_type2, shell_type3
        long ishell3_max
        double con_coeff, alpha0, alpha1, alpha2, alpha3
        double *r0, *r1, *r2, *r3
        double *scales0, *scales1, *scales2, *scales3
        long ibasis0, ibasis1, ibasis2, ibasis3

        # 'private' iterator fields
        long ishell0, ishell1, ishell2, ishell3
        long nprim0, nprim1, nprim2, nprim3
        long iprim0, iprim1, iprim2, iprim3
        long oprim0, oprim1, oprim2, oprim3
