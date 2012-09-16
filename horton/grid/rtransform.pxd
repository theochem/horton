# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


cdef extern from "rtransform.h":
    cdef cppclass BaseRTransform:
        BaseRTransform(int npoint) except +
        double radius(double t)
        double deriv(double t)
        double inv(double t)

        void radius_array(double* t, double* r, int n)
        void deriv_array(double* t, double* d, int n)
        void inv_array(double* r, double* t, int n)
        int get_npoint()


    cdef cppclass IdentityRTransform:
        IdentityRTransform(int npoint) except +


    cdef cppclass LinearRTransform:
        LinearRTransform(double rmin, double rmax, int npoint) except +

        double get_rmin()
        double get_rmax()
        double get_alpha()


    cdef cppclass LogRTransform:
        LogRTransform(double rmin, double rmax, int npoint) except +

        double get_rmin()
        double get_rmax()
        double get_alpha()
