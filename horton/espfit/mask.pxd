# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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

cimport cell
cimport uniform

cdef extern from "mask.h":
    void multiply_dens_mask(double* rho, double lnrho0, double sigma,
        double* weights, long npoint)

    void multiply_near_mask(double* center, uniform.UniformGrid* ugrid,
        double r0, double gamma, double* weights)

    void multiply_far_mask(double* centers, long ncenter,
        uniform.UniformGrid* ugrid, double r0, double gamma,
        double* weights)
