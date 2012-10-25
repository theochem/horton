# -*- coding: utf-8 -*-
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


import numpy as np
cimport numpy as np
np.import_array()

cimport nucpot

__all__ = [
    'compute_grid_nucpot',
]


def compute_grid_nucpot(np.ndarray[long, ndim=1] numbers not None,
                        np.ndarray[double, ndim=2] coordinates not None,
                        np.ndarray[double, ndim=2] points not None,
                        np.ndarray[double, ndim=1] output not None):
        assert numbers.flags['C_CONTIGUOUS']
        cdef long natom = numbers.shape[0]
        assert coordinates.flags['C_CONTIGUOUS']
        assert coordinates.shape[0] == natom
        assert coordinates.shape[1] == 3
        assert output.flags['C_CONTIGUOUS']
        cdef long npoint = output.shape[0]
        assert points.flags['C_CONTIGUOUS']
        assert points.shape[0] == npoint
        assert points.shape[1] == 3
        nucpot.compute_grid_nucpot(
            <long*>numbers.data, <double*>coordinates.data, natom,
            <double*>points.data, <double*>output.data, npoint)
