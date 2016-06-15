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
'''C++ extensions'''


import numpy as np
cimport numpy as np
np.import_array()

cimport electrostatics
cimport mask

cimport horton.cext
cimport horton.grid.cext

from horton.log import log


__all__ = [
    # electrostatics
    'pair_ewald',
    'setup_esp_cost_cube',
    'compute_esp_grid_cube',
    # mask
    'multiply_dens_mask', 'multiply_near_mask', 'multiply_far_mask',
]


#
# electrostatics
#


def setup_esp_cost_cube(horton.grid.cext.UniformGrid ugrid not None,
                        np.ndarray[double, ndim=3] vref not None,
                        np.ndarray[double, ndim=3] weights not None,
                        np.ndarray[double, ndim=2] centers not None,
                        np.ndarray[double, ndim=2] A not None,
                        np.ndarray[double, ndim=1] B not None,
                        np.ndarray[double, ndim=0] C not None,
                        double rcut, double alpha, double gcut):

    is0d = (ugrid.pbc == [0, 0, 0]).all()
    is3d = (ugrid.pbc == [1, 1, 1]).all()

    assert vref.flags['C_CONTIGUOUS']
    assert vref.shape[0] == ugrid.shape[0]
    assert vref.shape[1] == ugrid.shape[1]
    assert vref.shape[2] == ugrid.shape[2]
    assert weights.flags['C_CONTIGUOUS']
    assert weights.shape[0] == ugrid.shape[0]
    assert weights.shape[1] == ugrid.shape[1]
    assert weights.shape[2] == ugrid.shape[2]
    assert centers.flags['C_CONTIGUOUS']
    ncenter = centers.shape[0]
    assert ncenter > 0
    assert centers.shape[1] == 3
    assert A.flags['C_CONTIGUOUS']
    assert A.shape[0] == ncenter+is3d
    assert A.shape[1] == ncenter+is3d
    assert B.flags['C_CONTIGUOUS']
    assert B.shape[0] == ncenter+is3d
    assert C.flags['C_CONTIGUOUS']
    if is0d:
        assert rcut == 0.0
        assert alpha == 0.0
        assert gcut == 0.0
    elif is3d:
        assert rcut > 0
        assert alpha > 0
        assert gcut > 0
    else:
        raise NotImplementedError

    electrostatics.setup_esp_cost_cube(ugrid._this, &vref[0, 0, 0],
        &weights[0, 0, 0], &centers[0, 0], &A[0, 0],
        &B[0], <double*>np.PyArray_DATA(C), ncenter, rcut, alpha, gcut)


def compute_esp_grid_cube(horton.grid.cext.UniformGrid ugrid not None,
                          np.ndarray[double, ndim=3] esp not None,
                          np.ndarray[double, ndim=2] centers not None,
                          np.ndarray[double, ndim=1] charges not None,
                          double rcut, double alpha, double gcut):
    assert centers.flags['C_CONTIGUOUS']
    ncenter = centers.shape[0]
    assert ncenter > 0
    assert centers.shape[1] == 3
    assert charges.flags['C_CONTIGUOUS']
    assert charges.shape[0] == ncenter
    assert esp.flags['C_CONTIGUOUS']
    assert esp.shape[0] == ugrid.shape[0]
    assert esp.shape[1] == ugrid.shape[1]
    assert esp.shape[2] == ugrid.shape[2]
    assert rcut > 0
    assert alpha > 0
    assert gcut > 0

    if ugrid.pbc.sum() in [0,3]:
        electrostatics.compute_esp_cube(ugrid._this, &esp[0, 0, 0],
            &centers[0, 0], &charges[0], ncenter, rcut, alpha,
            gcut)
    else:
        raise NotImplementedError


def pair_ewald(np.ndarray[double, ndim=1] delta not None,
               horton.cext.Cell cell not None,
               double rcut, double alpha, double gcut):

    assert delta.flags['C_CONTIGUOUS']
    assert delta.shape[0] == 3
    assert rcut > 0
    assert alpha > 0
    assert gcut > 0

    if cell.nvec == 3:
        return electrostatics.pair_ewald3d(&delta[0], cell._this, rcut, alpha, gcut)
    else:
        raise NotImplementedError


#
# mask
#


def multiply_dens_mask(np.ndarray[double, ndim=3] rho not None,
                       double lnrho0, double sigma,
                       np.ndarray[double, ndim=3] weights not None):
    assert rho.flags['C_CONTIGUOUS']
    assert weights.flags['C_CONTIGUOUS']
    assert weights.shape[0] == rho.shape[0]
    assert weights.shape[1] == rho.shape[1]
    assert weights.shape[2] == rho.shape[2]
    cdef long npoint = weights.size

    log.cite('hu2007', 'the density-based weight function for ESP fitting')
    mask.multiply_dens_mask(&rho[0, 0, 0], lnrho0, sigma, &weights[0, 0, 0], npoint)


def multiply_near_mask(np.ndarray[double, ndim=1] center not None,
                       horton.grid.cext.UniformGrid ugrid not None,
                       double r0, double gamma,
                       np.ndarray[double, ndim=3] weights not None):
    assert center.flags['C_CONTIGUOUS']
    assert center.shape[0] == 3
    assert r0 > 0
    assert gamma > 0
    assert weights.flags['C_CONTIGUOUS']
    assert weights.shape[0] == ugrid.shape[0]
    assert weights.shape[1] == ugrid.shape[1]
    assert weights.shape[2] == ugrid.shape[2]

    mask.multiply_near_mask(&center[0], ugrid._this, r0, gamma,
        &weights[0, 0, 0])


def multiply_far_mask(np.ndarray[double, ndim=2] centers not None,
                      horton.grid.cext.UniformGrid ugrid not None,
                      double r0, double gamma,
                      np.ndarray[double, ndim=3] weights not None):
    assert centers.flags['C_CONTIGUOUS']
    cdef long ncenter = centers.shape[0]
    assert ncenter > 0
    assert centers.shape[1] == 3
    assert r0 > 0
    assert gamma > 0
    assert weights.flags['C_CONTIGUOUS']
    assert weights.shape[0] == ugrid.shape[0]
    assert weights.shape[1] == ugrid.shape[1]
    assert weights.shape[2] == ugrid.shape[2]

    mask.multiply_far_mask(&centers[0, 0], ncenter, ugrid._this,
        r0, gamma, &weights[0, 0, 0])
