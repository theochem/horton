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


import numpy as np
cimport numpy as np
np.import_array()

cimport cell
cimport nucpot

__all__ = [
    'Cell', 'compute_grid_nucpot',
]


cdef class Cell:
    '''Representation of periodic boundary conditions.

       0, 1, 2 and 3 dimensional systems are supported. The cell vectors need
       not to be orthogonal.
    '''
    cdef cell.Cell* _this

    def __cinit__(self, *args, **kwargs):
        self._this = new cell.Cell()

    def __dealloc__(self):
        del self._this

    def __init__(self, np.ndarray[double, ndim=2] rvecs=None):
        '''
           **Arguments:**

           rvecs
                A numpy array with at most three cell vectors, layed out as
                rows in a rank-2 matrix. For non-periodic systems, this array
                must have shape (0,3).
        '''
        self.update_rvecs(rvecs)

    def update_rvecs(self, np.ndarray[double, ndim=2] rvecs):
        '''update_rvecs(rvecs)

           Change the cell vectors and recompute the reciprocal cell vectors.

           rvecs
                A numpy array with at most three cell vectors, layed out as
                rows in a rank-2 matrix. For non-periodic systems, this array
                must have shape (0,3).
        '''
        cdef np.ndarray[double, ndim=2] gvecs
        cdef int nvec
        if rvecs is None or rvecs.size == 0:
            gvecs = np.identity(3, float)
            nvec = 0
        else:
            if not rvecs.ndim==2 or not rvecs.flags['C_CONTIGUOUS'] or rvecs.shape[0] > 3 or rvecs.shape[1] != 3:
                raise TypeError('rvecs must be a C contiguous array with three columns and at most three rows.')
            nvec = len(rvecs)
            Up, Sp, Vt = np.linalg.svd(rvecs, full_matrices=True)
            S = np.ones(3, float)
            S[:nvec] = Sp
            U = np.identity(3, float)
            U[:nvec,:nvec] = Up
            gvecs = np.dot(U/S, Vt)
        self._this.update(<double*>rvecs.data, <double*>gvecs.data, nvec)

    def _get_nvec(self):
        '''The number of cell vectors'''
        return self._this.get_nvec()

    nvec = property(_get_nvec)

    def _get_volume(self):
        '''The generalized volume of the unit cell (length, area or volume)'''
        return self._this.get_volume()

    volume = property(_get_volume)

    def _get_rvecs(self):
        '''The real-space cell vectors, layed out as rows.'''
        cdef np.ndarray[double, ndim=2] result
        result = np.zeros((self.nvec, 3), float)
        self._this.copy_rvecs(<double*>result.data)
        result.setflags(write=False)
        return result

    rvecs = property(_get_rvecs)

    def _get_gvecs(self):
        '''The reciporcal-space cell vectors, layed out as rows.'''
        cdef np.ndarray[double, ndim=2] result
        result = np.zeros((self.nvec, 3), float)
        self._this.copy_gvecs(<double*>result.data)
        result.setflags(write=False)
        return result

    gvecs = property(_get_gvecs)

    def _get_rspacings(self):
        '''The (orthogonal) spacing between opposite sides of the real-space unit cell.'''
        cdef np.ndarray[double, ndim=1] result
        result = np.zeros(self.nvec, float)
        self._this.copy_rspacings(<double*>result.data)
        result.setflags(write=False)
        return result

    rspacings = property(_get_rspacings)

    def _get_gspacings(self, full=False):
        '''The (orthogonal) spacing between opposite sides of the reciprocal-space unit cell.'''
        cdef np.ndarray[double, ndim=1] result
        result = np.zeros(self.nvec, float)
        self._this.copy_gspacings(<double*>result.data)
        result.setflags(write=False)
        return result

    gspacings = property(_get_gspacings)

    def _get_parameters(self):
        '''The cell parameters (lengths and angles)'''
        rvecs = self.rvecs
        tmp = np.dot(rvecs, rvecs.T)
        lengths = np.sqrt(np.diag(tmp))
        tmp /= lengths
        tmp /= lengths.reshape((-1,1))
        if len(rvecs) < 2:
            cosines = np.arrays([])
        elif len(rvecs) == 2:
            cosines = np.array([tmp[0,1]])
        else:
            cosines = np.array([tmp[1,2], tmp[2,0], tmp[0,1]])
        angles = np.arccos(np.clip(cosines, -1, 1))
        return lengths, angles

    parameters = property(_get_parameters)

    def mic(self, np.ndarray[double, ndim=1] delta):
        """mic(delta)

           Apply the minimum image convention to delta in-place
        """
        assert delta.size == 3
        self._this.mic(<double*> delta.data)

    def to_center(self, np.ndarray[double, ndim=1] pos):
        '''to_center(pos)

           Return the corresponding position in the central cell
        '''
        assert pos.size == 3
        cdef np.ndarray[long, ndim=1] result
        result = np.zeros(self.nvec, int)
        self._this.to_center(<double*> pos.data, <long*> result.data)
        return result

    def add_vec(self, np.ndarray[double, ndim=1] delta, np.ndarray[long, ndim=1] r):
        """add_vec(delta, r)

           Add a linear combination of cell vectors, ``r``, to ``delta`` in-place
        """
        assert delta.size == 3
        assert r.size == self.nvec
        self._this.add_vec(<double*> delta.data, <long*> r.data)



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
