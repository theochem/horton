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
'''C++ extensions'''

import numpy as np
cimport numpy as np
np.import_array()

cimport cell
cimport moments
cimport nucpot

from horton.utils import typecheck_geo

__all__ = [
    # cell.cpp
    'Cell', 'smart_wrap',
    # moments.cpp
    'fill_cartesian_polynomials', 'fill_pure_polynomials', 'fill_radial_polynomials',
    # nucpot.cpp
    'compute_grid_nucpot', 'compute_nucnuc',
]


#
# cell.cpp
#


cdef class Cell:
    '''Representation of periodic boundary conditions.

       0, 1, 2 and 3 dimensional systems are supported. The cell vectors don't
       need to be orthogonal.
    '''

    def __cinit__(self, np.ndarray[double, ndim=2] rvecs=None, initvoid=False):
        if initvoid:
            self._this = NULL
        elif rvecs is None:
            self._this = new cell.Cell(NULL, 0)
        else:
            assert rvecs.flags['C_CONTIGUOUS']
            assert rvecs.shape[0] <= 3
            assert rvecs.shape[1] == 3
            nvec = rvecs.shape[0]
            self._this = new cell.Cell(<double*>np.PyArray_DATA(rvecs), nvec)

    def __init__(self, np.ndarray[double, ndim=2] rvecs=None):
        pass

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct a Cell object from data in an HDF5 group'''
        if grp['rvecs'].size > 0:
            rvecs = np.array(grp['rvecs'])
            return cls(rvecs)
        else:
            return cls(None)

    def to_hdf5(self, grp):
        '''Write the cell object to an HDF5 group'''
        grp.create_dataset('rvecs', data=self.rvecs, maxshape=(None,None))

    @classmethod
    def from_parameters(cls, lengths, angles):
        """Construct a cell with the given parameters

           The a vector is always parallel with the x-axis and they point in the
           same direction. The b vector is always in the xy plane and points
           towards the positive y-direction. The c vector points towards the
           positive z-direction.

           The number of elements in the lengths and angles arrays determines
           the number of cell vectors. There are four cases:

           * len(lengths) == 0 and len(angles) == 0: 0 rvecs

           * len(lengths) == 1 and len(angles) == 0: 1 rvecs

           * len(lengths) == 2 and len(angles) == 1: 2 rvecs

           * len(lengths) == 3 and len(angles) == 3: 3 rvecs
        """
        if len(lengths) == 0 and len(angles) != 0:
            raise TypeError('When no lengths are given, no angles are expected.')
        elif len(lengths) == 1 and len(angles) != 0:
            raise TypeError('When one length is given, no angles are expected.')
        elif len(lengths) == 2 and len(angles) != 1:
            raise TypeError('When two lengths are given, one angle is expected.')
        elif len(lengths) == 3 and len(angles) != 3:
            raise TypeError('When three lengths are given, three angles are expected.')
        elif len(lengths) > 3:
            raise ValueError('More than three lengths are given.')

        for length in lengths:
            if length <= 0:
                raise ValueError("The length parameters must be strictly positive.")
        for angle in angles:
            if angle <= 0 or angle >= np.pi:
                raise ValueError("The angle parameters must lie in the range ]0 deg, 180 deg[.")

        if len(lengths) == 0:
            return Cell(None)

        rvecs = np.zeros((len(lengths), 3), float)

        if len(lengths) > 0:
            # first cell vector along x-axis
            rvecs[0, 0] = lengths[0]

        if len(lengths) > 1:
            # second cell vector in x-y plane
            if len(lengths) == 2:
                angle = angles[0]
            else:
                angle = angles[2]
            rvecs[1, 0] = np.cos(angle)*lengths[1]
            rvecs[1, 1] = np.sin(angle)*lengths[1]

        if len(lengths) > 2:
            # Finding the third cell vector is slightly more difficult. :-)
            # It works like this:
            # The dot products of a with c, b with c and c with c are known. the
            # vector a has only an x component, b has no z component. This results
            # in the following equations:
            u_a = lengths[0]*lengths[2]*np.cos(angles[1])
            u_b = lengths[1]*lengths[2]*np.cos(angles[0])
            rvecs[2, 0] = u_a/rvecs[0, 0]
            rvecs[2, 1] = (u_b - rvecs[1, 0]*rvecs[2, 0])/rvecs[1, 1]
            u_c = lengths[2]**2 - rvecs[2, 0]**2 - rvecs[2, 1]**2
            if u_c < 0:
                raise ValueError("The given cell parameters do not correspond to a unit cell.")
            rvecs[2, 2] = np.sqrt(u_c)

        return cls(rvecs)

    property nvec:
        '''The number of cell vectors'''
        def __get__(self):
            return self._this.get_nvec()

    property volume:
        '''The generalized volume of the unit cell (length, area or volume)'''
        def __get__(self):
            return self._this.get_volume()

    property rvecs:
        '''The real-space cell vectors, layed out as rows.'''
        def __get__(self):
            cdef np.ndarray[double, ndim=2] result
            result = np.zeros((self.nvec, 3), float)
            self._this.copy_rvecs(<double*>np.PyArray_DATA(result))
            result.setflags(write=False)
            return result

    property gvecs:
        '''The reciporcal-space cell vectors, layed out as rows.'''
        def __get__(self):
            cdef np.ndarray[double, ndim=2] result
            result = np.zeros((self.nvec, 3), float)
            self._this.copy_gvecs(<double*>np.PyArray_DATA(result))
            result.setflags(write=False)
            return result

    property rlengths:
        '''The lengths of the real-space vectors.'''
        def __get__(self):
            cdef np.ndarray[double, ndim=1] result
            result = np.zeros(self.nvec, float)
            self._this.copy_rlengths(<double*>np.PyArray_DATA(result))
            result.setflags(write=False)
            return result

    property glengths:
        '''The lengths of the reciprocal-space vectors.'''
        def __get__(self):
            cdef np.ndarray[double, ndim=1] result
            result = np.zeros(self.nvec, float)
            self._this.copy_glengths(<double*>np.PyArray_DATA(result))
            result.setflags(write=False)
            return result

    property rspacings:
        '''The (orthogonal) spacing between opposite sides of the real-space unit cell.'''
        def __get__(self):
            cdef np.ndarray[double, ndim=1] result
            result = np.zeros(self.nvec, float)
            self._this.copy_rspacings(<double*>np.PyArray_DATA(result))
            result.setflags(write=False)
            return result

    property gspacings:
        '''The (orthogonal) spacing between opposite sides of the reciprocal-space unit cell.'''
        def __get__(self):
            cdef np.ndarray[double, ndim=1] result
            result = np.zeros(self.nvec, float)
            self._this.copy_gspacings(<double*>np.PyArray_DATA(result))
            result.setflags(write=False)
            return result

    def get_rlength(self, int i):
        '''Get the length of the i-the real-space cell vector.'''
        return self._this.get_rlength(i);

    def get_glength(self, int i):
        '''Get the length of the i-the reciprocal cell vector.'''
        return self._this.get_glength(i);

    def get_rspacing(self, int i):
        '''Get the spacing between the i-the real-space cell planes.'''
        return self._this.get_rspacing(i);

    def get_gspacing(self, int i):
        '''Get the spacing between the i-the reciprocal cell planes.'''
        return self._this.get_gspacing(i);

    property parameters:
        '''The cell parameters (lengths and angles)'''
        def __get__(self):
            lengths = self.rlengths
            rvecs = self.rvecs
            tmp = np.dot(rvecs, rvecs.T)
            tmp /= lengths
            tmp /= lengths.reshape((-1,1))
            if len(rvecs) < 2:
                cosines = np.array([])
            elif len(rvecs) == 2:
                cosines = np.array([tmp[0,1]])
            else:
                cosines = np.array([tmp[1,2], tmp[2,0], tmp[0,1]])
            angles = np.arccos(np.clip(cosines, -1, 1))
            return lengths, angles

    def mic(self, np.ndarray[double, ndim=1] delta not None):
        '''Apply the minimum image convention to delta in-place'''
        assert delta.flags['C_CONTIGUOUS']
        assert delta.size == 3
        self._this.mic(&delta[0])

    def to_frac(self, np.ndarray[double, ndim=1] cart not None):
        '''Return the corresponding fractional coordinates'''
        assert cart.flags['C_CONTIGUOUS']
        assert cart.size == 3
        cdef np.ndarray[double, ndim=1] result
        result = np.zeros(3, float)
        self._this.to_frac(&cart[0], &result[0])
        return result

    def to_cart(self, np.ndarray[double, ndim=1] frac not None):
        '''Return the corresponding Cartesian coordinates'''
        assert frac.flags['C_CONTIGUOUS']
        assert frac.size == 3
        cdef np.ndarray[double, ndim=1] result
        result = np.zeros(3, float)
        self._this.to_cart(&frac[0], &result[0])
        return result

    def g_lincomb(self, np.ndarray[double, ndim=1] coeffs not None):
        '''Return a linear combination of reciprocal cell vectors'''
        assert coeffs.flags['C_CONTIGUOUS']
        assert coeffs.size == 3
        cdef np.ndarray[double, ndim=1] result
        result = np.zeros(3, float)
        self._this.g_lincomb(&coeffs[0], &result[0])
        return result

    def dot_rvecs(self, np.ndarray[double, ndim=1] cart not None):
        '''Return the corresponding dot product with the rvecs'''
        assert cart.flags['C_CONTIGUOUS']
        assert cart.size == 3
        cdef np.ndarray[double, ndim=1] result
        result = np.zeros(3, float)
        self._this.dot_rvecs(&cart[0], &result[0])
        return result

    def add_rvec(self, np.ndarray[double, ndim=1] delta not None,
                 np.ndarray[long, ndim=1] r not None):
        """Add a linear combination of real cell vectors, ``r``, to ``delta`` in-place"""
        assert delta.flags['C_CONTIGUOUS']
        assert delta.size == 3
        assert r.flags['C_CONTIGUOUS']
        assert r.size == self.nvec
        self._this.add_rvec(&delta[0], <long*>np.PyArray_DATA(r))

    def get_ranges_rcut(self, np.ndarray[double, ndim=1] center not None, double rcut):
        '''Return the integer ranges for linear combinations of cell vectors.

           **Arguments:**

           center
                The origin of the cutoff sphere

           rcut
                A cutoff radius

           The returned ranges span the linear combination of cell vectors that
           can be added to delta to obtain all periodic images within the cutoff
           sphere.
        '''
        assert center.flags['C_CONTIGUOUS']
        assert center.size == 3
        assert rcut >= 0

        cdef np.ndarray[long, ndim=1] ranges_begin = np.zeros(self.nvec, int)
        cdef np.ndarray[long, ndim=1] ranges_end = np.zeros(self.nvec, int)
        self._this.set_ranges_rcut(
            &center[0], rcut,  <long*>np.PyArray_DATA(ranges_begin),
            <long*>np.PyArray_DATA(ranges_end))
        return ranges_begin, ranges_end

def smart_wrap(long i, long shape, long pbc):
    '''Returned a standardize modulo operation i % shape if pbc is nonzero, -1 otherwise.'''
    return cell.smart_wrap(i, shape, pbc)



#
# moments.cpp
#


def fill_cartesian_polynomials(np.ndarray[double, ndim=1] output not None, long lmax):
    '''Fill the output vector with cartesian polynomials

       **Arguments:**

       output
            A double precision numpy array where the first three values are
            x, y and z coordinates.

       lmax
            The maximum angular momentum to compute.

       The polynomials are stored according to the conventions set in
       ``get_cartesian_powers``.

       **Returns:**

       The index of the first element of the array that contains the polynomials
       of the outermost shell.
    '''
    assert output.flags['C_CONTIGUOUS']
    if output.shape[0] < ((lmax+1)*(lmax+2)*(lmax+3))/6-1:
        raise ValueError('The size of the output array is not sufficient to store the polynomials.')
    return moments.fill_cartesian_polynomials(&output[0], lmax)


def fill_pure_polynomials(np.ndarray output not None, long lmax):
    '''Fill the output vector with pure polynomials

       **Arguments:**

       output
            This can either be a double precission Numpy vector or 2D array. In
            the first case, the first three values are z, x, and y coordinates.
            In the second case, the first three columns contain z, x and y
            coordinates.

       lmax
            The maximum angular momentum to compute.

       **Returns:**

       The index of the first element of the array that contains the polynomials
       of the outermost shell.
    '''
    assert output.flags['C_CONTIGUOUS']
    if output.ndim == 1:
        if output.shape[0] < (lmax+1)**2-1:
            raise ValueError('The size of the output array is not sufficient to store the polynomials.')
        return moments.fill_pure_polynomials(<double*>np.PyArray_DATA(output), lmax)
    elif output.ndim == 2:
        if output.shape[1] < (lmax+1)**2-1:
            raise ValueError('The size of the output array is not sufficient to store the polynomials.')
        return moments.fill_pure_polynomials_array(<double*>np.PyArray_DATA(output), lmax, output.shape[0], output.shape[1])
    else:
        raise NotImplementedError


def fill_radial_polynomials(np.ndarray[double, ndim=1] output not None, long lmax):
    '''Fill the output vector with radial polynomials

       **Arguments:**

       output
            A double precision numpy array where the first element is the radius

       lmax
            The maximum angular momentum to compute.

       All elements after the first will be filled up with increasing powers of
       the first element, up to lmax.
    '''
    assert output.flags['C_CONTIGUOUS']
    if output.shape[0] < lmax:
        raise ValueError('The size of the output array is not sufficient to store the polynomials.')
    moments.fill_radial_polynomials(&output[0], lmax)



#
# nucpot.cpp
#


def compute_grid_nucpot(np.ndarray[double, ndim=2] coordinates not None,
                        np.ndarray[double, ndim=1] charges not None,
                        np.ndarray[double, ndim=2] points not None,
                        np.ndarray[double, ndim=1] output not None):
    '''Compute the potential due to a set of (nuclear) point charges

       coordinates
            A (N, 3) float numpy array with Cartesian coordinates of the
            atoms.

       charges
            A (N,) numpy vector with the atomic charges.

       points
            An (M, 3) array with grid points where the potential must be
            computed.

       output
            An (M,) output array in which the potential is stored.
    '''
    # type checking
    assert coordinates.flags['C_CONTIGUOUS']
    assert charges.flags['C_CONTIGUOUS']
    ncharge, coordinates, charges = typecheck_geo(coordinates, None, charges, need_numbers=False)
    assert output.flags['C_CONTIGUOUS']
    cdef long npoint = output.shape[0]
    assert points.flags['C_CONTIGUOUS']
    assert points.shape[0] == npoint
    assert points.shape[1] == 3
    # actual computation
    nucpot.compute_grid_nucpot(
        &coordinates[0,0], &charges[0], ncharge,
        &points[0,0], &output[0], npoint)


def compute_nucnuc(np.ndarray[double, ndim=2] coordinates not None,
                   np.ndarray[double, ndim=1] charges not None):
    '''Compute interaction energy of the nuclei

       **Arguments:**

       coordinates
            A (N, 3) float numpy array with Cartesian coordinates of the
            atoms.

       charges
            A (N,) numpy vector with the atomic charges.
    '''
    # type checking
    assert coordinates.flags['C_CONTIGUOUS']
    assert charges.flags['C_CONTIGUOUS']
    ncharge, coordinates, charges = typecheck_geo(coordinates, None, charges, need_numbers=False)
    # actual computation
    result = 0.0
    natom = len(charges)
    for i in xrange(ncharge):
        for j in xrange(i):
            distance = np.linalg.norm(coordinates[i]-coordinates[j])
            result += charges[i]*charges[j]/distance
    return result
