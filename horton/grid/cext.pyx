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

from libc.stdlib cimport malloc, free

cimport lebedev_laikov
cimport becke
cimport cubic_spline
cimport rtransform
cimport utils


__all__ = [
    # lebedev_laikov
    'lebedev_laikov_npoint', 'lebedev_laikov_sphere', 'lebedev_laikov_npoints',
    # becke
    'becke_helper_atom',
    # cubic_spline
    'tridiag_solve', 'tridiagsym_solve', 'CubicSpline',
    # rtransform
    'BaseRTransform', 'IdentityRTransform', 'LinearRTransform', 'LogRTransform',
    # utils
    'dot_multi',
]


#
# lebedev_laikov
#


lebedev_laikov_npoints = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230,
                          266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730,
                          2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294,
                          5810]


def lebedev_laikov_npoint(int lvalue):
    '''lebedev_laikov_npoint(lvalue)

       Return the number of Lebedev-Laikov grid points for a given angular
       momentum.
    '''
    return lebedev_laikov.lebedev_laikov_npoint(lvalue)


def lebedev_laikov_sphere(np.ndarray[double, ndim=2] points,
                          np.ndarray[double, ndim=1] weights):
    '''lebedev_laikov_sphere(grid)

       Fill the grid with a Lebedev Laikov grid points of a given size.

       **Arguments:**

       points
            The output array for the grid points, shape (npoint,3).

       weights
            The output array for the grid weights, shape (npoint,).
    '''
    assert points.flags['C_CONTIGUOUS']
    assert points.shape[1] == 3
    npoint = points.shape[0]
    assert weights.flags['C_CONTIGUOUS']
    assert weights.shape[0] == npoint
    lebedev_laikov.lebedev_laikov_sphere(npoint, <double*>points.data,
                                         <double*>weights.data)


#
# becke
#


def becke_helper_atom(np.ndarray[double, ndim=2] points,
                      np.ndarray[double, ndim=1] weights,
                      np.ndarray[double, ndim=1] radii,
                      np.ndarray[double, ndim=2] centers,
                      int select, int order):
    '''beck_helper_atom(points, weights, radii, centers, i, k)

       Compute the Becke weights for a given atom an a grid.

       **Arguments:**

       points
            The Cartesian coordinates of the grid points. Numpy array with
            shape (npoint, 3)

       weights
            The output array where the Becke partitioning weights are written.
            Numpy array with shape (npoint,)

       radii
            The covalent radii used to shrink/enlarge basins in the Becke
            scheme.

       centers
            The positions of the nuclei.

       select
            The selected atom for which the weights should be created.

       order
            The order of the switching functions. (That is k in Becke's paper.)

       See Becke's paper for the details: http://dx.doi.org/10.1063/1.454033
    '''
    assert points.flags['C_CONTIGUOUS']
    assert points.shape[1] == 3
    npoint = points.shape[0]
    assert weights.flags['C_CONTIGUOUS']
    assert weights.shape[0] == npoint
    assert radii.flags['C_CONTIGUOUS']
    natom = radii.shape[0]
    assert centers.flags['C_CONTIGUOUS']
    assert centers.shape[0] == natom
    assert centers.shape[1] == 3
    assert select >= 0 and select < natom
    assert order > 0

    becke.becke_helper_atom(points.shape[0], <double*>points.data,
                            <double*>weights.data, natom, <double*>radii.data,
                            <double*>centers.data, select, order)


#
# cubic_spline
#

def tridiag_solve(np.ndarray[double, ndim=1] diag_low,
                  np.ndarray[double, ndim=1] diag_mid,
                  np.ndarray[double, ndim=1] diag_up,
                  np.ndarray[double, ndim=1] right,
                  np.ndarray[double, ndim=1] solution):
    assert diag_mid.flags['C_CONTIGUOUS']
    n = diag_mid.shape[0]
    assert n > 1
    assert diag_low.flags['C_CONTIGUOUS']
    assert diag_low.shape[0] == n-1
    assert diag_up.flags['C_CONTIGUOUS']
    assert diag_up.shape[0] == n-1
    assert right.flags['C_CONTIGUOUS']
    assert right.shape[0] == n
    assert solution.flags['C_CONTIGUOUS']
    assert solution.shape[0] == n
    cubic_spline.tridiag_solve(<double*>diag_low.data, <double*>diag_mid.data,
                               <double*>diag_up.data, <double*>right.data,
                               <double*>solution.data, n)


def tridiagsym_solve(np.ndarray[double, ndim=1] diag_mid,
                     np.ndarray[double, ndim=1] diag_up,
                     np.ndarray[double, ndim=1] right,
                     np.ndarray[double, ndim=1] solution):
    assert diag_mid.flags['C_CONTIGUOUS']
    n = diag_mid.shape[0]
    assert n > 1
    assert diag_up.flags['C_CONTIGUOUS']
    assert diag_up.shape[0] == n-1
    assert right.flags['C_CONTIGUOUS']
    assert right.shape[0] == n
    assert solution.flags['C_CONTIGUOUS']
    assert solution.shape[0] == n
    cubic_spline.tridiagsym_solve(<double*>diag_mid.data, <double*>diag_up.data,
                                  <double*>right.data, <double*>solution.data,
                                  n)


cdef class CubicSpline(object):
    cdef cubic_spline.CubicSpline* _this
    cdef cubic_spline.Extrapolation* _ep

    def __cinit__(self, np.ndarray[double, ndim=1] y not None, np.ndarray[double, ndim=1] d=None):
        cdef double* ddata
        assert y.flags['C_CONTIGUOUS']
        n = y.shape[0]
        if d is None:
            ddata = <double*>NULL
        else:
            assert d.flags['C_CONTIGUOUS']
            assert d.shape[0] == n
            ddata = <double*>d.data
        # Only exponential extrapolation is needed for now.
        self._ep = <cubic_spline.Extrapolation*>(new cubic_spline.ExponentialExtrapolation())
        self._this = new cubic_spline.CubicSpline(
            <double*>y.data, ddata, self._ep, n
        )

    def __dealloc__(self):
        del self._this
        del self._ep

    def copy_y(self):
        cdef np.npy_intp shape[1]
        shape[0] = self._this.n
        tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, <void*> self._this.y)
        return tmp.copy()

    def copy_d(self):
        cdef np.npy_intp shape[1]
        shape[0] = self._this.n
        tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, <void*> self._this.d)
        return tmp.copy()

    def __call__(self, np.ndarray[double, ndim=1] new_x not None, np.ndarray[double, ndim=1] new_y=None):
        assert new_x.flags['C_CONTIGUOUS']
        new_n = new_x.shape[0]
        if new_y is None:
            new_y = np.zeros(new_n, float)
        else:
            assert new_y.flags['C_CONTIGUOUS']
            assert new_y.shape[0] == new_n
        self._this.eval(<double*>new_x.data, <double*>new_y.data, new_n)
        return new_y

    def deriv(self, np.ndarray[double, ndim=1] new_x not None, np.ndarray[double, ndim=1] new_d=None):
        assert new_x.flags['C_CONTIGUOUS']
        new_n = new_x.shape[0]
        if new_d is None:
            new_d = np.zeros(new_n, float)
        else:
            assert new_d.flags['C_CONTIGUOUS']
            assert new_d.shape[0] == new_n
        self._this.eval_deriv(<double*>new_x.data, <double*>new_d.data, new_n)
        return new_d

    def deriv2(self, np.ndarray[double, ndim=1] new_x not None, np.ndarray[double, ndim=1] new_d2=None):
        assert new_x.flags['C_CONTIGUOUS']
        new_n = new_x.shape[0]
        if new_d2 is None:
            new_d2 = np.zeros(new_n, float)
        else:
            assert new_d2.flags['C_CONTIGUOUS']
            assert new_d2.shape[0] == new_n
        self._this.eval_deriv2(<double*>new_x.data, <double*>new_d2.data, new_n)
        return new_d2

    def integrate(self):
        return self._this.integrate()


#
# rtransform
#


cdef class BaseRTransform(object):
    cdef rtransform.BaseRTransform* _this

    property npoint:
        def __get__(self):
            return self._this.get_npoint()

    def radius(self, double t):
        return self._this.radius(t)

    def deriv(self, double t):
        return self._this.deriv(t)

    def inv(self, double r):
        return self._this.inv(r)

    def radius_array(self, np.ndarray[double, ndim=1] t not None,
                           np.ndarray[double, ndim=1] r not None):
        assert t.flags['C_CONTIGUOUS']
        cdef int n = t.shape[0]
        assert r.flags['C_CONTIGUOUS']
        assert r.shape[0] == n
        self._this.radius_array(<double*>t.data, <double*>r.data, n)

    def deriv_array(self, np.ndarray[double, ndim=1] t not None,
                          np.ndarray[double, ndim=1] d not None):
        assert t.flags['C_CONTIGUOUS']
        cdef int n = t.shape[0]
        assert d.flags['C_CONTIGUOUS']
        assert d.shape[0] == n
        self._this.deriv_array(<double*>t.data, <double*>d.data, n)

    def inv_array(self, np.ndarray[double, ndim=1] r not None,
                        np.ndarray[double, ndim=1] t not None):
        assert r.flags['C_CONTIGUOUS']
        cdef int n = r.shape[0]
        assert t.flags['C_CONTIGUOUS']
        assert t.shape[0] == n
        self._this.inv_array(<double*>r.data, <double*>t.data, n)

    def get_radii(self):
        '''Return an array with radii'''
        result = np.arange(self.npoint, dtype=float)
        self.radius_array(result, result)
        return result

    def get_volume_elements(self):
        '''Return an array with volume elements associated with the transform'''
        result = np.arange(self.npoint, dtype=float)
        self.deriv_array(result, result)
        return result

    @classmethod
    def from_string(cls, s):
        '''Construct a BaseRTransform subclass from a string.'''
        words = s.split()
        clsname = words[0]
        args = words[1:]
        if clsname == 'IdentityRTransform':
            if len(args) != 1:
                raise ValueError('The IdentityRTransform needs one argument, got %i.' % len(words))
            npoint = int(args[0])
            return IdentityRTransform(npoint)
        if clsname == 'LinearRTransform':
            if len(args) != 3:
                raise ValueError('The LinearRTransform needs three arguments, got %i.' % len(words))
            rmin = float(args[0])
            rmax = float(args[1])
            npoint = int(args[2])
            return LinearRTransform(rmin, rmax, npoint)
        if clsname == 'LogRTransform':
            if len(args) != 3:
                raise ValueError('The LogRTransform needs three arguments, got %i.' % len(words))
            rmin = float(args[0])
            rmax = float(args[1])
            npoint = int(args[2])
            return LogRTransform(rmin, rmax, npoint)
        else:
            raise TypeError('Unkown BaseRTransform subclass: %s' % clsname)

    def to_string(self):
        raise NotImplementedError


cdef class IdentityRTransform(BaseRTransform):
    '''For testing only'''
    def __cinit__(self, int npoint):
        self._this = <rtransform.BaseRTransform*>(new rtransform.IdentityRTransform(npoint))

    def to_string(self):
        return ' '.join(['IdentityRTransform', repr(self.npoint)])


cdef class LinearRTransform(BaseRTransform):
    '''A linear grid.

       The grid points are distributed as follows:

       .. math:: r_i = \\alpha i + r_0

       with

       .. math:: \\alpha = (r_{N-1} -r_0)/(N-1).
    '''
    def __cinit__(self, double rmin, double rmax, int npoint):
        self._this = <rtransform.BaseRTransform*>(new rtransform.LinearRTransform(rmin, rmax, npoint))

    property rmin:
        def __get__(self):
            return (<rtransform.LinearRTransform*>self._this).get_rmin()

    property rmax:
        def __get__(self):
            return (<rtransform.LinearRTransform*>self._this).get_rmax()

    property alpha:
        def __get__(self):
            return (<rtransform.LinearRTransform*>self._this).get_alpha()

    def to_string(self):
        return ' '.join(['LinearRTransform', repr(self.rmin), repr(self.rmax), repr(self.npoint)])


cdef class LogRTransform(BaseRTransform):
    '''A logarithmic grid.

       The grid points are distributed as follows:

       .. math:: r_i = r_0 \\alpha^i

       with

       .. math:: \\alpha = \log(r_{N-1}/r_0)/(N-1).
    '''
    def __cinit__(self, double rmin, double rmax, int npoint):
        self._this = <rtransform.BaseRTransform*>(new rtransform.LogRTransform(rmin, rmax, npoint))

    property rmin:
        def __get__(self):
            return (<rtransform.LogRTransform*>self._this).get_rmin()

    property rmax:
        def __get__(self):
            return (<rtransform.LogRTransform*>self._this).get_rmax()

    property alpha:
        def __get__(self):
            return (<rtransform.LogRTransform*>self._this).get_alpha()

    def to_string(self):
        return ' '.join(['LogRTransform', repr(self.rmin), repr(self.rmax), repr(self.npoint)])

#
# utils
#


def dot_multi(*args):
    if len(args) == 0:
        return 0.0

    cdef double** pointers = <double **>malloc(len(args)*sizeof(double*))
    if pointers == NULL:
        raise MemoryError()
    cdef np.ndarray[double, ndim=1] arg

    npoint = None
    try:
        for i in xrange(len(args)):
            arg = args[i]
            assert arg.flags['C_CONTIGUOUS']
            if npoint is None:
                npoint = arg.shape[0]
            else:
                assert npoint == arg.shape[0]
            pointers[i] = <double*>arg.data
        result = utils.dot_multi(npoint, len(args), pointers)
    finally:
        free(pointers)
    return result
