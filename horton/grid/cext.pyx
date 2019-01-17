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
from horton.log import log
from horton.grid.utils import parse_args_integrate
import numbers


cimport numpy as np
np.import_array()

from libc.stdlib cimport malloc, free
cimport libcpp

cimport lebedev_laikov
cimport becke
cimport cubic_spline
cimport evaluate
cimport ode2
cimport rtransform
cimport uniform
cimport utils

cimport horton.cext


__all__ = [
    # lebedev_laikov
    'lebedev_laikov_npoints', 'lebedev_laikov_lmaxs', 'lebedev_laikov_sphere',
    # becke
    'becke_helper_atom',
    # cubic_spline
    'Extrapolation', 'ZeroExtrapolation', 'CuspExtrapolation',
    'PowerExtrapolation', 'PotentialExtrapolation', 'tridiagsym_solve', 'CubicSpline',
    'compute_cubic_spline_int_weights',
    # evaluate
    'eval_spline_grid', 'eval_decomposition_grid',
    # ode2
    'hermite_overlap2', 'hermite_overlap3', 'hermite_node', 'hermite_product2',
    'build_ode2',
    # rtransform
    'RTransform', 'IdentityRTransform', 'LinearRTransform', 'ExpRTransform',
    'PowerRTransform', 'HyperbolicRTransform',
    # UniformGrid
    'UniformGrid', 'index_wrap',
    # utils
    'dot_multi', 'dot_multi_moments_cube', 'dot_multi_moments',
]


#
# lebedev_laikov
#


lebedev_laikov_npoints = {
    6: 3, 14: 5, 26: 7, 38: 9, 50: 11, 74: 13, 86: 15, 110: 17, 146: 19, 170:
    21, 194: 23, 230: 25, 266: 27, 302: 29, 350: 31, 434: 35, 590: 41, 770: 47,
    974: 53, 1202: 59, 1454: 65, 1730: 71, 2030: 77, 2354: 83, 2702: 89, 3074:
    95, 3470: 101, 3890: 107, 4334: 113, 4802: 119, 5294: 125, 5810: 131,
}

lebedev_laikov_lmaxs = {
    0: 6, 1: 6, 2: 6, 3: 6, 4: 14, 5: 14, 6: 26, 7: 26, 8: 38, 9: 38, 10: 50,
    11: 50, 12: 74, 13: 74, 14: 86, 15: 86, 16: 110, 17: 110, 18: 146, 19: 146,
    20: 170, 21: 170, 22: 194, 23: 194, 24: 230, 25: 230, 26: 266, 27: 266, 28:
    302, 29: 302, 30: 350, 31: 350, 32: 434, 33: 434, 34: 434, 35: 434, 36: 590,
    37: 590, 38: 590, 39: 590, 40: 590, 41: 590, 42: 770, 43: 770, 44: 770, 45:
    770, 46: 770, 47: 770, 48: 974, 49: 974, 50: 974, 51: 974, 52: 974, 53: 974,
    54: 1202, 55: 1202, 56: 1202, 57: 1202, 58: 1202, 59: 1202, 60: 1454, 61:
    1454, 62: 1454, 63: 1454, 64: 1454, 65: 1454, 66: 1730, 67: 1730, 68: 1730,
    69: 1730, 70: 1730, 71: 1730, 72: 2030, 73: 2030, 74: 2030, 75: 2030, 76:
    2030, 77: 2030, 78: 2354, 79: 2354, 80: 2354, 81: 2354, 82: 2354, 83: 2354,
    84: 2702, 85: 2702, 86: 2702, 87: 2702, 88: 2702, 89: 2702, 90: 3074, 91:
    3074, 92: 3074, 93: 3074, 94: 3074, 95: 3074, 96: 3470, 97: 3470, 98: 3470,
    99: 3470, 100: 3470, 101: 3470, 102: 3890, 103: 3890, 104: 3890, 105: 3890,
    106: 3890, 107: 3890, 108: 4334, 109: 4334, 110: 4334, 111: 4334, 112: 4334,
    113: 4334, 114: 4802, 115: 4802, 116: 4802, 117: 4802, 118: 4802, 119: 4802,
    120: 5294, 121: 5294, 122: 5294, 123: 5294, 124: 5294, 125: 5294, 126: 5810,
    127: 5810, 128: 5810, 129: 5810, 130: 5810, 131: 5810,
}


def lebedev_laikov_sphere(np.ndarray[double, ndim=2] points not None,
                          np.ndarray[double, ndim=1] weights not None):
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
    lebedev_laikov.lebedev_laikov_sphere(npoint, &points[0, 0], &weights[0])


#
# becke
#


def becke_helper_atom(np.ndarray[double, ndim=2] points not None,
                      np.ndarray[double, ndim=1] weights not None,
                      np.ndarray[double, ndim=1] radii not None,
                      np.ndarray[double, ndim=2] centers not None,
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

    becke.becke_helper_atom(points.shape[0], &points[0, 0], &weights[0], natom,
                            &radii[0], &centers[0, 0], select, order)


#
# cubic_spline
#

def tridiagsym_solve(np.ndarray[double, ndim=1] diag_mid not None,
                     np.ndarray[double, ndim=1] diag_up not None,
                     np.ndarray[double, ndim=1] right not None,
                     np.ndarray[double, ndim=1] solution not None):
    assert diag_mid.flags['C_CONTIGUOUS']
    n = diag_mid.shape[0]
    assert n > 1
    assert diag_up.flags['C_CONTIGUOUS']
    assert diag_up.shape[0] == n-1
    assert right.flags['C_CONTIGUOUS']
    assert right.shape[0] == n
    assert solution.flags['C_CONTIGUOUS']
    assert solution.shape[0] == n
    cubic_spline.tridiagsym_solve(&diag_mid[0], &diag_up[0],
                                  &right[0], &solution[0],
                                  n)


cdef class Extrapolation(object):
    def __dealloc__(self):
        del self._this

    def eval_left(self, double x):
        '''Evaluate the extrapolation function at the left of the cubic spline interval'''
        return self._this.eval_left(x)

    def eval_right(self, double x):
        '''Evaluate the extrapolation function at the right of the cubic spline interval'''
        return self._this.eval_right(x)

    def deriv_left(self, double x):
        '''Evaluate the extrapolation function derivative at the left of the cubic spline interval'''
        return self._this.deriv_left(x)

    def deriv_right(self, double x):
        '''Evaluate the extrapolation function derivative at the right of the cubic spline interval'''
        return self._this.deriv_right(x)

    def to_string(self):
        '''Return an extrapolation object in string respresentation'''
        return self.__class__.__name__

    @classmethod
    def from_string(cls, s):
        '''Create an extrpolation object from a string description'''
        words = s.split()
        if words[0] == 'ZeroExtrapolation':
            if len(words) != 1:
                raise ValueError('ZeroExtrapolation takes not arguments.')
            return ZeroExtrapolation()
        elif words[0] == 'CuspExtrapolation':
            if len(words) != 1:
                raise ValueError('CuspExtrapolation takes not arguments.')
            return CuspExtrapolation()
        elif words[0] == 'PowerExtrapolation':
            if len(words) != 2:
                raise ValueError('PowerExtrapolation takes one argument.')
            power = float(words[1])
            return PowerExtrapolation(power)
        elif words[0] == 'PotentialExtrapolation':
            if len(words) != 2:
                raise ValueError('PotentialExtrapolation takes one argument.')
            l = int(words[1])
            return PotentialExtrapolation(l)
        else:
            raise NotImplementedError


cdef class ZeroExtrapolation(Extrapolation):
    '''Zero left and right of the cubic spline interval'''
    def __cinit__(self):
        self._this = <cubic_spline.Extrapolation*>(new cubic_spline.ZeroExtrapolation())


cdef class CuspExtrapolation(Extrapolation):
    '''Exponential extrapolation at the left side, zero at the right size'''
    def __cinit__(self):
        self._this = <cubic_spline.Extrapolation*>(new cubic_spline.CuspExtrapolation())


cdef class PowerExtrapolation(Extrapolation):
    '''Zero at the right side, power law at the left side'''
    def __cinit__(self, double power):
        self._this = <cubic_spline.Extrapolation*>(new cubic_spline.PowerExtrapolation(power))

    property power:
        '''The power parameter'''
        def __get__(self):
            return (<cubic_spline.PowerExtrapolation*>self._this).get_power()

    def to_string(self):
        return 'PowerExtrapolation %s' % repr(self.power)


cdef class PotentialExtrapolation(Extrapolation):
    """Zero at the right side, power law at the left side"""
    def __cinit__(self, long l):
        self._this = <cubic_spline.Extrapolation*>(new cubic_spline.PotentialExtrapolation(l))

    property l:
        """The angular momentum parameter"""
        def __get__(self):
            return (<cubic_spline.PotentialExtrapolation*>self._this).get_l()

    property amp_left:
        """The amplitude of the polynomial on the left side"""
        def __get__(self):
            return (<cubic_spline.PotentialExtrapolation*>self._this).get_amp_left()

    property amp_right:
        """The amplitude of the polynomial on the right side"""
        def __get__(self):
            return (<cubic_spline.PotentialExtrapolation*>self._this).get_amp_right()

    def to_string(self):
        return 'PotentialExtrapolation %s' % repr(self.l)


cdef class CubicSpline(object):
    '''A cubic spline object

       **Arguments:**

       y
            The function values at the 1D grid.

       **Optional arguments:**

       dx
            The derivative of the function values at the 1D grid. If not given,
            they are determined such that the second derivative of the cubic
            spline is continuous at the grid points.

       rtransform
            The transformation object that specifies the 1D grid. If not given,
            an identity transform is used

       extrapolation
            The extrapolation object that specifies the spline function outside
            the interval determined by the 1D grid. By default,
            CuspExtrapolation() is used.
    '''

    def __cinit__(self,
                  np.ndarray[double, ndim=1] y not None,
                  np.ndarray[double, ndim=1] dx=None,
                  RTransform rtransform=None,
                  Extrapolation extrapolation=None):
        assert y.flags['C_CONTIGUOUS']
        self._y = y
        n = y.shape[0]

        # Set the rtransform
        if rtransform is None:
            self._rtransform = IdentityRTransform(n)
        else:
            self._rtransform = rtransform
            assert rtransform.npoint == n

        # use given derivatives or construct new ones.
        v = self._rtransform.get_deriv()
        if dx is None:
            self._dt = np.zeros(n, float)
            cubic_spline.solve_cubic_spline_system(&y[0], <double*>self._dt.data, n)
            self._dx = self._dt/v
        else:
            assert dx.flags['C_CONTIGUOUS']
            assert dx.shape[0] == n
            self._dx = dx
            self._dt = dx*v

        # Only exponential extrapolation is needed for now
        if extrapolation is None:
            self._extrapolation = CuspExtrapolation()
        else:
            self._extrapolation = extrapolation

        self._this = new cubic_spline.CubicSpline(
            <double*>self._y.data, <double*>self._dt.data, self._extrapolation._this,
            self._rtransform._this, n
        )

    def __dealloc__(self):
        del self._this

    @classmethod
    def from_hdf5(cls, grp):
        return cls(
            grp['y'][:],
            grp['d'][:],
            RTransform.from_string(grp.attrs['rtransform']),
            Extrapolation.from_string(grp.attrs['extrapolation']),
        )

    def to_hdf5(self, grp):
        grp['y'] = self.y
        grp['d'] = self.dx
        grp.attrs['rtransform'] = self.rtransform.to_string()
        grp.attrs['extrapolation'] = self.extrapolation.to_string()

    property rtransform:
        '''The RTransform object used for this spline'''
        def __get__(self):
            return self._rtransform

    property extrapolation:
        '''The extrapolation object used for this spline'''
        def __get__(self):
            return self._extrapolation

    property y:
        '''Array with function values at the grid points'''
        def __get__(self):
            return self._y.view()

    property dx:
        '''Array with derivatives (towards x) at the grid points'''
        def __get__(self):
            return self._dx.view()

    property dt:
        '''Array with derivatives (towards t) at the grid points'''
        def __get__(self):
            return self._dt.view()

    def __call__(self, np.ndarray[double, ndim=1] new_x not None,
                 np.ndarray[double, ndim=1] new_y=None):
        '''evaluate the spline on a grid

           **Arguments:**

           new_x
                A numpy array with the x-values at which the spline must be
                evaluated.

           **Optional arguments:**

           new_y
                When given, it is used as output argument. This array must have
                the same size of new_x.

           **Returns:** new_y
        '''
        assert new_x.flags['C_CONTIGUOUS']
        new_n = new_x.shape[0]
        if new_y is None:
            new_y = np.zeros(new_n, float)
        else:
            assert new_y.flags['C_CONTIGUOUS']
            assert new_y.shape[0] == new_n
        self._this.eval(&new_x[0], &new_y[0], new_n)
        return new_y

    def deriv(self, np.ndarray[double, ndim=1] new_x not None,
              np.ndarray[double, ndim=1] new_dx=None):
        '''Evaluate the derivative of the spline (towards x) on a grid

           **Arguments:**

           new_x
                A numpy array with the x-values at which the spline must be
                evaluated.

           **Optional arguments:**

           new_dx
                When given, it is used as output argument. This array must have
                the same size of new_x.

           **Returns:** new_dx
        '''
        assert new_x.flags['C_CONTIGUOUS']
        new_n = new_x.shape[0]
        if new_dx is None:
            new_dx = np.zeros(new_n, float)
        else:
            assert new_dx.flags['C_CONTIGUOUS']
            assert new_dx.shape[0] == new_n
        self._this.eval_deriv(&new_x[0], &new_dx[0], new_n)
        return new_dx


def compute_cubic_spline_int_weights(np.ndarray[double, ndim=1] weights not None):
    assert weights.flags['C_CONTIGUOUS']
    npoint = weights.shape[0]
    cubic_spline.compute_cubic_spline_int_weights(&weights[0], npoint)


#
# evaluate
#


def eval_spline_grid(CubicSpline spline not None,
                     np.ndarray[double, ndim=1] center not None,
                     np.ndarray[double, ndim=1] output not None,
                     np.ndarray[double, ndim=2] points not None,
                     horton.cext.Cell cell not None):
    '''Evaluate a spherically symmetric function on a general grid

       **Arguments:**

       spline
            The cubic spline that contains the radial dependence of the
            spherically symmetric function.

       center
            The center of the spherically symmetric function.

       points
            An array with grid points, with shape (N, 3)

       output
            The output array in which the result is stored.

       cell
            A specification of the periodic boundary conditions.
    '''
    assert center.flags['C_CONTIGUOUS']
    assert center.shape[0] == 3
    assert output.flags['C_CONTIGUOUS']
    assert points.flags['C_CONTIGUOUS']
    assert points.shape[1] == 3
    assert points.shape[0] == output.shape[0]

    evaluate.eval_spline_grid(spline._this, &center[0], &output[0],
                              &points[0, 0], cell._this, output.shape[0])


def eval_decomposition_grid(splines not None,
                     np.ndarray[double, ndim=1] center not None,
                     np.ndarray[double, ndim=1] output not None,
                     np.ndarray[double, ndim=2] points not None,
                     horton.cext.Cell cell not None):
    '''Evaluate a sphericall decomposition on a general grid

       **Arguments:**

       splines
            The splines with the spherical decomposition. These are usually
            generated with AtomicGrid.get_spherical_decomposition.

       center
            The center of the spherically symmetric function.

       points
            An array with grid points, with shape (N, 3)

       output
            The output array in which the result is stored.

       cell
            A specification of the periodic boundary conditions.
    '''

    # parse the splines argument and construct an array of c++ cubic spline objects
    assert len(splines) > 0
    cdef CubicSpline spline
    cdef cubic_spline.CubicSpline** cpp_splines = <cubic_spline.CubicSpline**>malloc(len(splines)*sizeof(cubic_spline.CubicSpline*))
    if cpp_splines == NULL:
        raise MemoryError()

    try:
        for i in xrange(len(splines)):
            spline = splines[i]
            cpp_splines[i] = spline._this

        assert center.flags['C_CONTIGUOUS']
        assert center.shape[0] == 3
        assert output.flags['C_CONTIGUOUS']
        assert points.flags['C_CONTIGUOUS']
        assert points.shape[1] == 3
        assert points.shape[0] == output.shape[0]

        evaluate.eval_decomposition_grid(cpp_splines, <double*>center.data,
            <double*>output.data, <double*>points.data, cell._this, len(splines),
            output.shape[0])
    finally:
        free(cpp_splines)


#
# ode2
#

# Except for build_ode2, all functions wrapped below are only used for testing
# purposes.

def hermite_overlap2(long xmax, long i0, libcpp.bool deriv0, long i1, libcpp.bool deriv1):
    return ode2.hermite_overlap2(xmax, i0, deriv0, i1, deriv1)

def hermite_overlap3(long xmax, long i0, libcpp.bool deriv0, long i1, libcpp.bool deriv1, long i2, libcpp.bool deriv2):
    return ode2.hermite_overlap3(xmax, i0, deriv0, i1, deriv1, i2, deriv2)

def hermite_node(long x, long center, libcpp.bool kind, libcpp.bool deriv):
    return ode2.hermite_node(x, center, kind, deriv)

def hermite_product2(long x, long i0, libcpp.bool deriv0, long i1, libcpp.bool deriv1):
    return ode2.hermite_product2(x, i0, deriv0, i1, deriv1)


def build_ode2(np.ndarray[double, ndim=1] by not None,
               np.ndarray[double, ndim=1] bd not None,
               np.ndarray[double, ndim=1] ay not None,
               np.ndarray[double, ndim=1] ad not None,
               np.ndarray[double, ndim=1] fy not None,
               np.ndarray[double, ndim=1] fd not None,
               bcs):
    '''Build set of equations for a second order ODE problem

       The ODE has the following form:

       .. math::
           u''(x) + b(x) u'(x) + a(x) u(x) = f(x)

       A linear system is constructed to approximate the solution for this
       equation on an equidistant grid with spacing 1. The function values
       and the derivatives of the known functions must be given, together with
       the boundary conditions, i.e. the function value of u(x) at the first and
       last grid point.

       **Arguments:**

       by
            An array with function values of b(x) at the grid points.

       bd
            An array with the first derivatives of b(x) at the grid points.

       ay
            An array with function values of a(x) at the grid points.

       ad
            An array with the first derivatives of a(x) at the grid points.

       fy
            An array with function values of f(x) at the grid points.

       fd
            An array with the first derivatives of f(x) at the grid points.

       bcs
            A four-tuple with boundary condition specifications: (uyfirst,
            udfirst, uylast, ydlast). Excatly two of these four values must
            be None, with the two other ones fixing the boundary conditions.

       The arrays by, bd, ay, ad, fy and fd must have the same length.
    '''

    def merge(np.ndarray[double, ndim=1] y not None, np.ndarray[double, ndim=1] d not None):
        '''Put y and d in one vector'''
        return np.array([y, d]).T.ravel()

    # parse array arguments
    npoint = by.shape[0]
    assert bd.shape[0] == npoint
    assert ay.shape[0] == npoint
    assert ad.shape[0] == npoint
    assert fy.shape[0] == npoint
    assert fd.shape[0] == npoint
    cdef np.ndarray[double] b = merge(by, bd)
    cdef np.ndarray[double] a = merge(ay, ad)
    cdef np.ndarray[double] f = merge(fy, fd)

    # parse boundary conditions argument
    if len(bcs) != 4:
        raise ValueError('bcs must have four elements.')
    if sum([bc is None for bc in bcs]) != 2:
        raise ValueError('bcs must contain two None elements.')

    cdef double c_bcs[4]
    for i in xrange(4):
        c_bcs[i] = 0.0 if bcs[i] is None else bcs[i]

    cdef double* p_bcs[4]
    for i in xrange(4):
        if bcs[i] is None:
            p_bcs[i] = NULL
        else:
            p_bcs[i] = &c_bcs[i]

    # prepare output
    cdef np.ndarray[double, ndim=2] coeffs = np.zeros((2*npoint, 2*npoint), float)
    cdef np.ndarray[double] rhs = np.zeros(2*npoint, float)

    # call c routine
    ode2.build_ode2(&b[0], &a[0], &f[0], p_bcs, &coeffs[0,0], &rhs[0], npoint)

    # done
    return coeffs, rhs

#
# rtransform
#


# The following class contains a lot of duplicate code because Cython can not do this (yet):
# ctypedef double (rtransform.RTransform::*fn_ptr)(double)
# ctypedef void (rtransform.RTransform::*fn_array_ptr)(double*, double*,int)

cdef class RTransform(object):
    '''A definition of (radial) grid points by means of a transformation.

       The definition starts from a uniform 1D grid with spacing 1 and starting
       point 0: 0, 1, 2, 3, ... npoint-1. These values are defined on the
       so-called t-axis. The transformation is a function r=f(t) that defines
       the actual grid points on the r-axis: f(0), f(1), f(2), ... f(npoint-1).
       Different implementation for the function f are available.
    '''

    property npoint:
        def __get__(self):
            return self._this.get_npoint()

    def radius(self, t, output=None):
        '''Return the 1D grid points for the given index(es)

           **Arguments:**

           t
                A number or an array of numbers for the indexes. t may be
                fractional or integer.
        '''
        if isinstance(t, numbers.Number):
            assert output is None
            return self._this.radius(t)
        elif isinstance(t, np.ndarray):
            return self._radius_array(t.astype(float), output)
        else:
            raise NotImplementedError

    def _radius_array(self, np.ndarray[double, ndim=1] t not None,
                            np.ndarray[double, ndim=1] output=None):
        assert t.flags['C_CONTIGUOUS']
        cdef int n = t.shape[0]
        if output is None:
            output = np.zeros(n, float)
        else:
            assert output.flags['C_CONTIGUOUS']
            assert output.shape[0] == n
        self._this.radius_array(&t[0], &output[0], n)
        return output

    def deriv(self, t, output=None):
        '''Return the derivative of the transformation for the given index(es)

           **Arguments:**

           t
                A number or an array of numbers for the indexes. t may be
                fractional or integer.
        '''
        if isinstance(t, numbers.Number):
            assert output is None
            return self._this.deriv(t.astype(float))
        elif isinstance(t, np.ndarray):
            return self._deriv_array(t, output)
        else:
            raise NotImplementedError

    def _deriv_array(self, np.ndarray[double, ndim=1] t not None,
                            np.ndarray[double, ndim=1] output=None):
        assert t.flags['C_CONTIGUOUS']
        cdef int n = t.shape[0]
        if output is None:
            output = np.zeros(n, float)
        else:
            assert output.flags['C_CONTIGUOUS']
            assert output.shape[0] == n
        self._this.deriv_array(&t[0], &output[0], n)
        return output

    def deriv2(self, t, output=None):
        '''Return the second derivative of the transformation for the given index(es)

           **Arguments:**

           t
                A number or an array of numbers for the indexes. t may be
                fractional or integer.
        '''
        if isinstance(t, numbers.Number):
            assert output is None
            return self._this.deriv2(t.astype(float))
        elif isinstance(t, np.ndarray):
            return self._deriv2_array(t, output)
        else:
            raise NotImplementedError

    def _deriv2_array(self, np.ndarray[double, ndim=1] t not None,
                            np.ndarray[double, ndim=1] output=None):
        assert t.flags['C_CONTIGUOUS']
        cdef int n = t.shape[0]
        if output is None:
            output = np.zeros(n, float)
        else:
            assert output.flags['C_CONTIGUOUS']
            assert output.shape[0] == n
        self._this.deriv2_array(&t[0], &output[0], n)
        return output

    def deriv3(self, t, output=None):
        '''Return the third of the transformation for the given index(es)

           **Arguments:**

           t
                A number or an array of numbers for the indexes. t may be
                fractional or integer.
        '''
        if isinstance(t, numbers.Number):
            assert output is None
            return self._this.deriv3(t.astype(float))
        elif isinstance(t, np.ndarray):
            return self._deriv3_array(t, output)
        else:
            raise NotImplementedError

    def _deriv3_array(self, np.ndarray[double, ndim=1] t not None,
                            np.ndarray[double, ndim=1] output=None):
        assert t.flags['C_CONTIGUOUS']
        cdef int n = t.shape[0]
        if output is None:
            output = np.zeros(n, float)
        else:
            assert output.flags['C_CONTIGUOUS']
            assert output.shape[0] == n
        self._this.deriv3_array(&t[0], &output[0], n)
        return output

    def inv(self, r, output=None):
        '''Return the indexes for given radial grid points

           **Arguments:**

           r
                A number or an array of numbers for the radial grid points.
        '''
        if isinstance(r, numbers.Number):
            assert output is None
            return self._this.inv(r)
        elif isinstance(r, np.ndarray):
            return self._inv_array(r.astype(float), output)
        else:
            raise NotImplementedError

    def _inv_array(self, np.ndarray[double, ndim=1] r not None,
                         np.ndarray[double, ndim=1] output=None):
        assert r.flags['C_CONTIGUOUS']
        cdef int n = r.shape[0]
        if output is None:
            output = np.zeros(n, float)
        else:
            assert output.flags['C_CONTIGUOUS']
            assert output.shape[0] == n
        self._this.inv_array(&r[0], &output[0], n)
        return output

    def get_radii(self):
        '''Return an array with radii'''
        result = np.arange(self.npoint, dtype=float)
        self._radius_array(result, result)
        return result

    def get_deriv(self):
        '''Return an array with derivatives at the grid points'''
        result = np.arange(self.npoint, dtype=float)
        self._deriv_array(result, result)
        return result

    def get_deriv2(self):
        '''Return an array with second derivatives at the grid points'''
        result = np.arange(self.npoint, dtype=float)
        self._deriv2_array(result, result)
        return result

    def get_deriv3(self):
        '''Return an array with third derivatives at the grid points'''
        result = np.arange(self.npoint, dtype=float)
        self._deriv3_array(result, result)
        return result

    @classmethod
    def from_string(cls, s):
        '''Construct a RTransform subclass from a string.'''
        words = s.split()
        clsname = words[0]
        args = words[1:]
        if clsname == 'IdentityRTransform':
            if len(args) != 1:
                raise ValueError('The IdentityRTransform needs one argument, got %i.' % len(words))
            npoint = int(args[0])
            return IdentityRTransform(npoint)
        elif clsname == 'LinearRTransform':
            if len(args) != 3:
                raise ValueError('The LinearRTransform needs three arguments, got %i.' % len(words))
            rmin = float(args[0])
            rmax = float(args[1])
            npoint = int(args[2])
            return LinearRTransform(rmin, rmax, npoint)
        elif clsname == 'ExpRTransform':
            if len(args) != 3:
                raise ValueError('The ExpRTransform needs three arguments, got %i.' % len(words))
            rmin = float(args[0])
            rmax = float(args[1])
            npoint = int(args[2])
            return ExpRTransform(rmin, rmax, npoint)
        elif clsname == 'PowerRTransform':
            if len(args) != 3:
                raise ValueError('The PowerRTransform needs three arguments, got %i.' % len(words))
            rmin = float(args[0])
            rmax = float(args[1])
            npoint = int(args[2])
            return PowerRTransform(rmin, rmax, npoint)
        elif clsname == 'HyperbolicRTransform':
            if len(args) != 3:
                raise ValueError('The HyperbolicRTransform needs three arguments, got %i.' % len(words))
            a = float(args[0])
            b = float(args[1])
            npoint = int(args[2])
            return HyperbolicRTransform(a, b, npoint)
        else:
            raise TypeError('Unkown RTransform subclass: %s' % clsname)

    def to_string(self):
        '''Represent the rtransform object as a string'''
        raise NotImplementedError

    def chop(self, npoint):
        '''Return an rtransform with ``npoint`` number of grid points

           The remaining grid points are such that they coincide with those from
           the old rtransform.
        '''
        raise NotImplementedError

    def half(self):
        '''Return an rtransform with half the number of grid points

           The returned rtransform is such that old(2t+1) = new(t).
        '''
        raise NotImplementedError

    def get_default_int1d(self):
        '''Return the recommended 1D integrator for this rtransform'''
        from horton.grid.int1d import SimpsonIntegrator1D
        return SimpsonIntegrator1D()



cdef class IdentityRTransform(RTransform):
    '''For testing only'''
    def __cinit__(self, int npoint):
        self._this = <rtransform.RTransform*>(new rtransform.IdentityRTransform(npoint))

    def to_string(self):
        return ' '.join(['IdentityRTransform', repr(self.npoint)])

    def chop(self, npoint):
        return IdentityRTransform(npoint)


cdef class LinearRTransform(RTransform):
    '''A linear grid.

       The grid points are distributed as follows:

       .. math:: r_i = \\alpha i + r_0

       with

       .. math:: \\alpha = (r_{N-1} -r_0)/(N-1).
    '''
    def __cinit__(self, double rmin, double rmax, int npoint):
        self._this = <rtransform.RTransform*>(new rtransform.LinearRTransform(rmin, rmax, npoint))

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

    def chop(self, npoint):
        rmax = self.radius(npoint-1)
        return LinearRTransform(self.rmin, rmax, npoint)

    def half(self):
        if self.npoint %2 != 0:
            raise ValueError('Half method can only be called on a rtransform with an even number of points.')
        rmin = self.radius(1)
        return LinearRTransform(rmin, self.rmax, self.npoint/2)


cdef class ExpRTransform(RTransform):
    r'''An exponential grid.

       The grid points are distributed as follows:

       .. math:: r_i = r_0 \exp(\alpha i)

       with

       .. math:: \alpha = \log(r_{N-1}/r_0)/(N-1).
    '''
    def __cinit__(self, double rmin, double rmax, int npoint):
        self._this = <rtransform.RTransform*>(new rtransform.ExpRTransform(rmin, rmax, npoint))

    property rmin:
        def __get__(self):
            return (<rtransform.ExpRTransform*>self._this).get_rmin()

    property rmax:
        def __get__(self):
            return (<rtransform.ExpRTransform*>self._this).get_rmax()

    property alpha:
        def __get__(self):
            return (<rtransform.ExpRTransform*>self._this).get_alpha()

    def to_string(self):
        return ' '.join(['ExpRTransform', repr(self.rmin), repr(self.rmax), repr(self.npoint)])

    def chop(self, npoint):
        rmax = self.radius(npoint-1)
        return ExpRTransform(self.rmin, rmax, npoint)

    def half(self):
        if self.npoint %2 != 0:
            raise ValueError('Half method can only be called on a rtransform with an even number of points.')
        rmin = self.radius(1)
        return ExpRTransform(rmin, self.rmax, self.npoint/2)


cdef class PowerRTransform(RTransform):
    r'''A power grid.

       The grid points are distributed as follows:

       .. math:: r_i = r_0 i^{\alpha}

       with

       .. math::
            \alpha = \frac{\ln r_{N-1} - \ln r_0}{\ln N-1}
    '''
    def __cinit__(self, double rmin, double rmax, int npoint):
        self._this = <rtransform.RTransform*>(new rtransform.PowerRTransform(rmin, rmax, npoint))

    property rmin:
        def __get__(self):
            return (<rtransform.PowerRTransform*>self._this).get_rmin()

    property rmax:
        def __get__(self):
            return (<rtransform.PowerRTransform*>self._this).get_rmax()

    property power:
        def __get__(self):
            return (<rtransform.PowerRTransform*>self._this).get_power()

    def to_string(self):
        return ' '.join(['PowerRTransform', repr(self.rmin), repr(self.rmax), repr(self.npoint)])

    def chop(self, npoint):
        rmax = self.radius(npoint-1)
        return PowerRTransform(self.rmin, rmax, npoint)

    def half(self):
        if self.npoint %2 != 0:
            raise ValueError('Half method can only be called on a rtransform with an even number of points.')
        rmin = self.radius(1)
        return PowerRTransform(rmin, self.rmax, self.npoint/2)

    def get_default_int1d(self):
        from horton.grid.int1d import StubIntegrator1D
        return StubIntegrator1D()


cdef class HyperbolicRTransform(RTransform):
    r'''A Hyperbolic grid (as in the GPAW program).

       The grid points are distributed as follows:

       .. math:: r_i = \frac{ai}{1 - bi}
    '''
    def __cinit__(self, double a, double b, int npoint):
        self._this = <rtransform.RTransform*>(new rtransform.HyperbolicRTransform(a, b, npoint))

    property a:
        def __get__(self):
            return (<rtransform.HyperbolicRTransform*>self._this).get_a()

    property b:
        def __get__(self):
            return (<rtransform.HyperbolicRTransform*>self._this).get_b()

    def to_string(self):
        return ' '.join(['HyperbolicRTransform', repr(self.a), repr(self.b), repr(self.npoint)])

    def get_default_int1d(self):
        from horton.grid.int1d import StubIntegrator1D
        return StubIntegrator1D()


#
# uniform
#


cdef class UniformGrid(object):
    def __cinit__(self, np.ndarray[double, ndim=1] origin not None,
                  np.ndarray[double, ndim=2] grid_rvecs not None,
                  np.ndarray[long, ndim=1] shape not None,
                  np.ndarray[long, ndim=1] pbc not None):
        assert origin.flags['C_CONTIGUOUS']
        assert origin.shape[0] == 3
        assert grid_rvecs.flags['C_CONTIGUOUS']
        assert grid_rvecs.shape[0] == 3
        assert grid_rvecs.shape[1] == 3
        assert shape.flags['C_CONTIGUOUS']
        assert shape.shape[0] == 3
        assert pbc.flags['C_CONTIGUOUS']
        assert pbc.shape[0] == 3
        self._this = <uniform.UniformGrid*>(new uniform.UniformGrid(
            &origin[0], &grid_rvecs[0, 0], &shape[0], &pbc[0],
        ))

    def __init__(self, np.ndarray[double, ndim=1] origin not None,
                 np.ndarray[double, ndim=2] grid_rvecs not None,
                 np.ndarray[long, ndim=1] shape not None,
                 np.ndarray[long, ndim=1] pbc not None):
        """__init__(self, ndarray origin, ndarray grid_rvecs, ndarray shape, ndarray pbc)
        Initialize a UniformGrid instance.

        Parameters
        ----------
        origin : np.ndarray[double, ndim=1]
                 The origin of the uniform grid, where the first grid point is located.
                 shape=(3,)
        grid_rvecs : np.ndarray[double, ndim=2]
                     The rows are real-space basis vectors that define the spacings.
                     between the grids. shape=(3,3)
        shape : np.ndarray[long, ndim=1]
                The shape of the uniform grid, i.e. number of points along each basis.
        pbc : np.ndarray[long, ndim=1]
              Three flags (0 or 1) for periodicity along each row in grid_rvecs.
        """

    def __dealloc__(self):
        del self._this

    property origin:
        def __get__(self):
            cdef np.npy_intp dims[1]
            dims[0] = 3
            cdef np.ndarray result = np.PyArray_SimpleNewFromData(1, dims, np.NPY_DOUBLE, self._this.origin)
            np.set_array_base(result, self)
            return result

    property grid_rvecs:
        def __get__(self):
            cdef np.npy_intp dims[2]
            dims[0] = 3
            dims[1] = 3
            cdef np.ndarray result = np.PyArray_SimpleNewFromData(2, dims, np.NPY_DOUBLE, self._this.grid_rvecs)
            np.set_array_base(result, self)
            return result

    property shape:
        def __get__(self):
            cdef np.npy_intp dims[1]
            dims[0] = 3
            cdef np.ndarray result = np.PyArray_SimpleNewFromData(1, dims, np.NPY_LONG, self._this.shape)
            np.set_array_base(result, self)
            return result

    property pbc:
        def __get__(self):
            cdef np.npy_intp dims[1]
            dims[0] = 3
            cdef np.ndarray result = np.PyArray_SimpleNewFromData(1, dims, np.NPY_LONG, self._this.pbc)
            np.set_array_base(result, self)
            return result

    property size:
        def __get__(self):
            return np.product(self.shape)

    def get_cell(self):
        cdef horton.cext.Cell result = horton.cext.Cell.__new__(horton.cext.Cell, initvoid=True)
        result._this = self._this.get_cell()
        return result

    def get_grid_cell(self):
        cdef horton.cext.Cell result = horton.cext.Cell.__new__(horton.cext.Cell, initvoid=True)
        result._this = self._this.get_grid_cell()
        return result

    @classmethod
    def from_hdf5(cls, grp):
        return cls(
            grp['origin'][:],
            grp['grid_rvecs'][:],
            grp['shape'][:],
            grp['pbc'][:],
        )

    def to_hdf5(self, grp):
        grp['grid_rvecs'] = self.grid_rvecs
        grp['origin'] = self.origin
        grp['shape'] = self.shape
        grp['pbc'] = self.pbc

    def zeros(self):
        return np.zeros(self.shape, float)

    def integrate(self, *args):
        '''Integrate the product of all arguments

           **Arguments:**

           data1, data2, ...
                All arguments must be arrays with the same size as the number
                of grid points. The arrays contain the functions, evaluated
                at the grid points, that must be multiplied and integrated.

        '''
        # This is often convenient for cube grid data:
        args = [arg.ravel() for arg in args if arg is not None]
        # Similar to conventional integration routine:
        return dot_multi(*args)*abs(np.linalg.det(self.grid_rvecs))

    def get_ranges_rcut(self, np.ndarray[double, ndim=1] center not None, double rcut):
        '''Return the ranges if indexes that lie within the cutoff sphere.

           **Arguments:**

           center
                The center of the cutoff sphere

           rcut
                The radius of the cutoff sphere

           The ranges are trimmed to avoid points that fall of non-periodic
           boundaries of the grid.
        '''
        assert center.flags['C_CONTIGUOUS']
        assert center.size == 3
        assert rcut >= 0

        cdef np.ndarray[long, ndim=1] ranges_begin = np.zeros(3, int)
        cdef np.ndarray[long, ndim=1] ranges_end = np.zeros(3, int)
        self._this.set_ranges_rcut(
            &center[0], rcut, &ranges_begin[0],
            &ranges_end[0])
        return ranges_begin, ranges_end

    def dist_grid_point(self, np.ndarray[double, ndim=1] center not None,
                        np.ndarray[long, ndim=1] indexes not None):
        '''Return the distance between a center and a grid point

           **Arguments:**

           center
                The center

           indexes
                The integer indexes of the grid point (may fall outside of shape)
        '''
        assert center.flags['C_CONTIGUOUS']
        assert center.size == 3
        assert indexes.flags['C_CONTIGUOUS']
        assert indexes.size == 3
        return self._this.dist_grid_point(&center[0], &indexes[0])

    def delta_grid_point(self, np.ndarray[double, ndim=1] center not None,
                         np.ndarray[long, ndim=1] indexes not None):
        '''Return the vector **from** a center **to** a grid point

           **Arguments:**

           center
                The center

           indexes
                The integer indexes of the grid point (may fall outside of shape)
        '''
        assert center.flags['C_CONTIGUOUS']
        assert center.size == 3
        assert indexes.flags['C_CONTIGUOUS']
        assert indexes.size == 3
        cdef np.ndarray[double, ndim=1] result = center.copy()
        self._this.delta_grid_point(&result[0], &indexes[0])
        return result


def index_wrap(long i, long high):
    return uniform.index_wrap(i, high)


#
# utils
#


def _check_integranda(integranda, npoint=None):
    assert len(integranda) > 0
    for integrandum in integranda:
        assert integrandum.flags['C_CONTIGUOUS']
        if npoint is None:
            npoint = integrandum.shape[0]
        else:
            assert npoint == integrandum.shape[0]
    return npoint


cdef double** _parse_integranda(integranda):
    cdef double** result = <double **>malloc(len(integranda)*sizeof(double*))
    if result == NULL:
        raise MemoryError()
    cdef np.ndarray[double, ndim=1] integrandum
    for i in xrange(len(integranda)):
        integrandum = integranda[i]
        result[i] = &integrandum[0]
    return result


def _parse_segments(segments, npoint):
    if segments is None:
        segments = np.array([npoint])
    assert segments.flags['C_CONTIGUOUS']
    nsegment = segments.shape[0]
    return segments, nsegment


def dot_multi(*integranda, np.ndarray[long, ndim=1] segments=None):
    '''Multiply the arguments piecewise and sum up the products

       **Arguments:**

       data1, data2, ...
            Arrays of the same size, whose elements will be multiplied piecewise
            and then added.

       **Optional arguments:**

       segments
            An array with segment sizes (integer). If given, the summation is
            carried out in segments of the given sizes and the return value is
            an array with the same size as segments.
    '''
    cdef long npoint = _check_integranda(integranda)
    segments, nsegment = _parse_segments(segments, npoint)
    cdef np.ndarray[double, ndim=1] output = np.zeros(nsegment)
    cdef double** pointers = _parse_integranda(integranda)
    try:
        utils.dot_multi(npoint, len(integranda), pointers, &segments[0], &output[0])
    finally:
        free(pointers)
    if nsegment == 1:
        return output[0]
    else:
        return output


def _get_nmoment(long lmax, long mtype):
    if mtype==1:
        # cartesian moments
        return ((lmax+1)*(lmax+2)*(lmax+3))/6
    elif mtype==2:
        # pure moments
        return (lmax+1)**2
    elif mtype==3:
        # radial moments
        return lmax+1
    elif mtype==4:
        # surface moments
        return (lmax+1)**2
    else:
        raise ValueError('Unsupported mtype.')



def dot_multi_moments_cube(integranda, UniformGrid ugrid not None,
                           np.ndarray[double, ndim=1] center not None,
                           long lmax, long mtype):
    '''Multiply the arguments piecewise, including one of a series of multipole functions at a time, and sum up the products.

       **Arguments:**

       data1, data2, ...
            Arrays of the same size, whose elements will be multiplied piecewise
            and then added.

       center
            The origin for the multipole functions

       lmax
            The maximum angular momentum for the moments

       mtype
            The type of moments: 1=``cartesian``, 2=``pure``,
            3=``radial``, 4=``surface``.

       **Returns:** an array where the number of elements matches the number of
       multipole moments for the given combiantion of lmax and mtype.
    '''


    cdef long npoint = _check_integranda(integranda, ugrid.size)
    # Only non-periodic grids are supported to guarantee an unambiguous definition of the polynomial.
    assert ugrid.pbc[0] == 0
    assert ugrid.pbc[1] == 0
    assert ugrid.pbc[2] == 0
    #
    assert center.flags['C_CONTIGUOUS']
    assert center.shape[0] == 3

    cdef long nmoment = _get_nmoment(lmax, mtype)
    cdef np.ndarray[double, ndim=1] output = np.zeros(nmoment)
    cdef double** pointers = _parse_integranda(integranda)
    try:
        utils.dot_multi_moments_cube(len(integranda), pointers, ugrid._this, &center[0], lmax, mtype, &output[0], nmoment)
    finally:
        free(pointers)
    return output


def dot_multi_moments(integranda,
                      np.ndarray[double, ndim=2] points not None,
                      np.ndarray[double, ndim=1] center not None,
                      long lmax, long mtype,
                      np.ndarray[long, ndim=1] segments):
    '''Multiply the arguments piecewise, including one of a series of multipole functions at a time, and sum up the products.

       **Arguments:**

       data1, data2, ...
            Arrays of the same size, whose elements will be multiplied piecewise
            and then added.

       center
            The origin for the multipole functions

       lmax
            The maximum angular momentum for the moments

       mtype
            The type of moments: 1=``cartesian``, 2=``pure``,
            3=``radial``, 4=``surface``.

       segments
            An array with segment sizes (integer). If given, the summation is
            carried out in segments of the given sizes.

       **Returns:** an array where the number of elements matches the number of
       multipole moments for the given combiantion of lmax and mtype. If
       ``segments`` is given, the return array has two indices, the first one
       running over the segments and the second one running over the moments.
    '''

    assert points.flags['C_CONTIGUOUS']
    assert points.shape[1] == 3
    cdef long npoint = _check_integranda(integranda, points.shape[0])
    #
    assert center.flags['C_CONTIGUOUS']
    assert center.shape[0] == 3

    cdef long nmoment = _get_nmoment(lmax, mtype)
    segments, nsegment = _parse_segments(segments, npoint)
    cdef np.ndarray[double, ndim=2] output = np.zeros((nsegment, nmoment))
    cdef double** pointers = _parse_integranda(integranda)
    try:
        utils.dot_multi_moments(npoint, len(integranda), pointers,
            &points[0, 0], &center[0], lmax, mtype,
            &segments[0], &output[0, 0], nmoment)
    finally:
        free(pointers)
    if nsegment == 1:
        return output[0]
    else:
        return output
