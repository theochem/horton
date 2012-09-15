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


__all__ = ['BaseRTransform', 'LogRTransform']


class BaseRTransform(object):
    def __init__(self, int1d):
        '''
           **Arguments:**

           int1d
                A one-dimensional integration algorithm, instance of a subclass
                of BaseIntegrator1D.
        '''
        self.int1d = int1d

    def get_radii(self, npoint, shift=0):
        '''Return an array of radii with the given size

           **Optional argument:**

           shift
                Displacement to be applied to the uniform reference grid.
        '''
        raise NotImplementedError

    def get_int_weights(self, npoint):
        '''Return an array with radial integration weights'''
        return self.int1d.get_weights(npoint)*self.get_volume_elements(npoint)

    def get_volume_elements(self, npoint):
        '''Return an array with volume elements associated with the transform'''
        raise NotImplementedError


class LogRTransform(BaseRTransform):
    '''A plain logarithmic grid.

       The grid points are distributed as follows:

       .. math:: r_i = r_0 \\alpha^i

       with

       .. math:: alpha = \log(r_{N-1}/r_0)/(N-1).
    '''

    def __init__(self, int1d, rmin, rmax, npoint):
        '''
           **Arguments:**

           int1d
                A one-dimensional integration algorithm, instance of a subclass
                of BaseIntegrator1D.

           rmin
                The first radial grid point.

           rmax
                The last radial grid point.

           npoint
                The number of radial points.

           Note that ``rmax`` and ``npoint`` are not strictly enforced. They are
           is used as convenient parameters to determine ``alpha``.
        '''
        if npoint < 2:
            raise ValueError('At least two points are needed')
        self._rmin = rmin
        self._rmax = rmax
        self._npoint = npoint
        self._alpha = np.log(rmax/rmin)/(npoint-1)
        BaseRTransform.__init__(self, int1d)

    def _get_rmin(self):
        '''The first grid point.'''
        return self._rmin

    rmin = property(_get_rmin)

    def _get_rmax(self):
        '''The `last` grid point.'''
        return self._rmax

    rmax = property(_get_rmax)

    def _get_npoint(self):
        '''The `number of grid points`.'''
        return self._npoint

    npoint = property(_get_npoint)

    def _get_alpha(self):
        '''The alpha parameter of the grid.'''
        return self._alpha

    alpha = property(_get_alpha)

    def get_radii(self, npoint, shift=0):
        '''Return an array of radii with the given size

           **Optional argument:**

           shift
                Displacement to be applied to the uniform reference grid.
        '''
        return self._rmin*np.exp((np.arange(npoint)+shift)*self._alpha)

    def get_volume_elements(self, npoint):
        '''Return an array with volume elements associated with the transform'''
        return self._rmin*np.exp(np.arange(npoint)*self._alpha)*self._alpha
