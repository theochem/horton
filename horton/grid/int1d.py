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
'''1D integration algorithms'''


import numpy as np
from horton.grid.cext import compute_cubic_spline_int_weights


__all__ = [
    'Integrator1D', 'TrapezoidIntegrator1D', 'CubicIntegrator1D', 'SimpsonIntegrator1D'
]


class Integrator1D(object):
    npoint_min = None
    '''Base class for integration algorithms'''
    def get_weights(self, npoint):
        '''Return integration weights for linear grid.'''
        raise NotImplementedError


class TrapezoidIntegrator1D(Integrator1D):
    npoint_min = 2
    '''Trapezoid rule integration algorithm'''
    def get_weights(self, npoint):
        '''Return integration weights for linear grid.'''
        assert npoint >= self.npoint_min
        result = np.ones(npoint, dtype=float)
        result[0] = 0.5
        result[-1] = 0.5
        return result


class CubicIntegrator1D(Integrator1D):
    npoint_min = 2
    '''Cubic spline integration algorithm'''
    def get_weights(self, npoint):
        '''Return integration weights for linear grid.'''
        assert npoint >= self.npoint_min
        result = np.ones(npoint, dtype=float)
        compute_cubic_spline_int_weights(result)
        return result


class SimpsonIntegrator1D(Integrator1D):
    npoint_min = 8
    '''Composite Simpson's rule integration algorithm'''
    def get_weights(self, npoint):
        '''Return integration weights for linear grid.'''
        assert npoint >= self.npoint_min
        result = np.ones(npoint, dtype=float)
        result[0] = 17.0/48.0
        result[-1] = 17.0/48.0
        result[1] = 59.0/48.0
        result[-2] = 59.0/48.0
        result[2] = 43.0/48.0
        result[-3] = 43.0/48.0
        result[3] = 49.0/48.0
        result[-4] = 49.0/48.0
        return result
