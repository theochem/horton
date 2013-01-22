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
from horton.grid.cext import compute_cubic_spline_int_weights


__all__ = [
    'Integrator1D', 'TrapezoidIntegrator1D', 'CubicIntegrator1D', 'SimpsonIntegrator1D'
]


class Integrator1D(object):
    '''Base class for integration algorithms'''
    def get_weights(self, npoint):
        '''Return integration weights for linear grid.'''
        raise NotImplementedError


class TrapezoidIntegrator1D(Integrator1D):
    '''Trapezoid rule integration algorithm'''
    def get_weights(self, npoint):
        '''Return integration weights for linear grid.'''
        assert npoint >= 2
        result = np.ones(npoint, dtype=float)
        result[0] = 0.5
        result[-1] = 0.5
        return result


class CubicIntegrator1D(Integrator1D):
    '''Cubic spline integration algorithm'''
    def get_weights(self, npoint):
        '''Return integration weights for linear grid.'''
        assert npoint >= 2
        result = np.ones(npoint, dtype=float)
        compute_cubic_spline_int_weights(result)
        return result


class SimpsonIntegrator1D(Integrator1D):
    '''Composite Simpson's rule integration algorithm'''
    def get_weights(self, npoint):
        '''Return integration weights for linear grid.'''
        assert npoint >= 8
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
