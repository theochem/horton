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

from horton.log import log
from horton.grid.cext import CubicSpline
from horton.part.base import WPart, CPart

__all__ = [
    'StockholderWPart', 'StockholderCPart',
]


class StockHolderMixin(object):
    def get_proatom_rho(self, index, *args, **kwargs):
        raise NotImplementedError

    def fix_proatom_rho(self, index, rho):
        '''Check if the radial density for the proatom is correct and fix as needed.

           **Arguments:**

           index
                The atom for which this proatom rho is created.

           rho
                The radial density
        '''
        # TODO: should not refer to proatomdb
        number = self.system.numbers[index]
        weights = self.proatomdb.get_radial_weights(number)

        # Check for negative parts
        original = np.dot(rho, weights)
        if rho.min() < 0:
            rho[rho<0] = 0.0
            error = np.dot(rho, weights) - original
            if log.do_medium:
                log('                    Pro-atom not positive everywhere. Lost %.1e electrons' % error)

        return rho

    def get_proatom_spline(self, index, *args, **kwargs):
        # TODO: should not refer to proatomdb
        # Get the radial density
        number = self.system.numbers[index]
        rho = self.get_proatom_rho(index, *args, **kwargs)

        # Double check and fix if needed
        rho = self.fix_proatom_rho(index, rho)

        # Make a spline
        return CubicSpline(rho, rtf=self.proatomdb.get_rtransform(number))

    def eval_spline(self, index, spline, output, grid=None, label='noname'):
        center = self.system.coordinates[index]
        if grid is None:
            grid = self.get_grid(index)
        if log.do_medium:
            log('  Evaluating spline (%s) for atom %i on %i grid points' % (label, index, grid.size))
        grid.eval_spline(spline, center, output)

    def eval_proatom(self, index, output, grid=None):
        spline = self.get_proatom_spline(index)
        output[:] = 0.0
        self.eval_spline(index, spline, output, grid, label='proatom')
        output += 1e-100

    def update_at_weights(self):
        # This will reconstruct the promolecular density and atomic weights
        # based on the current proatomic splines.
        promoldens = self.cache.load('promoldens', alloc=self.grid.shape)[0]
        promoldens[:] = 0

        # update the promolecule density and store the proatoms in the at_weights
        # arrays for later.
        for index in xrange(self.system.natom):
            grid = self.get_grid(index)
            at_weights = self.cache.load('at_weights', index, alloc=grid.shape)[0]
            self.update_pro(index, at_weights, promoldens)

        # Compute the atomic weights by taking the ratios between proatoms and
        # promolecules.
        for index in xrange(self.system.natom):
            at_weights = self.cache.load('at_weights', index)
            at_weights /= self.to_atomic_grid(index, promoldens)
            #np.clip(at_weights, 0, 1, at_weights) # TODO: would this help?

    def update_pro(self, index, proatdens, promoldens):
        raise NotImplementedError


class StockholderWPart(StockHolderMixin, WPart):
    def update_pro(self, index, proatdens, promoldens):
        work = self.grid.zeros()
        self.eval_proatom(index, work, self.grid)
        promoldens += work
        proatdens[:] = self.to_atomic_grid(index, work)


class StockholderCPart(StockHolderMixin, CPart):
    def update_pro(self, index, proatdens, promoldens):
        self.eval_proatom(index, proatdens)
        promoldens += self.to_sys_grid(index, proatdens)

    def to_sys_grid(self, index, data):
        if self.local:
            result = self.grid.zeros()
            grid = self.get_grid(index)
            grid.wrap(data, result)
            return result
        else:
            return data
