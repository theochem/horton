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
            if log.do_high:
                log('                    Pro-atom not positive everywhere. Lost %.5f electrons' % error)

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


class StockholderWPart(StockHolderMixin, WPart):
    def compute_at_weights(self, i0, at_weights):
        grid = self.get_grid(i0)
        promoldens, new = self.cache.load('promoldens', grid.size, alloc=grid.size)
        if new or self.local:
            # In case of local grids, the pro-molecule must always be recomputed.
            promoldens[:] = 0.0 # needed if not new and local.
            for i1 in xrange(self.system.natom):
                spline = self.get_proatom_spline(i1)
                if i1 == i0:
                    at_weights[:] = 0.0
                    grid.eval_spline(spline, self.system.coordinates[i1], at_weights)
                    promoldens += at_weights
                else:
                    grid.eval_spline(spline, self.system.coordinates[i1], promoldens)
            # The following seems worse than it is. It does nothing to the
            # relevant numbers. It just avoids troubles in the division.
            promoldens[:] += 1e-100
        else:
            # In case of a global grid and when the pro-molecule is up to date,
            # only the pro-atom needs to be recomputed.
            spline = self.get_proatom_spline(i0)
            at_weights[:] = 0.0
            grid.eval_spline(spline, self.system.coordinates[i0], at_weights)
        # Finally compute the ratio
        at_weights[:] /= promoldens


class StockholderCPart(StockHolderMixin, CPart):
    def __init__(self, system, grid, moldens, store, wcor_numbers, wcor_rcut_max=2.0, wcor_rcond=0.1):
        '''
           See CPart base class for the description of the arguments.
        '''
        CPart.__init__(self, system, grid, moldens, store, wcor_numbers, wcor_rcut_max, wcor_rcond)
        assert self._cache.has('promoldens')

    def compute_spline(self, index, spline, output, name='noname', window=None):
        if log.do_medium:
            if window is None:
                shape = self.grid.shape
            else:
                shape = window.shape
            log('Computing spline (%s) for atom  %i/%i on grid [%i %i %i]' % (name, index, self.system.natom-1, shape[0], shape[1], shape[2]))

        center = self._system.coordinates[index]
        output[:] = 0.0
        if window is None:
            self.grid.eval_spline(spline, center, output)
        else:
            window.eval_spline(spline, center, output)

    def compute_proatom(self, index, output, window=None):
        # This routine must compute the pro-atom, not from scratch, but based
        # on the info available after the partitioning. The default behavior
        # is to load the pro-atom from the store, but in case the store is fake,
        # this implementation computes the pro-atom on the fly.
        spline = self.get_proatom_spline(index)
        self.compute_spline(index, spline, output, 'proatom', window)
        output += 1e-100

    def compute_at_weights(self, index, output, window=None):
        self.compute_proatom(index, output, window)
        if window is None:
            output /= self._cache.load('promoldens')
        else:
            promoldens = window.zeros()
            window.extend(self._cache.load('promoldens'), promoldens)
            output /= promoldens
