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

    def compute_at_weights(self, i0, at_weights):
        promoldens = self.get_promoldens(i0)
        spline = self.get_proatom_spline(i0)
        at_weights[:] = 0.0
        grid = self.get_grid(i0)
        grid.eval_spline(spline, self.system.coordinates[i0], at_weights)
        at_weights[:] /= promoldens
        #np.clip(at_weights, 0, 1, at_weights) # TODO: would this help?

    def get_promoldens(self, index):
        # Get/Construct the promolecule on the entire grid if needed
        promoldens, new = self.cache.load('promoldens', alloc=self.grid.shape)
        if new:
            for i in xrange(self.system.natom):
                spline = self.get_proatom_spline(i)
                self.grid.eval_spline(spline, self.system.coordinates[i], promoldens)
            # The following seems worse than it is. It does nothing to the
            # relevant numbers. It just avoids troubles in the division.
            promoldens[:] += 1e-100

        # Take the relevant part
        return self.to_atomic_grid(index, promoldens)


class StockholderWPart(StockHolderMixin, WPart):
    pass


class StockholderCPart(StockHolderMixin, CPart):
    pass
