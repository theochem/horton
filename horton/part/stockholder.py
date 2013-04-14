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
from horton.part.base import DPart, CPart

__all__ = [
    'StockholderDPart', 'StockholderCPart',
]


# TODO: reduce duplicate code


def check_proatom(spline, expected_pop):
    # TODO: complete this routine and use it whenever a pro-atom spline is generateds
    proradrho = spline.copy_y()

    target_pseudo_population = self.system.pseudo_numbers[index] - target_charge
    error = np.dot(proradrho, weights) - target_pseudo_population
    assert abs(error) < 1e-5

    # Check for negative parts
    if proradrho.min() < 0:
        proradrho[proradrho<0] = 0.0
        error = np.dot(proradrho, weights) - target_pseudo_population
        if log.do_medium:
            log('                    Pro-atom not positive everywhere. Lost %.5f electrons' % error)



class StockholderDPart(DPart):
    def compute_at_weights(self, i0, grid, at_weights):
        promol, new = self.cache.load('promol', grid.size, alloc=grid.size)
        if new or self.local:
            # In case of local grids, the pro-molecule must always be recomputed.
            promol[:] = 0.0 # needed if not new and local.
            for i1 in xrange(self.system.natom):
                spline = self.get_proatom_spline(i1)
                if i1 == i0:
                    at_weights[:] = 0.0
                    grid.eval_spline(spline, self.system.coordinates[i1], at_weights)
                    promol += at_weights
                else:
                    grid.eval_spline(spline, self.system.coordinates[i1], promol)
            # The following seems worse than it is. It does nothing to the
            # relevant numbers. It just avoids troubles in the division.
            promol[:] += 1e-100
        else:
            # In case of a global grid and when the pro-molecule is up to date,
            # only the pro-atom needs to be recomputed.
            spline = self.get_proatom_spline(i0)
            at_weights[:] = 0.0
            grid.eval_spline(spline, self.system.coordinates[i0], at_weights)
        # Finally compute the ratio
        at_weights[:] /= promol

    def get_proatom_spline(self, index):
        raise NotImplementedError


class StockholderCPart(CPart):
    def __init__(self, system, ui_grid, moldens, store, smooth=False):
        '''
           See CPart base class for the description of the arguments.
        '''
        CPart.__init__(self, system, ui_grid, moldens, store, smooth)
        assert self._cache.has('promoldens')

    def compute_spline(self, index, spline, output, name='noname', window=None):
        if log.do_medium:
            if window is None:
                shape = self._ui_grid.shape
            else:
                shape = window.shape
            log('Computing spline (%s) for atom  %i/%i on grid [%i %i %i]' % (name, index, self.system.natom-1, shape[0], shape[1], shape[2]))

        center = self._system.coordinates[index]
        output[:] = 0.0
        if window is None:
            self._ui_grid.eval_spline(spline, center, output)
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

    def get_proatom_spline(self, index):
        raise NotImplementedError
