# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
'''Base classes for all stockholder partitioning schemes'''


import numpy as np

from horton.log import log
from horton.grid.cext import CubicSpline
from horton.part.base import WPart
from horton.grid.poisson import solve_poisson_becke


__all__ = [
    'StockholderWPart',
]


class StockHolderMixin(object):
    def get_rgrid(self, index):
        raise NotImplementedError

    def get_proatom_rho(self, index, *args, **kwargs):
        raise NotImplementedError

    def fix_proatom_rho(self, index, rho, deriv):
        '''Check if the radial density for the proatom is correct and fix as needed.

           **Arguments:**

           index
                The atom for which this proatom rho is created.

           rho
                The radial density

           deriv
                the derivative of the radial density or None.
        '''
        rgrid = self.get_rgrid(index)

        # Check for negative parts
        original = rgrid.integrate(rho)
        if rho.min() < 0:
            rho[rho<0] = 0.0
            deriv = None
            error = rgrid.integrate(rho) - original
            if log.do_medium:
                log('                    Pro-atom not positive everywhere. Lost %.1e electrons' % error)

        return rho, deriv

    def get_proatom_spline(self, index, *args, **kwargs):
        # Get the radial density
        rho, deriv = self.get_proatom_rho(index, *args, **kwargs)

        # Double check and fix if needed
        rho, deriv = self.fix_proatom_rho(index, rho, deriv)

        # Make a spline
        rtf = self.get_rgrid(index).rtransform
        return CubicSpline(rho, deriv, rtf)

    def eval_spline(self, index, spline, output, grid, label='noname'):
        center = self.coordinates[index]
        if log.do_debug:
            number = self.numbers[index]
            log('  Evaluating spline (%s) for atom %i (n=%i) on %i grid points' % (label, index, number, grid.size))
        grid.eval_spline(spline, center, output)

    def eval_proatom(self, index, output, grid):
        spline = self.get_proatom_spline(index)
        output[:] = 0.0
        self.eval_spline(index, spline, output, grid, label='proatom')
        output += 1e-100
        assert np.isfinite(output).all()

    def update_at_weights(self):
        # This will reconstruct the promolecular density and atomic weights
        # based on the current proatomic splines.
        promoldens = self.cache.load('promoldens', alloc=self.grid.shape)[0]
        promoldens[:] = 0

        # update the promolecule density and store the proatoms in the at_weights
        # arrays for later.
        for index in xrange(self.natom):
            grid = self.get_grid(index)
            at_weights = self.cache.load('at_weights', index, alloc=grid.shape)[0]
            self.update_pro(index, at_weights, promoldens)

        # Compute the atomic weights by taking the ratios between proatoms and
        # promolecules.
        for index in xrange(self.natom):
            at_weights = self.cache.load('at_weights', index)
            at_weights /= self.to_atomic_grid(index, promoldens)
            np.clip(at_weights, 0, 1, out=at_weights)

    def update_pro(self, index, proatdens, promoldens):
        raise NotImplementedError

    def do_prosplines(self):
        for index in xrange(self.natom):
            # density
            key = ('spline_prodensity', index)
            if key not in self.cache:
                if log.medium:
                    log('Storing proatom density spline for atom %i.' % index)
                spline = self.get_proatom_spline(index)
                self.cache.dump(key, spline, tags='o')
            # hartree potential
            key = ('spline_prohartree', index)
            if key not in self.cache:
                if log.medium:
                    log('Computing proatom hartree potential spline for atom %i.' % index)
                rho_spline = self.cache.load('spline_prodensity', index)
                v_spline = solve_poisson_becke([rho_spline])[0]
                self.cache.dump(key, v_spline, tags='o')


class StockholderWPart(StockHolderMixin, WPart):
    def update_pro(self, index, proatdens, promoldens):
        work = self.grid.zeros()
        self.eval_proatom(index, work, self.grid)
        promoldens += work
        proatdens[:] = self.to_atomic_grid(index, work)
