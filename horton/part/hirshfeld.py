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

from horton.cache import just_once
from horton.log import log
from horton.part.stockholder import StockholderCPart
from horton.part.base import DPart


__all__ = ['HirshfeldDPart', 'HirshfeldCPart']


# TODO: isolate common code in mixin class
# TODO: proofread and add tests for pseudo densities


class HirshfeldDPart(DPart):
    name = 'h'
    options = ['local']

    '''Base class for Hirshfeld partitioning'''
    def __init__(self, molgrid, proatomdb, local=True):
        self._proatomdb = proatomdb
        DPart.__init__(self, molgrid, local)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def _init_log(self):
        DPart._init_log(self)
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld'),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('hirshfeld1977', 'the use of Hirshfeld partitioning')

    def _compute_at_weights(self, i0, grid, at_weights):
        promol, new = self.cache.load('promol', grid.size, alloc=grid.size)
        if new or self.local:
            # In case of local grids, the pro-molecule must always be recomputed.
            promol[:] = 0.0 # needed if not new and local.
            for i1 in xrange(self.system.natom):
                proatom_fn = self.cache.load('proatom_fn', i1)
                if i1 == i0:
                    at_weights[:] = 0.0
                    grid.eval_spline(proatom_fn, self.system.coordinates[i1], at_weights)
                    promol += at_weights
                else:
                    grid.eval_spline(proatom_fn, self.system.coordinates[i1], promol)
            # The following seems worse than it is. It does nothing to the
            # relevant numbers. It just avoids troubles in the division.
            promol[:] += 1e-100
        else:
            # In case of a global grid and when the pro-molecule is up to date,
            # only the pro-atom needs to be recomputed.
            proatom_fn = self.cache.load('proatom_fn', i0)
            at_weights[:] = 0.0
            grid.eval_spline(proatom_fn, self.system.coordinates[i0], at_weights)
        # Finally compute the ratio
        at_weights[:] /= promol

    def _at_weights_cleanup(self):
        # Get rid of cached work arrays
        for i, grid in self.iter_grids():
            self.cache.discard('at_weights_work', id(grid.size))

    @just_once
    def _init_at_weights(self):
        # TODO. Splines should not be stored in the cache. Use same structure as
        # in CPart: get_proatom_spline and _*_propars methods.
        # Create the splines and store them in the cache
        for i in xrange(self.system.natom):
            proatom_fn = self.cache.load('proatom_fn', i, default=None)
            if proatom_fn is None:
                n = self.system.numbers[i]
                proatom_fn = self._proatomdb.get_spline(n)
                self.cache.dump('proatom_fn', i, proatom_fn)

        # Compute the weights
        for i0, grid in self.iter_grids():
            at_weights, new = self.cache.load('at_weights', i0, alloc=grid.size)
            if new:
                self._compute_at_weights(i0, grid, at_weights)

        self._at_weights_cleanup()


class HirshfeldCPart(StockholderCPart):
    name = 'h'

    def __init__(self, system, ui_grid, moldens, proatomdb, store, smooth=False):
        '''
           See CPart base class for the description of the arguments.
        '''
        self._proatomdb = proatomdb
        StockholderCPart.__init__(self, system, ui_grid, moldens, store, smooth)

        # Check of the same (or similar) psuedo potentials were used for the
        # system and the proatoms.
        for i in xrange(self._system.natom):
            number = self._system.numbers[i]
            pseudo_number = self._system.pseudo_numbers[i]
            pseudo_number_expected = self._proatomdb.get_record(number, 0).pseudo_number
            if pseudo_number_expected != pseudo_number:
                raise ValueError('The pseudo number of atom %i does not match with the proatom database (%i!=%i)' % (
                    i, pseudo_number, pseudo_number_expected))

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def _init_weight_corrections(self):
        funcs = []
        for i in xrange(self._system.natom):
            number = self._system.numbers[i]
            funcs.append((
                self._system.coordinates[i],
                [self._proatomdb.get_spline(number)],
            ))
        wcor = self._ui_grid.compute_weight_corrections(funcs)
        self._cache.dump('wcor', wcor)

    def _init_partitioning(self):
        work = self._ui_grid.zeros()
        self._update_promolecule(work)
        self._store_at_weights(work)

    def _update_promolecule(self, work):
        promoldens = self._cache.load('promoldens', alloc=self._ui_grid.shape)[0]
        promoldens[:] = 0.0
        for index in xrange(self._system.natom):
            self.compute_proatom(index, work)
            promoldens += work
            # Merely for efficiency:
            self._store.dump(work, 'at_weights', index)

    def _store_at_weights(self, work):
        # Compute the atomic weight functions if this is useful. This is merely
        # a matter of efficiency.
        promoldens = self._cache.load('promoldens')
        for index in xrange(self._system.natom):
            if ('at_weights', index) in self._store:
                self._store.load(work, 'at_weights', index)
                work /= promoldens
                self._store.dump(work, 'at_weights', index)

    def get_cutoff_radius(self, index):
        '''The radius at which the weight function goes to zero'''
        rtf = self._proatomdb.get_rtransform(self._system.numbers[index])
        return rtf.radius(rtf.npoint-1)

    def get_proatom_spline(self, index):
        return self._proatomdb.get_spline(self._system.numbers[index])

    def do_all(self):
        names = StockholderCPart.do_all(self)
        self.do_dispersion()
        return names + ['volumes', 'volume_ratios', 'c6s']

    @just_once
    def do_dispersion(self):
        # Method by Alexandre Tkatchenko and Matthias Scheffler:
        #   PhysRevLett-v102-p073005-y2009.pdf
        #   http://dx.doi.org/10.1103/PhysRevLett.102.073005
        # Reference C6 values taken from X. Chu and A. Dalgarno:
        #   JChemPhys-v121-p4083-y2004.pdf
        #   http://dx.doi.org/10.1063/1.1779576
        #   (Corrected values from tables VII to XII)
        # C6 value for Hydrogen taken from Zong-Chao Yan, James F. Babb, A. Dalgarno and G. W. F. Drake
        #   PhysRevA-v54-p2824-y1996.pdf
        #   http://dx.doi.org/10.1103/PhysRevA.54.2824
        ref_c6s = { # reference C6 values in atomic units
            1: 6.499, 2: 1.42, 3: 1392.0, 4: 227.0, 5: 99.5, 6: 46.6, 7: 24.2,
            8: 15.6, 9: 9.52, 10: 6.20, 11: 1518.0, 12: 626.0, 13: 528.0, 14:
            305.0, 15: 185.0, 16: 134.0, 17: 94.6, 18: 64.2, 19: 3923.0, 20:
            2163.0, 21: 1383.0, 22: 1044.0, 23: 832.0, 24: 602.0, 25: 552.0, 26:
            482.0, 27: 408.0, 28: 373.0, 29: 253.0, 30: 284.0, 31: 498.0, 32:
            354.0, 33: 246.0, 34: 210.0, 35: 162.0, 36: 130.0, 37: 4769.0, 38:
            3175.0, 49: 779.0, 50: 659.0, 51: 492.0, 52: 445.0, 53: 385.0,
        }

        volumes, new_volumes = self._cache.load('volumes', alloc=self.system.natom)
        volume_ratios, new_volume_ratios = self._cache.load('volume_ratios', alloc=self.system.natom)
        c6s, new_c6s = self._cache.load('c6s', alloc=self.system.natom)

        if new_volumes or new_volume_ratios or new_c6s:
            self.do_populations()
            self.do_moments()
            k3 = 2
            radial_powers = self._cache.load('radial_powers')
            assert radial_powers[k3] == 3
            radial_moments = self._cache.load('radial_moments')
            populations = self._cache.load('populations')

            if log.do_medium:
                log('Computing atomic dispersion coefficients.')

            for i in xrange(self.system.natom):
                n = self.system.numbers[i]
                if n not in ref_c6s:
                    raise NotImplementedError('No reference C6 value available for atom number %i.' % n)
                volumes[i] = radial_moments[i,k3]/populations[i]
                ref_volume = self.proatomdb.get_record(n, 0).get_moment(3)/n
                volume_ratios[i] = volumes[i]/ref_volume
                c6s[i] = (volume_ratios[i])**2*ref_c6s[n]
