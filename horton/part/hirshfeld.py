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
from horton.part.stockholder import StockholderDPart, StockholderCPart


__all__ = ['HirshfeldDPart', 'HirshfeldCPart']


def check_proatomdb(system, proatomdb):
    # Check if the same pseudo numbers (effective core charges) are used for the
    # system and the proatoms.
    for i in xrange(system.natom):
        number = system.numbers[i]
        pseudo_number = system.pseudo_numbers[i]
        pseudo_number_expected = proatomdb.get_record(number, 0).pseudo_number
        if pseudo_number_expected != pseudo_number:
            raise ValueError('The pseudo number of atom %i does not match with the proatom database (%i!=%i)' % (
                i, pseudo_number, pseudo_number_expected))


class HirshfeldMixin(object):
    def __init__(self, proatomdb):
        self._proatomdb = proatomdb

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def get_proatom_rho(self, index):
        return self.proatomdb.get_rho(self._system.numbers[index])

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


class HirshfeldDPart(HirshfeldMixin, StockholderDPart):
    name = 'h'
    options = ['local']

    '''Base class for Hirshfeld partitioning'''
    def __init__(self, system, grid, proatomdb, local=True):
        check_proatomdb(system, proatomdb)
        HirshfeldMixin. __init__(self, proatomdb)
        StockholderDPart.__init__(self, system, grid, local)

    def _init_log(self):
        StockholderDPart._init_log(self)
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld'),
                ('Proatomic DB',  self.proatomdb),
            ])
            log.cite('hirshfeld1977', 'the use of Hirshfeld partitioning')

    @just_once
    def _init_partitioning(self):
        # Compute the weights
        for i0 in xrange(self.natom):
            grid = self.get_grid(i0)
            at_weights, new = self.cache.load('at_weights', i0, alloc=grid.size)
            if new:
                self.compute_at_weights(i0, at_weights)

    def do_all(self):
        names = StockholderDPart.do_all(self)
        self.do_dispersion()
        return names + ['volumes', 'volume_ratios', 'c6s']


class HirshfeldCPart(HirshfeldMixin, StockholderCPart):
    name = 'h'

    def __init__(self, system, grid, moldens, proatomdb, store, smooth=False):
        '''
           See CPart base class for the description of the arguments.
        '''
        check_proatomdb(system, proatomdb)
        HirshfeldMixin. __init__(self, proatomdb)
        StockholderCPart.__init__(self, system, grid, moldens, store, smooth)

    def _init_weight_corrections(self):
        funcs = []
        for i in xrange(self._system.natom):
            number = self._system.numbers[i]
            funcs.append((
                self._system.coordinates[i],
                [self.proatomdb.get_spline(number)],
            ))
        wcor = self.grid.compute_weight_corrections(funcs)
        self._cache.dump('wcor', wcor)

    def _init_partitioning(self):
        self._update_promolecule()
        self._store_at_weights()
        self._work_cleanup()

    def _update_promolecule(self):
        promoldens = self.cache.load('promoldens', alloc=self.grid.shape)[0]
        promoldens[:] = 0.0
        work = self.cache.load('work0', alloc=self.grid.shape)[0]
        for index in xrange(self._system.natom):
            self.compute_proatom(index, work)
            promoldens += work
            # Merely for efficiency:
            self._store.dump(work, 'at_weights', index)

    def _store_at_weights(self):
        # Compute the atomic weight functions if this is useful. This is merely
        # a matter of efficiency.
        promoldens = self.cache.load('promoldens')
        work = self.cache.load('work0', alloc=self.grid.shape)[0]
        for index in xrange(self._system.natom):
            if ('at_weights', index) in self._store:
                self._store.load(work, 'at_weights', index)
                work /= promoldens
                self._store.dump(work, 'at_weights', index)

    def _work_cleanup(self):
        self.cache.discard('work0')

    def get_cutoff_radius(self, index):
        '''The radius at which the weight function goes to zero'''
        rtf = self.proatomdb.get_rtransform(self._system.numbers[index])
        return rtf.radius(rtf.npoint-1)

    def do_all(self):
        names = StockholderCPart.do_all(self)
        self.do_dispersion()
        return names + ['volumes', 'volume_ratios', 'c6s']
