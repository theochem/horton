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
'''Hirshfeld Partitioning'''


from horton.cache import just_once
from horton.log import log
from horton.part.stockholder import StockholderWPart, StockholderCPart


__all__ = ['HirshfeldWPart', 'HirshfeldCPart']


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
    name = 'h'
    options = []
    linear = True

    def __init__(self, proatomdb):
        self._proatomdb = proatomdb

    def _init_log_scheme(self):
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld'),
                ('Proatomic DB',  self.proatomdb),
            ])
            log.cite('hirshfeld1977', 'the use of Hirshfeld partitioning')

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def get_rgrid(self, index):
        number = self.system.numbers[index]
        return self.proatomdb.get_rgrid(number)

    def get_proatom_rho(self, index):
        return self.proatomdb.get_rho(self._system.numbers[index])

    @just_once
    def do_dispersion(self):

        if log.do_medium:
            log.cite('tkatchenko2009', 'the method to evaluate atoms-in-molecules C6 parameters')
            log.cite('chu2004', 'the reference C6 parameters of isolated atoms')
            log.cite('yan1996', 'the isolated hydrogen C6 parameter')

        ref_c6s = { # reference C6 values in atomic units
            1: 6.499, 2: 1.42, 3: 1392.0, 4: 227.0, 5: 99.5, 6: 46.6, 7: 24.2,
            8: 15.6, 9: 9.52, 10: 6.20, 11: 1518.0, 12: 626.0, 13: 528.0, 14:
            305.0, 15: 185.0, 16: 134.0, 17: 94.6, 18: 64.2, 19: 3923.0, 20:
            2163.0, 21: 1383.0, 22: 1044.0, 23: 832.0, 24: 602.0, 25: 552.0, 26:
            482.0, 27: 408.0, 28: 373.0, 29: 253.0, 30: 284.0, 31: 498.0, 32:
            354.0, 33: 246.0, 34: 210.0, 35: 162.0, 36: 130.0, 37: 4769.0, 38:
            3175.0, 49: 779.0, 50: 659.0, 51: 492.0, 52: 445.0, 53: 385.0,
        }

        volumes, new_volumes = self._cache.load('volumes', alloc=self.system.natom, tags='o')
        volume_ratios, new_volume_ratios = self._cache.load('volume_ratios', alloc=self.system.natom, tags='o')
        c6s, new_c6s = self._cache.load('c6s', alloc=self.system.natom, tags='o')

        if new_volumes or new_volume_ratios or new_c6s:
            self.do_populations()
            self.do_moments()
            radial_moments = self._cache.load('radial_moments')
            populations = self._cache.load('populations')

            if log.do_medium:
                log('Computing atomic dispersion coefficients.')

            for i in xrange(self.system.natom):
                n = self.system.numbers[i]
                volumes[i] = radial_moments[i,2]/populations[i]
                ref_volume = self.proatomdb.get_record(n, 0).get_moment(3)/n
                volume_ratios[i] = volumes[i]/ref_volume
                if n in ref_c6s:
                    c6s[i] = (volume_ratios[i])**2*ref_c6s[n]
                else:
                    c6s[i] = -1 # This is just to indicate that no value is available.


class HirshfeldWPart(HirshfeldMixin, StockholderWPart):
    options = HirshfeldMixin.options + ['epsilon']

    def __init__(self, system, grid, proatomdb, local=True, epsilon=0):
        check_proatomdb(system, proatomdb)
        HirshfeldMixin. __init__(self, proatomdb)
        StockholderWPart.__init__(self, system, grid, local, epsilon)


class HirshfeldCPart(HirshfeldMixin, StockholderCPart):
    def __init__(self, system, grid, local, moldens, proatomdb, wcor_numbers, wcor_rcut_max=2.0, wcor_rcond=0.1):
        '''
           See CPart base class for the description of the arguments.
        '''
        check_proatomdb(system, proatomdb)
        HirshfeldMixin. __init__(self, proatomdb)
        StockholderCPart.__init__(self, system, grid, local, moldens, wcor_numbers, wcor_rcut_max, wcor_rcond)

    def get_cutoff_radius(self, index):
        '''The radius at which the weight function goes to zero'''
        rtf = self.get_rgrid(index).rtransform
        return rtf.radius(rtf.npoint-1)

    def get_wcor_funcs(self, index):
        number = self.system.numbers[index]
        if number in self.wcor_numbers:
            return [(self._system.coordinates[index], [self._proatomdb.get_spline(number)])]
        else:
            return []
