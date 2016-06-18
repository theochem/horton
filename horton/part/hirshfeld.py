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
'''Hirshfeld partitioning'''


from horton.cache import just_once
from horton.log import log
from horton.part.stockholder import StockholderWPart


__all__ = ['HirshfeldWPart']


def check_proatomdb(numbers, pseudo_numbers, proatomdb):
    # Check if the same pseudo numbers (effective core charges) are used for the
    # molecule and the proatoms.
    for i in xrange(len(numbers)):
        number = numbers[i]
        pseudo_number = pseudo_numbers[i]
        pseudo_number_expected = proatomdb.get_record(number, 0).pseudo_number
        if pseudo_number_expected != pseudo_number:
            raise ValueError('The pseudo number of atom %i does not match with the proatom database (%i!=%i)' % (
                i, pseudo_number, pseudo_number_expected))


class HirshfeldMixin(object):
    name = 'h'
    options = ['lmax']
    linear = True

    def __init__(self, numbers, pseudo_numbers, proatomdb):
        check_proatomdb(numbers, pseudo_numbers, proatomdb)
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
        number = self.numbers[index]
        return self.proatomdb.get_rgrid(number)

    def get_proatom_rho(self, index):
        return self.proatomdb.get_rho(self.numbers[index], do_deriv=True)

    @just_once
    def do_dispersion(self):
        if self.lmax < 3:
            if log.do_warning:
                log.warn('Skipping the computation of dispersion coefficients because lmax=%i<3' % self.lmax)
            return

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

        volumes, new_volumes = self._cache.load('volumes', alloc=self.natom, tags='o')
        volume_ratios, new_volume_ratios = self._cache.load('volume_ratios', alloc=self.natom, tags='o')
        c6s, new_c6s = self._cache.load('c6s', alloc=self.natom, tags='o')

        if new_volumes or new_volume_ratios or new_c6s:
            self.do_moments()
            radial_moments = self._cache.load('radial_moments')

            if log.do_medium:
                log('Computing atomic dispersion coefficients.')

            for i in xrange(self.natom):
                n = self.numbers[i]
                volumes[i] = radial_moments[i, 3]
                ref_volume = self.proatomdb.get_record(n, 0).get_moment(3)
                volume_ratios[i] = volumes[i]/ref_volume
                if n in ref_c6s:
                    c6s[i] = (volume_ratios[i])**2*ref_c6s[n]
                else:
                    c6s[i] = -1 # This is just to indicate that no value is available.


class HirshfeldWPart(HirshfeldMixin, StockholderWPart):
    '''Hirshfeld partitioning with Becke-Lebedev grids'''

    def __init__(self, coordinates, numbers, pseudo_numbers, grid, moldens,
                 proatomdb, spindens=None, local=True, lmax=3):
        '''
           **Arguments:** (that are not defined in ``WPart``)

           proatomdb
                In instance of ProAtomDB that contains all the reference atomic
                densities.
        '''
        HirshfeldMixin. __init__(self, numbers, pseudo_numbers, proatomdb)
        StockholderWPart.__init__(self, coordinates, numbers, pseudo_numbers,
                                  grid, moldens, spindens, local, lmax)
