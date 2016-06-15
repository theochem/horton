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
'''Minimal Basis Iterative Stockholder (MBIS) partitioning'''


import numpy as np
from horton.log import log
from horton.part.stockholder import StockholderWPart
from horton.part.iterstock import IterativeProatomMixin


__all__ = ['MBISWPart']


def _get_nshell(number):
    noble = np.array([2, 10, 18, 36, 54, 86, 118])
    return noble.searchsorted(number) + 1


def _get_initial_mbis_propars(number):
    nshell = _get_nshell(number)
    propars = np.zeros(2*nshell, float)
    S0 = 2.0*number
    if nshell > 1:
        S1 = 2.0
        alpha = (S1/S0)**(1.0/(nshell-1))
    else:
        alpha = 1.0
    nel_in_shell = np.array([2.0, 8.0, 8.0, 18.0, 18.0, 32.0, 32.0])
    for ishell in xrange(nshell):
        propars[2*ishell] = nel_in_shell[ishell]
        propars[2*ishell+1] = S0*alpha**ishell
    propars[-2] = number - propars[:-2:2].sum()
    return propars


def _opt_mbis_propars(rho, propars, rgrid, threshold):
    assert len(propars)%2 == 0
    nshell = len(propars)/2
    r = rgrid.radii
    terms = np.zeros((nshell, len(r)), float)
    oldpro = None
    for irep in xrange(1000):
        # compute the contributions to the pro-atom
        for ishell in xrange(nshell):
            N = propars[2*ishell]
            S = propars[2*ishell+1]
            terms[ishell] = (N*S**3/(8*np.pi))*np.exp(-S*r)
        pro = terms.sum(axis=0)
        # transform to partitions
        terms *= rho/pro
        # the partitions and the updated parameters
        for ishell in xrange(nshell):
            m0 = rgrid.integrate(terms[ishell])
            m1 = rgrid.integrate(terms[ishell], r)
            propars[2*ishell] = m0
            propars[2*ishell+1] = 3*m0/m1
        # check for convergence
        if oldpro is None:
            change = 1e100
        else:
            error = oldpro - pro
            change = np.sqrt(rgrid.integrate(error, error))
        if change < threshold:
            return propars
        oldpro = pro
    assert False


class MBISWPart(IterativeProatomMixin, StockholderWPart):
    '''Iterative Stockholder Partitioning with Becke-Lebedev grids'''
    name = 'mbis'
    options = ['lmax', 'threshold', 'maxiter']
    linear = False

    def __init__(self, coordinates, numbers, pseudo_numbers, grid, moldens,
                 spindens=None, lmax=3, threshold=1e-6, maxiter=500):
        '''
           **Optional arguments:** (that are not defined in ``WPart``)

           threshold
                The procedure is considered to be converged when the maximum
                change of the charges between two iterations drops below this
                threshold.

           maxiter
                The maximum number of iterations. If no convergence is reached
                in the end, no warning is given.
                Reduce the CPU cost at the expense of more memory consumption.
        '''
        self._threshold = threshold
        self._maxiter = maxiter
        StockholderWPart.__init__(self, coordinates, numbers, pseudo_numbers,
                                  grid, moldens, spindens, True, lmax)

    def _init_log_scheme(self):
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Minimal Basis Iterative Stockholder (MBIS)'),
                ('Convergence threshold', '%.1e' % self._threshold),
                ('Maximum iterations', self._maxiter),
            ])

    def get_rgrid(self, iatom):
        return self.get_grid(iatom).rgrid

    def get_proatom_rho(self, iatom, propars=None):
        if propars is None:
            propars = self.cache.load('propars')
        rgrid = self.get_rgrid(iatom)
        r = rgrid.radii
        y = np.zeros(len(r), float)
        d = np.zeros(len(r), float)
        my_propars = propars[self._ranges[iatom]:self._ranges[iatom+1]]
        for ishell in xrange(self._nshells[iatom]):
            N, S = my_propars[2*ishell:2*ishell+2]
            f = (N*S**3/(8*np.pi))*np.exp(-S*r)
            y += f
            d += -S*f
        return y, d

    def _init_propars(self):
        IterativeProatomMixin._init_propars(self)
        self._ranges = [0]
        self._nshells = []
        for iatom in xrange(self.natom):
            nshell = _get_nshell(self.numbers[iatom])
            self._ranges.append(self._ranges[-1]+2*nshell)
            self._nshells.append(nshell)
        ntotal = self._ranges[-1]
        propars = self.cache.load('propars', alloc=ntotal, tags='o')[0]
        for iatom in xrange(self.natom):
            propars[self._ranges[iatom]:self._ranges[iatom+1]] = _get_initial_mbis_propars(self.numbers[iatom])
        return propars

    def _update_propars_atom(self, iatom):
        # compute spherical average
        atgrid = self.get_grid(iatom)
        rgrid = atgrid.rgrid
        dens = self.get_moldens(iatom)
        at_weights = self.cache.load('at_weights', iatom)
        spherical_average = np.clip(atgrid.get_spherical_average(at_weights, dens), 1e-100, np.inf)

        # assign as new propars
        my_propars = self.cache.load('propars')[self._ranges[iatom]:self._ranges[iatom+1]]
        my_propars[:] = _opt_mbis_propars(spherical_average, my_propars.copy(), rgrid, self._threshold)

        # compute the new charge
        pseudo_population = rgrid.integrate(spherical_average)
        charges = self.cache.load('charges', alloc=self.natom, tags='o')[0]
        charges[iatom] = self.pseudo_numbers[iatom] - pseudo_population

    def _finalize_propars(self):
        IterativeProatomMixin._finalize_propars(self)
        propars = self.cache.load('propars')
        core_charges = []
        valence_charges = []
        valence_widths = []
        for iatom in xrange(self.natom):
            my_propars = propars[self._ranges[iatom]:self._ranges[iatom+1]]
            valence_charges.append(-my_propars[-2])
            valence_widths.append(1.0/my_propars[-1])
        valence_charges = np.array(valence_charges)
        valence_widths = np.array(valence_widths)
        core_charges = self._cache.load('charges') - valence_charges
        self.cache.dump('core_charges', core_charges, tags='o')
        self.cache.dump('valence_charges', valence_charges, tags='o')
        self.cache.dump('valence_widths', valence_widths, tags='o')
