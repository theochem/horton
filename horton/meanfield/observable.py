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
'''Base classes for energy terms and other observables of the wavefunction'''


from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = [
    'Observable', 'OneBodyTerm', 'DirectTerm', 'ExchangeTerm',
]


class Observable(object):
    def __init__(self, label):
        self.label = label
        self._hamiltonian = None

    def set_hamiltonian(self, hamiltonian):
        if not self._hamiltonian is None:
            raise ValueError('This term is already assigned to a Hamiltonian.')
        self._hamiltonian = hamiltonian

    # The following four properties are added for convenience:

    def _get_cache(self):
        '''The cache of the hamiltonian object, cleared after a change in wfn.'''
        return self._hamiltonian.cache #FIXME: remove the hamiltonian dependency

    cache = property(_get_cache)

    def compute(self):
        raise NotImplementedError

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        '''Add contributions to alpha (and beta) Fock matrix(es).

           **Arguments:**

           fock_alpha
                A One-Body operator output argument for the alpha fock matrix.

           fock_alpha
                A One-Body operator output argument for the beta fock matrix.

           **Optional arguments:**

           scale
                A scale factor for this contribution

           In the case of a closed-shell computation, the argument fock_beta is
           ``None``.
        '''
        raise NotImplementedError


class OneBodyTerm(Observable):
    '''Class for all observables that are linear in the density matrix

       A restriction of this implementation is that the fock operator for the
       alpha and the beta electrons are the same.
    '''
    def __init__(self, operator, wfn, label):
        self._operator = operator
        self._wfn = wfn
        Observable.__init__(self, label)

    def get_operator(self):
        # subclasses should override this method such that it
        # returns the operator.
        raise NotImplementedError

    def compute(self):
        '''See Observable.compute'''
        if isinstance(self._wfn, RestrictedWFN):
            return 2 * self._operator.expectation_value(self._wfn.dm_alpha)
        else:
            return self._operator.expectation_value(self._wfn.dm_full)

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        '''See Observable.add_fock_matrix'''
        for fock in fock_alpha, fock_beta:
            if fock is not None:
                fock.iadd(self._operator, scale)


class DirectTerm(Observable):
    def __init__(self, operator, wfn, label):
        self._operator = operator
        self._wfn = wfn
        Observable.__init__(self, label)

    def _update_direct(self):
        '''Recompute the direct operator if it has become invalid'''
        if isinstance(self._wfn, RestrictedWFN):
            direct, new = self.cache.load('op_%s' % self.label,
                                          alloc=self._wfn.dm_alpha.new)
            if new:
                self._operator.apply_direct(self._wfn.dm_alpha, direct)
                direct.iscale(2)
        else:
            direct, new = self.cache.load('op_%s' % self.label,
                                          alloc=self._wfn.dm_full.new)
            if new:
                self._operator.apply_direct(self._wfn.dm_full, direct)

    def compute(self):
        self._update_direct()
        direct = self.cache.load('op_%s' % self.label)
        if isinstance(self._wfn, RestrictedWFN):
            return direct.expectation_value(self._wfn.dm_alpha)
        else:
            return 0.5 * direct.expectation_value(self._wfn.dm_full)

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        self._update_direct()
        direct = self.cache.load('op_%s' % self.label)
        fock_alpha.iadd(direct, scale)
        if isinstance(self._wfn, UnrestrictedWFN):
            fock_beta.iadd(direct, scale)


class ExchangeTerm(Observable):
    def __init__(self, operator, wfn, label, fraction_exchange=1.0):
        self._operator = operator
        self._wfn = wfn
        self.fraction_exchange = fraction_exchange
        Observable.__init__(self, label)

    def _update_exchange(self):
        '''Recompute the Exchange operator(s) if invalid'''
        def helper(select):
            dm = self._wfn.get_dm(select)
            exchange, new = self.cache.load('op_%s_%s' % (self.label, select),
                                            alloc=dm.new)
            if new:
                self._operator.apply_exchange(dm, exchange)

        helper('alpha')
        if isinstance(self._wfn, UnrestrictedWFN):
            helper('beta')

    def compute(self):
        self._update_exchange()
        exchange_alpha = self.cache.load('op_%s_alpha' % self.label)
        if isinstance(self._wfn, RestrictedWFN):
            return -self.fraction_exchange * exchange_alpha.expectation_value(self._wfn.dm_alpha)
        else:
            exchange_beta = self.cache.load('op_%s_beta' % self.label)
            return -0.5 * self.fraction_exchange * exchange_alpha.expectation_value(self._wfn.dm_alpha) \
                   -0.5 * self.fraction_exchange * exchange_beta.expectation_value(self._wfn.dm_beta)

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        self._update_exchange()
        fock_alpha.iadd(self.cache.load('op_%s_alpha' % self.label),
                        -self.fraction_exchange * scale)
        if fock_beta is not None:
            fock_beta.iadd(self.cache.load('op_%s_beta' % self.label),
                           -self.fraction_exchange * scale)
