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
'''Mean-field DFT/HF Hamiltonian data structures'''


from horton.log import log
from horton.cache import Cache
from horton.meanfield.wfn import UnrestrictedWFN


__all__ = [
    'Hamiltonian',
]


class Hamiltonian(object):
    def __init__(self, terms, external=None):
        '''
           **Arguments:**

           terms
                The terms in the Hamiltonian.

           **Optional arguments:**

           external
                A dictionary with external energy contributions that do not
                depend on the wavefunction, e.g. nuclear-nuclear interactions
                or QM/MM mechanical embedding terms. Use ``nn`` as key for the
                nuclear-nuclear term.
        '''
        # check arguments:
        if len(terms) == 0:
            raise ValueError('At least one term must be present in the Hamiltonian.')

        # Assign attributes
        self.terms = list(terms)
        self.external = {} if external is None else external

        # Create a cache for shared intermediate results. This cache should only
        # be used for derived quantities that depend on the wavefunction and
        # need to be updated at each SCF cycle.
        self.cache = Cache()

        # bind the terms to this hamiltonian such that certain shared
        # intermediated results can be reused for the sake of efficiency.
        for term in self.terms:
            term.set_hamiltonian(self)

    def add_term(self, term):
        '''Add a new term to the hamiltonian'''
        self.terms.append(term)
        term.set_hamiltonian(self)

    def clear(self):
        '''Mark the properties derived from the wfn as outdated.

           This method does not recompute anything, but just marks operators
           as outdated. They are recomputed as they are needed.
        '''
        self.cache.clear()

    def compute(self):
        '''Compute the energy.

           **Returns:**

           The total energy, including nuclear-nuclear repulsion.
        '''
        total = 0.0
        for term in self.terms:
            energy = term.compute()
            self.cache['energy_%s' % term.label] = energy
            total += energy
        for key, energy in self.external.iteritems():
            self.cache['energy_%s' % key] = energy
            total += energy
        self.cache['energy'] = total
        return total

    def log_energy(self):
        '''Write an overview of the last energy computation on screen'''
        log('Contributions to the energy:')
        log.hline()
        log('                                       Energy term                 Value')
        log.hline()
        for term in self.terms:
            energy = self.cache['energy_%s' % term.label]
            log('%50s  %20.12f' % (term.label, energy))
        for key, energy in self.external.iteritems():
            log('%50s  %20.12f' % (key, energy))
        log('%50s  %20.12f' % ('total', self.cache['energy']))
        log.hline()
        log.blank()

    def compute_fock(self, fock_alpha, fock_beta):
        '''Compute alpha (and beta) Fock matrix(es).

           **Arguments:**

           fock_alpha
                A One-Body operator output argument for the alpha fock matrix.

           fock_alpha
                A One-Body operator output argument for the beta fock matrix.

           In the case of a closed-shell computation, the argument fock_beta is
           ``None``.
        '''
        # Loop over all terms and add contributions to the Fock matrix.
        for term in self.terms:
            term.add_fock_matrix(fock_alpha, fock_beta)
