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


import numpy as np

from horton.log import log, timer
from horton.cache import Cache
from horton.meanfield.core import KineticEnergy, ExternalPotential
from horton.meanfield.builtin import Hartree


__all__ = [
    'Hamiltonian',
]


class Hamiltonian(object):
    def __init__(self, system, terms, grid=None, auto_complete=True):
        '''
           **Arguments:**

           system
                The System object for which the energy must be computed.

           terms
                The terms in the Hamiltonian. Kinetic energy, external
                potential (nuclei), and Hartree are added automatically.

           **Optional arguments:**

           grid
                The integration grid, in case some terms need one.

           auto_complete
                When set to False, the kinetic energy, external potential and
                Hartree terms are not added automatically.
        '''
        # check arguments:
        if len(terms) == 0:
            raise ValueError('At least one term must be present in the Hamiltonian.')
        for term in terms:
            if term.require_grid and grid is None:
                raise TypeError('The term %s requires a grid, but not grid is given.' % term)

        # Assign attributes
        self.system = system
        self.terms = list(terms)
        self.grid = grid

        if auto_complete:
            # Add standard terms if missing
            #  1) Kinetic energy
            if sum(isinstance(term, KineticEnergy) for term in terms) == 0:
                self.terms.append(KineticEnergy())
            #  2) Hartree (or HatreeFock, which is a subclass of Hartree)
            if sum(isinstance(term, Hartree) for term in terms) == 0:
                self.terms.append(Hartree())
            #  3) External Potential
            if sum(isinstance(term, ExternalPotential) for term in terms) == 0:
                self.terms.append(ExternalPotential())

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

    def clear(self, dealloc=False):
        '''Mark the properties derived from the wfn as outdated.

           **Optional arguments:**
           
           dealloc
                When True, arrays are really deallocated instead of keeping the
                memory for later use.

           This method does not recompute anything, but just marks operators
           as outdated. They are recomputed as they are needed.
        '''
        self.cache.clear(dealloc=dealloc)

    def compute(self):
        '''Compute the energy.

           **Returns:**

           The total energy, including nuclear-nuclear repulsion.
        '''
        if log.do_high:
            log('Computing the energy of the system.')
            log.hline()
            log('                   Energy term                 Value')
            log.hline()

        total = 0.0
        for term in self.terms:
            energy = term.compute()
            if log.do_high:
                log('%30s  %20.10f' % (term.label, energy))
            self.system.extra['energy_%s' % term.label] = energy
            total += energy

        energy = self.system.compute_nucnuc()
        total += energy
        # Store result in chk file
        self.system.extra['energy'] = total
        self.system.update_chk('extra')

        if log.do_high:
            log('%30s  %20.10f' % ('nn', energy))
            log('%30s  %20.10f' % ('total', total))
            log.hline()

        return total

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
        for term in self.terms:
            term.add_fock_matrix(fock_alpha, fock_beta)
