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
'''Base class for energy terms and other observables of the wavefunction'''


__all__ = [
    'Observable',
]


class Observable(object):
    # The following are needed for the idiot-proof option of the Hamiltonian
    # class:
    kinetic = False    # Set to True for kinetic energys terms.
    hartree = False    # Set to True for hartree terms.
    exchange = False   # Set to True for exhange terms.
    external = False   # Set to True for an external potential.

    def __init__(self, lf, label):
        self.label = label
        self._hamiltonian = None
        self._lf = lf

    def set_hamiltonian(self, hamiltonian):
        if not self._hamiltonian is None:
            raise ValueError('This term is already assigned to a Hamiltonian.')
        self._hamiltonian = hamiltonian

    # The following four properties are added for convenience:

    def _get_lf(self):
        return self._lf

    lf = property(_get_lf)

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
