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
'''Linear energy terms (Core Hamiltonian)'''


import numpy as np

from horton.log import log, timer
from horton.cache import Cache
from horton.meanfield.observable import Observable
from horton.meanfield.wfn import RestrictedWFN


__all__ = [
    'LinearObservable', 'KineticEnergy', 'ExternalPotential',
    'CustomLinearObservable',
]


class LinearObservable(Observable):
    '''Base class for all observables that are linear in the density matrix

       This is (technically) a special class because the Fock operator does not
       have to be recomputed when the density matrix changes.
    '''
    def get_operator(self):
        # subclasses should return the operator
        raise NotImplementedError

    def compute(self):
        operator = self.get_operator()
        if isinstance(self.system.wfn, RestrictedWFN):
            return 2*operator.expectation_value(self.system.wfn.dm_alpha)
        else:
            return operator.expectation_value(self.system.wfn.dm_full)

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        operator = self.get_operator()
        for fock in fock_alpha, fock_beta:
            if fock is not None:
                fock.iadd(operator, scale)


class KineticEnergy(LinearObservable):
    def __init__(self, label='kin'):
        LinearObservable.__init__(self, label)

    def get_operator(self):
        return self.system.get_kinetic()


class ExternalPotential(LinearObservable):
    def __init__(self, label='ne'):
        LinearObservable.__init__(self, label)

    def get_operator(self):
        return self.system.get_nuclear_attraction()


class CustomLinearObservable(LinearObservable):
    '''This is a user-defined term that is linear in the density matrix

       This term can be used to implemented perturbations by finite fields.
    '''
    def __init__(self, label, get_operator):
        self.get_operator = get_operator
        LinearObservable.__init__(self, label)
