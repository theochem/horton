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
from horton.meanfield.observable import Observable


__all__ = [
    'LinearObservable', 'KineticEnergy', 'ExternalPotential',
]


class LinearObservable(Observable):
    '''Base class for all terms that are linear in the density matrix

       This is (technically) a special class because the Fock operator does not
       have to be recomputed when the density matrix changes.
    '''
    def get_operator(self, system):
        # subclasses should return the operator and a suffix.
        raise NotImplementedError

    def prepare_system(self, system, cache, grid):
        Observable.prepare_system(self, system, cache, grid)
        self.operator = self.get_operator(system)

    def compute(self):
        if self.system.wfn.closed_shell:
            return 2*self.operator.expectation_value(self.system.wfn.dm_alpha)
        else:
            return self.operator.expectation_value(self.system.wfn.dm_full)

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        for fock in fock_alpha, fock_beta:
            if fock is not None:
                fock.iadd(self.operator, scale)


class KineticEnergy(LinearObservable):
    def __init__(self, label='kin'):
        LinearObservable.__init__(self, label)

    def get_operator(self, system):
        return system.get_kinetic()


class ExternalPotential(LinearObservable):
    def __init__(self, label='ne'):
        LinearObservable.__init__(self, label)

    def get_operator(self, system):
        tmp = system.get_nuclear_attraction().copy() # take copy because of next line
        tmp.iscale(-1)
        return tmp
