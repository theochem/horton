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


from horton.meanfield.linear import LinearObservable


__all__ = [
    'KineticEnergy', 'ExternalPotential',
]


class KineticEnergy(LinearObservable):
    kinetic = True

    def __init__(self, obasis, lf, wfn, label='kin'):
        LinearObservable.__init__(self, obasis, lf, wfn, label)

    def get_operator(self):
        kinetic, new = self.cache.load('kin', alloc=self._lf.create_one_body,
                                        tags='o')
        if new:
            self._obasis.compute_kinetic(kinetic)
        return kinetic


class ExternalPotential(LinearObservable):
    external = True

    def __init__(self, obasis, lf, wfn, numbers, coordinates, label='ne'):
        self.numbers = numbers
        self.coordinates = coordinates
        LinearObservable.__init__(self, obasis, lf, wfn, label)

    def get_operator(self):
        nuclear_attraction, new = self.cache.load('na', alloc=self._lf.create_one_body,
                                                  tags='o')
        if new:
            # TODO: ghost atoms and extra charges
            self._obasis.compute_nuclear_attraction(self.numbers.astype(float),
                                                   self.coordinates,
                                                   nuclear_attraction)
        return nuclear_attraction
