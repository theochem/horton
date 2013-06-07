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


from horton.meanfield.observable import Observable
from horton.meanfield.core import LinearObservable


__all__ = ['CustomLinearObservable', 'CustomGridLinearObservable']


class CustomLinearObservable(LinearObservable):
    '''This is a user-defined term that is linear in the density matrix

       This term can be used to implemented perturbations by finite fields.
    '''
    def __init__(self, label, operator):
        self.operator = operator
        LinearObservable.__init__(self, label)

    def prepare_system(self, system, cache, grid):
        # override because assignment of self.operator is not needed.
        Observable.prepare_system(self, system, cache, grid)


class CustomGridLinearObservable(LinearObservable):
    '''This is a user-defined term, through a grid, that is linear in the density matrix

       This term can be used to implemented perturbations by finite potentials
       defined on a real-space grid.
    '''
    def __init__(self, label, custom_grid, potential):
        LinearObservable.__init__(self, label)
        self.custom_grid = custom_grid
        self.potential = potential

    def get_operator(self, system):
        operator = system.lf.create_one_body()
        self.system.compute_grid_density_fock(self.custom_grid.points, self.custom_grid.weights, self.potential, operator)
        return operator
