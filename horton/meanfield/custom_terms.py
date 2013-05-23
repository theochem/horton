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


from horton.meanfield.term import HamiltonianTerm
from horton.meanfield.core_terms import FixedTerm


__all__ = ['CustomFixedTerm', 'CustomGridFixedTerm']


class CustomFixedTerm(FixedTerm):
    '''This is a user-defined term that is linear in the density matrix

       This term can be used to implemented perturbations by finite fields.
    '''
    def __init__(self, operator, suffix):
        self.operator = operator
        self.suffix = suffix

    def prepare_system(self, system, cache, grid):
        # override because assignment of self.operator and self.suffix are not needed.
        HamiltonianTerm.prepare_system(self, system, cache, grid)


class CustomGridFixedTerm(FixedTerm):
    '''This is a user-defined term, through a grid, that is linear in the density matrix

       This term can be used to implemented perturbations by finite potentials
       defined on a real-space grid.
    '''
    def __init__(self, custom_grid, potential, suffix):
        self.custom_grid = custom_grid
        self.potential = potential
        self.suffix = suffix

    def get_operator(self, system):
        operator = system.lf.create_one_body()
        self.system.compute_grid_density_fock(self.custom_grid.points, self.custom_grid.weights, self.potential, operator)
        return operator, self.suffix
