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

from horton.gbasis.cext import get_shell_nbasis


__all__ = ['assign_mulliken_operator', 'get_mulliken_operators']


def assign_mulliken_operator(operator, system, index):
    '''Fill in the mulliken operator in the first argument

       **Arguments:**

       operator
            A One body operator for the output

       system
            The system for which the Mulliken operator is to be constructed

       index
            The index of the atom (center) for which the Mulliken operator
            needs to be constructed

       This routine implies that the first ``natom`` centers in the object
       system.obasis corresponds to the atoms in the system object.
    '''
    mask = np.zeros(system.obasis.nbasis, dtype=bool)
    begin = 0
    for ishell in xrange(system.obasis.nshell):
        end = begin + get_shell_nbasis(system.obasis.shell_types[ishell])
        if system.obasis.shell_map[ishell] != index:
            mask[begin:end] = True
        begin = end
    olp = system.get_overlap()
    operator.assign(olp)
    operator._array[mask] = 0.0
    operator._array[:] = 0.5*(operator._array + operator._array.T)


def get_mulliken_operators(system):
    '''Return a list of Mulliken operators for the given system.'''
    operators = []
    for icenter in xrange(system.obasis.ncenter):
        operator = system.lf.create_one_body()
        assign_mulliken_operator(operator, system, icenter)
        operators.append(operator)
    return operators
