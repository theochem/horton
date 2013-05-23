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


__all__ = ['get_mulliken_operators']


def get_mulliken_operators(sys):
    '''Return a list of mulliken operators for the given system.'''
    operators = []
    for icenter in xrange(sys.obasis.ncenter):
        mask = np.zeros(sys.obasis.nbasis, dtype=bool)
        begin = 0
        for ishell in xrange(sys.obasis.nshell):
            end = begin + get_shell_nbasis(sys.obasis.shell_types[ishell])
            if sys.obasis.shell_map[ishell] != icenter:
                mask[begin:end] = True
            begin = end
        pop = sys.get_overlap().copy()
        pop._array[mask] = 0.0
        pop._array[:] = 0.5*(pop._array + pop._array.T)
        operators.append(pop)
    return operators
