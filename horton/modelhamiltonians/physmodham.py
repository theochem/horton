# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
'''Physical Model Hamiltonian'''


__all__ = [
    'PhysModHam', 'Hubbard'
]


class PhysModHam(object):
    '''Base class for the Physical Model Hamiltonians: 1-D Hubbard,
    PPP, Ising, etc '''
    def __init__(self, pbc=True):

        '''
           **Attributes:**

            pdb
                Periodic boundary conditions. Default, pdb=true

        '''
        self._pbc = pbc

    def _get_pbc(self):
        '''The periodic boundary conditions'''
        return self._pbc

    pbc = property(_get_pbc)

    def compute_kinetic(self, lf, tparam):
        '''Calculate the one-body term of the 1D Hubbard Hamiltonian'''
        raise NotImplementedError

    def compute_er(self, lf, uparam):
        '''Calculate the the-body term of the 1D Hubbard Hamiltonian'''
        raise NotImplementedError

    def compute_overlap(self, lf):
        '''Calculate overlap of the 1D Hubbard Hamiltonian, (identity matrix)'''
        raise NotImplementedError


class Hubbard(PhysModHam):
    '''The 1-D Hubbard class Hamiltonian'''
    def compute_kinetic(self, lf, tparam):
        '''Calculate the one-body term of the 1D Hubbard Hamiltonian'''
        result = lf.create_two_index(lf.default_nbasis)
        for i in range (lf.default_nbasis-1):
            result.set_element(i, i+1, tparam)
        if  self.pbc == True:
            result.set_element(lf.default_nbasis-1, 0, tparam)
        return result

    def compute_er(self, lf, uparam):
        '''Calculate the the-body term of the 1D Hubbard Hamiltonian'''
        result = lf.create_four_index(lf.default_nbasis)
        for i in range (lf.default_nbasis):
            result.set_element(i, i, i, i, uparam)
        return result

    def compute_overlap(self, lf):
        '''Calculate overlap of the 1D Hubbard Hamiltonian, (identity matrix)'''
        result = lf.create_two_index(lf.default_nbasis)
        result.assign_diagonal(1.0)
        return result
