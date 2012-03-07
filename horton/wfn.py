# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>, ...
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
"""Wavefunction implementations

   The essential part of the wavefunction consists of the expansion coefficients
   in a certain basis. If all the relevant matrix elements with respect to the
   basis are already known, it is not critical to have a detailed specification
   of the basis set.
"""


import numpy as np

from horton.matrix import DenseExpansion


class ClosedShellWFN(object):
    def __init__(self, nep, basis=None, nbasis=None):
        """
           **Arguments:**

           nep
                The number of electron pairs in the wave function.

           **Optional arguments:**

           basis
                A specification of the basis set.

           nbasis
                The number of basis functions.

           Either basis or nbasis must be given.
        """
        self._nep = nep
        self._basis = basis
        self._nbasis = nbasis

        if self.nep <= 0:
            raise ValueError('At least one pair of electrons is required.')
        if self.nbasis < self.nep:
            raise ValueError('The number of spatial basis functions must not be lower than the number of electron pairs.')

        norb = self.nbasis
        self._expansion = DenseExpansion(norb, self.nbasis)
        self._epsilons = np.zeros(norb, float)

    def get_nep(self):
        return self._nep

    nep = property(get_nep)

    def get_basis(self):
        return self._basis

    basis = property(get_basis)

    def get_nbasis(self):
        if self.basis is not None:
            return self._basis.size
        else:
            return self._nbasis

    nbasis = property(get_nbasis)

    def get_density_matrix(self):
        return self._expansion.get_density_matrix(self.nep)

    def apply_two_body(self, two_body, direct, exchange):
        """Compute the direct and exchange one-body operators

           **Arguments:**

           two_body
                The two_body operators, e.g. an instance of the Dense4 class

           direct, exchange
                The two output arguments, which must be Dense2 instances.

           Returns the one-body operator for the direct and the exchange
           contribution to the derivative of the expectation value of the
           two-body operator.
        """
        dm = self.get_density_matrix()
        two_body.apply_direct(dm, direct)
        two_body.apply_exchange(dm, exchange)
        return direct, exchange
