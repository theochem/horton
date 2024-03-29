# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2022 The HORTON Development Team
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
"""Physical Model Hamiltonian."""


import numpy as np


__all__ = [
    'PhysModHam', 'Hubbard'
]


class PhysModHam(object):
    """Base class for the Physical Model Hamiltonians: 1-D Hubbard, PPP, Ising, etc."""

    def __init__(self, nbasis, pbc=True):
        """Initialize a PhysModHam object.

        Parameters
        ----------
        nbasis : int
            The number of sites.
        pdb : bool
            Periodic boundary conditions. Default, pdb=true
        """
        self._nbasis = nbasis
        self._pbc = pbc

    @property
    def pbc(self):
        """The periodic boundary condition."""
        return self._pbc

    @property
    def nbasis(self):
        """The number of sites."""
        return self._nbasis

    def compute_kinetic(self, tparam):
        """Calculate the one-body term of the 1D Hubbard Hamiltonian."""
        raise NotImplementedError

    def compute_er(self, uparam):
        """Calculate the the-body term of the 1D Hubbard Hamiltonian."""
        raise NotImplementedError

    def compute_overlap(self):
        """Calculate overlap of the 1D Hubbard Hamiltonian, (identity matrix)."""
        raise NotImplementedError


class Hubbard(PhysModHam):
    """The 1-D Hubbard class Hamiltonian."""

    def compute_kinetic(self, tparam):
        """Calculate the one-body term of the 1D Hubbard Hamiltonian."""
        result = np.zeros((self.nbasis, self.nbasis))
        for i in range(self.nbasis - 1):
            result[i, i + 1] = tparam
            result[i + 1, i] = tparam
        if  self.pbc == True:
            result[self.nbasis - 1, 0] = tparam
            result[0, self.nbasis - 1] = tparam
        return result

    def compute_er(self, uparam):
        """Calculate the the-body term of the 1D Hubbard Hamiltonian."""
        result = np.zeros((self.nbasis, self.nbasis, self.nbasis, self.nbasis))
        for i in range(self.nbasis):
            result[i, i, i, i] = uparam
        return result

    def compute_overlap(self):
        """Calculate overlap of the 1D Hubbard Hamiltonian, (identity matrix)."""
        return np.identity(self.nbasis)
