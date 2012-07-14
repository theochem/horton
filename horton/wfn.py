# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


__all__ = ['ClosedShellWFN', 'OpenShellWFN']


class ClosedShellWFN(object):
    def __init__(self, nep, lf, basis=None, nbasis=None, norb=None):
        """
           **Arguments:**

           nep
                The number of electron pairs in the wave function.

           lf
                A LinalgFactor instance

           **Optional arguments:**

           basis
                A specification of the basis set.

           nbasis
                The number of basis functions.

           norb
               the number of orbitals (occupied + virtual). When not given,
               it is set to nbasis.

           Either basis or nbasis must be given.
        """
        self._nep = nep
        self._lf = lf
        self._basis = basis
        self._nbasis = nbasis
        self._norb = norb

        if self.nep <= 0:
            raise ValueError('At least one pair of electrons is required.')
        if self.nbasis < self.nep:
            raise ValueError('The number of spatial basis functions must not be lower than the number of electron pairs.')

        self._expansion = lf.create_expansion(self.nbasis, self.norb, do_energies=True)

    def get_nep(self):
        return self._nep

    nep = property(get_nep)

    def get_basis(self):
        return self._basis

    basis = property(get_basis)

    def get_nbasis(self):
        if self.basis is not None:
            return self._basis.nbasis
        else:
            return self._nbasis

    nbasis = property(get_nbasis)

    def get_norb(self):
        if self._norb is None:
            return self.nbasis
        else:
            return self._norb

    norb = property(get_norb)

    def get_expansion(self):
        return self._expansion

    expansion = property(get_expansion)

    def compute_density_matrix(self, dm):
        """Compute the density matrix

           **Arguments:**

           dm
                An output density matrix. This must be an instance of the
                One-body operator class of the linalg factory, self.lf.
        """
        self._expansion.compute_density_matrix(self.nep, dm)


class OpenShellWFN(object):
    def __init__(self, nalpha, nbeta, lf, basis=None, nbasis=None, norb=None):
        """
           An unrestricted open-shell wavefunction.

           **Arguments:**

           nalpha
                The number of alpha electrons in the wave function.

           nbeta
                The number of beta electrons in the wave function.

           lf
                A LinalgFactor instance

           **Optional arguments:**

           basis
                A specification of the basis set.

           nbasis
                The number of basis functions.

           norb
               the number of orbitals (occupied + virtual). When not given,
               it is set to nbasis.

           Either basis or nbasis must be given.
        """
        self._nalpha = nalpha
        self._nbeta = nbeta
        self._lf = lf
        self._basis = basis
        self._nbasis = nbasis
        self._norb = norb

        if self.nalpha + self.nbeta <= 0:
            raise ValueError('At least one alpha or beta electron is required.')
        if self.nbasis < self.nalpha or self.nbasis < self.nbeta:
            raise ValueError('The number of spatial basis functions must not be lower than the number of alpha or beta electrons.')

        self._alpha_expansion = lf.create_expansion(self.nbasis, self.norb, do_energies=True)
        self._beta_expansion = lf.create_expansion(self.nbasis, self.norb, do_energies=True)

    def get_nalpha(self):
        return self._nalpha

    nalpha = property(get_nalpha)

    def get_nbeta(self):
        return self._nbeta

    nbeta = property(get_nbeta)

    def get_basis(self):
        return self._basis

    basis = property(get_basis)

    def get_nbasis(self):
        if self.basis is not None:
            return self._basis.nbasis
        else:
            return self._nbasis

    nbasis = property(get_nbasis)

    def get_norb(self):
        if self._norb is None:
            return self.nbasis
        else:
            return self._norb

    norb = property(get_norb)

    def get_alpha_expansion(self):
        return self._alpha_expansion

    alpha_expansion = property(get_alpha_expansion)

    def get_beta_expansion(self):
        return self._beta_expansion

    beta_expansion = property(get_beta_expansion)

    def compute_density_matrix(self, dm):
        """Compute the density matrix

           **Arguments:**

           dm
                An output density matrix. This must be an instance of the
                One-body operator class of the linalg factory, self.lf.
        """
        self._alpha_expansion.compute_density_matrix(self.nalpha, dm)
        self._beta_expansion.compute_density_matrix(self.nbeta, dm, factor=1)

    def compute_spin_density_matrix(self, dm):
        """Compute the spin density matrix

           **Arguments:**

           dm
                An output density matrix. This must be an instance of the
                One-body operator class of the linalg factory, self.lf.
        """
        self._alpha_expansion.compute_density_matrix(self.nalpha, dm)
        self._beta_expansion.compute_density_matrix(self.nbeta, dm, factor=-1)

    def compute_alpha_density_matrix(self, dm):
        """Compute the alpha density matrix

           **Arguments:**

           dm
                An output density matrix. This must be an instance of the
                One-body operator class of the linalg factory, self.lf.
        """
        self._alpha_expansion.compute_density_matrix(self.nalpha, dm)

    def compute_beta_density_matrix(self, dm):
        """Compute the beta density matrix

           **Arguments:**

           dm
                An output density matrix. This must be an instance of the
                One-body operator class of the linalg factory, self.lf.
        """
        self._beta_expansion.compute_density_matrix(self.nbeta, dm)
