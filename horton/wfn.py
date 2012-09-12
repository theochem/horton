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
   in a certain basis.
"""


import numpy as np


__all__ = ['BaseWFN', 'ClosedShellWFN', 'OpenShellWFN']


# TODO: move common code to Base
class BaseWFN(object):
    @classmethod
    def from_hdf5(cls, grp, lf):
        clsname = grp.attrs['class']
        if clsname == 'ClosedShellWFN':
            return ClosedShellWFN.from_hdf5(grp, lf)
        elif clsname == 'OpenShellWFN':
            return OpenShellWFN.from_hdf5(grp, lf)
        else:
            raise TypeError('Can not load wavefunction of class %s from HDF5.' % clsname)


# TODO: underscore convention get methods
class ClosedShellWFN(BaseWFN):
    def __init__(self, nep, lf, nbasis, norb=None):
        """
           **Arguments:**

           nep
                The number of electron pairs in the wave function.

           lf
                A LinalgFactor instance

           nbasis
                The number of basis functions.

           **Optional arguments:**

           norb
               the number of orbitals (occupied + virtual). When not given,
               it is set to nbasis.
        """
        self._nep = nep
        self._lf = lf
        self._nbasis = nbasis
        if norb is None:
            self._norb = nbasis
        else:
            self._norb = norb

        if self.nep <= 0:
            raise ValueError('At least one pair of electrons is required.')
        if self.nbasis < self.nep:
            raise ValueError('The number of spatial basis functions must not be lower than the number of electron pairs.')

        self._expansion = lf.create_expansion(self.nbasis, self.norb, do_energies=True)

    @classmethod
    def from_hdf5(cls, grp, lf):
        result = ClosedShellWFN(grp['nep'][()], lf, nbasis=grp['nbasis'][()], norb=grp['norb'][()])
        result._expansion.read_from_hdf5(grp['expansion'])
        return result

    def to_hdf5(self, grp):
        grp['nep'] = self._nep
        grp['nbasis'] = self._nbasis
        grp['norb'] = self._norb
        grp_expansion = grp.create_group('expansion')
        self._expansion.to_hdf5(grp_expansion)

    def _get_nel(self):
        '''The number of electrons'''
        return self._nep*2

    nel = property(_get_nel)

    def get_nep(self):
        return self._nep

    nep = property(get_nep)

    def get_nbasis(self):
        return self._nbasis

    nbasis = property(get_nbasis)

    def get_norb(self):
        return self._norb

    norb = property(get_norb)

    def get_expansion(self):
        return self._expansion

    expansion = property(get_expansion)

    def iter_expansions(self):
        yield self._expansion, self._nep, 2

    def compute_alpha_density_matrix(self, dm):
        """Compute the alpha density matrix

           **Arguments:**

           dm
                An output density matrix. This must be an instance of the
                One-body operator class of the linalg factory, self.lf.
        """
        self._expansion.compute_density_matrix(self.nep, dm)

    def compute_density_matrix(self, dm):
        """Compute the density matrix

           **Arguments:**

           dm
                An output density matrix. This must be an instance of the
                One-body operator class of the linalg factory, self.lf.
        """
        self.compute_alpha_density_matrix(dm)
        dm.iscale(2)

    def apply_basis_permutation(self, permutation):
        """Reorder the expansion coefficients"""
        self._expansion.apply_basis_permutation(permutation)


class OpenShellWFN(BaseWFN):
    def __init__(self, nalpha, nbeta, lf, nbasis, norb=None):
        """
           An unrestricted open-shell wavefunction.

           **Arguments:**

           nalpha
                The number of alpha electrons in the wave function.

           nbeta
                The number of beta electrons in the wave function.

           lf
                A LinalgFactor instance

           nbasis
                The number of basis functions.

           **Optional arguments:**

           norb
               the number of orbitals (occupied + virtual). When not given,
               it is set to nbasis.
        """
        self._nalpha = nalpha
        self._nbeta = nbeta
        self._lf = lf
        self._nbasis = nbasis
        if norb is None:
            self._norb = nbasis
        else:
            self._norb = norb

        if self.nalpha + self.nbeta <= 0:
            raise ValueError('At least one alpha or beta electron is required.')
        if self.nbasis < self.nalpha or self.nbasis < self.nbeta:
            raise ValueError('The number of spatial basis functions must not be lower than the number of alpha or beta electrons.')

        self._alpha_expansion = lf.create_expansion(self.nbasis, self.norb, do_energies=True)
        self._beta_expansion = lf.create_expansion(self.nbasis, self.norb, do_energies=True)

    @classmethod
    def from_hdf5(cls, grp, lf):
        result = OpenShellWFN(grp['nalpha'][()], grp['nbeta'][()], lf, nbasis=grp['nbasis'][()], norb=grp['norb'][()])
        result._alpha_expansion.read_from_hdf5(grp['alpha_expansion'])
        result._beta_expansion.read_from_hdf5(grp['beta_expansion'])
        return result

    def to_hdf5(self, grp):
        grp['nalpha'] = self._nalpha
        grp['nbeta'] = self._nbeta
        grp['nbasis'] = self._nbasis
        grp['norb'] = self._norb
        grp_alpha_expansion = grp.create_group('alpha_expansion')
        self._alpha_expansion.to_hdf5(grp_alpha_expansion)
        grp_beta_expansion = grp.create_group('beta_expansion')
        self._beta_expansion.to_hdf5(grp_beta_expansion)

    def _get_nel(self):
        '''The number of electrons'''
        return self._nalpha + self._nbeta

    nel = property(_get_nel)

    def get_nalpha(self):
        return self._nalpha

    nalpha = property(get_nalpha)

    def get_nbeta(self):
        return self._nbeta

    nbeta = property(get_nbeta)

    def get_nbasis(self):
        return self._nbasis

    nbasis = property(get_nbasis)

    def get_norb(self):
        return self._norb

    norb = property(get_norb)

    def get_alpha_expansion(self):
        return self._alpha_expansion

    alpha_expansion = property(get_alpha_expansion)

    def get_beta_expansion(self):
        return self._beta_expansion

    beta_expansion = property(get_beta_expansion)

    def iter_expansions(self):
        yield self._alpha_expansion, self._nalpha, 1
        yield self._beta_expansion, self._nbeta, 1

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


    def apply_basis_permutation(self, permutation):
        """Reorder the expansion coefficients"""
        self._alpha_expansion.apply_basis_permutation(permutation)
        self._beta_expansion.apply_basis_permutation(permutation)
