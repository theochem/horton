# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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
"""Orbital class."""

import numpy as np
from scipy.linalg import eigh

from horton.utils import check_type
from horton.log import log

__all__ = ['Orbitals']


class Orbitals(object):
    """Orbital coefficient, energies and occupation numbers (single spin channel)."""

    #
    # Constructor and destructor
    #

    def __init__(self, nbasis, nfn=None):
        """Initialize an Orbitals object.

        Parameters
        ----------
        nbasis : int
            The number of basis functions.
        nfn : int
            The number of functions to store. Defaults to nbasis.
        """
        if nfn is None:
            nfn = nbasis
        self._coeffs = np.zeros((nbasis, nfn))
        self._energies = np.zeros(nfn)
        self._occupations = np.zeros(nfn)

    #
    # Properties
    #

    def _get_nbasis(self):
        """The number of basis functions."""
        return self._coeffs.shape[0]

    nbasis = property(_get_nbasis)

    def _get_nfn(self):
        """The number of orbitals (or functions in general)."""
        return self._coeffs.shape[1]

    nfn = property(_get_nfn)

    def _get_coeffs(self):
        """The matrix with the expansion coefficients."""
        return self._coeffs.view()

    coeffs = property(_get_coeffs)

    def _get_energies(self):
        """The orbital energies."""
        return self._energies.view()

    energies = property(_get_energies)

    def _get_occupations(self):
        """The orbital occupations."""
        return self._occupations.view()

    occupations = property(_get_occupations)

    #
    # Methods
    #

    def __eq__(self, other):
        """Compare self with other."""
        return isinstance(other, Orbitals) and \
               other.nbasis == self.nbasis and \
               other.nfn == self.nfn and \
               (other._coeffs == self._coeffs).all() and \
               (other._energies == self._energies).all() and \
               (other._occupations == self._occupations).all()

    @classmethod
    def from_hdf5(cls, grp):
        """Construct an instance from data previously stored in an h5py.Group.

        Parameters
        ----------
        grp : h5py.Group
            Group used to take all data to initialize an Orbitals object
        """
        if grp.attrs['class'] != cls.__name__:
            raise TypeError('The class of the expansion in the HDF5 file does not match.')
        nbasis, nfn = grp['coeffs'].shape
        result = cls(nbasis, nfn)
        grp['coeffs'].read_direct(result._coeffs)
        grp['energies'].read_direct(result._energies)
        grp['occupations'].read_direct(result._occupations)
        return result

    def to_hdf5(self, grp):
        """Dump this object in an h5py.Group.

        Parameters
        ----------
        grp : h5py.Group
            Destination where the data are stored.
        """
        grp.attrs['class'] = self.__class__.__name__
        grp['coeffs'] = self._coeffs
        grp['energies'] = self._energies
        grp['occupations'] = self._occupations

    def clear(self):
        """Reset all elements to zero."""
        self._coeffs[:] = 0.0
        self._energies[:] = 0.0
        self._occupations[:] = 0.0

    def copy(self):
        """Return a copy of the object."""
        result = Orbitals(self.nbasis, self.nfn)
        result._coeffs[:] = self._coeffs
        result._energies[:] = self._energies
        result._occupations[:] = self._occupations
        return result

    def assign(self, other):
        """Assign with the contents of another object.

        Parameters
        ----------
        other : Orbitals
            Another Orbitals object
        """
        check_type('other', other, Orbitals)
        self._coeffs[:] = other._coeffs
        self._energies[:] = other._energies
        self._occupations[:] = other._occupations

    def randomize(self):
        """Fill with random normal data."""
        self._coeffs[:] = np.random.normal(0, 1, self._coeffs.shape)
        self._energies[:] = np.random.normal(0, 1, self._energies.shape)
        self._occupations[:] = np.random.normal(0, 1, self._occupations.shape)

    def permute_basis(self, permutation):
        """Reorder the coefficients for a given permutation of basis functions (rows).

        Parameters
        ----------
        permutation : np.ndarray, dtype=int, shape=(nbasis,)
            An array that defines the new order of the basis functions.
        """
        self._coeffs[:] = self.coeffs[permutation]

    def permute_orbitals(self, permutation):
        """Reorder the coefficients for a given permutation of orbitals (columns).

        Parameters
        ----------
        permutation : np.ndarray, dtype=int, shape=(nbasis,)
            An array that defines the new order of the orbitals.
        """
        self._coeffs[:] = self.coeffs[:, permutation]

    def change_basis_signs(self, signs):
        """Correct for different sign conventions of the basis functions.

        Parameters
        ----------
        signs : np.ndarray, dtype=int, shape=(nbasis,)
            An array with sign changes indicated by +1 and -1.
        """
        self._coeffs *= signs.reshape(-1, 1)

    def check_normalization(self, overlap, eps=1e-4):
        """Check that the occupied orbitals are normalized.

        When the orbitals are not normalized, an AssertionError is raised.

        Parameters
        ----------
        overlap : np.ndarray, shape=(nbasis, nbasis)
            The overlap matrix.
        eps : float
            The allowed deviation from unity, very loose by default.
        """
        for i in range(self.nfn):
            if self.occupations[i] == 0:
                continue
            norm = np.dot(self._coeffs[:, i], np.dot(overlap, self._coeffs[:, i]))
            # print i, norm
            assert abs(norm - 1) < eps, 'The orbitals are not normalized!'

    def check_orthonormality(self, overlap, eps=1e-4):
        """Check that the occupied orbitals are orthogonal and normalized.

        When the orbitals are not orthonormal, an AssertionError is raised.

        Parameters
        ----------
        overlap : np.ndarray, shape=(nbasis, nbasis)
            The overlap matrix.
        eps : float
            The allowed deviation from unity, very loose by default.
        """
        for i0 in range(self.nfn):
            if self.occupations[i0] == 0:
                continue
            for i1 in range(i0 + 1):
                if self.occupations[i1] == 0:
                    continue
                dot = np.dot(self._coeffs[:, i0], np.dot(overlap, self._coeffs[:, i1]))
                if i0 == i1:
                    assert abs(dot - 1) < eps
                else:
                    assert abs(dot) < eps

    def error_eigen(self, fock, overlap):
        """Compute the error of the orbitals with respect to the eigenproblem.

        Parameters
        ----------
        fock : np.ndarray, shape=(nbasis, nbasis)
            The fock matrix.
        overlap : np.ndarray, shape=(nbasis, nbasis)
            The overlap matrix.

        Returns
        -------
        error : float
            The RMSD error on the orbital energies.
        """
        errors = np.dot(fock, (self.coeffs)) \
                 - self.energies * np.dot(overlap, (self.coeffs))
        return np.sqrt((abs(errors) ** 2).mean())

    def from_fock(self, fock, overlap):
        """Diagonalize a Fock matrix to obtain orbitals and energies.

        This method updated the attributes ``coeffs`` and ``energies`` in-place.

        Parameters
        ----------
        fock : np.ndarray, shape=(nbasis, nbasis)
            The fock matrix.
        overlap : np.ndarray, shape=(nbasis, nbasis)
            The overlap matrix.
        """
        evals, evecs = eigh(fock, overlap)
        self._energies[:] = evals[:self.nfn]
        self._coeffs[:] = evecs[:, :self.nfn]

    def from_fock_and_dm(self, fock, dm, overlap, epstol=1e-8):
        """Combined diagonalization of a Fock and a density matrix.

        This routine first diagonalizes the Fock matrix to obtain orbitals and orbital
        energies. Then, using first order (degenerate) perturbation theory, the occupation
        numbers are computed and, if needed, the the degeneracies of the Fock orbitals are
        lifted. It is assumed that the Fock and the density matrices commute. This method
        updated the attributes ``coeffs``, ``energies`` and ``occupations`` in-place.

        Parameters
        ----------
        fock : np.ndarray, shape=(nbasis, nbasis)
            The fock matrix.
        dm : np.ndarray, shape=(nbasis, nbasis)
            The density matrix.
        overlap : np.ndarray, shape=(nbasis, nbasis)
            The overlap matrix.
        epstol : float
            The threshold for recognizing degenerate energy levels. When two subsequent
            energy levels are separated by an orbital energy less than ``epstol``, they
            are considered to be degenerate. When a series of energy levels have an
            orbital energy spacing between subsequent levels that is smaller than
            ``epstol``, they are all considered to be part of the same degenerate group.
            For every degenerate set of orbitals, the density matrix is used to (try to)
            lift the degeneracy.
        """
        # Diagonalize the Fock Matrix
        self.from_fock(fock, overlap)

        # Build clusters of degenerate orbitals. Rely on the fact that the
        # energy levels are sorted (one way or the other).
        clusters = []
        begin = 0
        for ifn in range(1, self.nfn):
            if abs(self.energies[ifn] - self.energies[ifn - 1]) > epstol:
                end = ifn
                clusters.append([begin, end])
                begin = ifn
        end = self.nfn
        clusters.append([begin, end])

        # Lift degeneracies using the density matrix
        sds = np.dot(overlap.T, np.dot(dm, overlap))
        for begin, end in clusters:
            if end - begin == 1:
                self.occupations[begin] = np.dot(self.coeffs[:, begin], np.dot(sds, self.coeffs[:, begin]))
            else:
                # Build matrix
                mat = np.dot(self.coeffs[:, begin:end].T, np.dot(sds, self.coeffs[:, begin:end]))
                # Diagonalize and reverse order
                evals, evecs = np.linalg.eigh(mat)
                evals = evals[::-1]
                evecs = evecs[:, ::-1]
                # Rotate the orbitals
                self.coeffs[:, begin:end] = np.dot(self.coeffs[:, begin:end], evecs)
                # Compute expectation values
                self.occupations[begin:end] = evals
                for i0 in range(end - begin):
                    self.energies[begin + i0] = np.dot(self.coeffs[:, begin + i0],
                                                       np.dot(fock, self.coeffs[:, begin + i0]))

    def derive_naturals(self, dm, overlap):
        """Derive natural orbitals from a given density matrix and assign the result to self.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis)
            The density matrix.
        overlap : np.ndarray, shape=(nbasis, nbasis)
            The overlap matrix
        """
        # Transform density matrix to Fock-like form
        sds = np.dot(overlap.T, np.dot(dm, overlap))
        # Diagonalize and compute eigenvalues
        evals, evecs = eigh(sds, overlap)
        self._coeffs[:] = evecs[:, :self.nfn]
        self._occupations[:] = evals
        self._energies[:] = 0.0

    def get_homo_index(self, offset=0):
        """Return the index of a HOMO orbital.

        Parameters
        ----------
        offset : int
            By default, the (highest) homo energy is returned. When this index is above
            zero, the corresponding lower homo energy is returned.
        """
        if offset < 0:
            raise ValueError('Offset must be zero or positive.')
        homo_indexes = self.occupations.nonzero()[0]
        if len(homo_indexes) > offset:
            return homo_indexes[len(homo_indexes) - offset - 1]

    homo_index = property(get_homo_index)

    def get_homo_energy(self, offset=0):
        """Return a homo energy.

        Parameters
        ----------
        offset : int
            By default, the (highest) homo energy is returned. When this index is above
            zero, the corresponding lower homo energy is returned.
        """
        index = self.get_homo_index(offset)
        if index is not None:
            return self.energies[index]

    homo_energy = property(get_homo_energy)

    def get_lumo_index(self, offset=0):
        """Return the index of a LUMO orbital.

        Parameters
        ----------
        offset : int
            By default, the (lowest) lumo energy is returned. When this index is above
            zero, the corresponding higher homo energy is returned.
        """
        if offset < 0:
            raise ValueError('Offset must be zero or positive.')
        lumo_indexes = (self.occupations == 0.0).nonzero()[0]
        if len(lumo_indexes) > offset:
            return lumo_indexes[offset]

    lumo_index = property(get_lumo_index)

    def get_lumo_energy(self, offset=0):
        """Return a lumo energy.

        Parameters
        ----------
        offset : int
            By default, the (lowest) lumo energy is returned. When this index is above
            zero, the corresponding higher homo energy is returned.
        """
        index = self.get_lumo_index(offset)
        if index is not None:
            return self.energies[index]

    lumo_energy = property(get_lumo_energy)

    def to_dm(self, other=None):
        """Compute the density matrix.

        Parameters
        ----------
        other : Orbitals
            Another Orbitals object to construct a transfer-density matrix.

        Returns
        -------
        dm : np.ndarray, shape=(nbasis, nbasis)
            The density matrix.
        """
        if other is None:
            return np.dot(self._coeffs * self.occupations, self._coeffs.T)
        else:
            return np.dot(self._coeffs * (self.occupations * other.occupations) ** 0.5, other._coeffs.T)

    def rotate_random(self):
        """Apply random unitary transformation distributed with Haar measure.

        The attributes ``energies`` and ``occupations`` are not altered.
        """
        z = np.random.normal(0, 1, (self.nfn, self.nfn))
        q, r = np.linalg.qr(z)
        self.coeffs[:] = np.dot(self.coeffs, q)

    def rotate_2orbitals(self, angle=0.7853981633974483, index0=None, index1=None):
        """Rotate two orbitals.

        Parameters
        ----------
        angle : float
            The rotation angle, defaults to 45 deg.
        index0, index1 : int
            The orbitals to rotate, defaults to HOMO and LUMO,

        The attributes ``energies`` and ``occupations`` are not altered.
        """
        if index0 == None:
            index0 = self.homo_index
        if index1 == None:
            index1 = self.lumo_index
        old0 = self.coeffs[:, index0].copy()
        old1 = self.coeffs[:, index1].copy()
        self.coeffs[:, index0] = np.cos(angle) * old0 - np.sin(angle) * old1
        self.coeffs[:, index1] = np.sin(angle) * old0 + np.cos(angle) * old1

    def swap_orbitals(self, swaps):
        """Change the order of the orbitals using pair-exchange.

        Parameters
        ----------
        swaps : np.ndarray, shape=(m, 2), dtype=int
            An integer numpy array with two columns where every row corresponds to one
            swap.

        The attributes ``energies`` and ``occupations`` are also reordered.
        """
        if not (swaps.shape[1] == 2 and swaps.ndim == 2 and np.issubdtype(swaps.dtype, np.int)):
            raise TypeError('The argument swaps has the wrong shape/type.')
        for iswap in range(len(swaps)):
            index0, index1 = swaps[iswap]
            if log.do_medium:
                log('  Swapping orbitals %i and %i' % (index0, index1))
            tmp = self.coeffs[:, index0].copy()
            self.coeffs[:, index0] = self.coeffs[:, index1]
            self.coeffs[:, index1] = tmp
            self.energies[index0], self.energies[index1] = \
                self.energies[index1], self.energies[index0]
            self.occupations[index0], self.occupations[index1] = \
                self.occupations[index1], self.occupations[index0]
