# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2015 The HORTON Development Team
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
"""One- and two-orbital entanglement measures

   Abbreviations used in this module:

   * dm1 = a 1-RDM
   * dm2 = a 2-RDM
   * dm3 = a 3-RDM
   * dm4 = a 4-RDM
   * odm1 = one-orbital reduced density matrix
   * odm2 = two-orbital reduced density matrix
   * soentropy = single-orbital entropy
   * toentropy = two-orbital entropy
   * mutualinfo = mutual information
"""


import numpy as np
from horton.log import log
from horton.utils import check_options


__all__ = [
    'OrbitalEntanglement',
    'OrbitalEntanglementAp1rog',
]



class OrbitalEntanglement(object):
    def __init__(self, lf, one_dm, two_dm, three_dm=None, four_dm=None):
        '''
            **Arguments:**

            lf
                Instance of :py:class:`horton.matrix.dense.DenseLinalgFactory` or :py:class:`horton.matrix.cholesky.CholeskyLinalgFactory`

            one_dm
                Instance of :py:class:`horton.matrix.base.OneIndex`

                A 1-RDM

            two_dm
                List of instances of :py:class:`horton.matrix.base.TwoIndex`

                List of 2-RDM's. e.g. [:math:`\Gamma_{pp}^{qq}`, :math:`\Gamma_{pq}^{pq}`]

           **Optional arguments:**

           three_dm
                Instance of :py:class:`horton.matrix.base.ThreeIndex`

                A 3-RDM.

           four_dm
                Instance of :py:class:`horton.matrix.base.FourIndex`

                A 4-RDM.
        '''
        self._lf = lf
        self._nbasis = lf.default_nbasis
        self._dm1 = one_dm
        self._dm2 = two_dm
        self._dm3 = three_dm
        self._dm4 = four_dm
        self._odm1 = []
        self._odm2 = []
        self._soentropy = lf.create_one_index()
        self._toentropy = lf.create_two_index()
        self._mutualinfo = lf.create_two_index()

    def __call__(self):
        """dumps single-orbital entropy and orbital-pair mututal information

        see :py:meth:`horton.orbital_entanglement.orbital_entanglement.OrbitalEntanglement.dump_output` for more info
        """

        if log.do_medium:
            log(' ')
            log('Calculating orbital entanglement measures')
            log(' ')

        #
        # Compute single-orbital entropy and mutual information
        #
        if log.do_medium:
            log('  Computing s(1) and I_ij')
        self.compute_single_orbital_entropy()
        self.compute_two_orbital_entropy()
        self.compute_mutual_information()

        #
        # Dump output to file:
        #
        if log.do_medium:
            log('  Dumping output files')
            log(' ')
            log.hline('=')
        self.dump_output()

    def _get_lf(self):
        '''The LinalgFactory.'''
        return self._lf

    lf = property(_get_lf)

    def _get_nbasis(self):
        '''The number of basis functions.'''
        return self._nbasis

    nbasis = property(_get_nbasis)

    def _get_dm1(self):
        '''Some input 1-RDM'''
        return self._dm1

    dm1 = property(_get_dm1)

    def _get_dm2(self):
        '''Some input 2-RDM'''
        return self._dm2

    dm2 = property(_get_dm2)

    def _get_dm3(self):
        '''Some input 3-RDM'''
        return self._dm3

    dm3 = property(_get_dm3)

    def _get_dm4(self):
        '''Some input 4-RDM'''
        return self._dm4

    dm4 = property(_get_dm4)

    def _get_odm1(self):
        '''The 1-ORDM'''
        return self._odm1

    odm1 = property(_get_odm1)

    def _get_odm2(self):
        '''The 2-ORDM'''
        return self._odm2

    odm2 = property(_get_odm2)

    def _get_soentropy(self):
        '''The single-orbital entropy'''
        return self._soentropy

    so_entropy = property(_get_soentropy)

    def _get_toentropy(self):
        '''The two-orbital entropy'''
        return self._toentropy

    to_entropy = property(_get_toentropy)

    def _get_mutualinfo(self):
        '''The mutual information'''
        return self._mutualinfo

    mutual_info = property(_get_mutualinfo)

    def _get_cache(self):
        '''The cache object in which the main attributes are stored'''
        return self._cache

    cache = property(_get_cache)

    def __clear__(self):
        self.clear()

    def clear(self):
        '''Clear all information'''
        self._cache.clear()

    def clear_dm(self):
        '''Clear the orbital density matrices'''
        self._odm1 = []
        self._odm2 = []

    def append_odm1(self, index, matrix):
        '''Append index and one-orbital-reduced density matrix of orbital
           'index' to list.

           **Arguments:**

           index
                Orbital index.

           matrix
                The OneIndex object containing the ODMs.

           **Optional arguments:**

        '''
        self._odm1.append((index, matrix))

    def append_odm2(self, index1, index2, matrix):
        '''Append indices and two-orbital-reduced density matrix of orbital
           pair 'index1/index2' to list.

           **Arguments:**

           index1
                First orbital index.

           index2
                Second orbital index.

           matrix
                The OneIndex object containing the ODMs.

           **Optional arguments:**

        '''
        self._odm2.append((index1, index2, matrix))

    def compute_odm1(self, index1):
        raise NotImplementedError

    def compute_odm2(self, index1, index2):
        raise NotImplementedError

    def calculate_entropy_term(self, val, select='vonNeumann'):
        '''Calculate entropic term

           **Arguements**

           val
                Used to determine entropy

           **Optional arguments:**

           select
                Select entropy function. Default: von Neumann.
        '''
        check_options('select', select, 'vonNeumann')
        if select=='vonNeumann':
            if val > 0.0:
                return np.log(val)*val
            else:
                if abs(val) > 1e-6:
                    log('Neglecting negative value %f in entropy function' % val)
                return 0.0

    def compute_single_orbital_entropy(self, select='vonNeumann'):
        '''Compute single-orbital entropy for each orbital in the active space.
           Currently, only the von Neumann entropy is supported.

           The 1-ODM is assumed to be diagonalized.

           **Optional arguments:**

           select
                Select entropy function. Default: von Neumann.
        '''
        check_options('select', select, 'vonNeumann')
        for index in range(self.nbasis):
            self.compute_odm1(index)
        for item in self.odm1:
            mat = item[1]
            term = 0.0
            for ind in range(mat.shape[0]):
                term -= self.calculate_entropy_term(mat.get_element(ind), select)
            self.so_entropy.set_element(item[0], term)

    def compute_two_orbital_entropy(self, select='vonNeumann'):
        '''Compute two-orbital entropy for each orbital in the active space.
           Currently, only the von Neumann entropy is supported.

           The 1-ODM and 2-ODM are assumed to be diagonalized.

           **Optional arguments:**

           select
                Select entropy function. Default: von Neumann.
        '''
        check_options('select', select, 'vonNeumann')
        for index1 in range(self.nbasis):
            for index2 in range(self.nbasis):
                if index2 is not index1:
                    self.compute_odm2(index1, index2)
        for item in self.odm2:
            mat = item[2]
            term = 0.0
            for ind in range(mat.shape[0]):
                term -= self.calculate_entropy_term(mat.get_element(ind), select)
            self.to_entropy.set_element(item[0], item[1], term)

    def compute_mutual_information(self):
        '''Compute mutual information using the single-orbital entropy and the
           two-orbital entropy.

           **Arguments:**

           one_entropy
                Single-orbital entropy.

           two_entropy
                Two-orbital entropy.

           **Optional arguments:**

        '''
        self.mutual_info.assign(self.to_entropy)
        self.mutual_info.iadd(self.so_entropy, -1.0)
        self.mutual_info.iadd_t(self.so_entropy, -1.0)
        self.mutual_info.iscale(-0.5)
        self.mutual_info.assign_diagonal(0.0)
        assert self.mutual_info.is_symmetric()

    def diagonalize(self, mat):
        '''Returns eigenvalues of a TwoIndex instance. Only real eigenvalues
           are returned, any imaginary part will be ignored

           **Arguments:**

           mat
                An TwoIndex instance to be diagonalized.
        '''
        evalue = mat.diagonalize()
        #
        # Delete negative eigenvalues.
        #
        if np.amin(np.real(evalue)) < 0.0:
            for i in evalue:
                if i < 0.0 and abs(i) > 1e-5:
                    log('Warning, negative eigenvalue of %f' %i)
            evalue[np.where(np.real(evalue) < 0.0)] = 0.0

        out = self.lf.create_one_index(mat.nbasis)
        out.assign(np.real(evalue))
        return out

    def dump_output(self, file1='s1.dat', file2='i12.dat'):
        ''' Dump entanglement output files for postprocessing.
            Output files can be visualized using the
            build_orbital_entanglement_diagrams.sh script, which uses gnuplot.

            **Optional arguments:**

            file1
                string

                Filename for storing single-orbital entropy

                Default: s1.dat

                Format: First column is orbital index. Second column is single-orbital entropy

            file2
                string

                Filename for storing orbital-pair mutual information

                Default: i12.dat

                Format: First two columns are orbital indices. Second column is corresponding orbital-pair mutual information
        '''
        #
        # Write single-orbital entropy to file:
        #
        f = open(file1, 'w')
        for i in range(self.nbasis):
            f.write('%3i  %1.12f \n' %(i+1, self.so_entropy.get_element(i)))
        f.close()
        #
        # Write mutual information to file:
        #
        f2 = open(file2, 'w')
        for i in range(self.nbasis):
            for j in range(i+1, self.nbasis):
                f2.write('%3i  %3i  %1.12f \n' %(i+1, j+1, self.mutual_info.get_element(i,j)))
        f2.close()


class OrbitalEntanglementAp1rog(OrbitalEntanglement):
    '''Orbital entanglement class for AP1roG. Note that the response 2-RDM
       is a list with elements 'dm2_ppqq' [0], and 'dm2_pqpq' [1].
    '''
    def compute_odm1(self, index1):
        '''Compute ODM for orbital index1.

           **Arguments:**

           index1
                First orbital index.
        '''
        odmat = self.lf.create_one_index(2)
        term = self.dm1.get_element(index1)
        odmat.set_element(0, (1-term))
        odmat.set_element(1, (term))
        self.append_odm1(index1, odmat)

    def compute_odm2(self, index1, index2):
        '''Compute 2-ODM for orbital indices 'index1/index2'.

           **Arguments:**

           index1
                First orbital index.

           index2
                Second orbital index.
        '''
        mat = self.lf.create_two_index(4, 4)
        term = 1.0-self.dm1.get_element(index1) \
                  -self.dm1.get_element(index2) \
                  +self.dm2[1].get_element(index1, index2)
        mat.set_element(0, 0, term)
        mat.set_element(1, 1, self.dm2[1].get_element(index1, index2))
        mat.set_element(2, 3, -self.dm2[0].get_element(index2, index1))
        mat.set_element(3, 2, -self.dm2[0].get_element(index1, index2))
        term = self.dm1.get_element(index2)-self.dm2[1].get_element(index1, index2)
        mat.set_element(2, 2, term)
        term = self.dm1.get_element(index1)-self.dm2[1].get_element(index1, index2)
        mat.set_element(3, 3, term)

        sol = self.diagonalize(mat)
        self.append_odm2(index1, index2, (sol))
