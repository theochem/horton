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
'''Correlated wavefunction implementations

   This is a general geminals class.
'''

import numpy as np
import math as math

from horton.cache import Cache
from horton.matrix import DenseTwoIndex, TwoIndex, DenseExpansion, Expansion
from horton.log import log, timer
from horton.utils import check_options, check_type

from itertools import permutations
from operator import mul

__all__ = [
    'Geminal',
]


class PropertyHelper(object):
    '''Auxiliary class to set up x_dm_y attributes of Geminal class.'''
    def __init__(self, method, arg, doc):
        self.method = method
        self.arg = arg
        self.__doc__ = doc

    def __get__(self, obj, objtype):
        #
        # For doc strings:
        #
        if obj is None:
            return self
        #
        # For actual use
        #
        try:
            return self.method(obj, self.arg)
        except KeyError, e:
            raise AttributeError('The requested attribute (%s) is not available.' % e.args[0])


class Geminal(object):
    '''A collection of geminals and optimization routines.

       This is just a base class that serves as a template for
       specific implementations.
    '''
    def __init__(self, lf, occ_model, npairs=None, nvirt=None):
        '''
           **Arguments:**

           lf
                A LinalgFactory instance.

           occ_model
                Occupation model

           **Optional arguments:**

           npairs
                Number of electron pairs, if not specified,
                npairs = number of occupied orbitals

           nvirt
                Number of virtual orbitals, if not specified,
                nvirt = (nbasis-npairs)
        '''
        check_type('pairs', npairs, int, type(None))
        check_type('virtuals', nvirt, int, type(None))
        self._lf = lf
        self._nocc = occ_model.noccs[0]
        self._nbasis = lf.default_nbasis
        if npairs is None:
            npairs = occ_model.noccs[0]
        elif npairs >= lf.default_nbasis:
            raise ValueError('Number of electron pairs (%i) larger than number of basis functions (%i)' %(npairs, self.nbasis))
        if nvirt is None:
            nvirt = (lf.default_nbasis-npairs)
        elif nvirt >= lf.default_nbasis:
            raise ValueError('Number of virtuals (%i) larger than number of basis functions (%i)' %(nvirt, self.nbasis))
        self._npairs = npairs
        self._nvirt = nvirt
        self._cache = Cache()
        self._ecore = 0
        self._geminal = lf.create_two_index(npairs, nvirt)
        self._lagrange = lf.create_two_index(npairs, nvirt)

    def __call__(self, one, two, core, orb, olp, scf, **kwargs):
        '''Optimize geminal coefficients and---if required---find
           optimal set of orbitals.

           **Arguments:**

           one, two
                One- and two-body integrals (some Hamiltonian matrix elements)
                expressed in the AO basis.

           core
                The core energy (not included in 'one' and 'two').

           orb
                An expansion instance. It contains the AO/MO coefficients
                (orbitals).

           olp
                The AO overlap matrix. A TwoIndex instance.

           scf
                A boolean. If True: Initializes orbital optimization.

           **Keywords:**
                See :py:meth:`RAp1rog.solve`
                and :py:meth:`RAp1rog.solve_scf`
        '''
        if scf:
            return self.solve_scf(one, two, core, orb, olp, **kwargs)
        else:
            return self.solve(one, two, core, orb, olp, **kwargs)

    def solve(self, one, two, core, orb, olp, **kwargs):
        raise NotImplementedError

    def solve_scf(self, one, two, core, orb, olp, **kwargs):
        raise NotImplementedError

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._nbasis

    nbasis = property(_get_nbasis)

    def _get_nocc(self):
        '''The number of occupied orbitals'''
        return self._nocc

    nocc = property(_get_nocc)

    def _get_nvirt(self):
        '''The number of virtual orbitals'''
        return self._nvirt

    nvirt = property(_get_nvirt)

    def _get_npairs(self):
        '''The number of electron pairs'''
        return self._npairs

    npairs = property(_get_npairs)

    def _get_lf(self):
        '''The LinalgFactory instance'''
        return self._lf

    lf = property(_get_lf)

    def _get_ecore(self):
        '''The core energy'''
        return self._ecore

    ecore = property(_get_ecore)

    def _get_dimension(self):
        '''The number of unknowns (i.e. the number of geminal coefficients)'''
        return self._npairs*self._nvirt

    dimension = property(_get_dimension)

    def _get_geminal(self):
        '''The geminal coefficients'''
        return self._geminal

    geminal = property(_get_geminal)

    def _get_lagrange(self):
        '''The Lagrange multipliers'''
        return self._lagrange

    lagrange = property(_get_lagrange)

    def __clear__(self):
        self.clear()

    def clear(self):
        '''Clear all wavefunction information'''
        self._cache.clear()

    def clear_dm(self):
        '''Clear RDM information'''
        self._cache.clear(tags='d', dealloc=True)

    def clear_geminal(self):
        '''Clear geminal information'''
        self._geminal.clear()

    def clear_lagrange(self):
        '''Clear lagrange information'''
        self._lagrange.clear()

    def update_ecore(self, new):
        '''Update core energy'''
        self._ecore = new

    def update_geminal(self, geminal=None):
        '''Update geminal matrix

           **Optional arguments:**

           geminal
                When provided, this geminal matrix is stored.
        '''
        if geminal is None:
            raise NotImplementedError
        else:
            self._geminal.assign(geminal)
        return geminal

    def update_lagrange(self, lagrange=None, dim1=None, dim2=None):
        '''Update Lagragne multipliers

           **Optional arguments:**

           lagrange
                When provided, this set of Lagrange multipliers is stored.
        '''
        if lagrange is None:
            raise NotImplementedError
        else:
            self.lagrange.assign(lagrange)
        return lagrange

    def update_auxmatrix(self, select, two_mo, one_mo=None):
        '''Update auxiliary matrices'''
        raise NotImplementedError

    def get_auxmatrix(self, select):
        '''Get auxiliary matrices'''
        raise NotImplementedError

    def init_one_dm(self, select):
        '''Initialize 1-RDM as OneIndex object

           The 1-RDM expressed in the natural orbital basis is diagonal and
           only the diagonal elements are stored.

           **Arguments**

           select
                'ps2' or 'response'.
        '''
        check_options('onedm', select, 'ps2', 'response')
        dm, new = self._cache.load('one_dm_%s' % select, alloc=(self._lf.create_one_index, self.nbasis), tags='d')
        if not new:
            raise RuntimeError('The density matrix one_dm_%s already exists. Call one_dm_%s.clear prior to updating the 1DM.' % select)
        return dm

    def init_two_dm(self, select):
        '''Initialize 2-RDM as TwoIndex object

           Only the symmetry-unique elements of the (response) 2-RDM are
           stored. These are matrix elements of type
                Gamma_{p\bar{q}p\bar{q}} (spin-up and spin-down (bar-sign)),
           or
                Gamma_{p\bar{p}q\bar{q}}
           and are stored as elements {pq} of two_dm_pqpq, and two_dm_ppqq.

           **Arguments**

           select
                '(r(esponse))ppqq', or '(r(esponse))pqpq'.
        '''
        check_options('twodm', select, 'ppqq', 'pqpq', 'rppqq', 'rpqpq')
        dm, new = self._cache.load('two_dm_%s' % select, alloc=(self._lf.create_two_index, self.nbasis), tags='d')
        if not new:
            raise RuntimeError('The density matrix two_dm_%s already exists. Call two_dm_%s.clear prior to updating the 2DM.' % select)
        return dm

    def init_three_dm(self, select):
        '''Initialize 3-RDM

           **Arguments**

           select
        '''
        raise NotImplementedError

    def init_four_dm(self, select):
        '''Initialize 4-RDM

           **Arguments**

           select
        '''
        raise NotImplementedError

    def get_one_dm(self, select):
        '''Get a density matrix (1-RDM). If not available, it will be created
           (if possible)

           **Arguments:**

           select
                'ps2', or 'response'.
        '''
        if not 'one_dm_%s' % select in self._cache:
            self.update_one_dm(select)
        return self._cache.load('one_dm_%s' % select)

    def get_two_dm(self, select):
        '''Get a density matrix (2-RDM). If not available, it will be created
           (if possible)

           **Arguments:**

           select
                '(r(esponse))ppqq', or '(r(esponse))pqpq'.
        '''
        if not 'two_dm_%s' % select in self._cache:
            self.update_two_dm(select)
        return self._cache.load('two_dm_%s' % select)

    def get_three_dm(self, select):
        '''Get a density matrix (3-RDM). If not available, it will be created
           (if possible)

           **Arguments:**

           select
        '''
        raise NotImplementedError

    def get_four_dm(self, select):
        '''Get a density matrix (4-RDM). If not available, it will be created
           (if possible)

           **Arguments:**

           select
        '''
        raise NotImplementedError

    one_dm_ps2 = PropertyHelper(get_one_dm, 'ps2', 'Alpha 1-RDM')
    one_dm_response = PropertyHelper(get_one_dm, 'response', 'Alpha 1-RDM')
    two_dm_ppqq = PropertyHelper(get_two_dm, 'ppqq', 'Alpha-beta (ppqq) 2-RDM')
    two_dm_pqpq = PropertyHelper(get_two_dm, 'pqpq', 'Alpha-beta (pqpq) 2-RDM')
    two_dm_rppqq = PropertyHelper(get_two_dm, 'rppqq', 'Alpha-beta (ppqq) 2-RDM')
    two_dm_rpqpq = PropertyHelper(get_two_dm, 'rpqpq', 'Alpha-beta (pqpq) 2-RDM')

    def update_one_dm(self, one_dm=None):
        '''Update 1-RDM

           **Optional arguments:**

           one_dm
                When provided, this 1-RDM is stored.
        '''
        raise NotImplementedError

    def update_two_dm(self, two_dm=None):
        '''Update 2-RDM

           **Optional arguments:**

           two_dm
                When provided, this 2-RDM is stored.
        '''
        raise NotImplementedError

    def update_three_dm(self, three_dm=None):
        '''Update 3-RDM

           **Optional arguments:**

           three_dm
                When provided, this 3-RDM is stored.
        '''
        raise NotImplementedError

    def update_four_dm(self, four_dm=None):
        '''Update 2-RDM

           **Optional arguments:**

           four_dm
                When provided, this 4-RDM is stored.
        '''
        raise NotImplementedError

    # Initial guess generators:
    def generate_guess(self, guess, dim=None):
        '''Generate a guess of type 'guess'.

           **Arguments:**

           guess
               A dictionary, containing the type of guess.

           **Optional arguments:**

           dim
               Length of guess.
        '''
        if dim is None:
            dim = self.dimension
        if guess['type'] == 'random':
            return np.random.random(dim)*guess['factor']
        elif guess['type'] == 'const':
            return np.ones(dim)*guess['factor']
        else:
            raise ValueError('Guess not supported.')

    def compute_rotation_matrix(self, coeff):
        '''Compute orbital rotation matrix'''
        raise NotImplementedError

    # Check convergence:
    def check_convergence(self, e0, e1, gradient, thresh):
        '''Check convergence.

           **Arguements:**

           e0, e1
                Used to calculate energy difference e0-e1

           gradient
                The gradient, a OneIndex instance

           thresh
                Dictionary containing threshold parameters ('energy', 'gradientmax',
                'gradientnorm')

           **Returns:**
                True if energy difference, norm of orbital gradient, largest
                element of orbital gradient are smaller than some threshold
                values.
        '''
        return math.fabs(e0-e1) < thresh['energy'] and \
               gradient.get_max() < thresh['gradientmax'] and \
               gradient.norm() < thresh['gradientnorm']

    def check_stepsearch(self, linesearch):
        '''Check trustradius. Abort calculation if trustradius is smaller than
           1e-8
        '''
        return linesearch.method == 'trust-region' and \
               linesearch.trustradius < 1e-8

    def prod(self, lst):
        return reduce(mul, lst)

    def perm(self, a):
        '''Calculate the permament of a matrix

           **Arguements**

           a
                A np array
        '''
        check_type('matrix', a, np.ndarray)
        n = len(a)
        r = range(n)
        s = permutations(r)
        return math.fsum(self.prod(a[i][sigma[i]] for i in r) for sigma in s)
