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

   Abbreviations used in this module:

   * wfn = wavefunction
'''

import numpy as np
import math as math

from horton.cache import Cache
from horton.matrix import DenseTwoIndex, Expansion
from horton.log import timer


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
        # For doc strings:
        if obj is None:
            return self
        # For actual use
        try:
            return self.method(obj, self.arg)
        except KeyError, e:
            raise AttributeError('The requested attribute (%s) is not available.' % e.args[0])


class Geminal(object):
    def __init__(self, lf, occ_model, nvirt=None):
        """
           **Arguments:**

           lf
                A LinalgFactory instance.

           occ_model
                Occupation model

           **Optional arguments:**

        """
        self._lf = lf
        self._nocc = occ_model.noccs[0]
        self._nbasis = lf.default_nbasis
        if nvirt is None:
            nvirt = (lf.default_nbasis-occ_model.noccs[0])
        self._nvirt = nvirt
        self._npairs = occ_model.noccs[0]
        self._cache = Cache()
        self._ecore = 0
        self._geminal = lf.create_two_index(self._npairs, nvirt)
        self._lagrange = lf.create_two_index(self._npairs, nvirt)

    def __call__(self, one, two, core, exps, orb, **kwargs):
        raise NotImplementedError

    def _get_nbasis(self):
        return self._nbasis

    nbasis = property(_get_nbasis)

    def _get_nocc(self):
        return self._nocc

    nocc = property(_get_nocc)

    def _get_nvirt(self):
        return self._nvirt

    nvirt = property(_get_nvirt)

    def _get_npairs(self):
        return self._npairs

    npairs = property(_get_npairs)

    def _get_nel(self):
        return 2*self._npairs

    nel = property(_get_nel)

    def _get_lf(self):
        return self._lf

    lf = property(_get_lf)

    def _get_ecore(self):
        return self._ecore

    ecore = property(_get_ecore)

    def _get_dimension(self):
        return self._npairs*self._nvirt

    dimension = property(_get_dimension)

    def _get_geminal(self):
        return self._geminal

    geminal = property(_get_geminal)

    def _get_lagrange(self):
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
        self.geminal.clear()

    def clear_lagrange(self):
        self.lagrange.clear()

    def update_ecore(self, new):
        '''Update core energy
        '''
        self._ecore = new

    def update_geminal(self, geminal=None, dim1=None, dim2=None):
        '''Update geminal matrix

           **Optional arguments:**

           geminal
                When provided, this geminal matrix is stored.
        '''
        if geminal is None:
            raise NotImplementedError
        else:
            if isinstance(geminal, DenseTwoIndex):
                self._geminal.assign(geminal)
            else:
                self._geminal.assign_array(geminal, self.nocc, self.nvirt)
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
            if isinstance(lagrange, DenseTwoIndex):
                self._lagrange.assign(lagrange)
            else:
                self._lagrange.assign_array(lagrange, dim1, dim2)
        return lagrange

    def update_matrix(self, select, two_mo, one_mo=None):
        raise NotImplementedError

    def get_matrix(self, select):
        raise NotImplementedError

    def init_one_dm(self, select):
        '''Initialize 1-RDM as OneIndex object

           The 1-RDM expressed in the natural orbital basis is diagonal and
           only the diagonal elements are stored.

           **Arguments**

           select
                'ps2' or 'response'.
        '''
        if select not in ['ps2', 'response']:
            raise ValueError('The select argument must be one of ps2 or response.')
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
                '(r)ppqq', or '(r)pqpq'.
        '''
        if select not in ['ppqq', 'pqpq', 'rppqq', 'rpqpq']:
            raise ValueError('The select argument must be one of ppqq, pqpq, rppqq, or rpqpq.')
        dm, new = self._cache.load('two_dm_%s' % select, alloc=(self._lf.create_two_index, self.nbasis), tags='d')
        if not new:
            raise RuntimeError('The density matrix two_dm_%s already exists. Call two_dm_%s.clear prior to updating the 2DM.' % select)
        return dm

    def init_three_dm(self, select):
        '''Initialize 3-RDM as ThreeIndex object

           **Arguments**

           select
                'uuu', 'uud', 'uudoff', 'udd', 'uddoff'.
        '''
        if select not in ['uuu', 'uudoff', 'uud', 'udd', 'uddoff']:
            raise ValueError('The select argument must be one of uuu, uud, udd, uudoff, uddoff.')
        dm, new = self._cache.load('three_dm_%s' % select, alloc=(self._lf.create_three_index, self.nbasis), tags='d')
        if not new:
            raise RuntimeError('The density matrix three_dm_%s already exists. Call three_dm_%s.clear prior to updating the 3DM.' % select)
        return dm

    def init_four_dm(self, select):
        '''Initialize 4-RDM as SomeIndex object. Currently, only one block is
           supported (stored as TwoIndex).

           **Arguments**

           select
                'udud'.
        '''
        if select not in ['udud']:
            raise ValueError('The select argument must be one of udud.')
        dm, new = self._cache.load('four_dm_%s' % select, alloc=(self._lf.create_two_index, self.nbasis), tags='d')
        if not new:
            raise RuntimeError('The density matrix four_dm_%s already exists. Call four_dm_%s.clear prior to updating the 4DM.' % select)
        return dm

    def get_one_dm(self, select):
        '''Get a density matrix (1-RDM). If not available, it will be created (if possible)

           **Arguments:**

           select
                'ps2', or 'response'.
        '''
        if not 'one_dm_%s' % select in self._cache:
            self.update_one_dm(select)
        return self._cache.load('one_dm_%s' % select)

    def get_two_dm(self, select):
        '''Get a density matrix (2-RDM). If not available, it will be created (if possible)

           **Arguments:**

           select
                '(r)ppqq', or '(r)pqpq'.
        '''
        if not 'two_dm_%s' % select in self._cache:
            self.update_two_dm(select)
        return self._cache.load('two_dm_%s' % select)

    def get_three_dm(self, select):
        '''Get a density matrix (3-RDM). If not available, it will be created (if possible)

           **Arguments:**

           select
                'uuu', 'uudoff', 'uud', 'udd', 'uddoff'.
        '''
        if not 'three_dm_%s' % select in self._cache:
            self.update_three_dm(select)
        return self._cache.load('three_dm_%s' % select)

    def get_four_dm(self, select):
        '''Get a density matrix (3-RDM). If not available, it will be created (if possible)

           **Arguments:**

           select
                'udud'.
        '''
        if not 'four_dm_%s' % select in self._cache:
            self.update_four_dm(select)
        return self._cache.load('four_dm_%s' % select)

    one_dm_ps2 = PropertyHelper(get_one_dm, 'ps2', 'Alpha 1-RDM')
    one_dm_response = PropertyHelper(get_one_dm, 'response', 'Alpha 1-RDM')
    two_dm_ppqq = PropertyHelper(get_two_dm, 'ppqq', 'Alpha-beta (ppqq) 2-RDM')
    two_dm_pqpq = PropertyHelper(get_two_dm, 'pqpq', 'Alpha-beta (pqpq) 2-RDM')
    two_dm_rppqq = PropertyHelper(get_two_dm, 'rppqq', 'Alpha-beta (ppqq) 2-RDM')
    two_dm_rpqpq = PropertyHelper(get_two_dm, 'rpqpq', 'Alpha-beta (pqpq) 2-RDM')
    three_dm_uuu = PropertyHelper(get_three_dm, 'uuu', 'Alpha-alpha-alpha 3-RDM')
    three_dm_uud = PropertyHelper(get_three_dm,'uud', 'Alpha-alpha-beta 3-RDM')
    three_dm_uudoff = PropertyHelper(get_three_dm,'uudoff', 'Alpha-alpha-beta 3-RDM (off-diagonal)')
    three_dm_udd = PropertyHelper(get_three_dm, 'udd', 'Alpha-beta-beta 3-RDM')
    three_dm_uddoff = PropertyHelper(get_three_dm, 'uddoff', 'Alpha-beta-beta 3-RDM (off-diagonal)')
    four_dm_udud = PropertyHelper(get_four_dm, 'udud', 'Alpha-beta-alpha-beta 4-RDM')

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

    @timer.with_section('4-Index Trans')
    def apply_4index_trans(self, two_ao, orb, orb2, orb3,
                           orb4, indextrans='tensordot'):
        '''Do four-index transformation. A subroutine that calls one function
           of the matrix class.

           **Arguments:**

           two_ao
               Two-electron integrals in the AO basis.

           orb
               AO/MO coefficient matrix.

           **Optional arguments:**

           indextrans
               Choice of 4-index transformation. Default 'tensordot'.
        '''
        two_mo = two_ao.copy()
        if indextrans == 'einstein':
            two_mo.apply_four_index_transform_einsum(two_ao, orb, orb2, orb3, orb4)
        elif indextrans == 'tensordot':
            # this is the prefered one! It is fast, but uses a lot of memory...
            two_mo.apply_four_index_transform_tensordot(two_ao, orb, orb2, orb3, orb4)
        return two_mo

    def update_mo_integrals(self, one, two, indextrans='tensordot', *exps):
        '''Update MO integrals. Returns list of transformed 1- and 2-electron
           integrals according to a list of expansion coefficients.

           **Arguments:**

           one
               One-electron integrals in the AO basis.

           two
               Two-electron integrals in the AO basis.

           **Optional arguments:**

           indextrans
               Choice of 4-index transformation. Default 'tensordot'.

           args
               The expansion coefficients.

        '''
        exp = []
        two_mo = []
        one_mo = []
        if not all([isinstance(i, Expansion) for i in exps]):
            raise TypeError('argument of unsupported type: %s' %i)
        for arg in exps:
            exp.append(arg)

        for i in xrange(len(exps)):
            for j in xrange(i, len(exps)):
                #Transform AO integrals into MO basis
                two_mo.append(self.apply_4index_trans(two, exp[i], exp[j], exp[i], exp[j], indextrans))

        # Set up one-electron part of the Hamiltonian and transform it to MO basis
            tmp = self.lf.create_two_index()
            tmp.apply_2index_trans(one, exps[i])
            one_mo.append(tmp)
        return one_mo, two_mo

    # Initial guess generators:
    def generate_guess(self, guess, dim=None, nguess=1):
        '''Generates a guess of type 'guess'.

           **Arguments:**

           guess
               Type of guess.

           **Optional arguments:**

           dim
               Length of guess.

           nguess
               Number of guesses generated.
        '''
        guessv = []
        if dim is None:
            dim = self.dimension
        if guess['guess']=='random':
            for i in range(nguess):
                guessv.append(np.random.random(dim)*guess['factor'])
        elif guess['guess']=='const':
            for i in range(nguess):
                guessv.append(np.ones(dim)*guess['factor'])
        else:
            raise ValueError('Guess not supported.')
        return guessv

    def get_rotation_matrix(self, coeff):
        raise NotImplementedError

    # Rotate MO coefficient vector:
    def rotate_orbitals(self, *args):
        '''Apply rotation to MO coefficient matrix.
        '''
        exps = []
        rotmat = []
        for arg in args:
            if isinstance(arg, np.ndarray):
                rotmat.append(arg)
            elif isinstance(arg, Expansion):
                exps.append(arg)
            else:
                raise TypeError('argument of unsupported type: %s' % arg)
        for i in xrange(len(exps)):
            output = self.lf.create_expansion()
            output.dotarray(exps[i], rotmat[i])
            exps[i].assign(output)

    # Check convergence:
    def check_convergence(self, energy, energy_old, gradient, ethresh=1e-8,
                          gradthresh=1e-5, gradnormthresh=1e-8):
        '''Check convergence.
           Return True, otherwise False in terms of energy difference,
           norm of orbital gradient, largest element of orbital
           gradient.
        '''
        cenergy =  False
        cgradient = False
        cgradientnorm = False
        if (math.fabs(energy-energy_old) < ethresh):
            cenergy = True
        if (np.max(np.absolute(gradient)) < gradthresh):
            cgradient = True
        if (np.dot(gradient,gradient) < gradnormthresh):
            cgradientnorm = True
        if cenergy and cgradient and cgradientnorm:
            return True
        else:
            return False

    def check_stepsearch(self, linesearch):
        '''Check convergence.
        '''
        if linesearch.method in ['trust-region']:
            if (linesearch.trustradius < 1e-8):
                return True
        else:
            return False
