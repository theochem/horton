# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2015 The Horton Development Team
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
'''Perturbation theory module

   Variables used in this module:
    :nocc:       number of occupied orbitals in the principle configuration
    :nvirt:      number of virtual orbitals in the principle configuration
    :nbasis:     total number of basis functions
    :energy:     the energy correction, list that can contain different
                 contributions
    :amplitudes: the optimized amplitudes, list that can contain different
                 contributions

    Indexing convention:
     :i,j,k,..: occupied orbitals of principle configuration
     :a,b,c,..: virtual orbitals of principle configuration
     :p,q,r,..: general indices (occupied, virtual)
'''

import numpy as np
import math
from scipy import optimize as opt
from horton.log import log, timer
from horton.cache import Cache
from horton.utils import check_type, check_options
from horton.orbital_utils import transform_integrals
from horton.matrix.base import LinalgFactory, LinalgObject, OneIndex, \
    Expansion, TwoIndex, ThreeIndex, FourIndex


__all__ = [
    'Perturbation',
    'RMP2',
    'PTa',
    'PTb',
]




class Perturbation(object):
    '''Perturbation class

       Purpose:
       Optimize amplitudes and determine energy correction to some
       reference wavefunction.

       Currently supported wavefunction models:
        * RHF
        * RAP1roG

       Currently supported Perturbation Theory models:
        * MP2 (if Psi_0 = RHF)
        * PTa/PTb (if Psi_0 = RAP1roG)
    '''

    def __init__(self, lf, occ_model):
        '''
           **Arguments:**

           lf
                A LinalgFactory instance.

           occ_model
                Aufbau model

           **Optional arguments:**

        '''
        self._lf = lf
        self._nocc = occ_model.noccs[0]
        self._nbasis = lf.default_nbasis
        self._nvirt = (lf.default_nbasis-occ_model.noccs[0])
        self._cache = Cache()
        self._energy = []
        self._amplitudes = []

    def _get_lf(self):
        '''The linalg factory'''
        return self._lf

    lf = property(_get_lf)

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

    def _get_energy(self):
        '''The PT energy'''
        return self._energy

    energy = property(_get_energy)

    def _get_amplitudes(self):
        '''The PT amplitudes'''
        return self._amplitudes

    amplitudes = property(_get_amplitudes)

    def update_energy(self, args):
        '''Update PT energy'''
        for arg in args:
            self.energy.append(arg)

    def update_amplitudes(self, new):
        raise NotImplementedError

    def vfunction(self, coeffs, matrix):
        raise NotImplementedError

    def calculate_energy(self, amplitudes, args):
        raise NotImplementedError

    def get_guess(self):
        raise NotImplementedError

    def clear(self):
        '''Clear all wavefunction information'''
        self._cache.clear()

    def clear_aux_matrix(self):
        '''Clear the auxiliary matrices'''
        self._cache.clear(tags='m', dealloc=True)

    def get_aux_matrix(self, select):
        '''Get an auxiliary matrix.

           **Arguments:**

           select
                Suffix of auxiliary matrix. See :py:meth:`RMP2.init_aux_matrix`,
                :py:meth:`PTa.init_aux_matrix`, and
                :py:meth:`PTb.init_aux_matrix` for possible choices
        '''
        if not '%s' % select in self._cache:
            raise ValueError("The auxmatrix %s not found in cache. Did you use init_aux_matrix?" %select)
        return self._cache.load('%s' % select)

    @timer.with_section('Perturbation')
    def __call__(self, one, two, *args, **kwargs):
        '''Performs a perturbation theory calculation.

           **Arguments:**

           one, two
               One- (TwoIndex) and two-body (FourIndex) integrals (some
               Hamiltonian matrix elements)

           args
               If Psi_0 = RHF, first argument is the MO coefficient matrix
               (Expansion instance), if Psi_0 = AP1roG, first argument is again
               the MO coefficient matrix the second argument is the geminal
               coefficient matrix (TwoIndex).

           **Keywords:**
               Contains reference energy and solver specific input parameters:
                * eref: (float) reference energy (default float('nan'))
                * ecore: (float) core energy (default float('nan'))
                * threshold: (float) tolerance for amplitudes (default 1e-6)
                * maxiter: (int) maximum number of iterations (default 200)
                * guess: (np.array) initial guess (default None)
                * indextrans: (str) 4-index transformation. One of ``einsum``,
                              ``tensordot`` (default ``tensordot``)

           **Returns**
                List of energy contributions (total energy, seniority-0,
                seniority-2, and seniority-4) and PT amplitudes (doubles)

        '''
        log.hline('=')
        log(' ')
        log('Entering perturbation theory module')
        log(' ')
        log.hline('~')

        #
        # Print method specific information
        #
        self.print_info(**kwargs)

        #
        # Check input parameters
        #
        self.check_input(**kwargs)

        #
        # Append arguments, used as arguments in root finding:
        #
        fargs = []
        for arg in args:
            fargs.append(arg)

        #
        # Transform integrals:
        #
        indextrans = kwargs.get('indextrans', 'tensordot')
        mo1, mo2 = transform_integrals(one, two, indextrans, fargs[0])
        for int1 in mo1:
            fargs.append(int1)
        for int2 in mo2:
            fargs.append(int2)

        #
        # Construct auxiliary matrices (checks also type of arguments):
        #
        matrix = self.calculate_aux_matrix(*fargs)

        #
        # Append arguments, used as arguments in root finding:
        #
        for mat in matrix:
            fargs.append(mat)

        #
        # Solve for energy and amplitudes:
        #
        energy, amplitudes = self.solve(*fargs, **kwargs)
        self.update_energy(energy)
        self.update_amplitudes(amplitudes)

        #
        # Print some output information for user:
        #
        self.print_energy(**kwargs)

        #
        # Make sanity checks. If somethings is wrong, abort calculation:
        #
        self.check_result(**kwargs)

        if log.do_medium:
            log(' ')
            log.hline('=')

        return self.energy, self.amplitudes




class RMP2(Perturbation):
    '''Moller-Plesset Perturbation Theory of second order

       Purpose:
       Optimize amplitudes and determine energy correction to Hartree-Fock
       reference wavefunction.
    '''

    @timer.with_section('MP2Solver')
    def solve(self, *args, **kwargs):
        '''Solve for energy and amplitudes

           **Arguments:**

           args
                Contains one- and two-electron integrals in the MO basis:
                    * [0]:  wfn expansion coefficients
                    * [1]:  1-el MO integrals
                    * [2]:  2-el MO integrals

           **Keywords**
                :eref: reference energy
                :threshold: threshold for symmetry check of MP2 amplitudes
        '''
        check_type('args[0]', args[0], Expansion)
        check_type('args[1]', args[1], TwoIndex)
        check_type('args[2]', args[2], FourIndex)
        #
        # Calculate excitation matrix <jb|kc>
        #
        ex_matrix = self.calculate_ex_matrix(args[2])
        #
        # Calculate MP2 amplitudes t_jbkc
        #
        amplitudes = self.calculate_amplitudes(ex_matrix)
        #
        # Calculate MP2 energy correction
        #
        energy = amplitudes.contract_four('abcd,abcd',ex_matrix, 2.0)
        energy-= amplitudes.contract_four('abcd,adcb',ex_matrix)

        return [energy], amplitudes

    @timer.with_section('ExMatrix')
    def calculate_ex_matrix(self, mo2):
        '''Calculates excitation matrix with elements (jbkc)

           **Arguments:**

           mo2
               2-el MO integrals
        '''
        # B_jk|bc
        out = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        mo2.slice_to_four('abcd->acbd', out, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        return out

    @timer.with_section('MP2Amplitudes')
    def calculate_amplitudes(self, matrix):
        '''Calculates MP2 amplitudes

           **Arguments:**

           matrix
               Sliced 2-el MO integrals
        '''
        fock = self.get_aux_matrix('fock')
        amplitudes = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        for i in range(self.nocc):
            for j in range(self.nocc):
                for a in range(self.nvirt):
                    for b in range(self.nvirt):
                        val = 1.0/(fock.get_element(i)+fock.get_element(j) \
                                  -fock.get_element(self.nocc+a) \
                                  -fock.get_element(self.nocc+b))
                        amplitudes.set_element(i,a,j,b,val, symmetry=1)
        amplitudes.imul(matrix)
        return amplitudes

    def calculate_aux_matrix(self, *args):
        '''Compute auxiliary matrices

           **Arguments:**

           args
                One- and two-electron integrals (some Hamiltonian matrix
                elements) in the MO basis.
        '''
        check_type('args[1]', args[1], TwoIndex)
        check_type('args[2]', args[2], FourIndex)
        self.clear_aux_matrix()
        return self.update_aux_matrix(args[1], args[2])

    def init_aux_matrix(self, select):
        '''Initialize auxiliary matrices

           **Arguments:**

           select
                One of 'fock'
        '''
        check_options('select', select, 'fock')
        if select=='fock':
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_one_index, self.lf.default_nbasis), tags='m')
        if not new:
            raise RuntimeError('The matrix %s already exists. Call clear prior \
                                to updating the wfn.' % select)
        return matrix

    @timer.with_section('UpdatingAuxMat')
    def update_aux_matrix(self, mo1, mo2):
        '''Derive all auxiliary matrices.
           fock_pp:     one_pp + sum_m(2<pm|pm> - <pm|mp>),

           **Arguments:**

           mo1, mo2
                one- and two-electron integrals.
        '''
        auxmat1 = self.init_aux_matrix('fock')
        #
        # 2-el part of inactive Fock (diagonal)
        #
        tmp = self.lf.create_two_index()
        mo2.slice_to_two('abab->ab', tmp, 2.0, True)
        mo2.slice_to_two('abba->ab', tmp,-1.0, False)
        tmp.contract_to_one('ab->a', auxmat1, 1.0, True, 0, self.nbasis, 0, self.nocc)

        #
        # 1-el part of inactive Fock (diagonal)
        #
        tmp1 = self.lf.create_one_index()
        mo1.copy_diagonal(tmp1)
        auxmat1.iadd(tmp1)
        del tmp1, tmp

        return [auxmat1]

    def print_energy(self, **kwargs):
        if log.do_medium:
            log('E_ref:     %16.8f a.u.' %(kwargs.get('eref')))
            log('E_MP2:     %16.8f a.u.' %self.energy[0])
            log.hline('-')
            log('E_tot:     %16.8f a.u.' %(self.energy[0]+kwargs.get('eref')))

    def update_amplitudes(self, new):
        '''Update MP2 amplitudes

           **Arguments:**

           new
                PT amplitudes. A FourIndex instance
        '''
        if self.amplitudes:
            raise ValueError('Warning: List of MP2 amplitudes not empty!')
        else:
            self.amplitudes.append(new)

    def print_info(self, **kwargs):
        if log.do_medium:
            log('MP2 perturbation module')
            log(' ')
            log('OPTIMIZATION PARAMETERS:')
            log('Reference Function           %s' %('RHF'))
            log('Number of pairs:             %i' %self.nocc)
            log('Number of virtuals:          %i' %self.nvirt)
            log.hline()

    def check_input(self, **kwargs):
        '''Check input parameters'''
        for name, value in kwargs.items():
            check_options('name', name, 'eref', 'threshold')
        eref = kwargs.get('eref', float('nan'))
        if math.isnan(eref):
            raise ValueError('Warning: Cannot find reference energy in MP2 module!')

    def check_result(self, **kwargs):
        '''Check if amplitudes are symmetric (within a given threshold).'''
        thresh = kwargs.get('threshold', 1e-6)

        if not self.amplitudes[0].is_symmetric('cdab', thresh):
            raise ValueError('Warning: Cluster amplitudes not symmetric!')




class PTa(Perturbation):

    @timer.with_section('PTaSolver')
    def solve(self, *args, **kwargs):
        '''Solve for energy and amplitudes

           **Arguments:**

           args
                Contains geminal coefficients, 1- and 2-el integrals, and
                auxiliary matrices

           **Keywords**
                :eref:      reference energy
                :ecore:     core energy
                :threshold: threshold when checking symmetry of amplitudes
        '''
        #
        # Get "excitation" matrices w.r.t. |AP1roG> (psi0)
        #
        exjbkc = self.vfunction_psi0(*args, **kwargs)
        #
        # Calculate amplitudes
        #
        amplitudes = self.calculate_amplitudes(exjbkc)
        #
        # Calculate energy contributions of different seniority sectors
        # and total energy
        #
        energy = self.calculate_energy(amplitudes, *args)

        return energy, amplitudes

    @timer.with_section('VecFctPTa')
    def vfunction_psi0(self, *args, **kwargs):
        '''Elements of <bcjk|H|AP1roG>.

           **Arguments:**

           args
                All function arguments needed to calculate the vector
                function:

                * [0]:  wfn expansion coefficients
                * [1]:  geminal coefficients
                * [2]:  1-el MO integrals
                * [3]:  2-el MO integrals
                * [4]:  inactive Fock matrix
                * [5]:  ocjbc auxilary matrix
                * [6]:  vcjkb auxilary matrix
                * [7]:  vcjbc auxilary matrix
                * [8]:  ocjkb auxilary matrix
                * [9]:  dcjb auxilary matrix
        '''
        check_type('args[0]', args[0], Expansion)
        check_type('args[1]', args[1], TwoIndex)
        check_type('args[2]', args[2], TwoIndex)
        check_type('args[3]', args[3], FourIndex)
        check_type('args[4]', args[4], TwoIndex)
        check_type('args[5]', args[5], ThreeIndex)
        check_type('args[6]', args[6], ThreeIndex)
        check_type('args[7]', args[7], ThreeIndex)
        check_type('args[8]', args[8], ThreeIndex)
        check_type('args[9]', args[9], TwoIndex)

        e0 = kwargs.get('eref')-kwargs.get('ecore')

        #
        # output
        #
        out = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        #
        # temporary storage
        #
        tmp4ind = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        tmp3ind0 = self.lf.create_three_index(self.nocc, self.nvirt, self.nvirt)
        tmp3ind1 = self.lf.create_three_index(self.nocc, self.nocc, self.nvirt)
        onesum = args[2].trace(0, self.nocc, 0, self.nocc)
        fockdiagsum = args[4].trace(0, self.nocc, 0, self.nocc)

        #
        # <jk|bc>
        #
        args[3].slice_to_four('abcd->acbd', out, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # c_kc <jc||bk>
        #
        args[3].contract_two_to_four('abcd,bd->cabd', args[1], out, 1.0, False, self.nocc, self.nbasis, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis)
        args[3].contract_two_to_four('abcd,bc->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jb <jc||bk>
        #
        args[3].contract_two_to_four('abcd,ca->cabd', args[1], out, 1.0, False, self.nocc, self.nbasis, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis)
        args[3].contract_two_to_four('abcd,da->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jc <bk|cj>
        #
        args[3].contract_two_to_four('abcd,dc->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_kb <bk|cj>
        #
        args[3].contract_two_to_four('abcd,ba->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jkbc <bk|jc>
        #
        args[3].contract_two_to_four('abcd,ac->acbd', args[1], tmp4ind, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp4ind.contract_two_to_four('abcd,cd->abcd', args[1], out, 1.0, False)
        args[3].contract_two_to_four('abcd,bc->acbd', args[1], tmp4ind, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp4ind.contract_two_to_four('abcd,ad->abcd', args[1], out, 1.0, False)
        #
        # delta_jk [ c_jc F_bc ]
        #
        tmp3ind0.iadd_expand_two_two('ac,bc->abc', args[1], args[4], 1.0, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # delta_jk [ c_jb F_bc ]
        #
        tmp3ind0.iadd_expand_two_two('ab,bc->abc', args[1], args[4], 1.0, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # delta_jk [ oc_jbc ]
        #
        tmp3ind0.iadd(args[5], -1.0)
        #
        # delta_jk [ vc_jbc ]
        #
        tmp3ind0.iadd(args[7], 1.0)
        #
        # Add delta_jk-contribution
        #
        out.iadd_expand_three_to_four('1-3-1-2', tmp3ind0, 1.0)
        #
        # delta_bc [ c_jb F_jk ]
        #
        tmp3ind1.iadd_expand_two_two('ab,ac->acb', args[1], args[4],-1.0, 0, self.nocc, 0, self.nocc)
        #
        # delta_bc [ c_kb F_jk ]
        #
        tmp3ind1.iadd_expand_two_two('cb,ac->acb', args[1], args[4],-1.0, 0, self.nocc, 0, self.nocc)
        #
        # delta_bc [ vc_jkb ]
        #
        tmp3ind1.iadd(args[6],-1.0)
        #
        # delta_bc [ oc_jkb ]
        #
        tmp3ind1.iadd(args[8], 1.0)
        #
        # Add delta_bc-contribution
        #
        out.iadd_expand_three_to_four('0-2-0-1', tmp3ind1, 1.0)
        #
        # delta_bc,jk [ (sum_m h_mm+F_mm)*c_jb ]
        #
        out.iadd_expand_two_to_four('diag', args[1],(-e0+onesum+fockdiagsum))
        #
        # delta_bc,jk [ (sum_md <mm|dd> c_jmbd) ]
        #
        out.iadd_expand_two_to_four('diag', args[9], 1.0)

        return out

    def vfunction_0(self, *args):
        '''Elements of <bcjk|H|0>.

           **Arguments:**

           args
                All function arguments needed to calculate the vector
                function:

                * [0]:  wfn expansion coefficients
                * [1]:  geminal coefficients
                * [2]:  1-el MO integrals
                * [3]:  2-el MO integrals
                * [4]:  inactive Fock matrix
                * [5]:  ocjbc auxilary matrix
                * [6]:  vcjkb auxilary matrix
                * [7]:  vcjbc auxilary matrix
                * [8]:  ocjkb auxilary matrix
        '''
        #
        # B_jb|kc
        #
        out = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        args[3].slice_to_four('abcd->cadb', out, 1.0, True, self.nocc, self.nbasis, self.nocc, self.nbasis, 0, self.nocc, 0, self.nocc)
        return out

    @timer.with_section('PTaAmplitudes')
    def calculate_amplitudes(self, matrix):
        '''Calculate amplitudes

           **Arguments:**

           matrix
                A FourIndex instance.
        '''
        fock = self.get_aux_matrix('fock')
        amplitudes = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        for i in range(self.nocc):
            for j in range(self.nocc):
                for a in range(self.nvirt):
                    aa = a+self.nocc
                    for b in range(self.nvirt):
                        bb = b+self.nocc
                        val = 1.0/(fock.get_element(i,i)+fock.get_element(j,j) \
                                  -fock.get_element(aa,aa) \
                                  -fock.get_element(bb,bb))
                        amplitudes.set_element(i,a,j,b,val, symmetry=1)
        amplitudes.imul(matrix)
        return amplitudes

    def calculate_aux_matrix(self, *args):
        '''Compute auxiliary matrices

           **Arguments:**

           args
                List of arguements. Only geminal coefficients [1], one- [2] and
                two-body [3] integrals are used.
        '''
        check_type('args[1]', args[1], TwoIndex)
        check_type('args[2]', args[2], TwoIndex)
        check_type('args[3]', args[3], FourIndex)
        self.clear_aux_matrix()
        return self.update_aux_matrix(args[2], args[3], args[1])

    def init_aux_matrix(self, select, dim=None, dim2=None):
        '''Initialize auxiliary matrices

           **Arguments:**

           select
                One of ``fock``, ``ocjbc``, ``vcjkb``, ``ocjkb``,
                ``vcjbc``, ``dcjb``

           **Optional arguments:**

           dim, dim2
                The dimension of the auxiliary matrices for specific axes.

        '''
        check_options('select', select, 'fock', 'ocjbc', 'vcjkb', 'ocjkb',
            'vcjbc', 'dcjb')
        if select=='fock':
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_two_index, self.lf.default_nbasis), tags='m')
        if select=='dcjb':
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_two_index, dim, dim2), tags='m')
        if select in ['ocjbc', 'vcjbc']:
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_three_index, dim, dim2, dim2), tags='m')
        if select in ['vcjkb', 'ocjkb']:
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_three_index, dim, dim, dim2), tags='m')
        if not new:
            raise RuntimeError('The matrix %s already exists. Call clear prior \
                                to updating the wfn.' % select)
        return matrix

    @timer.with_section('UpdatingAuxMat')
    def update_aux_matrix(self, mo1, mo2, cia):
        '''Derive all matrices.
           fock_pq:     one_pq + sum_m(2<pm|qm> - <pm|mq>),
           oc_jbc:      sum_m(<mm|bc> c_jm^bc),
           vc_jkb:      sum_d(<dd|jk> c_jk^bd),
           vc_jbc:      sum_d(<dd|bc> c_j^d),
           oc_jkb:      sum_m(<mm|jk> c_m^b),
           dc_jb:       sum_md(<mm|dd> c_jm^bd),

           **Arguments:**

           mo1, mo2
                one- and two-electron integrals to be sorted.

           cia
                The geminal coefficients. A TwoIndex instance
        '''
        #
        # Inactive Fock matrix
        #
        auxmat1 = self.init_aux_matrix('fock')
        mo2.contract_to_two('abcb->ac', auxmat1, 2.0, True, 0, self.nbasis, 0, self.nocc, 0, self.nbasis, 0, self.nocc)
        mo2.contract_to_two('abbc->ac', auxmat1,-1.0, False, 0, self.nbasis, 0, self.nocc, 0, self.nocc, 0, self.nbasis)
        auxmat1.iadd(mo1)
        #
        # oc_jbc = sum_m <mm|bc> c_jm^bc
        #
        tmp = self.lf.create_four_index(self.nocc, self.nocc, self.nvirt, self.nvirt)
        auxmat2 = self.init_aux_matrix('ocjbc', cia.nbasis, cia.nfn)
        mo2.contract_two_to_four('aabc,db->adbc', cia, tmp, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp.contract_two_to_three('abcd,ad->bcd', cia, auxmat2, 1.0, True)
        mo2.contract_two_to_four('aabc,dc->adbc', cia, tmp, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp.contract_two_to_three('abcd,ac->bcd', cia, auxmat2, 1.0, False)
        #
        # vc_jkb = sum_d <dd|jk> c_jk^bd
        #
        auxmat3 = self.init_aux_matrix('vcjkb', cia.nbasis, cia.nfn)
        mo2.contract_two_to_four('abcc,ad->abcd', cia, tmp, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp.contract_two_to_three('abcd,bc->abd', cia, auxmat3, 1.0, True)
        mo2.contract_two_to_four('abcc,bd->abcd', cia, tmp, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp.contract_two_to_three('abcd,ac->abd', cia, auxmat3, 1.0, False)
        #
        # vc_jbc = sum_d <bc|dd> c_j^d
        #
        auxmat4 = self.init_aux_matrix('vcjbc', cia.nbasis, cia.nfn)
        mo2.contract_two_to_three('abcc,dc->dab', cia, auxmat4, 1.0, True, self.nocc, self.nbasis, self.nocc, self.nbasis, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # oc_jkb = sum_m <mm|jk> c_m^b
        #
        auxmat5 = self.init_aux_matrix('ocjkb', cia.nbasis, cia.nfn)
        mo2.contract_two_to_three('aabc,ad->bcd', cia, auxmat5, 1.0, True, 0, self.nocc, 0, self.nocc, 0, self.nocc, 0, self.nocc)
        #
        # dc_jb = sum_md <mm|dd> c_jm^bd
        #
        auxmat6= self.init_aux_matrix('dcjb', cia.nbasis, cia.nfn)
        tmp = self.lf.create_two_index(self.nbasis, self.nbasis)
        tmp2 = self.lf.create_two_index(self.nocc, self.nocc)
        # There is a bug in np.einsum that forces us to slice first...
        mo2.slice_to_two('aabb->ab', tmp, 1.0, True)
        factor = cia.contract_two('ab,ab', tmp, 0, self.nocc, self.nocc, self.nbasis)
        auxmat6.iadd(cia, factor)
        cia.contract_two_to_two('ab,cb->ac', tmp, tmp2, 1.0, True, 0, self.nocc, self.nocc, self.nbasis)
        cia.contract_two_to_two('ab,ca->cb', tmp2, auxmat6, 1.0, False)

        return [auxmat1, auxmat2, auxmat3, auxmat4, auxmat5, auxmat6]

    def update_amplitudes(self, new):
        '''Update cluster amplitudes

           **Arguments:**

           new
                PT amplitudes. A FourIndex instance
        '''
        if self.amplitudes:
            raise ValueError('Warning: List of PT amplitudes not empty!')
        else:
            self.amplitudes.append(new)

    def calculate_energy(self, amplitudes, *args):
        '''Calculate PT energy and energy contribution of seniority sectors

           **Arguments:**

           amplitudes
                PT amplitudes. A FourIndex instance

           args
                List containing the geminal coefficients, 1- and 2-el
                integrals, and auxiliary matrices
        '''
        #
        # Get "excitation" matrices w.r.t. |HF> (0)
        #
        exjbkc0 = self.vfunction_0(*args)
        #
        # Seniority-0 sector
        #
        e_seniority_0 = amplitudes.contract_four('abab,abab', exjbkc0, 1.0)
        #
        # Seniority-2 sector
        #
        e_seniority_2 = amplitudes.contract_four('abad,abad', exjbkc0, 1.0)
        e_seniority_2+= amplitudes.contract_four('abdb,abdb', exjbkc0, 1.0)
        energy = amplitudes.contract_four('abcd,abcd', exjbkc0, 2.0)
        energy-= amplitudes.contract_four('abcd,adcb', exjbkc0, 1.0)
        #
        # Seniority-4 sector
        #
        e_seniority_4 = energy-e_seniority_2-e_seniority_0

        return [energy, e_seniority_0, e_seniority_2, e_seniority_4]

    def print_energy(self, **kwargs):
        if log.do_medium:
            log('E_PTa(Seniority 0):     %18.12f a.u.' %(self.energy[1]))
            log('E_PTa(Seniority 2):     %18.12f a.u.' %(self.energy[2]))
            log('E_PTa(Seniority 4):     %18.12f a.u.' %(self.energy[3]))
            log.hline('-')
            log('E_ref:                  %18.12f a.u.' %(kwargs.get('eref')))
            log('E_PTa:                  %18.12f a.u.' %self.energy[0])
            log.hline('-')
            log('E_tot:                  %18.12f a.u.' %(self.energy[0]+kwargs.get('eref')))

    def print_info(self, **kwargs):
        if log.do_medium:
            log('PTa perturbation module')
            log(' ')
            log('OPTIMIZATION PARAMETERS:')
            log('Reference Function           %s' %('AP1roG'))
            log('Number of pairs:             %i' %self.nocc)
            log('Number of virtuals:          %i' %self.nvirt)
            log.hline()

    def check_input(self, **kwargs):
        '''Check input parameters.'''
        for name, value in kwargs.items():
            check_options(name, name, 'ecore', 'eref', 'threshold', 'indextrans')
        eref = kwargs.get('eref', float('nan'))
        ecore = kwargs.get('ecore', float('nan'))
        if math.isnan(eref):
            raise ValueError('Warning: Cannot find reference energy in PTa module!')
        if math.isnan(ecore):
            raise ValueError('Warning: Cannot find core: energy in PTa module!')

    def check_result(self, **kwargs):
        '''Check if amplitudes are reasonable.'''
        thresh = kwargs.get('threshold', 1e-6)

        # check symmetry of amplitudes:
        if not self.amplitudes[0].is_symmetric('cdab', thresh):
            raise ValueError('Warning: Cluster amplitudes not symmetric. Aborting optimization!')

        # check if diagonal amplitudes are zero:
        tmp = self.amplitudes[0].slice_to_two('abab->ab')
        if tmp.sum() > thresh:
            raise ValueError('Warning: Diagonal cluster amplitudes not negligible. Aborting optimization!')




class PTb(Perturbation):

    @timer.with_section('PTbSolver')
    def solve(self, *args, **kwargs):
        '''Solve for energy and amplitudes

           **Arguments:**

           args
                Contains geminal coefficients, 1- and 2-el integrals, and
                auxiliary matrices

           **Keywords**
                :eref: reference energy
                :ecore: core energy
                :threshold: threshold when checking symmetry of amplitudes
                :guess: initial guess

           For more details, see :py:meth:`Pertubation.__call__`
        '''
        thresh = kwargs.get('threshold', 1e-6)
        maxiter = kwargs.get('maxiter', 200)
        eref = kwargs.get('eref', float('nan'))
        ecore = kwargs.get('ecore', float('nan'))
        initguess = kwargs.get('guess', None)

        #
        # Append el energy of principle det to function arguments:
        #
        args += ((eref-ecore),)

        #
        # Check argument type prior optimization:
        #
        check_type('args[0]', args[0], Expansion)
        check_type('args[1]', args[1], TwoIndex)
        check_type('args[2]', args[2], TwoIndex)
        check_type('args[3]', args[3], FourIndex)
        check_type('args[4]', args[4], TwoIndex)
        check_type('args[5]', args[5], ThreeIndex)
        check_type('args[6]', args[6], ThreeIndex)
        check_type('args[7]', args[7], ThreeIndex)
        check_type('args[8]', args[8], ThreeIndex)
        check_type('args[9]', args[9], TwoIndex)
        check_type('args[10]', args[10], float)
        #
        # Get initial guess:
        #
        if initguess is not None:
            guess = initguess
        else:
            guess = self.get_guess()

        #
        # Solve matrix equation using scipy.optimize.root routine:
        #
        amplitudes = opt.root(self.vfunction, guess, args=(args),
                              method='krylov',
                              options={'fatol': thresh, 'maxiter': maxiter})
        if not amplitudes.success:
            raise ValueError('ERROR: program terminated. Error in solving \
                              amplitude equations: %s' %amplitudes.message)
        if log.do_medium:
            log('Optimization of cluster amplitudes converged in %i iterations.' %(amplitudes.nit))
            log(' ')
            log('Calculating energy correction:')

        #
        # Store amplitudes
        #
        ptamplitudes = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        ptamplitudes.assign(amplitudes.x)

        #
        # Get correction to reference energy:
        #
        energy = self.calculate_energy(ptamplitudes, *args)

        return energy, ptamplitudes

    @timer.with_section('VecFctPTb')
    def vfunction(self, amplitudes, *args):
        '''.

           **Arguments:**

           amplitudes
                Cluster amplitudes. Need to be determined.

           args
                All function arguments needed to calculated the vector
                function:

                * [0]:  wfn expansion coefficients
                * [1]:  geminal coefficients
                * [2]:  1-el MO integrals
                * [3]:  2-el MO integrals
                * [4]:  inactive Fock matrix
                * [5]:  ocjbc auxilary matrix
                * [6]:  vcjkb auxilary matrix
                * [7]:  vcjbc auxilary matrix
                * [8]:  ocjkb auxilary matrix
                * [9]:  dcjb auxilary matrix
                * [10]: eref-ecore
        '''
        e0 = args[10]

        #
        # output
        #
        out = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        #
        # temporary storage
        #
        tmp4ind = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        tmp3ind0 = self.lf.create_three_index(self.nocc, self.nvirt, self.nvirt)
        tmp3ind1 = self.lf.create_three_index(self.nocc, self.nocc, self.nvirt)
        onesum = args[2].trace(0, self.nocc, 0, self.nocc)
        fockdiagsum = args[4].trace(0, self.nocc, 0, self.nocc)
        #
        # PT amplitudes
        #
        ptamplitudes = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        ptamplitudes.assign(amplitudes)

        #
        # sum_aifo A^jbkc_iaof
        #
        # sum_c F_ac t_icjb
        #
        ptamplitudes.contract_two_to_four('abcd,eb->aecd', args[4], out, 0.5, True, begin4=self.nocc, end4=self.nbasis, begin5=self.nocc, end5=self.nbasis)
        # P_ijab
        ptamplitudes.contract_two_to_four('abcd,eb->cdae', args[4], out, 0.5, False, begin4=self.nocc, end4=self.nbasis, begin5=self.nocc, end5=self.nbasis)
        #
        # sum_c F_ac t_jbic
        #
        ptamplitudes.contract_two_to_four('abcd,ed->ceab', args[4], out, 0.5, False, begin4=self.nocc, end4=self.nbasis, begin5=self.nocc, end5=self.nbasis)
        # P_ijab
        ptamplitudes.contract_two_to_four('abcd,ed->abce', args[4], out, 0.5, False, begin4=self.nocc, end4=self.nbasis, begin5=self.nocc, end5=self.nbasis)
        #
        # sum_l F_lj t_lbia
        #
        ptamplitudes.contract_two_to_four('abcd,ae->cdeb', args[4], out,-0.5, False, begin4=0, end4=self.nocc, begin5=0, end5=self.nocc)
        # P_ijab
        ptamplitudes.contract_two_to_four('abcd,ae->ebcd', args[4], out,-0.5, False, begin4=0, end4=self.nocc, begin5=0, end5=self.nocc)
        #
        # sum_l F_lj t_ialb
        #
        ptamplitudes.contract_two_to_four('abcd,ce->edab', args[4], out,-0.5, False, begin4=0, end4=self.nocc, begin5=0, end5=self.nocc)
        # P_ijab
        ptamplitudes.contract_two_to_four('abcd,ce->abed', args[4], out,-0.5, False, begin4=0, end4=self.nocc, begin5=0, end5=self.nocc)

        #
        # <jk|bc>
        #
        args[3].slice_to_four('abcd->acbd', out, 1.0, False, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # c_kc <jc||bk>
        #
        args[3].contract_two_to_four('abcd,bd->cabd', args[1], out, 1.0, False, self.nocc, self.nbasis, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis)
        args[3].contract_two_to_four('abcd,bc->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jb <jc||bk>
        #
        args[3].contract_two_to_four('abcd,ca->cabd', args[1], out, 1.0, False, self.nocc, self.nbasis, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis)
        args[3].contract_two_to_four('abcd,da->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jc <bk|cj>
        #
        args[3].contract_two_to_four('abcd,dc->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_kb <bk|cj>
        #
        args[3].contract_two_to_four('abcd,ba->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jkbc <bk|jc>
        #
        args[3].contract_two_to_four('abcd,ac->acbd', args[1], tmp4ind, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp4ind.contract_two_to_four('abcd,cd->abcd', args[1], out, 1.0, False)
        args[3].contract_two_to_four('abcd,bc->acbd', args[1], tmp4ind, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp4ind.contract_two_to_four('abcd,ad->abcd', args[1], out, 1.0, False)
        #
        # delta_jk [ c_jc F_bc ]
        #
        tmp3ind0.iadd_expand_two_two('ac,bc->abc', args[1], args[4], 1.0, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # delta_jk [ c_jb F_bc ]
        #
        tmp3ind0.iadd_expand_two_two('ab,bc->abc', args[1], args[4], 1.0, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # delta_jk [ oc_jbc ]
        #
        tmp3ind0.iadd(args[5], -1.0)
        #
        # delta_jk [ vc_jbc ]
        #
        tmp3ind0.iadd(args[7], 1.0)
        #
        # Add delta_jk-contribution
        #
        out.iadd_expand_three_to_four('1-3-1-2', tmp3ind0, 1.0)
        #
        # delta_bc [ c_jb F_jk ]
        #
        tmp3ind1.iadd_expand_two_two('ab,ac->acb', args[1], args[4],-1.0, 0, self.nocc, 0, self.nocc)
        #
        # delta_bc [ c_kb F_jk ]
        #
        tmp3ind1.iadd_expand_two_two('cb,ac->acb', args[1], args[4],-1.0, 0, self.nocc, 0, self.nocc)
        #
        # delta_bc [ vc_jkb ]
        #
        tmp3ind1.iadd(args[6],-1.0)
        #
        # delta_bc [ oc_jkb ]
        #
        tmp3ind1.iadd(args[8], 1.0)
        #
        # Add delta_bc-contribution
        #
        out.iadd_expand_three_to_four('0-2-0-1', tmp3ind1, 1.0)
        #
        # delta_bc,jk [ (sum_m h_mm+F_mm)*c_jb ]
        #
        out.iadd_expand_two_to_four('diag', args[1],(-e0+onesum+fockdiagsum))
        #
        # delta_bc,jk [ (sum_md <mm|dd> c_jmbd) ]
        #
        out.iadd_expand_two_to_four('diag', args[9], 1.0)

        return out._array.ravel(order='C')

    def vfunction_psi0(self, *args, **kwargs):
        '''Elements of <bcjk|H|AP1roG>.

           **Arguments:**

           args
                All function arguments needed to calculate the vector
                function:

                * [0]:  wfn expansion coefficients
                * [1]:  geminal coefficients
                * [2]:  1-el MO integrals
                * [3]:  2-el MO integrals
                * [4]:  inactive Fock matrix
                * [5]:  ocjbc auxilary matrix
                * [6]:  vcjkb auxilary matrix
                * [7]:  vcjbc auxilary matrix
                * [8]:  ocjkb auxilary matrix
                * [9]:  dcjb auxilary matrix
                * [10]: eref-ecore
        '''
        e0 = args[10]

        #
        # output
        #
        out = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        #
        # temporary storage
        #
        tmp4ind = self.lf.create_four_index(self.nocc, self.nvirt, self.nocc, self.nvirt)
        tmp3ind0 = self.lf.create_three_index(self.nocc, self.nvirt, self.nvirt)
        tmp3ind1 = self.lf.create_three_index(self.nocc, self.nocc, self.nvirt)
        onesum = args[2].trace(0, self.nocc, 0, self.nocc)
        fockdiagsum = args[4].trace(0, self.nocc, 0, self.nocc)

        #
        # <jk|bc>
        #
        args[3].slice_to_four('abcd->acbd', out, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # c_kc <jc||bk>
        #
        args[3].contract_two_to_four('abcd,bd->cabd', args[1], out, 1.0, False, self.nocc, self.nbasis, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis)
        args[3].contract_two_to_four('abcd,bc->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jb <jc||bk>
        #
        args[3].contract_two_to_four('abcd,ca->cabd', args[1], out, 1.0, False, self.nocc, self.nbasis, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis)
        args[3].contract_two_to_four('abcd,da->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jc <bk|cj>
        #
        args[3].contract_two_to_four('abcd,dc->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_kb <bk|cj>
        #
        args[3].contract_two_to_four('abcd,ba->dabc', args[1], out,-1.0, False, self.nocc, self.nbasis, 0, self.nocc, self.nocc, self.nbasis, 0, self.nocc)
        #
        # c_jkbc <bk|jc>
        #
        args[3].contract_two_to_four('abcd,ac->acbd', args[1], tmp4ind, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp4ind.contract_two_to_four('abcd,cd->abcd', args[1], out, 1.0, False)
        args[3].contract_two_to_four('abcd,bc->acbd', args[1], tmp4ind, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp4ind.contract_two_to_four('abcd,ad->abcd', args[1], out, 1.0, False)
        #
        # delta_jk [ c_jc F_bc ]
        #
        tmp3ind0.iadd_expand_two_two('ac,bc->abc', args[1], args[4], 1.0, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # delta_jk [ c_jb F_bc ]
        #
        tmp3ind0.iadd_expand_two_two('ab,bc->abc', args[1], args[4], 1.0, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # delta_jk [ oc_jbc ]
        #
        tmp3ind0.iadd(args[5], -1.0)
        #
        # delta_jk [ vc_jbc ]
        #
        tmp3ind0.iadd(args[7], 1.0)
        #
        # Add delta_jk-contribution
        #
        out.iadd_expand_three_to_four('1-3-1-2', tmp3ind0, 1.0)
        #
        # delta_bc [ c_jb F_jk ]
        #
        tmp3ind1.iadd_expand_two_two('ab,ac->acb', args[1], args[4],-1.0, 0, self.nocc, 0, self.nocc)
        #
        # delta_bc [ c_kb F_jk ]
        #
        tmp3ind1.iadd_expand_two_two('cb,ac->acb', args[1], args[4],-1.0, 0, self.nocc, 0, self.nocc)
        #
        # delta_bc [ vc_jkb ]
        #
        tmp3ind1.iadd(args[6],-1.0)
        #
        # delta_bc [ oc_jkb ]
        #
        tmp3ind1.iadd(args[8], 1.0)
        #
        # Add delta_bc-contribution
        #
        out.iadd_expand_three_to_four('0-2-0-1', tmp3ind1, 1.0)
        #
        # delta_bc,jk [ (sum_m h_mm+F_mm)*c_jb ]
        #
        out.iadd_expand_two_to_four('diag', args[1],(-e0+onesum+fockdiagsum))
        #
        # delta_bc,jk [ (sum_md <mm|dd> c_jmbd) ]
        #
        out.iadd_expand_two_to_four('diag', args[9], 1.0)

        return out

    def get_guess(self):
        '''Generate initial guess for amplitudes'''
        tmp = np.random.rand(self.nocc*self.nvirt, self.nocc*self.nvirt)*0.01
        tmp = (tmp+tmp.T)/2
        np.fill_diagonal(tmp, 0.0)
        return tmp.ravel()

    def calculate_aux_matrix(self, *args):
        '''Compute auxiliary matrices

           **Arguments:**

           mo1, mo2
                One- and two-electron integrals (some Hamiltonian matrix
                elements) in the MO basis.

           args
                List of arguements. Only geminal coefficients are used.
        '''
        check_type('args[1]', args[1], TwoIndex)
        check_type('args[2]', args[2], TwoIndex)
        check_type('args[3]', args[3], FourIndex)
        self.clear_aux_matrix()
        return self.update_aux_matrix(args[2], args[3], args[1])

    def init_aux_matrix(self, select, dim=None, dim2=None):
        '''Initialize auxiliary matrices

           **Arguments:**

           select
                One of ``fock``, ``ocjbc``, ``vcjkb``, ``ocjkb``,
                ``vcjbc``, ``dcjb``

           **Optional arguments:**

           dim, dim2
                The dimension of the auxiliary matrices for specific axes.
        '''
        check_options('select', select, 'fock', 'ocjbc', 'vcjkb', 'ocjkb',
            'vcjbc', 'dcjb')
        if select=='fock':
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_two_index, self.lf.default_nbasis), tags='m')
        if select=='dcjb':
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_two_index, dim, dim2), tags='m')
        if select in ['ocjbc', 'vcjbc']:
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_three_index, dim, dim2, dim2), tags='m')
        if select in ['vcjkb', 'ocjkb']:
            matrix, new = self._cache.load('%s' %select, alloc=(self.lf.create_three_index, dim, dim, dim2), tags='m')
        if not new:
            raise RuntimeError('The matrix %s already exists. Call clear prior \
                                to updating the wfn.' % select)
        return matrix

    @timer.with_section('UpdatingAuxMat')
    def update_aux_matrix(self, mo1, mo2, cia):
        '''Derive all matrices.
           fock_pq:     one_pq + sum_m(2<pm|qm> - <pm|mq>),
           oc_jbc:      sum_m(<mm|bc> c_jm^bc),
           vc_jkb:      sum_d(<dd|jk> c_jk^bd),
           vc_jbc:      sum_d(<dd|bc> c_j^d),
           oc_jkb:      sum_m(<mm|jk> c_m^b),
           dc_jb:       sum_md(<mm|dd> c_jm^bd),

           **Arguments:**

           mo1, mo2
                one- and two-electron integrals to be sorted.

           cia
                The geminal coefficients. A TwoIndex instance
        '''
        #
        # Inactive Fock matrix
        #
        auxmat1 = self.init_aux_matrix('fock')
        mo2.contract_to_two('abcb->ac', auxmat1, 2.0, True, 0, self.nbasis, 0, self.nocc, 0, self.nbasis, 0, self.nocc)
        mo2.contract_to_two('abbc->ac', auxmat1,-1.0, False, 0, self.nbasis, 0, self.nocc, 0, self.nocc, 0, self.nbasis)
        auxmat1.iadd(mo1)
        #
        # oc_jbc = sum_m <mm|bc> c_jm^bc
        #
        tmp = self.lf.create_four_index(self.nocc, self.nocc, self.nvirt, self.nvirt)
        auxmat2 = self.init_aux_matrix('ocjbc', cia.nbasis, cia.nfn)
        mo2.contract_two_to_four('aabc,db->adbc', cia, tmp, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp.contract_two_to_three('abcd,ad->bcd', cia, auxmat2, 1.0, True)
        mo2.contract_two_to_four('aabc,dc->adbc', cia, tmp, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp.contract_two_to_three('abcd,ac->bcd', cia, auxmat2, 1.0, False)
        #
        # vc_jkb = sum_d <dd|jk> c_jk^bd
        #
        auxmat3 = self.init_aux_matrix('vcjkb', cia.nbasis, cia.nfn)
        mo2.contract_two_to_four('abcc,ad->abcd', cia, tmp, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp.contract_two_to_three('abcd,bc->abd', cia, auxmat3, 1.0, True)
        mo2.contract_two_to_four('abcc,bd->abcd', cia, tmp, 1.0, True, 0, self.nocc, 0, self.nocc, self.nocc, self.nbasis, self.nocc, self.nbasis)
        tmp.contract_two_to_three('abcd,ac->abd', cia, auxmat3, 1.0, False)
        #
        # vc_jbc = sum_d <bc|dd> c_j^d
        #
        auxmat4 = self.init_aux_matrix('vcjbc', cia.nbasis, cia.nfn)
        mo2.contract_two_to_three('abcc,dc->dab', cia, auxmat4, 1.0, True, self.nocc, self.nbasis, self.nocc, self.nbasis, self.nocc, self.nbasis, self.nocc, self.nbasis)
        #
        # oc_jkb = sum_m <mm|jk> c_m^b
        #
        auxmat5 = self.init_aux_matrix('ocjkb', cia.nbasis, cia.nfn)
        mo2.contract_two_to_three('aabc,ad->bcd', cia, auxmat5, 1.0, True, 0, self.nocc, 0, self.nocc, 0, self.nocc, 0, self.nocc)
        #
        # dc_jb = sum_md <mm|dd> c_jm^bd
        #
        auxmat6= self.init_aux_matrix('dcjb', cia.nbasis, cia.nfn)
        tmp = self.lf.create_two_index(self.nbasis, self.nbasis)
        tmp2 = self.lf.create_two_index(self.nocc, self.nocc)
        # There is a bug in np.einsum that forces us to slice first...
        mo2.slice_to_two('aabb->ab', tmp, 1.0, True)
        factor = cia.contract_two('ab,ab', tmp, 0, self.nocc, self.nocc, self.nbasis)
        auxmat6.iadd(cia, factor)
        cia.contract_two_to_two('ab,cb->ac', tmp, tmp2, 1.0, True, 0, self.nocc, self.nocc, self.nbasis)
        cia.contract_two_to_two('ab,ca->cb', tmp2, auxmat6, 1.0, False)

        return [auxmat1, auxmat2, auxmat3, auxmat4, auxmat5, auxmat6]

    def update_amplitudes(self, new):
        '''Update cluster amplitudes

           **Arguments:**

           new
                PT amplitudes. A FourIndex instance
        '''
        if self.amplitudes:
            raise ValueError('Warning: List of PT amplitudes not empty!')
        else:
            self.amplitudes.append(new)

    def calculate_energy(self, amplitudes, *args):
        '''Calculate PT energy and energy contribution of seniority sectors

           **Arguments:**

           amplitudes
                PT amplitudes. A FourIndex instance

           args
                List containing the geminal coefficients, 1- and 2-el
                integrals, and auxiliary matrices
        '''
        #
        # Get "excitation" matrices w.r.t. |AP1roG> (psi0)
        #
        exjbkc0 = self.vfunction_psi0(*args)
        #
        # Seniority-0 sector
        #
        e_seniority_0 = amplitudes.contract_four('abab,abab', exjbkc0, 1.0)
        #
        # Seniority-2 sector
        #
        e_seniority_2 = amplitudes.contract_four('abad,abad', exjbkc0, 1.0)
        e_seniority_2+= amplitudes.contract_four('abdb,abdb', exjbkc0, 1.0)
        energy = amplitudes.contract_four('abcd,abcd', exjbkc0, 2.0)
        energy-= amplitudes.contract_four('abcd,adcb', exjbkc0, 1.0)
        #
        # Seniority-4 sector
        #
        e_seniority_4 = energy-e_seniority_2-e_seniority_0

        return [energy, e_seniority_0, e_seniority_2, e_seniority_4]

    def print_energy(self, **kwargs):
        if log.do_medium:
            log('E_PTb(Seniority 0):     %18.12f a.u.' %(self.energy[1]))
            log('E_PTb(Seniority 2):     %18.12f a.u.' %(self.energy[2]))
            log('E_PTb(Seniority 4):     %18.12f a.u.' %(self.energy[3]))
            log.hline('-')
            log('E_ref:                  %18.12f a.u.' %(kwargs.get('eref')))
            log('E_PTb:                  %18.12f a.u.' %self.energy[0])
            log.hline('-')
            log('E_tot:                  %18.12f a.u.' %(self.energy[0]+kwargs.get('eref')))

    def print_info(self, **kwargs):
        thresh = kwargs.get('threshold', 1e-6)
        maxiter = kwargs.get('maxiter', 200)
        guess = kwargs.get('guess', None)
        if log.do_medium:
            log('PTb perturbation module')
            log(' ')
            log('OPTIMIZATION PARAMETERS:')
            log('Reference Function           %s' %('AP1roG'))
            log('Number of pairs:             %i' %self.nocc)
            log('Number of virtuals:          %i' %self.nvirt)
            if guess is None:
                log('Initial guess:               %s' %('random'))
            else:
                log('Initial guess:               %s' %('user'))
            log('Solvers:')
            log('  PTb amplitudes:            %s' %('krylov'))
            log('Optimization thresholds:')
            log('  PTb amplitudes:            %1.2e' %(thresh))
            log('  maxiter:                   %i' %(maxiter))
            log.hline()

    def check_input(self, **kwargs):
        '''Check input parameters.'''
        for name, value in kwargs.items():
            check_options(name, name, 'ecore', 'eref', 'threshold', 'maxiter',
                'guess')
        eref = kwargs.get('eref', float('nan'))
        ecore = kwargs.get('ecore', float('nan'))
        guess = kwargs.get('guess', None)
        if guess is not None:
            check_type('guess', guess, np.ndarray)
            if not len(guess) == self.nocc*self.nocc*self.nvirt*self.nvirt:
                raise ValueError('Length of guess array does not agree with number of unknowns')
        if math.isnan(eref):
            raise ValueError('Warning: Cannot find reference energy in PTb module!')
        if math.isnan(ecore):
            raise ValueError('Warning: Cannot find core: energy in PTb module!')

    def check_result(self, **kwargs):
        '''Check if amplitudes are reasonable.'''
        thresh = kwargs.get('threshold', 1e-6)

        # check symmetry of amplitudes:
        if not self.amplitudes[0].is_symmetric('cdab', thresh):
            raise ValueError('Warning: Cluster amplitudes not symmetric!')
