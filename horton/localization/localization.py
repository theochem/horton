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
#--
'''Orbital localization procedures

   Currently supported localization flavours:
    * Pipek-Mezey

   This is work in progress.
'''

import numpy as np

from horton.cache import Cache
from horton.matrix.base import Expansion
from horton.log import log, timer
from horton.orbital_utils import compute_unitary_matrix
from horton.correlatedwfn.stepsearch import RStepSearch
from horton.utils import check_options


__all__ = [
    'Localization',
    'PipekMezey',
]


class Localization(object):
    '''Base class for all localization methods'''
    def __init__(self, lf, occ_model, projector):
        '''Localize canonical HF orbitals.

           **Arguments:**

           lf
                A LinalgFactory instance.

           occ_model
                Occupation model.

           Projector
                Projectors for atomic basis function. A list of TwoIndex
                instances.

           **Optional arguments:**

        '''
        self._lf = lf
        self._proj = projector
        self._nocc = occ_model.noccs[0]
        self._nbasis = lf.default_nbasis
        self._nvirt = (lf.default_nbasis-occ_model.noccs[0])
        self._cache = Cache()
        self._locblock = None
        self._popmatrix = None

    @timer.with_section('Localization')
    def __call__(self, orb, select, **kwargs):
        '''Localizes the orbitals using a unitary transformation to rotate the
           AO/MO coefficient matrix. The orbitals are optimized by minimizing
           an objective function.

           This works only for restricted orbitals.

           **Arguments:**

           orb
                The AO/MO coefficients. An Expansion instance.

           select
                The orbital block to be localised (str). Any of ``occ`` (occupied
                orbitals), ``virt`` (virtual orbitals)

           **Keywords:**

           :maxiter: (int) maximum number of iterations for localization
                     (default 2000)
           :threshold: (float) localization threshold for objective function
                       (default 1e-6)
           :levelshift: level shift of Hessian (float) (default 1e-8)
           :stepsearch: step search options (dictionary) containing:

                        * method: step search method used (str). One of
                          ``trust-region`` (default), ``None``, ``backtracking``
                        * alpha: scaling factor for Newton step (float), used in
                          ``backtracking`` and ``None`` method (default 0.75)
                        * c1: parameter used in ``backtracking`` (float)
                          (default 1e-4)
                        * minalpha: minimum step length used in ``backracking``
                          (float) (default 1e-6)
                        * maxiterouter: maximum number of search steps (int)
                          (default 10)
                        * maxiterinner: maximum number of optimization
                          steps in each search step (int) (used only in ``pcg``,
                          default 500)
                        * maxeta: upper bound for estimated vs actual change in
                          ``trust-region`` (float) (default 0.75)
                        * mineta: lower bound for estimated vs actual change in
                          ``trust-region`` (float) (default 0.25)
                        * upscale: scaling factor to increase trustradius in
                          ``trust-region`` (float) (default 2.0)
                        * downscale: scaling factor to decrease trustradius in
                          ``trust-region`` (float) and scaling factor in
                          ``backtracking`` (default 0.25)
                        * trustradius: initial trustradius (float) (default
                          0.75)
                        * maxtrustradius: maximum trustradius (float) (default
                          0.75)
                        * threshold: trust-region optimization threshold, only
                          used in ``pcg`` (float) (default 1e-8)
                        * optimizer: optimizes step to boundary of trustradius
                          (str). One of ``pcg``, ``dogleg``, ``ddl`` (default
                          ddl)
        '''
        if log.do_medium:
            log('Performing localization of %s block' %(select))
        log.cite('pipek1989', 'the Pipek-Mezey localization scheme')
        #
        # Assign default keyword arguements
        #
        names = []
        def _helper(x,y):
            names.append(x)
            return kwargs.get(x,y)
        maxiter = _helper('maxiter', 2000)
        thresh = _helper('threshold', 1e-6)
        lshift = _helper('levelshift', 1e-8)
        stepsearch = _helper('stepsearch', dict({}))
        stepsearch.setdefault('method', 'trust-region')
        stepsearch.setdefault('minalpha', 1e-6)
        stepsearch.setdefault('alpha', 1.0)
        stepsearch.setdefault('c1', 0.0001)
        stepsearch.setdefault('maxiterouter', 10)
        stepsearch.setdefault('maxiterinner', 500)
        stepsearch.setdefault('maxeta', 0.75)
        stepsearch.setdefault('mineta', 0.25)
        stepsearch.setdefault('upscale', 2.0)
        stepsearch.setdefault('downscale', 0.25)
        stepsearch.setdefault('trustradius', 0.75)
        stepsearch.setdefault('maxtrustradius', 0.75)
        stepsearch.setdefault('threshold', 1e-8)
        stepsearch.setdefault('optimizer', 'ddl')

        for name, value in kwargs.items():
            if name not in names:
                raise ValueError("Unknown keyword argument %s" % name)
            if value < 0:
                raise ValueError('Illegal value for %s: %s' %(name, value))

        #
        # Update information about localization block
        #
        self.update_locblock(select)

        if log.do_medium:
            log('%3s  %12s  %10s' %('Iter', 'D(ObjectiveFunction)', 'Steplength'))
        #
        # Initialize step search
        #
        stepsearch_ = RStepSearch(self.lf, **stepsearch)
        #
        # Calculate initial objective function
        #
        self.solve_model(orb)
        objfct_ref = self.compute_objective_function()

        maxThresh = True
        maxIter = True
        it = 0
        while maxThresh and maxIter:
            #
            # Update population matrix for new orbitals
            #
            self.compute_population_matrix(orb)
            #
            # Calculate orbital gradient and diagonal approximation to the Hessian
            #
            kappa, gradient, hessian = self.orbital_rotation_step(lshift)
            #
            # Apply steps search to orbital rotation step 'kappa' and perform
            # orbital rotation
            #
            stepsearch_(self, None, None, orb,
                       **{'kappa': kappa, 'gradient': gradient, 'hessian': hessian
                         })
            #
            # update objective function
            #
            objfct = self.compute_objective_function()
            it += 1
            #
            # Print localization progress
            #
            if log.do_medium:
                log('%4i   %14.8f' %(it, abs(objfct-objfct_ref)))
            #
            # Check convergence
            #
            maxThresh = abs(objfct-objfct_ref)>thresh
            maxIter = it<maxiter
            #
            # Prepare for new iteration
            #
            objfct_ref = objfct
        if maxThresh and not maxIter:
            if log.do_medium:
                log(' ')
                log('Warning: Orbital localization not converged in %i iteration' %(it-1))
                log(' ')
        else:
            if log.do_medium:
                log(' ')
                log('Orbital localization converged in %i iteration' %(it-1))
                log(' ')

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

    def _get_lf(self):
        '''The LinalgFactory'''
        return self._lf

    lf = property(_get_lf)

    def _get_proj(self):
        '''The Projectors. A list of TwoIndex instances'''
        return self._proj

    proj = property(_get_proj)

    def _get_locblock(self):
        '''The orbital block to be localized'''
        return self._locblock

    locblock = property(_get_locblock)

    def _get_popmatrix(self):
        '''The population matrix. A list of TwoIndex instances'''
        return self._popmatrix

    popmatrix = property(_get_popmatrix)

    def update_locblock(self, new):
        '''Update localization block'''
        self._locblock = new

    def __clear__(self):
        self.clear()

    def clear(self):
        '''Clear all wavefunction information'''
        self._cache.clear()

    def compute_rotation_matrix(self, coeff):
        '''Determine orbital rotation matrix

           **Arguments:**

           coeff
                The non-reduntant orbital rotations, we need only values for
                p<q
        '''
        indl = np.tril_indices(self.nbasis, -1)
        kappa = self.lf.create_two_index(self.nbasis, self.nbasis)
        #
        # k_pq = -k_qp
        #
        kappa.assign(coeff, indl)
        kappa.iadd_t(kappa, -1.0)

        out = compute_unitary_matrix(kappa)
        return out

    def compute_population_matrix(self, exp):
        '''Determine population matrix

           **Arguments:**

           exp
                The current AO/MO coefficients. An Expansion instance
        '''
        #
        # Get orbital block to be localized, a OneIndex instance
        #
        block = self.assign_locblock()
        #
        # Calculate population matrices for orbital block
        #
        popmat = []
        for op in self.proj:
            pop = self.lf.create_two_index()
            expblock = exp.copy()
            expblock.imul(block)
            expblock.itranspose()
            pop.assign_dot(expblock, op)
            expblock.itranspose()
            pop.idot(expblock)

            popmat.append(pop)
        self._popmatrix = popmat


class PipekMezey(Localization):
    '''Perform Pipek-Mezey localization of occupied or virtual Hartree-Fock
       orbitals.
    '''

    def assign_locblock(self):
        '''Get localization block. A OneIndex instance'''
        check_options('block', self.locblock, 'occ', 'virt')
        block = self.lf.create_one_index()
        if self.locblock=='occ':
            block.assign(1.0, end0=self.nocc)
        elif self.locblock=='virt':
            block.assign(1.0, begin0=self.nocc)
        return block

    def grad(self):
        '''Gradient of objective function

           **Arguments:**

           popmatrix
                The population matrix. A list of TwoIndex instances
        '''
        grad = self.lf.create_two_index()
        for pop in self.popmatrix:
            #
            # 4 [ pop_kk - pop_ll ] * pop_kl
            #
            diag = pop.copy_diagonal()
            grad.iadd_contract_two_one('ab,a->ab', pop, diag, 4.0)
            grad.iadd_contract_two_one('ab,b->ab', pop, diag,-4.0)
        grad.iscale(-1)

        ind = np.tril_indices(self.nbasis, -1)
        return grad.copy_slice(ind)

    # TODO: use exact Hessian
    def hessian(self, lshift=1e-8):
        '''Diagonal Hessian of objective function

           **Arguments:**

           popmatrix
                The population matrix. A list of TwoIndex instances

           **Optinal arguments:**

           lshift
                Shift elements of Hessian (float)
        '''
        hessian = self.lf.create_two_index()
        for pop in self.popmatrix:
            #
            # pop_ll*pop_ll
            #
            diag2 = self.lf.create_one_index()
            diag = pop.copy_diagonal()
            diag.mult(diag, diag2)
            #
            # H_kl += -4 pop_ll*pop_ll - 4 pop_kk*pop_kk
            #
            hessian.iadd_t(diag2, 4.0)
            hessian.iadd(diag2, 4.0)
            #
            # H_kl += 8 pop_ll*pop_kk
            #
            hessian.iadd_one_mult(diag, diag,-8.0, transpose0=True)
            #
            # H_kl += 8 pop_kl*pop_kl
            #
            hessian.iadd_mult(pop, pop, 8.0)
        hessian.iadd_shift(lshift)

        ind = np.tril_indices(self.nbasis, -1)
        return hessian.copy_slice(ind)

    def compute_objective_function(self):
        '''Objective function of PM localization to be maximized. The current
           implementation minimizes -(objective_function).
        '''
        #
        # Get orbital block to be localized, a OneIndex instance
        #
        block = self.assign_locblock()
        result = 0.0
        #
        # sum_iA pop(A)_ii^2
        #
        for pop in self.popmatrix:
            diag2 = self.lf.create_one_index()
            diag = pop.copy_diagonal()
            diag.mult(diag, diag2)

            result += diag2.trace()
        return -result

    #
    # Don't change function name or implementation will break. The function
    # ``solve_model`` is used in the StepSearch module.
    #
    def solve_model(self, *args, **kwargs):
        '''Update population matrix used to calculate objective function'''
        for arg in args:
            if isinstance(arg, Expansion):
                self.compute_population_matrix(arg)

    def orbital_rotation_step(self, lshift):
        '''Calculate gradient, hessian, and orbital rotation step

           **Arguments:**

           lshift
                Shift elements of Hessian
        '''
        grad = self.grad()
        hess = self.hessian(lshift)

        kappa = grad.divide(hess, -1.0)

        return kappa, grad, hess
