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
'''Correlated wavefunction implementations

   This module contains restricted AP1roG

   Variables used in this module:
    :nocc:      number of (active) occupied orbitals in the principle configuration
    :pairs:     number of electron pairs
    :nvirt:     number of (active) virtual orbitals in the principle configuration
    :nbasis:    total number of basis functions

    Indexing convention:
     :i,j,k,..: occupied orbitals of principle configuration
     :a,b,c,..: virtual orbitals of principle configuration
     :p,q,r,..: general indices (occupied, virtual)

    Abbreviations used (if not mentioned in doc-strings):
     :L_pqrs: 2<pq|rs>-<pq|sr>
     :g_pqrs: <pq|rs>

    For more information see doc-strings.
'''

import numpy as np
import math as math
from scipy import optimize as opt
import warnings

from horton.cache import Cache
from horton.log import log, timer
from horton.orbital_utils import compute_unitary_matrix, transform_integrals
from horton.correlatedwfn.geminal import Geminal
from horton.correlatedwfn.stepsearch import RStepSearch
from horton.utils import check_type, check_options

from copy import copy


__all__ = [
    'RAp1rog',
]



class RAp1rog(Geminal):
    '''Restricted AP1roG wavefunction class'''

    @timer.with_section('AP1roG')
    def solve(self, one, two, core, orb, olp, **kwargs):
        '''Optimize AP1roG coefficients for some Hamiltonian.
           For restricted, closed-shell AP1roG.

           **Arguments:**

           one, two
                One- (TwoIndex instance) and two-body integrals (FourIndex or
                Cholesky instance) (some Hamiltonian matrix elements)

           core
                The core energy (not included in 'one' and 'two'); a float.

           orb
                An expansion instance. It contains the MO coefficients.

           olp
                The AO overlap matrix. A TwoIndex instance.

           **Keywords:**

            :indextrans: 4-index Transformation (str). Choice between
                         ``tensordot`` (default) and ``einsum``
            :warning: Print warnings (boolean) (default False)
            :guess: initial guess (dictionary) containing:

                    * type: guess type (str). One of ``random`` (random
                      numbers, default), ``const`` (a constant scaled
                      by a factor)
                    * factor: a scaling factor (float) (default -0.1)
                    * geminal: external guess for geminal coefficients
                      (1-d np.array); if provided, 'type' is
                      ignored (default None)

            :solver: wfn solver (dictionary) containing:

                     * wfn: wavefunction solver (str) (default ``krylov``)

            :maxiter: max number of iterations (dictionary) containing:

                      * wfniter: maximum number of iterations (int) for wfn
                        solver (default 200)

            :dumpci: dump ci coefficients (dictionary):

                     * amplitudestofile: write wfn amplitudes to file
                       (boolean) (default False)
                     * amplitudesfilename: (str) (default ap1rog_amplitudes.dat)

            :thresh: optimization thresholds (dictionary) containing:

                     * wfn: threshold for geminal coefficients (float)
                       (default 1e-12)

            :printoptions: print level (dictionary) containing:

                           * geminal: (boolean), if True, geminal matrix is
                             printed (default True)

                           * ci: threshold for CI coefficients (float) (requires
                             evaluation of a permanent), all coefficients
                             (for a given excitation order) larger than
                             'ci' are printed (default 0.01)
                           * excitationlevel: number of excited pairs w.r.t.
                             reference determinant for which wfn amplitudes
                             are reconstructed (int) (default 1)

            :swapa: swap orbitals (numpy 2-dim array), each row contains 2
                    orbital indices to be swapped (default np.array([[]]))
        '''
        if log.do_medium:
            log.hline('=')
            log('Entering AP1roG module')

        #
        # Assign keyword arguements
        #
        names = []
        def _helper(x,y):
            names.append(x)
            return kwargs.get(x,y)
        indextrans = _helper('indextrans', 'tensordot')
        warning = _helper('warning', False)
        swapa = _helper('swapa', np.array([[]]))
        guess = _helper('guess', {})
        guess.setdefault('type', 'random')
        guess.setdefault('factor', -0.1)
        guess.setdefault('geminal', None)
        solver = _helper('solver', {'wfn': 'krylov'})
        solver.setdefault('wfn', 'krylov')
        maxiter = _helper('maxiter', {'wfniter': 200})
        dumpci = _helper('dumpci', {})
        dumpci.setdefault('amplitudestofile', False)
        dumpci.setdefault('amplitudesfilename', "./ap1rog_amplitudes.dat")
        thresh = _helper('thresh', {'wfn':  1e-12})
        printoptions = _helper('printoptions', {})
        printoptions.setdefault('geminal', True)
        printoptions.setdefault('ci', 0.01)
        printoptions.setdefault('excitationlevel', 1)
        for name, value in kwargs.items():
            if name not in names:
                raise ValueError("Unknown keyword argument %s" % name)

        #
        # Check dictionaries in keyword arguments
        #
        self.check_keywords(guess, solver, maxiter, dumpci, thresh,
                            printoptions)
        check_options('warning', warning, False, True, 0, 1)

        if not warning:
            warnings.filterwarnings('ignore')

        #
        # Print optimization parameters
        #
        if log.do_medium:
            self.print_options(guess, solver, maxiter, thresh, printoptions,
                               indextrans)

        #
        # Do restart and swap orbitals if needed:
        #
        self.start_up(orb, swapa)

        #
        # Generate Guess
        #
        if guess['geminal'] is None:
            initial_guess = self.generate_guess(guess)
        else:
            check_type('guess.geminal', guess['geminal'], np.ndarray)
            if not guess['geminal'].shape[0] == self.dimension:
                raise ValueError('Size of geminal guess array does not match number of unknowns.')
            initial_guess = guess['geminal']

        #
        # Update core energy
        #
        self.update_ecore(core)

        #
        # Solve for wavefunction
        #
        self.solve_model(one, two, orb, **{'maxiter': maxiter, 'thresh': thresh,
                         'guess': initial_guess,
                         'solver': solver, 'indextrans': indextrans,
                         'orbitaloptimizer': False})

        #
        # Final print statements:
        #
        if log.do_medium:
            self.print_final()
            self.dump_final(orb, olp, printoptions, dumpci, -1)
        #
        # Sanity check for correlation energy:
        #
        if self.compute_correlation_energy() > 0:
            raise ValueError('Warning: Final correlation energy is positive! \
                              Improve initial guess!')

        return self.compute_total_energy(), self.geminal

    @timer.with_section('AP1roG-SCF')
    def solve_scf(self, one, two, core, orb, olp, **kwargs):
        '''Find Geminal expansion coefficient for some Hamiltonian.
           For restricted, closed-shell AP1roG.
           Perform orbital optimization.

           **Arguments:**

           one, two
                One- (TwoIndex instance) and two-body integrals (FourIndex or
                Cholesky instance) (some Hamiltonian matrix elements)

           core
                The core energy (not included in 'one' and 'two').

           orb
                An expansion instance. It contains the MO coefficients.

           olp
                The AO overlap matrix. A TwoIndex instance.

           **Keywords:**

            :indextrans: 4-index Transformation (str). Choice between
                         ``tensordot`` (default) and ``einsum``
            :warning: print warnings (boolean) (default False)
            :guess: initial guess (dictionary) containing:

                    * type: guess type (str). One of ``random`` (random
                      numbers, default), ``const`` (a constant scaled by a
                      factor)
                    * factor: a scaling factor (float) (default -0.1)
                    * geminal: external guess for geminal coefficients
                      (1-dim np.array) (default None)
                    * lagrange: external guess for Lagrange multipliers
                      (1-dim np.array) (default None)

            :solver: wfn/Lagrange solver (dictionary) containing:

                     * wfn: wavefunction solver (str) (default ``krylov``)
                     * lagrange: Lagrange multiplier solver (str)
                       (default krylov)

            :maxiter: max number of iterations (dictionary) containing:

                      * wfniter: maximum number of iterations for
                                 wfn/lagrange solver (int) (default 200)
                      * orbiter: maximum number of orbital optimization
                                 steps (int) (default 100)
            :dumpci: dump ci coefficient (dictionary) containing:

                     * amplitudestofile: write wfn amplitudes to file
                       (boolean) (default False)
                     * amplitudesfilename: (str) (default ap1rog_amplitudes.dat)

            :thresh: optimization thresholds (dictionary) containing:

                     * wfn: threshold for geminal coefficients and
                       Lagrange multipliers (float) (defaulf 1e-12)
                     * energy: convergence threshold for energy (float)
                       (default 1e-8)
                     * gradientnorm: convergence threshold for norm of
                       orbital gradient (float) (default 1e-4)
                     * gradientmax: threshold for maximum absolute value of
                       orbital gradient (float) (default 5e-5)

            :printoptions: print level; dictionary containing:

                           * geminal: (boolean) if True, geminal matrix is
                             printed (default True)
                           * ci: threshold for CI coefficients (requires
                             evaluation of a permanent) (float). All
                             coefficients (for a given excitation order)
                             larger than ci are printed (default 0.01)
                           * excitationlevel: number of excited pairs w.r.t.
                             reference determinant for which
                             wfn amplitudes are reconstructed
                             (int) (default 1)

            :stepsearch: step search options (dictionary) containing:

                         * method: step search method used (str). One of
                           ``trust-region`` (default), ``None``,
                           ``backtracking``
                         * alpha: scaling factor for Newton step (float),
                           used in ``backtracking`` and ``None`` method (default
                           1.00)
                         * c1: parameter used in ``backtracking`` (float)
                           (default 1e-4)
                         * minalpha: minimum scaling factor used in
                           ``backracking`` (float) (default 1e-6)
                         * maxiterouter: maximum number of search steps
                           (int) (default 10)
                         * maxiterinner: maximum number of optimization
                           step in each search step (int) (used only in ``pcg``,
                           default 500)
                         * maxeta: upper bound for estimated vs actual
                           change in ``trust-region`` (float) (default 0.75)
                         * mineta: lower bound for estimated vs actual change in
                           ``trust-region`` (float) (default 0.25)
                         * upscale: scaling factor to increase trustradius
                           in ``trust-region`` (float) (default 2.0)
                         * downscale: scaling factor to decrease trustradius
                           in ``trust-region`` and step length in
                           ``backtracking`` (float) (default 0.25)
                         * trustradius: initial trustradius (float) (default
                           0.75)
                         * maxtrustradius: maximum trustradius (float) (default
                           0.75)
                         * threshold: trust-region optimization threshold, only
                           used in ``pcg`` (float) (default 1e-8)
                         * optimizer: optimizes step to boundary of trustradius
                           (str). One of ``pcg``, ``dogleg``, ``ddl`` (default
                           ddl)

            :checkpoint: frequency of checkpointing (int). If > 0, writes
                         orbitals and overlap to a checkpont file (defatul 1)
            :checkpoint_fn: filename to use for the checkpoint file (default
                            "checkpoint.h5")
            :levelshift: level shift of Hessian (float) (default 1e-8)
            :absolute: (boolean), if True, take absolute values of Hessian
                       (default False)
            :sort: (boolean), if True, orbitals are sorted according to their
                   natural occupation numbers. This requires us to solve for
                   the wavefunction again. Works only if orbitaloptimizer
                   is set to ``variational``. (default True)
            :swapa: swap orbitals (numpy 2-dim array) each row contains 2
                    orbital indices to be swapped (default np.array([[]]))
            :givensrot: rotate two orbitals (numpy 2-dim array) each row
                        contains 2 orbital indices and rotation angle
                        (default np.array([[]]))
            :orbitaloptimizer: (str) switch between variational orbital
                               optimization (``variational``) and PS2c
                               (``ps2c``) (default ``variational``).
        '''
        if log.do_medium:
            log.hline('=')
            log('Entering orbital-optimized AP1roG module')

        #
        # Assign keyword arguements, checks only name of dictionary, keys and
        # values are checked below.
        #
        names = []
        def _helper(x,y):
            names.append(x)
            return kwargs.get(x,y)
        indextrans = _helper('indextrans', 'tensordot')
        warning = _helper('warning', False)
        checkpoint = _helper('checkpoint', 1)
        checkpoint_fn = _helper('checkpoint_fn', 'checkpoint.h5')
        lshift = _helper('levelshift', 1e-8)
        pos = _helper('absolute', False)
        givensrot = _helper('givensrot', np.array([[]]))
        swapa = _helper('swapa', np.array([[]]))
        sort = _helper('sort', True)
        guess = _helper('guess', {})
        guess.setdefault('type', 'random')
        guess.setdefault('factor', -0.1)
        guess.setdefault('geminal', None)
        guess.setdefault('lagrange', None)
        solver = _helper('solver', {})
        solver.setdefault('wfn', 'krylov')
        solver.setdefault('lagrange', 'krylov')
        maxiter = _helper('maxiter', {})
        maxiter.setdefault('wfniter', 200)
        maxiter.setdefault('orbiter', 100)
        dumpci = _helper('dumpci', {})
        dumpci.setdefault('amplitudestofile', False)
        dumpci.setdefault('amplitudesfilename', "./ap1rog_amplitudes.dat")
        thresh = _helper('thresh', {})
        thresh.setdefault('wfn',  1e-12)
        thresh.setdefault('energy', 1e-8)
        thresh.setdefault('gradientnorm', 1e-4)
        thresh.setdefault('gradientmax', 5e-5)
        printoptions = _helper('printoptions', {})
        printoptions.setdefault('geminal', True)
        printoptions.setdefault('ci', 0.01)
        printoptions.setdefault('excitationlevel', 1)
        stepsearch = _helper('stepsearch', {})
        stepsearch.setdefault('method', 'trust-region')
        stepsearch.setdefault('alpha', 1.0)
        stepsearch.setdefault('c1', 0.0001)
        stepsearch.setdefault('minalpha', 1e-6)
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
        orbitaloptimizer = _helper('orbitaloptimizer', 'variational')

        for name, value in kwargs.items():
            if name not in names:
                raise ValueError("Unknown keyword argument %s" % name)

        #
        # Check dictionaries in keyword arguments
        #
        self.check_keywords_scf(guess, solver, maxiter, dumpci, thresh,
                                printoptions, stepsearch)
        check_options('warning', warning, False, True, 0, 1)
        check_type('checkpoint', checkpoint, int)
        check_type('checkpoint_fn', checkpoint_fn, str)
        check_type('levelshift', lshift, float, int)
        check_options('absolute', pos, False, True, 0, 1)
        check_options('sort', sort, True, False, 0, 1)

        #
        # Set optimization parameters
        #
        if not warning:
            warnings.filterwarnings('ignore')
        if maxiter['orbiter'] < 0:
            raise ValueError('Number of iterations must be greater/equal 0!')

        #
        # Print optimization parameters
        #
        if log.do_medium:
            self.print_options_scf(guess, solver, maxiter, lshift, stepsearch,
                                   thresh, printoptions, checkpoint,
                                   checkpoint_fn, indextrans, orbitaloptimizer,
                                   sort)

        #
        # Generate Guess [geminal, lagrange]
        #
        if guess['geminal'] is None:
            initial_guess = [self.generate_guess(guess)]
        else:
            check_type('guess.geminal', guess['geminal'], np.ndarray)
            if not guess['geminal'].shape[0] == self.dimension:
                raise ValueError('Size of geminal guess array does not match number of unknowns.')
            initial_guess = [guess['geminal']]
        if guess['lagrange'] is None:
            initial_guess.append(self.generate_guess(guess))
        else:
            check_type('guess.lagrange', guess['lagrange'], np.ndarray)
            if not guess['lagrange'].shape[0] == self.dimension:
                raise ValueError('Size of Lagrange guess array does not match number of unknowns.')
            initial_guess.append(guess['lagrange'])

        #
        # Modify starting orbitals if needed:
        #
        self.start_up(orb, swapa, givensrot)

        #
        # Update core energy
        #
        self.update_ecore(core)

        #
        # First iteration:
        #
        self.solve_model(one, two, orb, **{'maxiter': maxiter, 'thresh': thresh,
                         'guess': initial_guess[0], 'guesslm': initial_guess[0],
                         'solver': solver, 'indextrans': indextrans,
                         'orbitaloptimizer': orbitaloptimizer})

        #
        # Update total/corr energy
        #
        etot = self.compute_total_energy()
        ecorr= self.compute_correlation_energy()

        if log.do_medium:
            log.hline('~')
            log('Initial step:')
            self.print_final('Initial')
            log('Entering orbital optimization of AP1roG...')
            log.hline(' ')
            if stepsearch['method']=='trust-region':
                log('%3s %10s %14s  %10s     %10s   %8s   %6s   %10s'
                     %('step', 'Etot', 'D(Etot)', 'Ecorr', 'D(Ecorr)', 'Max(Grad)',
                       '|Grad|', 'TrustRegion'))
            else:
                log('%3s %10s %14s  %10s     %10s   %8s   %6s     %4s' %('step',
                    'Etot', 'D(Etot)', 'Ecorr', 'D(Ecorr)', 'Max(Grad)', '|Grad|',
                    'Step'))

        i = 0

        #
        # Initialize step search
        #
        stepsearch_ = RStepSearch(self.lf, **stepsearch)
        while i < maxiter['orbiter']:
            #
            # Copy energies from previous iteration step
            #
            etot_old = copy(etot)
            ecorr_old = copy(ecorr)

            #
            # Calculate orbital gradient and diagonal approximation to the Hessian
            #
            kappa, gradient, hessian = self.orbital_rotation_step(lshift,
                                                                  pos,
                                                                  orbitaloptimizer)

            #
            # Apply step search to orbital rotation step 'kappa'
            #
            stepsearch_(self, one, two, orb,
                        **{'kappa': kappa, 'thresh': thresh, 'maxiter': maxiter,
                           'gradient': gradient, 'hessian': hessian,
                           'guess': initial_guess[0], 'guesslm': initial_guess[1],
                           'solver': solver, 'indextrans': indextrans,
                           'orbitaloptimizer': orbitaloptimizer})

            #
            # reorder orbitals according to natural occupation numbers
            # works only for varitiational orbital optimization
            #
            if sort and orbitaloptimizer == 'variational':
                self.sort_natural_orbitals(orb)
                # Recalculate WFN:
                # FIXME: this slows us down.
                self.solve_model(one, two, orb, **{'maxiter': maxiter, 'thresh': thresh,
                    'guess': initial_guess[0], 'guesslm': initial_guess[0],
                    'solver': solver, 'indextrans': indextrans,
                    'orbitaloptimizer': orbitaloptimizer})

            etot = self.compute_total_energy()
            ecorr = self.compute_correlation_energy()

            #
            # Print information of iteration step
            #
            if log.do_medium:
                if stepsearch['method']=='trust-region':
                    log('%3i  %14.8f  %11.8f  %10.8f  %11.8f   %6.5f   %5.2e   %1.2e' %(i+1, (etot), (etot-etot_old),
                         ecorr, (ecorr-ecorr_old),
                         gradient.get_max(),
                         gradient.norm(), stepsearch_.trustradius))
                else:
                    log('%3i  %14.8f  %11.8f  %10.8f  %11.8f   %6.5f   %5.2e   %1.2e' %(i+1, (etot), (etot-etot_old), ecorr,
                         (ecorr-ecorr_old),
                         gradient.get_max(),
                         gradient.norm(), stepsearch_.alpha))

            #
            # Checkpoint for orbitals
            #
            if (i+1)%checkpoint == 0 and checkpoint > 0 and (etot < etot_old):
                self.do_checkpoint(orb, olp, checkpoint_fn)

            #
            # Check convergence
            #
            if self.check_convergence(etot, etot_old, gradient, thresh):
                if log.do_medium:
                    log.hline(' ')
                    log('Orbital optimization converged in %i iterations' %(i+1))
                    log.hline(' ')
                break
            elif self.check_stepsearch(stepsearch_):
                if log.do_medium:
                    log.hline(' ')
                    log('Trustradius too small. Orbital optimization aborted!')
                    log.hline(' ')
                break
            i = i+1

        #
        # Check convergence if i = maxorbiter:
        # Don't raise annoying ValueError
        #
        if i >= maxiter['orbiter'] and i>0:
            if not self.check_convergence(etot, etot_old, gradient, thresh):
                if log.do_medium:
                    log.hline(' ')
                    log('WARNING: Orbital optimization NOT converged in %i iterations' %(i))
                    log.hline(' ')

        #
        # Print final information
        #
        if log.do_medium:
            self.print_final()
            self.dump_final(orb, olp, printoptions, dumpci, checkpoint, checkpoint_fn)

        return self.compute_total_energy(), self.geminal, self.lagrange

    def clear(self):
        '''Clear all wavefunction information'''
        self._cache.clear()

    def clear_auxmatrix(self):
        '''Clear auxiliary matrices'''
        self._cache.clear(tags='m', dealloc=True)

    #
    # Density matrices:
    #
    def update_one_dm(self, select, one_dm=None):
        '''Update 1-RDM

           **Arguments:**

           select
                One of ``ps2``, ``response``

           **Optional arguments:**

           one_dm
                When provided, this 1-RDM is stored. A OneIndex instance.
        '''
        cached_one_dm = self.init_one_dm(select)
        if one_dm is None:
            check_options('one_dm', select, 'ps2', 'response')
            if select == 'ps2':
                self.compute_1dm(cached_one_dm, self.geminal, self.geminal, 1, False)
            elif select == 'response':
                self.compute_1dm(cached_one_dm, self.geminal, self.lagrange, 1, True)
        else:
            cached_one_dm.assign(one_dm)

    def update_two_dm(self, select, two_dm=None):
        '''Update 2-RDM

           **Arguments:**

           select
                One of ``ps2``, ``response``

           **Optional arguments:**

           two_dm
                When provided, this 2-RDM is stored.
        '''
        if two_dm is not None:
            raise NotImplementedError
        check_options('two_dm', select, 'ps2', 'response')
        if select == 'response':
            cached_2dm1 = self.init_two_dm('rppqq')
            self.compute_2dm(cached_2dm1, self.one_dm_response, self.geminal, self.lagrange, 'ppqq', True)
            cached_2dm2 = self.init_two_dm('rpqpq')
            self.compute_2dm(cached_2dm2, self.one_dm_response, self.geminal, self.lagrange, 'pqpq', True)
        elif select == 'ps2':
            cached_2dm1 = self.init_two_dm('ppqq')
            self.compute_2dm(cached_2dm1, self.one_dm_ps2, self.geminal, self.geminal, 'ppqq', False)
            cached_2dm2 = self.init_two_dm('pqpq')
            self.compute_2dm(cached_2dm1, self.one_dm_ps2, self.geminal, self.geminal, 'pqpq', False)

    def compute_1dm(self, dmout, mat1, mat2, factor=1.0, response=True):
        '''Compute 1-RDM for AP1roG

           **Arguments:**

           mat1, mat2
                A DenseTwoIndex instance used to calculated 1-RDM. For response
                RDM, mat1 is the geminal matrix, mat2 the Lagrange multipliers

           **Optional arguments:**

           factor
                A scalar factor

           select
                Switch between response (True) and PS2 (False) 1-RDM
        '''
        summand = 1.0
        if not response:
            summand = 1+mat1.contract_two('ab,ab', mat2)

        #
        # Calculate occupied block
        #
        tmpocc = self.lf.create_one_index(self.npairs)
        tmpocc.assign(summand)
        mat1.contract_two_to_one('ab,ab->a', mat2, tmpocc, -1.0, False)

        #
        # Calculate virtual block
        #
        tmpvir = self.lf.create_one_index(self.nvirt)
        mat1.contract_two_to_one('ab,ab->b', mat2, tmpvir, 1.0, True)

        #
        # Combine both blocks and scale
        #
        dmout.assign(tmpocc, 0, self.npairs)
        dmout.assign(tmpvir, self.npairs, self.nbasis)
        dmout.iscale(factor)

    def compute_2dm(self, dmout, dm1, mat1, mat2, select, response=True):
        '''Compute response 2-RDM for AP1roG

           ** Arguments **

           one_dm
               A 1DM. A OneIndex instance.

           mat1, mat2
               TwoIndex instances used to calculate the 2-RDM. To get the
               response DM, mat1 is the geminal coefficient matrix, mat2 are the
               Lagrange multipliers

           select
               Either 'ppqq' or 'pqpq'. Note that the elements (iiii), i.e.,
               the 1DM, are stored in pqpq, while the elements (aaaa) are
               stored in ppqq.

           response
               If True, calculate response 2-RDM. Otherwise the PS2 2-RDM is
               calculated.
        '''
        check_options('select', select, 'ppqq','pqpq')
        lc = mat1.contract_two('ab,ab', mat2)
        factor1 = 1.0
        if response:
            factor1 = factor1-lc
        if select == 'ppqq':
            #
            # temporary storage
            #
            tmpvv = self.lf.create_two_index(self.nvirt, self.nvirt)
            tmpoo = self.lf.create_two_index(self.npairs, self.npairs)
            tmpvv2 = self.lf.create_two_index(self.nvirt, self.nvirt)
            tmpov = self.lf.create_two_index(self.npairs, self.nvirt)
            tmpo = self.lf.create_one_index(self.npairs)
            tmpv = self.lf.create_one_index(self.nvirt)

            #
            # o-o block
            #
            mat2.contract_two_to_two('ab,cb->ca', mat1, tmpoo, 1.0, True)
            tmpoo.assign_diagonal(0.0)
            dmout.iadd(tmpoo, 1.0, 0, self.npairs, 0, self.npairs)
            #
            # v-v block
            #
            mat2.contract_two_to_two('ab,ac->bc', mat1, tmpvv, 1.0, True)
            dmout.iadd(tmpvv, 1.0, self.npairs, self.nbasis, self.npairs, self.nbasis)
            #
            # v-o block
            #
            dmout.iadd_t(mat2, 1.0, self.npairs, self.nbasis, 0, self.npairs)
            #
            # o-v block
            #
            tmpvv2.iadd_tdot(mat2, mat1)
            tmpov.iadd_dot(mat1, tmpvv2)
            mat2.contract_two_to_one('ab,ab->a', mat1, tmpo, 1.0, True)
            tmpov.iadd_contract_two_one('ab,a->ab', mat1, tmpo, -2.0)
            mat2.contract_two_to_one('ab,ab->b', mat1, tmpv, 1.0, True)
            tmpov.iadd_contract_two_one('ab,b->ab', mat1, tmpv, -2.0)
            mat2_ = mat2.copy()
            mat2_.imul(mat1)
            mat2_.imul(mat1)
            tmpov.iadd(mat2_, 2.0)

            dmout.iadd(mat1, (factor1+lc), 0, self.npairs, self.npairs, self.nbasis)
            dmout.iadd(tmpov, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        elif select == 'pqpq':
            #
            # temporary storage
            #
            tmpo = self.lf.create_one_index(self.npairs)
            mat2.contract_two_to_one('ab,ab->a', mat1, tmpo, 1.0, True)
            dm1v = dm1.copy(self.npairs, self.nbasis)
            mat2_ = mat2.copy()
            mat2_.imul(mat1)
            for i in range(self.npairs):
                for j in range(i+1,self.npairs):
                    value = factor1+lc-tmpo.get_element(i)-tmpo.get_element(j)
                    dmout.set_element(i,j, value)
                    dmout.set_element(j,i, value)
                value = factor1+lc-tmpo.get_element(i)
                dmout.set_element(i,i, value)
            dmout.iadd_t(dm1v, 1.0, 0, self.npairs, self.npairs, self.nbasis)
            dmout.iadd(mat2_,-1.0, 0, self.npairs, self.npairs, self.nbasis)
            dmout.iadd(dm1v, 1.0, self.npairs, self.nbasis, 0, self.npairs)
            dmout.iadd_t(mat2_,-1.0, self.npairs, self.nbasis, 0, self.npairs)

    def update_three_dm(self, select, three_dm=None):
        '''Update 3-RDM

           **Optional arguments:**

           three_dm
                When provided, this 3-RDM is stored.
        '''
        raise NotImplementedError

    def update_four_dm(self, select, four_dm=None):
        '''Update 4-RDM

           **Optional arguments:**

           four_dm
                When provided, this 4-RDM is stored.
        '''
        raise NotImplementedError

    #
    # WFN and Lagrange solvers:
    #
    def solve_model(self, one, two, orb, **kwargs):
        '''Solve for geminal model.

           **Arguments:**

           one, two
                One- and two-body integrals (some Hamiltonian matrix elements).
                A TwoIndex and FourIndex/Cholesky instance

           orb
                An expansion instance which contains the MO coefficients.

           **Keywords:**
                guess: initial guess for wfn (1-dim np.array)
                guesslm: initial guess Lagrange multipliers (1-dim np.array)
                solver: wfn/Lagrange solver (dictionary)
                indextrans: 4-index Transformation (str).
                maxiter: maximum number of iterations (dictionary)
                thresh: thresholds (dictionary)
                orbitaloptimizer: orbital optimization method (str)

                For more details, see :py:meth:`RAp1rog.solve_scf`
        '''
        guess = kwargs.get('guess', None)
        guesslm = kwargs.get('guesslm', None)
        solver = kwargs.get('solver', None)
        indextrans = kwargs.get('indextrans', 'tensordot')
        maxiter = kwargs.get('maxiter', None)
        thresh = kwargs.get('thresh', None)
        orbitaloptimizer = kwargs.get('orbitaloptimizer', 'variational')

        #
        # Transform integrals into MO basis
        #
        one_mo, two_mo = transform_integrals(one, two, indextrans, orb)

        #
        # Generate auxiliary matrices needed for optimization
        #
        self.clear_auxmatrix()
        self.update_auxmatrix('scf', two_mo, one_mo)
        del one_mo, two_mo

        #
        # Optimize AP1roG wavefunction amplitudes:
        #
        coeff = self.solve_geminal(guess, solver, thresh['wfn'], maxiter['wfniter'])
        self.update_geminal(coeff)

        #
        # Optimize AP1roG Lagrange multipliers (lambda equations):
        #
        if orbitaloptimizer == 'variational':
            lcoeff = self.solve_lagrange(guesslm, solver, thresh['wfn'], maxiter['wfniter'])
            self.update_lagrange(lcoeff)

            #
            # Update occupation numbers (response 1-RDM)
            #
            self.clear_dm()
            self.update_one_dm('response')
            orb.assign_occupations(self.one_dm_response)


    @timer.with_section('ProjectedSEq')
    def solve_geminal(self, guess, solver, wfnthreshold, wfnmaxiter):
        '''Solves for geminal matrix

           **Arguments:**

           guess
                The initial guess. A 1-dim np array.

           solver
                The solver used. A dictionary with element 'wfn'.

           wfnthreshold
                The optimization threshold. A float.

           wfnmaxiter
                The maximum number of iterations. An integer.
        '''
        iiaa = self.get_auxmatrix('gppqq')
        iaia = self.get_auxmatrix('lpqpq')
        one = self.get_auxmatrix('t')
        fock = self.get_auxmatrix('fock')
        sol = opt.root(self.vector_function_geminal,
                       guess,
                       args=(iiaa, iaia, one, fock),
                       method=solver['wfn'],
                       options={'xtol': wfnthreshold, 'maxiter': wfnmaxiter},
                       callback=None)
        if not sol.success:
            raise ValueError('ERROR: program terminated. Error in solving \
                             geminal equations: %s' %sol.message)
        if log.do_high:
            log('Optimization of geminal coefficients converged in %i \
                 iterations.' %(sol.nit))
        return sol.x

    @timer.with_section('LagrangeMult')
    def solve_lagrange(self, guess, solver, wfnthreshold, wfnmaxiter):
        '''Solves for Lagrange multipliers

           **Arguments:**

           guess
                The initial guess. A 1-dim np array.

           solver
                The solver used. A dictionary with element 'lagrange'.

           wfnthreshold
                The optimization threshold. A float.

           wfnmaxiter
                The maximum number of iterations. An integer.
        '''
        iiaa = self.get_auxmatrix('gppqq')
        iaia = self.get_auxmatrix('lpqpq')
        one = self.get_auxmatrix('t')
        fock = self.get_auxmatrix('fock')
        sol = opt.root(self.vector_function_lagrange,
                       guess,
                       args=(self.geminal, iiaa, iaia, one, fock),
                       method=solver['lagrange'],
                       callback=None,
                       options={'xtol': wfnthreshold, 'maxiter': wfnmaxiter})
        if not sol.success:
            raise ValueError('ERROR: program terminated. Error in solving \
                              Lagrange multipliers: %s' %sol.message)
        if log.do_high:
            log('Optimization of Lagrange multipliers converged in %i \
                     iterations.' %(sol.nit))
        return sol.x

    #
    # Auxiliary matrices
    #
    def init_auxmatrix(self, select):
        '''Initialize auxiliary matrices

           **Arguments:**

           select
                One of ``t``, ``fock``, ``gppqq``, ``gpqpq``, ``gpqqp``,
                ``lpqpq``, ``lpqrq``, ``gpqrr``
        '''
        check_options('select', select, 't', 'fock', 'gppqq', 'gpqpq', 'gpqqp', 'lpqpq', 'lpqrq', 'gpqrr')
        if select in ['lpqrq', 'gpqrr']:
            matrix, new = self._cache.load('matrix_%s' % select, alloc=(self._lf.create_three_index, self.nbasis), tags='m')
        elif select in ['fock']:
            matrix, new = self._cache.load('matrix_%s' % select, alloc=(self._lf.create_one_index, self.nbasis), tags='m')
        else:
            matrix, new = self._cache.load('matrix_%s' % select, alloc=(self._lf.create_two_index, self.nbasis), tags='m')
        if not new:
            raise RuntimeError('The auxiliary matrix matrix_%s already exists. Call clear prior to updating the auxiliary matrices.' % select)
        return matrix

    @timer.with_section('IntegralSort')
    def update_auxmatrix(self, select, two_mo, one_mo):
        '''Derive all auxiliary matrices.
           gppqq:   <pp|qq>,
           gpqpq:   <pq|pq>,
           gpqqp:   <pq|qp>,
           lpqpq:   2<pq|pq>-<pq|qp>,
           lpqrq:   2<pq|rq>-<pq|qr>,
           gpqrr:   <pq|rr>,
           fock:    h_pp + sum_i(2<pi|pi>-<pi|ip>)

           **Arguments:**

           select
                ``scf``.

           two_mo
                two-electron integrals in the MO basis. A FourIndex instance.

           one_mo
                one-electron integrals in the MO basis. A TwoIndex instance.
        '''
        if select != 'scf':
            raise RuntimeError('The auxiliary matrix matrix_%s is not supported.' % select)
        matrix = self.init_auxmatrix('t')
        matrix.assign(one_mo[0])
        #
        # <pp|qq>
        #
        matrix1 = self.init_auxmatrix('gppqq')
        two_mo[0].slice_to_two('aabb->ab', matrix1)
        #
        # <pq|pq>
        #
        matrix2 = self.init_auxmatrix('gpqpq')
        two_mo[0].slice_to_two('abab->ab', matrix2)
        #
        # <pq|qp>
        #
        matrix3 = self.init_auxmatrix('gpqqp')
        two_mo[0].slice_to_two('abba->ab', matrix3)
        #
        # <pq||pq>+<pq|pq>
        #
        matrix4 = self.init_auxmatrix('lpqpq')
        two_mo[0].slice_to_two('abab->ab', matrix4, 2.0, True)
        two_mo[0].slice_to_two('abba->ab', matrix4,-1.0, False)
        #
        # <pq||rq>+<pq|rq>
        #
        matrix5 = self.init_auxmatrix('lpqrq')
        two_mo[0].slice_to_three('abcb->abc', matrix5, 2.0, True)
        two_mo[0].slice_to_three('abbc->abc', matrix5,-1.0, False)
        #
        # <pq|rr>
        #
        matrix6 = self.init_auxmatrix('gpqrr')
        two_mo[0].slice_to_three('abcc->abc', matrix6)
        #
        # inactive diagonal Fock = h_pp + sum_i (<pi||pi>+<pi|pi>)
        #
        matrix7 = self.init_auxmatrix('fock')
        tmp1 = self.lf.create_one_index()
        matrix.copy_diagonal(tmp1)
        matrix4.contract_to_one('ab->a', matrix7, 1.0, True, 0, self.nbasis, 0, self.npairs)
        matrix7.iadd(tmp1)

    def get_auxmatrix(self, select):
        '''Get an auxiliary matrix from cache.

           **Arguments:**

           select
                't', 'gpqpq', 'gpqqp', 'lpqpq', 'pqrq', 'gpqrr', 'fock'.
        '''
        if not 'matrix_%s' % select in self._cache:
            raise ValueError("The auxmatrix %s not found in cache. Did you use init_auxmatrix?" %select)
        return self._cache.load('matrix_%s' % select)

    #
    # Functions for energy evaluation:
    #
    def compute_correlation_energy(self, arg=None):
        '''Get correlation energy of restricted AP1roG

           **Optional arguments:**

           arg
                The geminal coefficient matrix (np.array or TwoIndex instance).
                If not provided, the correlation energy is calculated
                from self.geminal (default None)
        '''
        if arg is None:
            coeff = self.geminal
        else:
            coeff = self.lf.create_two_index(self.npairs, self.nvirt)
            coeff.assign(arg)
        kmat = self.get_auxmatrix('gppqq')
        #
        # Ecorr = sum_ia c_ia <ii|aa>
        #
        return coeff.contract_two('ab,ab', kmat, 0, self.npairs, self.npairs, self.nbasis)

    def compute_reference_energy(self):
        '''Get energy of reference determinant for restricted AP1roG'''
        one = self.get_auxmatrix('t')
        fock = self.get_auxmatrix('fock')
        #
        # Eref = sum_i (t_ii + F_ii)
        #
        energy = one.trace(0, self.npairs, 0, self.npairs)
        energy += fock.trace(0, self.npairs)
        return energy

    def compute_total_energy(self, coeff=None):
        '''Get total energy (reference+correlation) including nuclear-repulsion/
           core energy for restricted AP1roG.

           **Optional arguments:**

           arg
                The geminal coefficient matrix (np.array or TwoIndex instance).
                If not provided, the correlation energy is calculated
                from self.geminal (default None)
        '''
        return (self.compute_correlation_energy(coeff)+self.compute_reference_energy()+self.ecore)

    #
    # Functions for orbital optimization (except gradient/Hessian):
    #
    def compute_rotation_matrix(self, coeff):
        '''Determine orbital rotation matrix for (oo), (vo), and (vv) blocks

           **Arguments:**

           coeff
                The nonreduntant orbital rotations k_pq (1-dim np.array). The
                elements are sorted w.r.t. lower triangular indices (p>q).
        '''
        indl = np.tril_indices(self.nbasis, -1)
        kappa = self.lf.create_two_index(self.nbasis, self.nbasis)
        kappa.assign(coeff, indl)
        #
        # k_pq = -k_qp
        #
        kappa.iadd_t(kappa, -1.0)

        out = compute_unitary_matrix(kappa)
        return out

    def orbital_rotation_step(self, lshift=1e-8, pos=True,
                              optimizer='variational'):
        '''Get orbital rotation step (Newton--Raphson)

           **Arguments:**


           **Optional arguments:**

           lshift
               Level shift (float) for approximate Hessian added to small elements
               (default 1e-8)

           pos
               (boolean) Make all elements of Hessian positive if set to True
               (default True)

           optimizer
               Orbital otimization method (str) (default 'variational')

        '''
        check_options('orbitaloptimizer', optimizer, 'variational', 'ps2c')
        #
        # Switch between different orbital optimization schemes
        #
        ps2c = (optimizer == 'ps2c')
        #
        # Calculate orbital gradient and diagonal approximation to the Hessian
        #
        grad = self.compute_orbital_gradient(ps2c)
        #
        # We use a diagonal Hessian
        #
        hessian = self.compute_orbital_hessian(lshift, pos, ps2c)
        #
        # Orbital rotation step
        #
        kappa = grad.divide(hessian, -1.0)

        return kappa, grad, hessian

    #
    # Vector function for AP1roG:
    #
    @timer.with_section('VecFctGeminal')
    def vector_function_geminal(self, coeff, miiaa, miaia, one, diagfock):
        '''Construct vector function for optimization of geminal coefficients.

           **Arguments:**

           coeff
                The geminal coefficients (np.array)

           miiaa, miaia
                Auxiliary matrices (TwoIndex instances)

           one
                1-electron integrals (TwoIndex instance)

           diagfock
                Diagonal inactive Fock matrix (OneIndex instances)
        '''
        gmat = self.lf.create_two_index(self.npairs, self.nvirt)
        gmat.assign(coeff)

        #
        # correlation energy for current iteration step
        #
        ecorr = self.compute_correlation_energy(gmat)

        #
        # vectorFunction_ia
        #
        result = self.lf.create_two_index(self.npairs,self.nvirt)

        #
        # Add contributions to vectorFunction_ia:
        #
        # c_ia*ecorr
        #
        result.iadd(gmat, -ecorr)

        #
        # c_0*miiaa
        #
        result.iadd_slice(miiaa, 1.0, 0, self.npairs, self.npairs, self.nbasis)

        #
        # -2c_ia*f_ii
        #
        result.iadd_contract_two_one('ab,a->ab', gmat, diagfock, -2.0, begin2=0, end2=self.npairs)

        #
        # 2c_ia*f_aa
        #
        result.iadd_contract_two_one('ab,b->ab', gmat, diagfock, 2.0, begin2=self.npairs, end2=self.nbasis)

        #
        # -2c_ia*(<ia||ia>+<ia|ia>)
        #
        result.iadd_mult(gmat, miaia, -2.0, 0, self.npairs, self.npairs, self.nbasis)

        #
        # c_ja*(<ii|jj>)
        #
        result.iadd_dot(miiaa, gmat, 1.0, 0, self.npairs, 0, self.npairs)

        #
        # c_ib*(<bb|aa>)
        #
        result.iadd_dot(gmat, miiaa, 1.0, begin2=self.npairs, end2=self.nbasis, begin3=self.npairs, end3=self.nbasis)

        #
        # c_ia*(c_jb*<jj|bb>)
        #
        factor = gmat.contract_two('ab,ab', miiaa, 0, self.npairs, self.npairs, self.nbasis)
        result.iadd(gmat, factor)

        #
        # c_ib*<bb|jj>*c_ja
        #
        tmp = self.lf.create_two_index(self.npairs,self.npairs)
        tmp.iadd_dot(gmat, miiaa, 1.0, begin2=self.npairs, end2=self.nbasis, begin3=0, end3=self.npairs)
        result.iadd_dot(tmp, gmat, 1.0)

        #
        # -2c_ja*c_ia*<jj|aa>
        #
        tmp = gmat.contract_two_to_one('ab,ab->b', miiaa, None, 1.0, True, 0, self.npairs, self.npairs, self.nbasis)
        result.iadd_contract_two_one('ab,b->ab', gmat, tmp, -2.0)

        #
        # -2c_ib*c_ia*<ii|bb>
        #
        tmp = gmat.contract_two_to_one('ab,ab->a', miiaa, None, 1.0, True, 0, self.npairs, self.npairs, self.nbasis)
        result.iadd_contract_two_one('ab,a->ab', gmat, tmp, -2.0)

        #
        # +2c_ia*c_ia*<ii|aa>
        #
        tmp = gmat.copy()
        tmp.imul(gmat)
        result.iadd_mult(tmp, miiaa, 2.0, 0, self.npairs, self.npairs, self.nbasis)

        return result._array.ravel(order='C')

    #
    # Jacobian for AP1roG:
    #
    def jacobian_ap1rog(self, coeff, miiaa, miaia, one, fock):
        '''Construct Jacobian for optimization of geminal coefficients.
        '''
        raise NotImplementedError

    #
    # Vector function for Lagrange multipliers:
    #
    @timer.with_section('VecFctLagrange')
    def vector_function_lagrange(self, lagrange, gmat, miiaa, miaia, one, diagfock):
        '''Construct vector function for optimization of Lagrange multipliers.

           **Arguments:**

           lagrange
                The lagrange multipliers (np array)

           gmat
                The geminal coefficients (TwoIndex instance)

           miiaa, miaia
                Auxiliary matrices (TwoIndex instances)

           one
                1-electron integrals (TwoIndex instance)

           diagfock
                Diagonal inactive Fock matrix (OneIndex instance)
        '''
        lmat = self.lf.create_two_index(self.npairs, self.nvirt)
        lmat.assign(lagrange)

        #
        # Correlation energy
        #
        ecorr = self.compute_correlation_energy(gmat)
        #
        # intermediate variables:
        #  * lc = sum_ia c_ia l_ia
        #  * ciiaa = sum_ia c_ia <ii|aa>
        #  * diagint = <pp|pp>
        #
        lc = lmat.contract_two('ab,ab', gmat)
        ciiaa = gmat.contract_two('ab,ab', miiaa, 0, self.npairs, self.npairs, self.nbasis)
        diagint = miiaa.copy_diagonal()
        #
        # cgi = sum_b c_ib <ii|bb>
        # cga = sum_j c_ja <jj|aa>
        # lci = sum_b c_ib l_ib
        # lca = sum_j c_ja l_ja
        #
        cgi = self.lf.create_one_index(self.npairs)
        cga = self.lf.create_one_index(self.nvirt)
        lci = self.lf.create_one_index(self.npairs)
        lca = self.lf.create_one_index(self.nvirt)
        gmat.contract_two_to_one('ab,ab->a', miiaa, cgi, 1.0, True, 0, self.npairs, self.npairs, self.nbasis)
        gmat.contract_two_to_one('ab,ab->b', miiaa, cga, 1.0, True, 0, self.npairs, self.npairs, self.nbasis)
        lmat.contract_two_to_one('ab,ab->a', gmat, lci, 1.0, True)
        lmat.contract_two_to_one('ab,ab->b', gmat, lca, 1.0, True)

        #
        # vectorFunction_ia
        #
        result = self.lf.create_two_index(self.npairs,self.nvirt)

        #
        # miiaa
        #
        result.iadd_slice(miiaa, 1.0, 0, self.npairs, self.npairs, self.nbasis)

        #
        # -2l_ia*cgi
        #
        result.iadd_contract_two_one('ab,a->ab', lmat, cgi, -2.0)

        #
        # -2l_ia*cga
        #
        result.iadd_contract_two_one('ab,b->ab', lmat, cga, -2.0)

        #
        # -2miiaa*lci
        #
        result.iadd_contract_two_one('ab,a->ab', miiaa, lci, -2.0, 0, self.npairs, self.npairs, self.nbasis)

        #
        # -2miiaa*lca
        #
        result.iadd_contract_two_one('ab,b->ab', miiaa, lca, -2.0, 0, self.npairs, self.npairs, self.nbasis)

        #
        # -2l_ia*f_ii
        #
        result.iadd_contract_two_one('ab,a->ab', lmat, diagfock, -2.0, begin2=0, end2=self.npairs)

        #
        # 2l_ia*f_aa
        #
        result.iadd_contract_two_one('ab,b->ab', lmat, diagfock, 2.0, begin2=self.npairs, end2=self.nbasis)

        #
        # -2l_ia*(<ia||ia>+<ia|ia>)
        #
        result.iadd_mult(lmat, miaia, -2.0, 0, self.npairs, self.npairs, self.nbasis)

        #
        # l_ja*(<ii|jj>)
        #
        result.iadd_dot(miiaa, lmat, 1.0, 0, self.npairs, 0, self.npairs)

        #
        # l_ib*(<bb|aa>)
        #
        result.iadd_dot(lmat, miiaa, 1.0, begin2=self.npairs, end2=self.nbasis, begin3=self.npairs, end3=self.nbasis)

        #
        # 4l_ia*c_ia*<ii|aa>
        #
        tmp = lmat.copy()
        tmp.imul(gmat)
        result.iadd_mult(tmp, miiaa, 4.0, begin0=0, end0=self.npairs, begin1=self.npairs, end1=self.nbasis)

        #
        # <ii|bb>*c_jb*l_ja
        #
        tmp = self.lf.create_two_index(self.nvirt, self.nvirt)
        tmp.iadd_tdot(gmat, lmat)
        result.iadd_dot(miiaa, tmp, 1.0, 0, self.npairs, self.npairs, self.nbasis)

        #
        # l_ib*c_jb*<jj|aa>
        #
        tmp = self.lf.create_two_index(self.npairs, self.npairs)
        tmp.iadd_dott(lmat, gmat)
        result.iadd_dot(tmp, miiaa, 1.0, begin2=0, end2=self.npairs, begin3=self.npairs, end3=self.nbasis)

        return result._array.ravel(order='C')

    #
    # Jacobian for Lagrange multipliers of OAP1roG:
    #
    def jacobian_lambda(self, lagrange, gmat, miiaa, miaia, one, diagfock):
        '''Construct Jacobian for optimization of Lagrange multipliers for
           restricted AP1roG.
        '''
        raise NotImplementedError

    #
    # Calculate orbital gradient for OAP1roG (oo,vo,vv):
    #
    @timer.with_section('OOAP1roG grad')
    def compute_orbital_gradient(self, ps2c=False):
        '''Determine orbital gradient for all non-reduntant orbital rotations
           (oo,vo,vv).

           **Optional arguments:**

           ps2c
                (boolean) If True, switches to PS2c orbital optimization
                (default False)
        '''
        one = self.get_auxmatrix('t')
        vpqrq = self.get_auxmatrix('lpqrq')
        vpqrr = self.get_auxmatrix('gpqrr')

        #
        # Get 1- and 2-RDMs
        #
        self.clear_dm()
        if ps2c:
            self.update_one_dm('ps2')
            self.update_two_dm('ps2')
            onedm = self.one_dm_ps2
            twodmpqpq = self.two_dm_pqpq
            twodmppqq = self.two_dm_ppqq
        else:
            self.update_one_dm('response')
            self.update_two_dm('response')
            onedm = self.one_dm_response
            twodmpqpq = self.two_dm_rpqpq
            twodmppqq = self.two_dm_rppqq

        #
        # Orbital gradient g_ut
        #
        gradient = self.lf.create_two_index(self.nbasis, self.nbasis)

        #
        # Symmetrize G(amma)_ppqq 2-RDM. This reduces the number of operations
        #
        two_dm_av = twodmppqq.copy()
        two_dm_av.symmetrize()
        #
        # L_uptp*G_up
        #
        vpqrq.contract_two_to_two('abc,ab->ac', twodmpqpq, gradient, 4.0, False)
        #
        # <tu|pp>*G_up
        #
        vpqrr.contract_two_to_two('abc,bc->ba', two_dm_av, gradient, 4.0, False)
        #
        # -L_tpup*G_tp
        #
        vpqrq.contract_two_to_two('abc,ab->ca', twodmpqpq, gradient, -4.0, False)
        #
        # -<ut|pp>*G_tp
        #
        vpqrr.contract_two_to_two('abc,bc->ab', two_dm_av, gradient, -4.0, False)
        #
        # h_ut g_tt
        #
        gradient.iadd_contract_two_one('ab,a->ab', one, onedm, 4.0)
        #
        # h_ut g_uu
        #
        gradient.iadd_contract_two_one('ab,b->ab', one, onedm,-4.0)

        ind = np.tril_indices(self.nbasis, -1)

        #
        # return only lower triangle
        #
        return gradient.copy_slice(ind)

    #
    # Approximate diagonal Hessian for OAP1roG:
    #
    def compute_orbital_hessian(self, lshift=1e-8, pos=False, ps2c=False):
        '''Construct diagonal approximation to Hessian for orbital optimization
           of AP1roG.

           **Optional arguments:**

           lshift
                Level shift (float) (default 1e-8)

           pos
                Set all elements positive (boolean) (default False)

           ps2c
                (boolean) If True, switches to PS2c orbital optimization
                (default False)
        '''
        #
        # Get auxiliary matrices
        #
        one = self.get_auxmatrix('t')
        vpq = self.get_auxmatrix('lpqpq')
        vpp = self.get_auxmatrix('gppqq')
        vmat = self.get_auxmatrix('gpqpq')
        vmata = self.get_auxmatrix('gpqqp')

        #
        # Get 1- and 2-RDMs
        #
        self.clear_dm()
        if ps2c:
            self.update_one_dm('ps2')
            self.update_two_dm('ps2')
            onedm = self.one_dm_ps2
            twodmpqpq = self.two_dm_pqpq
            twodmppqq = self.two_dm_ppqq
        else:
            self.update_one_dm('response')
            self.update_two_dm('response')
            onedm = self.one_dm_response
            twodmpqpq = self.two_dm_rpqpq
            twodmppqq = self.two_dm_rppqq

        #
        # Symmetrize G(amma)_ppqq 2-RDM
        #
        two_dm_av = twodmppqq.copy()
        two_dm_av.symmetrize()
        #
        # Modify diagonal elements to simplify contractions
        #
        two_dm_av.assign_diagonal(0.0)
        twodmpqpq.assign_diagonal(0.0)

        #
        # Calculate additional auxiliary matrices
        #
        vmatdiag = vmat.copy_diagonal()
        one_diag = one.copy_diagonal()
        two_c = vpq.contract_two_to_one('ab,ab->a', twodmpqpq, None, 1.0)
        two_ca = vpp.contract_two_to_one('ab,ab->a', two_dm_av, None, 1.0)

        #
        # Diagonal orbital Hessian hessian_pq = hessian_(pq,pq)
        #
        hessian = self.lf.create_two_index(self.nbasis, self.nbasis)
        #
        # <qt> G_pt
        #
        hessian.iadd_dot(vpq, twodmpqpq, 4.0)
        #
        # <pt> G_qt
        #
        hessian.iadd_dot(twodmpqpq, vpq, 4.0)
        #
        # <qt> G_qt
        #
        hessian.iadd(two_c, -4.0)
        #
        # <pt> G_pt
        #
        hessian.iadd_t(two_c, -4.0)
        #
        # <qt> G_pt
        #
        hessian.iadd_dot(two_dm_av, vpp, 4.0)
        #
        # <pt> G_qt
        #
        hessian.iadd_dot(vpp, two_dm_av, 4.0)
        #
        # <qt> G_qt
        #
        hessian.iadd(two_ca, -4.0)
        #
        # <pt> G_pt
        #
        hessian.iadd_t(two_ca, -4.0)
        #
        # <pq> G_pq
        #
        hessian.iadd_mult(vmat, twodmpqpq, 8.0)
        hessian.iadd_mult(vmata, twodmpqpq,-8.0)
        hessian.iadd_mult(vpp, twodmpqpq,-16.0)
        hessian.iadd_mult(vmat, two_dm_av, -8.0)
        hessian.iadd_mult(vpp, two_dm_av, -8.0)
        #
        # <ppqq> G_pp
        #
        hessian.iadd_contract_two_one('ab,a->ab', vpp, onedm, 8.0)
        hessian.iadd_contract_two_one('ab,b->ab', vpp, onedm, 8.0)
        #
        # <qq> g_pp
        #
        hessian.iadd_one_mult(one_diag, onedm, 4.0, True, False)
        #
        # <pp> g_qq
        #
        hessian.iadd_one_mult(onedm, one_diag, 4.0, True, False)
        #
        # <qq> g_qq
        #
        hessian.iadd_one_mult(one_diag, onedm, -4.0, True, True)
        #
        # <pp> g_pp
        #
        hessian.iadd_one_mult(onedm, one_diag, -4.0, False, False)
        hessian.iadd_contract_two_one('ab,a->ab', vmat, onedm, 4.0)
        hessian.iadd_contract_two_one('ab,b->ab', vmat, onedm, 4.0)
        hessian.iadd_one_mult(vmatdiag, onedm, -4.0, True, True)
        hessian.iadd_one_mult(vmatdiag, onedm, -4.0, False, False)

        #
        # Make everything positive
        #
        if pos:
            hessian.iabs()
        #
        # Add levelshift:
        #
        if lshift:
            hessian.iadd_shift(lshift)

        ind = np.tril_indices(self.nbasis, -1)

        return hessian.copy_slice(ind)

    @timer.with_section('exact Hessian')
    def get_exact_hessian(self, mo1, mo2):
        '''Construct exact Hessian for orbital optimization of restricted OAP1roG.

           The usage of the exact orbital Hessian for the orbital optimization
           is currently not supported. The exact Hessian can only be evaluated
           a posteriori.

           **Arguements**

           mo1, mo2
                The 1- and 2-el integrals in the MO basis (TwoIndex and
                FourIndex instances)
        '''
        #
        # exact orbital Hessian output hessian_pq,rs
        #
        hessian = self.lf.create_four_index()
        #
        # calculate <pq|rs>-<pq|sr>
        #
        mo2ex = mo2.copy()
        mo2ex.iadd_exchange()

        #
        # Get response DMs
        #
        self.clear_dm()
        self.update_one_dm('response')
        self.update_two_dm('response')
        dm1 = self.one_dm_response
        dm2pqpq = self.two_dm_rpqpq
        dm2ppqq = self.two_dm_rppqq
        #
        # Symmetrize 2DM
        #
        dm2av = dm2ppqq.copy()
        dm2av.symmetrize()
        #
        # Reset diagonal elements of DMs
        #
        dm2pqpqex.assign_diagonal(0.0)
        dm2pqpq.assign_diagonal(dm1)
        dm2av.assign_diagonal(0.0)

        #
        # temporary storage
        #
        ind2 = self.lf.create_two_index()
        ind30 = self.lf.create_three_index()
        ind31 = self.lf.create_three_index()

        # (1)
        # aa
        mo2ex.contract_two_to_four('abcd,cd->acbd', dm2pqpqex, hessian, 4.0, True)
        mo2ex.contract_two_to_four('abcd,ab->acdb', dm2pqpqex, hessian,-4.0, False)
        mo2ex.contract_two_to_four('abcd,cd->acdb', dm2pqpqex, hessian,-4.0, False)
        mo2ex.contract_two_to_four('abcd,ab->acbd', dm2pqpqex, hessian, 4.0, False)
        # ab
        mo2.contract_two_to_four('abcd,cd->acbd', dm2pqpq, hessian, 4.0, False)
        mo2.contract_two_to_four('abcd,ab->acdb', dm2pqpq, hessian,-4.0, False)
        mo2.contract_two_to_four('abcd,cd->acdb', dm2pqpq, hessian,-4.0, False)
        mo2.contract_two_to_four('abcd,ab->acbd', dm2pqpq, hessian, 4.0, False)
        # (2)
        # aa
        mo2ex.contract_two_to_four('abcd,cb->acdb', dm2pqpqex, hessian, 4.0, False)
        mo2ex.contract_two_to_four('abcd,ad->acbd', dm2pqpqex, hessian,-4.0, False)
        mo2ex.contract_two_to_four('abcd,cb->acbd', dm2pqpqex, hessian,-4.0, False)
        mo2ex.contract_two_to_four('abcd,ad->acdb', dm2pqpqex, hessian, 4.0, False)
        # ab
        mo2.contract_two_to_four('abcd,cb->acdb', dm2pqpq, hessian, 4.0, False)
        mo2.contract_two_to_four('abcd,ad->acbd', dm2pqpq, hessian,-4.0, False)
        mo2.contract_two_to_four('abcd,cb->acbd', dm2pqpq, hessian,-4.0, False)
        mo2.contract_two_to_four('abcd,ad->acdb', dm2pqpq, hessian, 4.0, False)
        # (3)
        mo2.contract_two_to_four('abcd,bd->abcd', dm2av, hessian, 4.0, False)
        mo2.contract_two_to_four('abcd,ad->abcd', dm2av, hessian,-4.0, False)
        mo2.contract_two_to_four('abcd,bd->abdc', dm2av, hessian,-4.0, False)
        mo2.contract_two_to_four('abcd,ac->abcd', dm2av, hessian, 4.0, False)
        mo2.contract_two_to_four('abcd,bc->abdc', dm2av, hessian, 4.0, False)
        mo2.contract_two_to_four('abcd,ac->abdc', dm2av, hessian,-4.0, False)
        mo2.contract_two_to_four('abcd,bc->abcd', dm2av, hessian,-4.0, False)
        mo2.contract_two_to_four('abcd,ad->abdc', dm2av, hessian, 4.0, False)
        # Apq,qw (pw) (qv) (-qw) (-pv)
        ind2.iadd_contract_two_one('ab,b->ab', mo1, dm1, 2.0)
        ind2.iadd_contract_two_one('ab,a->ab', mo1, dm1, 2.0)
        # aa
        mo2ex.contract_two_to_two('abcb,cb->ac', dm2pqpqex, ind2, 2.0, False)
        mo2ex.contract_two_to_two('abcb,ab->ac', dm2pqpqex, ind2, 2.0, False)
        # ab
        mo2.contract_two_to_two('abcb,cb->ac', dm2pqpq, ind2, 2.0, False)
        mo2.contract_two_to_two('abcb,ab->ac', dm2pqpq, ind2, 2.0, False)
        mo2.contract_two_to_two('abcc,bc->ba', dm2av, ind2, 2.0, False)
        mo2.contract_two_to_two('abcc,bc->ab', dm2av, ind2, 2.0, False)
        # Apq,qw (pq,qw) (-pq,vq)
        ind30.iadd_expand_two_one('ab,c->cab', mo1, dm1, -4)
        # aa
        mo2ex.contract_two_to_three('abcb,db->adc', dm2pqpqex, ind31,-4.0, True)
        # ab
        mo2.contract_two_to_three('abcb,db->adc', dm2pqpq, ind31,-4.0, False)
        mo2.contract_two_to_three('abcc,dc->adb', dm2av, ind31,-4.0, False)
        # Apq,qw (pq,vp) (-pq,pw)
        ind31.iadd_expand_two_one('ab,c->acb', mo1, dm1, -4)
        # aa
        mo2.contract_two_to_three('abcb,db->dac', dm2pqpqex, ind30,-4.0, False)
        mo2.contract_two_to_three('abbc,db->dac', dm2pqpqex, ind30, 4.0, False)
        # ab
        mo2.contract_two_to_three('abcb,db->dac', dm2pqpq, ind30,-4.0, False)
        mo2.contract_two_to_three('abcc,dc->dab', dm2av, ind30,-4.0, False)

        #
        # collect 2- and 3-index terms
        #
        hessian.iadd_expand_two_to_four('0-2', ind2,-1.0)
        hessian.iadd_expand_two_to_four('0-3', ind2, 1.0)
        hessian.iadd_expand_two_to_four('1-2', ind2, 1.0)
        hessian.iadd_expand_two_to_four('1-3', ind2,-1.0)
        hessian.iadd_expand_three_to_four('1-3-1-2', ind30,-1.0)
        hessian.iadd_expand_three_to_four('1-2-1-2', ind30, 1.0)
        hessian.iadd_expand_three_to_four('0-3-0-2', ind31, 1.0)
        hessian.iadd_expand_three_to_four('0-2-0-2', ind31,-1.0)

        #
        # reorder elements (p,q,r,s)->(pq,rs)
        #
        dim = (self.nbasis*(self.nbasis-1))/2
        tril = np.tril_indices(self.nbasis, -1)
        out = self.lf.create_two_index(dim, dim)
        # FIXME: use matrix class
        out._array[:] = (hessian._array[:,:,tril[0],tril[1]])[tril]

        return out._array

    def sort_natural_orbitals(self, orb):
        '''Sort orbitals w.r.t. the natural occupation numbers

           **Arguments:**

           orb
                The AO/MO coefficients (natural orbitals) and the natural
                occupation numbers (Expansion instance)
        '''
        #
        # Take natural occupation numbers from orb.occupations
        #
        onedm = self.lf.create_one_index()
        onedm.assign(orb.occupations)
        order = onedm.sort_indices()
        orderref = np.arange(self.nbasis)
        if not (order == orderref).all():
            orb.permute_orbitals(order)

    def print_solution(self, cithresh=0.01, excitationlevel=2,
                       amplitudestofile=False, filename="./ap1rog_amplitudes"):
        '''Print coefficients :math:`{ci}` and Slater Determinant if :math:`|ci|` > threshold.
           Prints up to hextuply excited pairs.

           **Optional arguments:**

           cithresh
                Upper threshold (a float) for CI coefficients to be
                reconstructed. (Default 0.01)

           excitationslevel
                The maximum number of substitutions w.r.t. the reference
                determinant for which the CI coefficients are reconstructed.
                (Default 2)

           amplitudestofile
                A boolean. If True, the CI coefficients are stored in a separate
                file. (Default False)

           filename
                The file name for the wfn amplitudes.
                (Default ap1rog_amplitudes)
        '''
        it = []
        coeff = self.geminal._array.copy()
        matrix = np.identity(self.npairs)
        for i in range(excitationlevel):
            it.append(np.nditer(coeff,flags=['multi_index']))
        if amplitudestofile:
            with open(filename, 'w') as filea:
                filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant((-1,-1)), 1.0))
                filea.write("\n")
        else:
            log('%s %20.16f' %(self.get_slater_determinant((-1,-1)), 1.0))
        for i in range(excitationlevel):
            if i==0:
                for ci in it[0]:
                    if math.fabs(ci) >= (cithresh):
                        if amplitudestofile is True:
                            with open(filename, 'a') as filea:
                                filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index), float(ci)))
                                filea.write("\n")
                        else:
                            log('%s %20.16f' %(self.get_slater_determinant(it[0].multi_index), float(ci)))
            if self.npairs > 1:
                for index in range((i+1)):
                    it[index] = np.nditer(coeff,flags=['multi_index'])
                if i==1:
                    for ci in it[0]:
                        if math.fabs(ci) >= (cithresh/excitationlevel):
                            it[1] = it[0].copy()
                            for ci2 in it[1]:
                                if (it[1].multi_index[0]>it[0].multi_index[0]):
                                    if (it[1].multi_index[1]>it[0].multi_index[1]):
                                        matrix[it[0].multi_index[0],it[0].multi_index[0]] = float(ci)
                                        matrix[it[1].multi_index[0],it[1].multi_index[0]] = float(ci2)
                                        matrix[it[0].multi_index[0],it[1].multi_index[0]] = coeff[it[0].multi_index[0],it[1].multi_index[1]]
                                        matrix[it[1].multi_index[0],it[0].multi_index[0]] = coeff[it[1].multi_index[0],it[0].multi_index[1]]
                                        amplitude = self.perm(matrix)
                                        if (math.fabs(amplitude) >= cithresh):
                                            if amplitudestofile is True:
                                                with open(filename, 'a') as filea:
                                                    filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index), amplitude))
                                                    filea.write("\n")
                                            else:
                                                log('%s %20.16f' %(self.get_slater_determinant(it[0].multi_index,it[1].multi_index), amplitude))
                                        matrix = np.identity(self.npairs)
                if i==2:
                    for ci in it[0]:
                        if math.fabs(ci) >= (cithresh/excitationlevel):
                            it[1] = it[0].copy()
                            for ci2 in it[1]:
                                if (it[1].multi_index[0]>it[0].multi_index[0]):
                                    if (it[1].multi_index[1]>it[0].multi_index[1]):
                                        it[2] = it[1].copy()
                                        for ci3 in it[2]:
                                            if (it[2].multi_index[0]>it[1].multi_index[0]):
                                                if (it[2].multi_index[1]>it[1].multi_index[1]):
                                                    matrix[it[0].multi_index[0],it[0].multi_index[0]] = float(ci)
                                                    matrix[it[1].multi_index[0],it[1].multi_index[0]] = float(ci2)
                                                    matrix[it[2].multi_index[0],it[2].multi_index[0]] = float(ci3)

                                                    matrix[it[0].multi_index[0],it[1].multi_index[0]] = coeff[it[0].multi_index[0],it[1].multi_index[1]]
                                                    matrix[it[0].multi_index[0],it[2].multi_index[0]] = coeff[it[0].multi_index[0],it[2].multi_index[1]]

                                                    matrix[it[1].multi_index[0],it[0].multi_index[0]] = coeff[it[1].multi_index[0],it[0].multi_index[1]]
                                                    matrix[it[1].multi_index[0],it[2].multi_index[0]] = coeff[it[1].multi_index[0],it[2].multi_index[1]]

                                                    matrix[it[2].multi_index[0],it[0].multi_index[0]] = coeff[it[2].multi_index[0],it[0].multi_index[1]]
                                                    matrix[it[2].multi_index[0],it[1].multi_index[0]] = coeff[it[2].multi_index[0],it[1].multi_index[1]]
                                                    amplitude = self.perm(matrix)
                                                    if (math.fabs(amplitude) >= cithresh):
                                                        if amplitudestofile is True:
                                                            with open(filename, 'a') as filea:
                                                                filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                                   it[2].multi_index), amplitude))
                                                                filea.write("\n")
                                                        else:
                                                            log('%s %20.16f' %(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                           it[2].multi_index), amplitude))
                                                    matrix = np.identity(self.npairs)
                if i==3:
                    for ci in it[0]:
                        if math.fabs(ci) >= (cithresh/excitationlevel):
                            it[1] = it[0].copy()
                            for ci2 in it[1]:
                                if (it[1].multi_index[0]>it[0].multi_index[0]):
                                    if (it[1].multi_index[1]>it[0].multi_index[1]):
                                        it[2] = it[1].copy()
                                        for ci3 in it[2]:
                                            if (it[2].multi_index[0]>it[1].multi_index[0]):
                                                if (it[2].multi_index[1]>it[1].multi_index[1]):
                                                    it[3] = it[2].copy()
                                                    for ci4 in it[3]:
                                                        if (it[3].multi_index[0]>it[2].multi_index[0]):
                                                            if (it[3].multi_index[1]>it[2].multi_index[1]):
                                                                matrix[it[0].multi_index[0],it[0].multi_index[0]] = float(ci)
                                                                matrix[it[1].multi_index[0],it[1].multi_index[0]] = float(ci2)
                                                                matrix[it[2].multi_index[0],it[2].multi_index[0]] = float(ci3)
                                                                matrix[it[3].multi_index[0],it[3].multi_index[0]] = float(ci4)

                                                                matrix[it[0].multi_index[0],it[1].multi_index[0]] = coeff[it[0].multi_index[0],it[1].multi_index[1]]
                                                                matrix[it[0].multi_index[0],it[2].multi_index[0]] = coeff[it[0].multi_index[0],it[2].multi_index[1]]
                                                                matrix[it[0].multi_index[0],it[3].multi_index[0]] = coeff[it[0].multi_index[0],it[3].multi_index[1]]

                                                                matrix[it[1].multi_index[0],it[0].multi_index[0]] = coeff[it[1].multi_index[0],it[0].multi_index[1]]
                                                                matrix[it[1].multi_index[0],it[2].multi_index[0]] = coeff[it[1].multi_index[0],it[2].multi_index[1]]
                                                                matrix[it[1].multi_index[0],it[3].multi_index[0]] = coeff[it[1].multi_index[0],it[3].multi_index[1]]

                                                                matrix[it[2].multi_index[0],it[0].multi_index[0]] = coeff[it[2].multi_index[0],it[0].multi_index[1]]
                                                                matrix[it[2].multi_index[0],it[1].multi_index[0]] = coeff[it[2].multi_index[0],it[1].multi_index[1]]
                                                                matrix[it[2].multi_index[0],it[3].multi_index[0]] = coeff[it[2].multi_index[0],it[3].multi_index[1]]

                                                                matrix[it[3].multi_index[0],it[0].multi_index[0]] = coeff[it[3].multi_index[0],it[0].multi_index[1]]
                                                                matrix[it[3].multi_index[0],it[1].multi_index[0]] = coeff[it[3].multi_index[0],it[1].multi_index[1]]
                                                                matrix[it[3].multi_index[0],it[2].multi_index[0]] = coeff[it[3].multi_index[0],it[2].multi_index[1]]
                                                                amplitude = self.perm(matrix)
                                                                if (math.fabs(amplitude) >= cithresh):
                                                                    if amplitudestofile is True:
                                                                        with open(filename, 'a') as filea:
                                                                            filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                        it[2].multi_index,it[3].multi_index), amplitude))
                                                                            filea.write("\n")
                                                                    else:
                                                                        log('%s %20.16f' %(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                           it[2].multi_index,it[3].multi_index), amplitude))
                                                                matrix = np.identity(self.npairs)
                if i==4:
                    for ci in it[0]:
                        if math.fabs(ci) >= (cithresh/excitationlevel):
                            it[1] = it[0].copy()
                            for ci2 in it[1]:
                                if ((it[1].multi_index[0]>it[0].multi_index[0])
                                    and (it[1].multi_index[1]>it[0].multi_index[1])):
                                    it[2] = it[1].copy()
                                    for ci3 in it[2]:
                                        if ((it[2].multi_index[0]>it[1].multi_index[0])
                                            and (it[2].multi_index[1]>it[1].multi_index[1])):
                                            it[3] = it[2].copy()
                                            for ci4 in it[3]:
                                                if ((it[3].multi_index[0]>it[2].multi_index[0])
                                                    and (it[3].multi_index[1]>it[2].multi_index[1])):
                                                    it[4] = it[3].copy()
                                                    for ci5 in it[4]:
                                                        if ((it[4].multi_index[0]>it[3].multi_index[0])
                                                            and (it[4].multi_index[1]>it[3].multi_index[1])):
                                                            matrix[it[0].multi_index[0],it[0].multi_index[0]] = float(ci)
                                                            matrix[it[1].multi_index[0],it[1].multi_index[0]] = float(ci2)
                                                            matrix[it[2].multi_index[0],it[2].multi_index[0]] = float(ci3)
                                                            matrix[it[3].multi_index[0],it[3].multi_index[0]] = float(ci4)
                                                            matrix[it[4].multi_index[0],it[4].multi_index[0]] = float(ci5)

                                                            matrix[it[0].multi_index[0],it[1].multi_index[0]] = coeff[it[0].multi_index[0],it[1].multi_index[1]]
                                                            matrix[it[0].multi_index[0],it[2].multi_index[0]] = coeff[it[0].multi_index[0],it[2].multi_index[1]]
                                                            matrix[it[0].multi_index[0],it[3].multi_index[0]] = coeff[it[0].multi_index[0],it[3].multi_index[1]]
                                                            matrix[it[0].multi_index[0],it[4].multi_index[0]] = coeff[it[0].multi_index[0],it[4].multi_index[1]]

                                                            matrix[it[1].multi_index[0],it[0].multi_index[0]] = coeff[it[1].multi_index[0],it[0].multi_index[1]]
                                                            matrix[it[1].multi_index[0],it[2].multi_index[0]] = coeff[it[1].multi_index[0],it[2].multi_index[1]]
                                                            matrix[it[1].multi_index[0],it[3].multi_index[0]] = coeff[it[1].multi_index[0],it[3].multi_index[1]]
                                                            matrix[it[1].multi_index[0],it[4].multi_index[0]] = coeff[it[1].multi_index[0],it[4].multi_index[1]]

                                                            matrix[it[2].multi_index[0],it[0].multi_index[0]] = coeff[it[2].multi_index[0],it[0].multi_index[1]]
                                                            matrix[it[2].multi_index[0],it[1].multi_index[0]] = coeff[it[2].multi_index[0],it[1].multi_index[1]]
                                                            matrix[it[2].multi_index[0],it[3].multi_index[0]] = coeff[it[2].multi_index[0],it[3].multi_index[1]]
                                                            matrix[it[2].multi_index[0],it[4].multi_index[0]] = coeff[it[2].multi_index[0],it[4].multi_index[1]]

                                                            matrix[it[3].multi_index[0],it[0].multi_index[0]] = coeff[it[3].multi_index[0],it[0].multi_index[1]]
                                                            matrix[it[3].multi_index[0],it[1].multi_index[0]] = coeff[it[3].multi_index[0],it[1].multi_index[1]]
                                                            matrix[it[3].multi_index[0],it[2].multi_index[0]] = coeff[it[3].multi_index[0],it[2].multi_index[1]]
                                                            matrix[it[3].multi_index[0],it[4].multi_index[0]] = coeff[it[3].multi_index[0],it[4].multi_index[1]]

                                                            matrix[it[4].multi_index[0],it[0].multi_index[0]] = coeff[it[4].multi_index[0],it[0].multi_index[1]]
                                                            matrix[it[4].multi_index[0],it[1].multi_index[0]] = coeff[it[4].multi_index[0],it[1].multi_index[1]]
                                                            matrix[it[4].multi_index[0],it[2].multi_index[0]] = coeff[it[4].multi_index[0],it[2].multi_index[1]]
                                                            matrix[it[4].multi_index[0],it[3].multi_index[0]] = coeff[it[4].multi_index[0],it[3].multi_index[1]]
                                                            amplitude = self.perm(matrix)
                                                            if (math.fabs(amplitude) >= cithresh):
                                                                if amplitudestofile is True:
                                                                    with open(filename, 'a') as filea:
                                                                        filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                    it[2].multi_index,it[3].multi_index, it[4].multi_index), amplitude))
                                                                        filea.write("\n")
                                                                else:
                                                                    log('%s %20.16f' %(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                       it[2].multi_index,it[3].multi_index, it[4].multi_index), amplitude))
                                                            matrix = np.identity(self.npairs)
                if i==5:
                    for ci in it[0]:
                        if math.fabs(ci) >= (cithresh/excitationlevel):
                            it[1] = it[0].copy()
                            for ci2 in it[1]:
                                if ((it[1].multi_index[0]>it[0].multi_index[0])
                                    and (it[1].multi_index[1]>it[0].multi_index[1])):
                                    it[2] = it[1].copy()
                                    for ci3 in it[2]:
                                        if ((it[2].multi_index[0]>it[1].multi_index[0])
                                            and (it[2].multi_index[1]>it[1].multi_index[1])):
                                            it[3] = it[2].copy()
                                            for ci4 in it[3]:
                                                if ((it[3].multi_index[0]>it[2].multi_index[0])
                                                    and (it[3].multi_index[1]>it[2].multi_index[1])):
                                                    it[4] = it[3].copy()
                                                    for ci5 in it[4]:
                                                        if ((it[4].multi_index[0]>it[3].multi_index[0])
                                                            and (it[4].multi_index[1]>it[3].multi_index[1])):
                                                            it[5] = it[4].copy()
                                                            for ci6 in it[5]:
                                                                if ((it[5].multi_index[0]>it[4].multi_index[0])
                                                                    and (it[5].multi_index[1]>it[4].multi_index[1])):
                                                                    matrix[it[0].multi_index[0],it[0].multi_index[0]] = float(ci)
                                                                    matrix[it[1].multi_index[0],it[1].multi_index[0]] = float(ci2)
                                                                    matrix[it[2].multi_index[0],it[2].multi_index[0]] = float(ci3)
                                                                    matrix[it[3].multi_index[0],it[3].multi_index[0]] = float(ci4)
                                                                    matrix[it[4].multi_index[0],it[4].multi_index[0]] = float(ci5)
                                                                    matrix[it[5].multi_index[0],it[5].multi_index[0]] = float(ci6)

                                                                    matrix[it[0].multi_index[0],it[1].multi_index[0]] = coeff[it[0].multi_index[0],it[1].multi_index[1]]
                                                                    matrix[it[0].multi_index[0],it[2].multi_index[0]] = coeff[it[0].multi_index[0],it[2].multi_index[1]]
                                                                    matrix[it[0].multi_index[0],it[3].multi_index[0]] = coeff[it[0].multi_index[0],it[3].multi_index[1]]
                                                                    matrix[it[0].multi_index[0],it[4].multi_index[0]] = coeff[it[0].multi_index[0],it[4].multi_index[1]]
                                                                    matrix[it[0].multi_index[0],it[5].multi_index[0]] = coeff[it[0].multi_index[0],it[5].multi_index[1]]

                                                                    matrix[it[1].multi_index[0],it[0].multi_index[0]] = coeff[it[1].multi_index[0],it[0].multi_index[1]]
                                                                    matrix[it[1].multi_index[0],it[2].multi_index[0]] = coeff[it[1].multi_index[0],it[2].multi_index[1]]
                                                                    matrix[it[1].multi_index[0],it[3].multi_index[0]] = coeff[it[1].multi_index[0],it[3].multi_index[1]]
                                                                    matrix[it[1].multi_index[0],it[4].multi_index[0]] = coeff[it[1].multi_index[0],it[4].multi_index[1]]
                                                                    matrix[it[1].multi_index[0],it[5].multi_index[0]] = coeff[it[1].multi_index[0],it[5].multi_index[1]]

                                                                    matrix[it[2].multi_index[0],it[0].multi_index[0]] = coeff[it[2].multi_index[0],it[0].multi_index[1]]
                                                                    matrix[it[2].multi_index[0],it[1].multi_index[0]] = coeff[it[2].multi_index[0],it[1].multi_index[1]]
                                                                    matrix[it[2].multi_index[0],it[3].multi_index[0]] = coeff[it[2].multi_index[0],it[3].multi_index[1]]
                                                                    matrix[it[2].multi_index[0],it[4].multi_index[0]] = coeff[it[2].multi_index[0],it[4].multi_index[1]]
                                                                    matrix[it[2].multi_index[0],it[5].multi_index[0]] = coeff[it[2].multi_index[0],it[5].multi_index[1]]

                                                                    matrix[it[3].multi_index[0],it[0].multi_index[0]] = coeff[it[3].multi_index[0],it[0].multi_index[1]]
                                                                    matrix[it[3].multi_index[0],it[1].multi_index[0]] = coeff[it[3].multi_index[0],it[1].multi_index[1]]
                                                                    matrix[it[3].multi_index[0],it[2].multi_index[0]] = coeff[it[3].multi_index[0],it[2].multi_index[1]]
                                                                    matrix[it[3].multi_index[0],it[4].multi_index[0]] = coeff[it[3].multi_index[0],it[4].multi_index[1]]
                                                                    matrix[it[3].multi_index[0],it[5].multi_index[0]] = coeff[it[3].multi_index[0],it[5].multi_index[1]]

                                                                    matrix[it[4].multi_index[0],it[0].multi_index[0]] = coeff[it[4].multi_index[0],it[0].multi_index[1]]
                                                                    matrix[it[4].multi_index[0],it[1].multi_index[0]] = coeff[it[4].multi_index[0],it[1].multi_index[1]]
                                                                    matrix[it[4].multi_index[0],it[2].multi_index[0]] = coeff[it[4].multi_index[0],it[2].multi_index[1]]
                                                                    matrix[it[4].multi_index[0],it[3].multi_index[0]] = coeff[it[4].multi_index[0],it[3].multi_index[1]]
                                                                    matrix[it[4].multi_index[0],it[5].multi_index[0]] = coeff[it[4].multi_index[0],it[5].multi_index[1]]

                                                                    matrix[it[5].multi_index[0],it[0].multi_index[0]] = coeff[it[5].multi_index[0],it[0].multi_index[1]]
                                                                    matrix[it[5].multi_index[0],it[1].multi_index[0]] = coeff[it[5].multi_index[0],it[1].multi_index[1]]
                                                                    matrix[it[5].multi_index[0],it[2].multi_index[0]] = coeff[it[5].multi_index[0],it[2].multi_index[1]]
                                                                    matrix[it[5].multi_index[0],it[3].multi_index[0]] = coeff[it[5].multi_index[0],it[3].multi_index[1]]
                                                                    matrix[it[5].multi_index[0],it[4].multi_index[0]] = coeff[it[5].multi_index[0],it[4].multi_index[1]]

                                                                    amplitude = self.perm(matrix)
                                                                    if (math.fabs(amplitude) >= cithresh):
                                                                        if amplitudestofile is True:
                                                                            with open(filename, 'a') as filea:
                                                                                filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                            it[2].multi_index,it[3].multi_index, it[4].multi_index,it[5].multi_index),
                                                                                            amplitude))
                                                                                filea.write("\n")
                                                                        else:
                                                                            log('%s %20.16f' %(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                               it[2].multi_index,it[3].multi_index, it[4].multi_index,it[5].multi_index),
                                                                                               amplitude))
                                                                    matrix = np.identity(self.npairs)
                if i==6:
                    for ci in it[0]:
                        if math.fabs(ci) >= (cithresh/excitationlevel):
                            it[1] = it[0].copy()
                            for ci2 in it[1]:
                                if ((it[1].multi_index[0]>it[0].multi_index[0])
                                    and (it[1].multi_index[1]>it[0].multi_index[1])):
                                    it[2] = it[1].copy()
                                    for ci3 in it[2]:
                                        if ((it[2].multi_index[0]>it[1].multi_index[0])
                                            and (it[2].multi_index[1]>it[1].multi_index[1])):
                                            it[3] = it[2].copy()
                                            for ci4 in it[3]:
                                                if ((it[3].multi_index[0]>it[2].multi_index[0])
                                                    and (it[3].multi_index[1]>it[2].multi_index[1])):
                                                    it[4] = it[3].copy()
                                                    for ci5 in it[4]:
                                                        if ((it[4].multi_index[0]>it[3].multi_index[0])
                                                            and (it[4].multi_index[1]>it[3].multi_index[1])):
                                                            it[5] = it[4].copy()
                                                            for ci6 in it[5]:
                                                                if ((it[5].multi_index[0]>it[4].multi_index[0])
                                                                    and (it[5].multi_index[1]>it[4].multi_index[1])):
                                                                    it[6] = it[5].copy()
                                                                    for ci7 in it[6]:
                                                                        if ((it[6].multi_index[0]>it[5].multi_index[0])
                                                                            and (it[6].multi_index[1]>it[5].multi_index[1])):
                                                                            matrix[it[0].multi_index[0],it[0].multi_index[0]] = float(ci)
                                                                            matrix[it[1].multi_index[0],it[1].multi_index[0]] = float(ci2)
                                                                            matrix[it[2].multi_index[0],it[2].multi_index[0]] = float(ci3)
                                                                            matrix[it[3].multi_index[0],it[3].multi_index[0]] = float(ci4)
                                                                            matrix[it[4].multi_index[0],it[4].multi_index[0]] = float(ci5)
                                                                            matrix[it[5].multi_index[0],it[5].multi_index[0]] = float(ci6)
                                                                            matrix[it[6].multi_index[0],it[6].multi_index[0]] = float(ci7)

                                                                            matrix[it[0].multi_index[0],it[1].multi_index[0]] = coeff[it[0].multi_index[0],it[1].multi_index[1]]
                                                                            matrix[it[0].multi_index[0],it[2].multi_index[0]] = coeff[it[0].multi_index[0],it[2].multi_index[1]]
                                                                            matrix[it[0].multi_index[0],it[3].multi_index[0]] = coeff[it[0].multi_index[0],it[3].multi_index[1]]
                                                                            matrix[it[0].multi_index[0],it[4].multi_index[0]] = coeff[it[0].multi_index[0],it[4].multi_index[1]]
                                                                            matrix[it[0].multi_index[0],it[5].multi_index[0]] = coeff[it[0].multi_index[0],it[5].multi_index[1]]
                                                                            matrix[it[0].multi_index[0],it[6].multi_index[0]] = coeff[it[0].multi_index[0],it[6].multi_index[1]]

                                                                            matrix[it[1].multi_index[0],it[0].multi_index[0]] = coeff[it[1].multi_index[0],it[0].multi_index[1]]
                                                                            matrix[it[1].multi_index[0],it[2].multi_index[0]] = coeff[it[1].multi_index[0],it[2].multi_index[1]]
                                                                            matrix[it[1].multi_index[0],it[3].multi_index[0]] = coeff[it[1].multi_index[0],it[3].multi_index[1]]
                                                                            matrix[it[1].multi_index[0],it[4].multi_index[0]] = coeff[it[1].multi_index[0],it[4].multi_index[1]]
                                                                            matrix[it[1].multi_index[0],it[5].multi_index[0]] = coeff[it[1].multi_index[0],it[5].multi_index[1]]
                                                                            matrix[it[1].multi_index[0],it[6].multi_index[0]] = coeff[it[1].multi_index[0],it[6].multi_index[1]]

                                                                            matrix[it[2].multi_index[0],it[0].multi_index[0]] = coeff[it[2].multi_index[0],it[0].multi_index[1]]
                                                                            matrix[it[2].multi_index[0],it[1].multi_index[0]] = coeff[it[2].multi_index[0],it[1].multi_index[1]]
                                                                            matrix[it[2].multi_index[0],it[3].multi_index[0]] = coeff[it[2].multi_index[0],it[3].multi_index[1]]
                                                                            matrix[it[2].multi_index[0],it[4].multi_index[0]] = coeff[it[2].multi_index[0],it[4].multi_index[1]]
                                                                            matrix[it[2].multi_index[0],it[5].multi_index[0]] = coeff[it[2].multi_index[0],it[5].multi_index[1]]
                                                                            matrix[it[2].multi_index[0],it[6].multi_index[0]] = coeff[it[2].multi_index[0],it[6].multi_index[1]]

                                                                            matrix[it[3].multi_index[0],it[0].multi_index[0]] = coeff[it[3].multi_index[0],it[0].multi_index[1]]
                                                                            matrix[it[3].multi_index[0],it[1].multi_index[0]] = coeff[it[3].multi_index[0],it[1].multi_index[1]]
                                                                            matrix[it[3].multi_index[0],it[2].multi_index[0]] = coeff[it[3].multi_index[0],it[2].multi_index[1]]
                                                                            matrix[it[3].multi_index[0],it[4].multi_index[0]] = coeff[it[3].multi_index[0],it[4].multi_index[1]]
                                                                            matrix[it[3].multi_index[0],it[5].multi_index[0]] = coeff[it[3].multi_index[0],it[5].multi_index[1]]
                                                                            matrix[it[3].multi_index[0],it[6].multi_index[0]] = coeff[it[3].multi_index[0],it[6].multi_index[1]]

                                                                            matrix[it[4].multi_index[0],it[0].multi_index[0]] = coeff[it[4].multi_index[0],it[0].multi_index[1]]
                                                                            matrix[it[4].multi_index[0],it[1].multi_index[0]] = coeff[it[4].multi_index[0],it[1].multi_index[1]]
                                                                            matrix[it[4].multi_index[0],it[2].multi_index[0]] = coeff[it[4].multi_index[0],it[2].multi_index[1]]
                                                                            matrix[it[4].multi_index[0],it[3].multi_index[0]] = coeff[it[4].multi_index[0],it[3].multi_index[1]]
                                                                            matrix[it[4].multi_index[0],it[5].multi_index[0]] = coeff[it[4].multi_index[0],it[5].multi_index[1]]
                                                                            matrix[it[4].multi_index[0],it[6].multi_index[0]] = coeff[it[4].multi_index[0],it[6].multi_index[1]]

                                                                            matrix[it[5].multi_index[0],it[0].multi_index[0]] = coeff[it[5].multi_index[0],it[0].multi_index[1]]
                                                                            matrix[it[5].multi_index[0],it[1].multi_index[0]] = coeff[it[5].multi_index[0],it[1].multi_index[1]]
                                                                            matrix[it[5].multi_index[0],it[2].multi_index[0]] = coeff[it[5].multi_index[0],it[2].multi_index[1]]
                                                                            matrix[it[5].multi_index[0],it[3].multi_index[0]] = coeff[it[5].multi_index[0],it[3].multi_index[1]]
                                                                            matrix[it[5].multi_index[0],it[4].multi_index[0]] = coeff[it[5].multi_index[0],it[4].multi_index[1]]
                                                                            matrix[it[5].multi_index[0],it[6].multi_index[0]] = coeff[it[5].multi_index[0],it[6].multi_index[1]]

                                                                            matrix[it[6].multi_index[0],it[0].multi_index[0]] = coeff[it[6].multi_index[0],it[0].multi_index[1]]
                                                                            matrix[it[6].multi_index[0],it[1].multi_index[0]] = coeff[it[6].multi_index[0],it[1].multi_index[1]]
                                                                            matrix[it[6].multi_index[0],it[2].multi_index[0]] = coeff[it[6].multi_index[0],it[2].multi_index[1]]
                                                                            matrix[it[6].multi_index[0],it[3].multi_index[0]] = coeff[it[6].multi_index[0],it[3].multi_index[1]]
                                                                            matrix[it[6].multi_index[0],it[4].multi_index[0]] = coeff[it[6].multi_index[0],it[4].multi_index[1]]
                                                                            matrix[it[6].multi_index[0],it[5].multi_index[0]] = coeff[it[6].multi_index[0],it[5].multi_index[1]]
                                                                            amplitude = self.perm(matrix)
                                                                            if (math.fabs(amplitude) >= cithresh):
                                                                                if amplitudestofile is True:
                                                                                    with open(filename, 'a') as filea:
                                                                                        filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                    it[2].multi_index,it[3].multi_index, it[4].multi_index,it[5].multi_index,
                                                                                                    it[6].multi_index), amplitude))
                                                                                        filea.write("\n")
                                                                                else:
                                                                                    log('%s %20.16f' %(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                       it[2].multi_index,it[3].multi_index,
                                                                                                       it[4].multi_index,it[5].multi_index,
                                                                                                       it[6].multi_index), amplitude))
                                                                            matrix = np.identity(self.npairs)

    def get_slater_determinant(self, *indices):
        '''Return excited Slater Determinant.

           **Arguments:**

           indices
                A (list of) multi-index. First element contains occupied index
                w.r.t. reference determinant, second element contains virtual
                index w.r.t. reference determinant.
        '''
        orb_ref = []
        excited = []
        for ind in indices:
            orb_ref.append(ind[0])
            excited.append(ind[1]+self.npairs)

        sd = []
        for i in range(self.nbasis):
            if i in orb_ref:
                sd.append(0)
            elif i in excited:
                sd.append(2)
            else:
                if i < self.npairs:
                    sd.append(2)
                else:
                    sd.append(0)

        return str(sd).translate(None, ", ")

    def do_checkpoint(self, orb, olp, checkpoint_fn):
        '''Dump orbitals for restart

           **Arguments:**

           orb
                An expansion instance. AO/MO coefficients stored to disk

           olp
                A TwoIndex instance. The AO overlap matrix. Required for restart
                at different molecular geometry.

           checkpoint_fn
                The filename of the checkpoint file.
        '''
        # Use the horton.io interface to write the checkpoint file
        from horton.io.iodata import IOData
        mol = IOData(olp=olp, exp_alpha=orb)
        mol.to_file(checkpoint_fn)

    def print_final(self, s='Final'):
        '''Print energies

           **Optional arguments:**

           s
                A string.
        '''
        log.hline('-')
        log('%s reference Energy:     %16.12f' %(s,self.compute_reference_energy()+self.ecore))
        log('%s correlation Energy:   %16.12f' %(s,self.compute_correlation_energy()))
        log.hline('=')
        log('%s total Energy:         %16.12f' %(s, (self.compute_total_energy())))
        log.hline('=')

    def dump_final(self, orb, olp, printoptions, dumpci, checkpoint, checkpoint_fn='checkpoint.h5'):
        '''Dump final solution (orbitals, wavefunction amplitudes, geminal
           coefficients)

           **Arguments:**

           orb
                An expansion instance. AO/MO coefficients stored to disk

           olp
                A TwoIndex instance. The AO overlap matrix. Required for restart

           printoptions, dumpci
                A dictionary. See :py:meth:`RAp1rog.solve_scf`

           checkpoint:
                An integer. See :py:meth:`RAp1rog.solve_scf`

           checkpoint_fn
                The filename of the checkpoint file.
        '''
        if checkpoint > 0:
            if log.do_medium:
                log('Writing orbitals to file')
            self.do_checkpoint(orb, olp, checkpoint_fn)
            if log.do_medium:
                log.hline('-')
                log('Final solution for coefficients: (only printed if |c_i| > %f)' %printoptions['ci'])
        self.print_solution(printoptions['ci'],
                            printoptions['excitationlevel'],
                            amplitudestofile=dumpci['amplitudestofile'],
                            filename=dumpci['amplitudesfilename'])
        if printoptions['geminal']:
            log.hline('-')
            log('Geminal matrix: (rows=occupied, columns=virtual)')
            self.print_geminal_coefficients()
        if log.do_medium:
            log.hline('=')

    def print_geminal_coefficients(self):
        '''Print geminal coefficients
        '''
        log.hline('')
        s = ''.join(str('%11i' %(i+1+self.npairs)) for i in range(self.nvirt))
        log("  "+s)
        log.hline('')
        for line in range(self.npairs):
            s = str('%2i :     ' %(line+1))
            s2 = ''
            for row in range(self.nvirt):
                s2 += ' '+str('%10.6f' %(self.geminal.get_element(line, row)))
            log(s+s2)

    def check_keywords(self, guess, solver, maxiter, dumpci, thresh,
                       printoptions):
        '''Check dictionaries if they contain proper keys.

           **Arguments:**

           guess, solver, maxiter, dumpci, thresh, printoptions
                See :py:meth:`RAp1rog.solve`

        '''
        #
        # Check guess, values are checked separately
        #
        for key in guess:
            check_options('guess', key, 'type', 'factor', 'geminal')
        #
        # Check solver
        #
        for key in solver:
            check_options('solver', key, 'wfn')
        #
        # Check maxiter
        #
        for key, value in maxiter.items():
            check_options('maxiter', key, 'wfniter')
            check_type('maxiter', value, int)
        #
        # Check thresh
        #
        for key, value in thresh.items():
            check_options('thresh', key, 'wfn')
            check_type('thresh', value, float)
            if value < 0:
                raise ValueError('Negative convergence threshold for %s is not allowed!' %key)
        #
        # Check printoptions
        #
        for key, value in printoptions.items():
            check_options('printoptions', key, 'geminal', 'ci', 'excitationlevel')
            if key == 'geminal':
                check_options('printoptions.geminal', value, False, True, 0, 1)
            elif key == 'ci':
                check_type('printoptions.ci', value, float)
            elif key == 'excitationlevel':
                check_type('printoptions.excitationlevel', value, int)
        #
        # Check dumpci
        #
        for key, value in dumpci.items():
            check_options('dumpci', key, 'amplitudestofile', 'amplitudesfilename')
            if key == 'amplitudestofile':
                check_options('dumpci.amplitudestofile', value, False, True, 0, 1)
            if key == 'amplitudesfilename':
                check_type('dumpci.amplitudesfilename', value, str)

    def check_keywords_scf(self, guess, solver, maxiter, dumpci, thresh,
                           printoptions, stepsearch):
        '''Check dictionaries if they contain proper keys.

           **Arguments:**

           guess, solver, maxiter, dumpci, thresh, printoptions, stepseach
                See :py:meth:`RAp1rog.solve_scf`

        '''
        #
        # Check guess, values are checked separately
        #
        for key in guess:
            check_options('guess', key, 'type', 'factor', 'geminal', 'lagrange')
        #
        # Check solver
        #
        for key in solver:
            check_options('solver', key, 'wfn', 'lagrange')
        #
        # Check maxiter
        #
        for key, value in maxiter.items():
            check_options('maxiter', key, 'wfniter', 'orbiter')
            check_type('maxiter', value, int)
        #
        # Check thresh
        #
        for key, value in thresh.items():
            check_options('thresh', key, 'wfn', 'energy', 'gradientnorm',
                          'gradientmax')
            check_type('thresh', value, float)
            if value < 0:
                raise ValueError('Negative convergence threshold for %s is not allowed!' %key)
        #
        # Check printoptions
        #
        for key, value in printoptions.items():
            check_options('printoptions', key, 'geminal', 'ci', 'excitationlevel')
            if key == 'geminal':
                check_options('printoptions.geminal', value, False, True, 0, 1)
            elif key == 'ci':
                check_type('printoptions.ci', value, float)
            elif key == 'excitationlevel':
                check_type('printoptions.excitationlevel', value, int)
        #
        # Check dumpci
        #
        for key, value in dumpci.items():
            check_options('dumpci', key, 'amplitudestofile', 'amplitudesfilename')
            if key == 'amplitudestofile':
                check_options('dumpci.amplitudestofile', value, False, True, 0, 1)
            if key == 'amplitudesfilename':
                check_type('dumpci.amplitudesfilename', value, str)
        #
        # Check stepsearch, values are checked separately
        #
        for key in stepsearch:
            check_options('stepsearch', key, 'method', 'optimizer', 'alpha',
                'c1', 'minalpha', 'maxiterouter', 'maxiterinner', 'maxeta',
                'mineta', 'upscale', 'downscale', 'trustradius', 'maxtrustradius',
                'threshold')

    def print_options(self, guess, solver, maxiter, thresh,
                      printoptions, indextrans):
        '''Print optimization options.

           **Arguments:**

           guess, solver, maxiter, thresh, printoptions, indextrans
                See :py:meth:`RAp1rog.solve`

        '''
        log.hline()
        log(' ')
        log('Entering AP1roG optimization (restricted, closed-shell):')
        log(' ')
        log.hline()
        log('OPTIMIZATION PARAMETERS:')
        log('Number of pairs:             %i' %self.npairs)
        log('Number of virtuals:          %i' %self.nvirt)
        log('4-index transformation:      %s' %indextrans)
        log('Initial guess:')
        log('  type:                        %s' %guess['type'])
        log('  scaling factor:              %2.3f' %guess['factor'])
        log('Solvers:')
        log('  wavefunction amplitudes:   %s' %solver['wfn'])
        log('Optimization thresholds:')
        log('  wavefunction:              %1.2e' %thresh['wfn'])
        log('Printing options:')
        log('  c_i:                       %1.2e' %printoptions['ci'])
        log('  excitation level:          %i' %printoptions['excitationlevel'])

    def print_options_scf(self, guess, solver, maxiter, lshift, stepsearch,
                          thresh, printoptions, checkpoint, checkpoint_fn,
                          indextrans, orbitaloptimizer, sort):
        '''Print optimization options.

           **Arguments:**

           See :py:meth:`RAp1rog.solve_scf`

        '''
        if log.do_medium:
            log.hline()
            log(' ')
            log('Entering AP1roG optimization including orbital optimization (restricted, closed-shell):')
            log(' ')
            log.hline()
            log('OPTIMIZATION PARAMETERS:')
            log('Number of pairs:               %i' %self.npairs)
            log('Number of virtuals:            %i' %self.nvirt)
            log('Initial guess:')
            log('  type:                        %s' %guess['type'])
            log('  scaling factor:              %2.3f' %guess['factor'])
            log('Solvers:')
            log('  wavefunction amplitudes:     %s' %solver['wfn'])
            log('  Lagrange multipliers:        %s' %solver['lagrange'])
            log('  orbital optimization:        %s' %orbitaloptimizer)
            log('Number of iterations:          %i' %maxiter['orbiter'])
            log('Checkpointing:                 %i' %checkpoint)
            log('Checkpoint file:               %s' %checkpoint_fn)
            log('4-index transformation:        %s' %indextrans)
            log('Level shift:                   %3.3e' %lshift)
            log('Sorting natural orbitals:      %s' %sort)
            if stepsearch['method']=='trust-region':
                log('Apply trust region:')
                log('  initial trust radius:        %1.3f' %stepsearch['trustradius'])
                log('  maximum turst radius:        %2.3f' %stepsearch['maxtrustradius'])
                log('  upper trust ratio bound:     %1.2e' %stepsearch['maxeta'])
                log('  lower trust ratio bound:     %1.2e' %stepsearch['mineta'])
                log('  upper scaling factor:        %1.2e' %stepsearch['upscale'])
                log('  lower scaling factor:        %1.2e' %stepsearch['downscale'])
                log('  max number of iterations:    %i' %stepsearch['maxiterouter'])
                log('  max number of optimizations: %i' %stepsearch['maxiterinner'])
                log('  optimization threshold:      %1.2e' %stepsearch['threshold'])
                log('  optimizer:                   %s' %stepsearch['optimizer'])
            elif stepsearch['method']=='backtracking':
                log('Apply line search:')
                log('  line search method:          %s' %stepsearch['method'])
                log('  initial scaling factor:      %1.3f' %stepsearch['alpha'])
                log('  contraction factor:          %1.3f' %stepsearch['downscale'])
                log('  c1 factor:                   %1.3e' %stepsearch['c1'])
                log('  minimum scaling factor:      %1.3e' %stepsearch['minalpha'])
            else:
                log('No step search selected:')
                log('  scaling factor:              %1.3f' %stepsearch['alpha'])
            log('Optimization thresholds:')
            log('  wavefunction:                %1.2e' %thresh['wfn'])
            log('  energy:                      %1.2e' %thresh['energy'])
            log('  gradient:                    %1.2e' %thresh['gradientmax'])
            log('  gradient norm:               %1.2e' %thresh['gradientnorm'])
            log('Printing options:')
            log('  c_i:                         %1.2e' %printoptions['ci'])
            log('  excitation level:            %i' %printoptions['excitationlevel'])

    def start_up(self, orb, swapa, givensrot=np.array([[]])):
        '''Modify orbitals prior to optimization

           **Arguments:**

           swapa
                An integer numpy array with two columns where every row
                corresponds to one swap.

           **Optional arguments:**

           givensrot
               Givens rotation of orbital pairs. A numpy array where every
               row corresponds to one Givens rotation with elements
               (index1, index2, angle in deg).

        '''
        if swapa.any():
            if log.do_medium:
                log.hline('~')
                log('Swap orbitals:')
            orb.swap_orbitals(swapa)
        if givensrot.any():
            if log.do_medium:
                log.hline('~')
                log('Apply Givens rotation:')
            if not (givensrot.shape[1] == 3 and givensrot.ndim == 2 and issubclass(givensrot.dtype.type, int)):
                raise TypeError('The argument givensrot has the wrong shape/type.')
            for irot in range(len(givensrot)):
                index0, index1, angle = givensrot[irot]
                if log.do_medium:
                    log('Rotating orbitals %i and %i' %(index0, index1))
                orb.rotate_2orbitals(angle, index0, index1)

    def compute_objective_function(self, coeff=None):
        '''Objective function for line search optimization

           **Optional arguments:**

           coeff
                See :py:meth:`RAp1rog.compute_total_energy`
        '''
        return self.compute_total_energy(coeff)
