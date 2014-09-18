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
from scipy import optimize as opt
import warnings

from horton.log import log, timer
from horton.correlatedwfn.geminal import Geminal
from horton.correlatedwfn.stepsearch import RStepSearch
from horton.matrix import TwoIndex

from itertools import permutations
from operator import mul
from math import fsum
from copy import deepcopy
# from spermutations import spermutations

__all__ = [
    'RAp1rog',
]



class RAp1rog(Geminal):
    '''Restricted AP1roG wavefunction class'''

    @timer.with_section('AP1roG')
    def __call__(self, one, two, core, exps, olp, orb, **kwargs):
        '''Find Geminal expansion coefficient and---if required---find
           optimal set of orbitals.
           For restricted, closed-shell AP1roG.

           **Arguments:**

           one, two
                One- and four-index integrals (some Hamiltonian matrix elements)
                expressed in the AO basis (primitives).

           core
                The core energy (not included in 'one' and 'two').

           exps
                An expansion instance. It contains the AO/MO coefficients.

           olp
                The AO overlap matrix. A TwoIndex instance.

           orb
                A boolean. Initializes orbital optimization.

           **Keywords:**
        '''
        if orb:
            return self.solve_ap1rog_scf(one, two, core, exps, olp, **kwargs)
        else:
            return self.solve_ap1rog(one, two, core, exps, olp, **kwargs)


    def solve_ap1rog(self, one, two, core, exps, olp, **kwargs):
        '''Find Geminal expansion coefficient for some Hamiltonian.
           For restricted, closed-shell AP1roG.

           **Arguments:**

           one, two
                One- and four-index integrals (some Hamiltonian matrix elements)
                expressed in the AO basis (primitives).

           core
                The core energy (not included in 'one' and 'two').

           exps
                An expansion instance. It contains the AO/MO coefficients.

           olp
                The AO overlap matrix. A TwoIndex instance.

           **Keywords:**
                indextrans: 4-index Transformation (tensordot).
                warning: Print warnings (False).
                guess: initial guess; dictionary.
                solver: wfn solver; dictionary.
                maxiter: maximum number of iterations; dictionary.
                dumpci: dump ci coefficient; dictionary.
                thresh: thresholds; dictionary
                printoptions: print level; dictionary.
                swapa: swap orbitals; list.

        '''
        indextrans = kwargs.get('indextrans', 'tensordot')
        warning = kwargs.get('warning', False)
        swapa = kwargs.get('swapa', [])
        guess = kwargs.get('guess', dict({'guess': 'random', 'factor': -0.1,
                           'update': False}))
        guess.setdefault('guess', 'random')
        guess.setdefault('factor', -0.1)
        solver = kwargs.get('solver', dict({'wfn': 'krylov'}))
        solver.setdefault('wfn', 'krylov')
        maxiter = kwargs.get('maxiter', dict({'wfniter': 128}))
        maxiter.setdefault('wfniter', 128)
        dumpci = kwargs.get('dumpci', dict({'amplitudestofile': False,
                            'amplitudesfilename': "./ap1rog_amplitudes.dat"}))
        dumpci.setdefault('amplitudestofile', False)
        dumpci.setdefault('amplitudesfilename', "./ap1rog_amplitudes.dat")
        thresh = kwargs.get('thresh', dict({'wfn': 1e-12}))
        thresh.setdefault('wfn',  1e-12)
        printoptions = kwargs.get('printoptions', dict({'geminal': True, 'ci':
                                   0.01, 'excitationlevel': 1}))
        printoptions.setdefault('geminal', True)
        printoptions.setdefault('ci', 0.01)
        printoptions.setdefault('excitationlevel', 1)
        names = ["guess", "solver", "thresh", "dumpci",
                 "printoptions", "indextrans", "warning", "checkpoint",
                 "swapa", "maxiter"]
        for name, value in kwargs.items():
            if name not in names:
                raise ValueError("Unknown keyword argument %s" % name)

        if warning is False:
            warnings.filterwarnings('ignore')
        if log.do_medium:
            self.print_options(guess, solver, maxiter, thresh, printoptions,
                               indextrans)

        # Do restart and swap orbitals if needed:
        self.start_up(exps, swapa)

        # Do restart and swap orbitals if needed:
        self.start_up(exps, swapa)

        # Generate Guess
        initial_guess = self.generate_guess(guess)

        self.update_ecore(core)
        if log.do_medium:
            log.hline('=')
            log('Starting optimization of AP1roG...')
        # Solve for wave function
        self.solve_wfn(one, two, exps, **{'maxiter': maxiter, 'thresh': thresh,
                       'guess': initial_guess[0],
                       'solver': solver, 'indextrans': indextrans,
                       'orbitaloptimizer': None})

        # Final print statements:
        if log.do_medium:
            self.print_final()
            self.dump_final(exps, olp, printoptions, dumpci, -1)
        # Sanity check for correlation energy:
        if self.get_correlation_energy() > 0:
            raise ValueError('Warning: Final correlation energy is positive! \
                              Improve initial guess!')
        return self.get_total_energy(), self.geminal

    def solve_ap1rog_scf(self, one, two, core, exps, olp, **kwargs):
        '''Find Geminal expansion coefficient for some Hamiltonian.
           For restricted, closed-shell AP1roG.
           Perform orbital optimization.

           **Arguments:**

           one, two
                One- and four-index integrals (some Hamiltonian matrix elements)
                expressed in the AO basis (primitives).

           core
                The core energy (not included in 'one' and 'two').

           exps
                An expansion instance. It contains the AO/MO coefficients.

           olp
                The AO overlap matrix. A TwoIndex instance.

           **Keywords:**
                indextrans: 4-index Transformation (tensordot).
                warning: Print warnings (False).
                guess: initial guess; dictionary.
                solver: wfn/Lagrange solver; dictionary.
                maxiter: maximum number of iterations; dictionary.
                dumpci: dump ci coefficient; dictionary.
                thresh: thresholds; dictionary
                printoptions: print level; dictionary.
                stepsearch: line search options.
                checkpoint: frequency of checkpointing.
                lshift: level shift of Hessian.
                swapa: swap orbitals; list.
                givensrotationa: rotate orbitals; list.

        '''
        indextrans = kwargs.get('indextrans', 'tensordot')
        warning = kwargs.get('warning', False)
        checkpoint = kwargs.get('checkpoint', 1)
        lshift = kwargs.get('lshift', 1e-8)
        pos = kwargs.get('absolute', True)
        givensrotation = kwargs.get('givensrotationa', [])
        swapa = kwargs.get('swapa', [])
        guess = kwargs.get('guess', dict({'guess': 'random', 'factor': -0.1,
                           'update': False}))
        guess.setdefault('guess', 'random')
        guess.setdefault('factor', -0.1)
        solver = kwargs.get('solver', dict({'wfn': 'krylov', 'lagrange':
                            'krylov'}))
        solver.setdefault('wfn', 'krylov')
        solver.setdefault('lagrange', 'krylov')
        maxiter = kwargs.get('maxiter', dict({'wfniter': 128, 'orbiter': 500}))
        maxiter.setdefault('wfniter', 128)
        maxiter.setdefault('orbiter', 500)
        dumpci = kwargs.get('dumpci', dict({'amplitudestofile': False,
                            'amplitudesfilename': "./ap1rog_amplitudes.dat"}))
        dumpci.setdefault('amplitudestofile', False)
        dumpci.setdefault('amplitudesfilename', "./ap1rog_amplitudes.dat")
        thresh = kwargs.get('thresh', dict({'wfn': 1e-12, 'energy': 1e-8,
                            'gradientnorm': 1e-8, 'gradientmax': 5e-5}))
        thresh.setdefault('wfn',  1e-12)
        thresh.setdefault('energy', 1e-8)
        thresh.setdefault('gradientnorm', 1e-8)
        thresh.setdefault('gradientmax', 5e-5)
        printoptions = kwargs.get('printoptions', dict({'geminal': True, 'ci':
                                   0.01, 'excitationlevel': 1}))
        printoptions.setdefault('geminal', True)
        printoptions.setdefault('ci', 0.01)
        printoptions.setdefault('excitationlevel', 1)
        stepsearch = kwargs.get('stepsearch', dict({'method': 'trust-region',
                                'stepa': 1.0, 'c1': 0.0001, 'c2': 0.9,
                                'maxstep': 0.75, 'minstep': 1e-6,
                                'maxiterouter': 40, 'maxeta': 0.25, 'optimizer':
                                'dogleg', 'maxtrustradius': 0.75,
                                'mineta': 0.0, 'upscale': 1.5, 'downscale': 0.7,
                                'trustradius': 0.75, 'threshold': 1e-8,
                                'maxiterinner': 500}))
        stepsearch.setdefault('method', 'trust-region')
        stepsearch.setdefault('stepa', 1.0)
        stepsearch.setdefault('c1', 0.0001)
        stepsearch.setdefault('c2', 0.9)
        stepsearch.setdefault('maxstep', 0.75)
        stepsearch.setdefault('minstep', 1e-6)
        stepsearch.setdefault('maxiterouter', 40)
        stepsearch.setdefault('maxiterinner', 500)
        stepsearch.setdefault('maxeta', 0.25)
        stepsearch.setdefault('mineta', 0.0)
        stepsearch.setdefault('upscale', 1.5)
        stepsearch.setdefault('downscale', 0.7)
        stepsearch.setdefault('trustradius', 0.75)
        stepsearch.setdefault('maxtrustradius', 0.75)
        stepsearch.setdefault('threshold', 1e-8)
        stepsearch.setdefault('optimizer', 'dogleg')
        orbitaloptimizer = kwargs.get('orbitaloptimizer', 'variational')

        names = ["guess", "solver", "stepsearch", "thresh", "dumpci",
                 "printoptions", "indextrans", "warning", "checkpoint",
                 "lshift", "givensrotationa", "swapa", "maxiter",
                 "orbitaloptimizer"]
        for name, value in kwargs.items():
            if name not in names:
                raise ValueError("Unknown keyword argument %s" % name)

        # Set parameters for optimization
        if warning is False:
            warnings.filterwarnings('ignore')
        if maxiter['orbiter'] < 0:
            raise ValueError('Number of iterations must be larger than 0!')

        # Get parameters for optimization and do some sanity tests:
        if log.do_medium:
            self.print_options_scf(guess, solver, maxiter, lshift, stepsearch,
                              thresh, printoptions, checkpoint, indextrans,
                              orbitaloptimizer)

        # Generate Guess
        initial_guess = self.generate_guess(guess, nguess=2)

        # Do restart and swap orbitals if needed:
        self.start_up(exps, swapa, givensrotation)

        if log.do_medium:
            log.hline('=')
            log('Starting optimization of AP1roG (separate optimization of Lagrange multipliers and wavefunction amplitudes)...')

        self.update_ecore(core)

        # First iteration:
        # Read integrals for correlation calculation
        self.solve_wfn(one, two, exps, **{'maxiter': maxiter, 'thresh': thresh,
                       'guess': initial_guess[0], 'guesslm': initial_guess[0],
                       'solver': solver, 'indextrans': indextrans,
                       'orbitaloptimizer': orbitaloptimizer})

        # Get total, reference and correlation energies of AP1roG calculation (zero step):
        energy = self.get_total_energy()
        correlation = self.get_correlation_energy()

        if log.do_medium:
            log.hline('~')
            log('Initial step:')
            self.print_final('Initial')
            log('Entering orbital optimization of AP1roG...')
            log.hline(' ')
            if stepsearch['method']=='trust-region':
                log('%3s %10s %14s  %10s     %10s   %8s   %6s    %6s    %10s'
                     %('step', 'Etot', 'D(Etot)', 'Ecorr', 'D(Ecorr)', 'Max(Grad)',
                       '|Grad|', 'Step','TrustRegion'))
            else:
                log('%3s %10s %14s  %10s     %10s   %8s   %6s     %4s' %('step',
                    'Etot', 'D(Etot)', 'Ecorr', 'D(Ecorr)', 'Max(Grad)', '|Grad|',
                    'Step'))

        i = 0

        linesearch = RStepSearch(**stepsearch)
        while i < maxiter['orbiter']:
            # Store energies from previous iteration step
            energy_old = deepcopy(energy)
            correlation_old = deepcopy(correlation)

            # Calculate orbital gradient and diagonal approximation to the Hessian
            kappa, gradient, hessian = self.do_orbital_rotation_step(lshift,
                                                                     pos,
                                                                     orbitaloptimizer)

            # apply line search:
            linesearch(self, one, two, exps,
                       **{'kappa': kappa, 'thresh': thresh, 'maxiter': maxiter,
                          'gradient': gradient, 'hessian': hessian,
                          'guess': initial_guess[0], 'guesslm': initial_guess[1],
                          'solver': solver, 'indextrans': indextrans,
                          'orbitaloptimizer': orbitaloptimizer})

            energy = self.get_total_energy()
            correlation = self.get_correlation_energy()

            # Print information of iteration step
            if log.do_medium:
                if stepsearch['method']=='trust-region':
                    log('%3i  %14.8f  %11.8f  %10.8f  %11.8f   %6.5f   %5.2e   %1.2e   %1.2e' %(i+1, (energy), (energy-energy_old),
                         correlation, (correlation-correlation_old),
                         np.max(np.absolute(gradient)),
                         np.dot(gradient,gradient), linesearch.stepa,
                         linesearch.trustradius))
                else:
                    log('%3i  %14.8f  %11.8f  %10.8f  %11.8f   %6.5f   %5.2e   %1.2e' %(i+1, (energy), (energy-energy_old), correlation,
                         (correlation-correlation_old),
                         np.max(np.absolute(gradient)),
                         np.dot(gradient,gradient), linesearch.stepa))

            # Checkpoint for orbitals
            if (i+1)%checkpoint == 0 and checkpoint > 0 and (energy < energy_old):
                self.do_checkpoint(exps, olp)

            # Check convergence
            if (self.check_convergence(energy, energy_old, gradient, thresh['energy'], thresh['gradientmax'], thresh['gradientnorm'])):
                if log.do_medium:
                    log.hline(' ')
                    log('Orbital optimization converged in %i iterations' %(i+1))
                    log.hline(' ')
                break
            elif self.check_stepsearch(linesearch):
                if log.do_medium:
                    log.hline(' ')
                    log('Trustradius too small. Orbital optimization aborted!')
                    log.hline(' ')
                break
            else:
                i = i+1

        # Check convergence if i = maxorbiter:
        if i >= maxiter['orbiter'] and i>0:
            if not self.check_convergence(energy, energy_old, gradient, thresh['energy'], thresh['gradientmax'], thresh['gradientnorm']):
                if log.do_medium:
                    log.hline(' ')
                    log('!!!WARNING: Orbital optimization NOT converged in %i iterations' %(i))
                    log.hline(' ')

        if log.do_medium:
            self.print_final()
            self.dump_final(exps, olp, printoptions, dumpci, checkpoint)
        return self.get_total_energy(), self.geminal, self.lagrange

    def update_one_dm(self, select, one_dm=None):
        '''Update 1-RDM

           **Optional arguments:**

           one_dm
                When provided, this 1-RDM is stored.
        '''
        cached_one_dm = self.init_one_dm(select)
        if one_dm is None:
            if select == 'ps2':
                cached_one_dm.compute_ps2_one_dm_ap1rog(self.geminal, self.geminal, factor=1)
            elif select == 'response':
                cached_one_dm.compute_response_one_dm_ap1rog(self.geminal, self.lagrange, factor=1)
            else:
                raise NotImplementedError
        else:
            cached_one_dm.assign(one_dm)
        return cached_one_dm

    def update_two_dm(self, select, two_dm=None):
        '''Update 2-RDM

           **Optional arguments:**

           two_dm
                When provided, this 2-RDM is stored.
        '''
        if two_dm is not None:
            raise NotImplementedError
        if select == 'rppqq':
            cached_two_dm = self.init_two_dm(select)
            cached_two_dm.compute_response_two_dm_ap1rog(self.one_dm_response, self.geminal, self.lagrange, 'ppqq')
        elif select == 'rpqpq':
            cached_two_dm = self.init_two_dm(select)
            cached_two_dm.compute_response_two_dm_ap1rog(self.one_dm_response, self.geminal, self.lagrange, 'pqpq')
        elif select == 'response':
            cached_two_dm1 = self.init_two_dm('rppqq')
            cached_two_dm1.compute_response_two_dm_ap1rog(self.one_dm_response, self.geminal, self.lagrange, 'ppqq')
            cached_two_dm2 = self.init_two_dm('rpqpq')
            cached_two_dm2.compute_response_two_dm_ap1rog(self.one_dm_response, self.geminal, self.lagrange, 'pqpq')
            return cached_two_dm1, cached_two_dm2
        elif select == 'ppqq':
            cached_two_dm = self.init_two_dm(select)
            cached_two_dm.compute_ps2_two_dm_ap1rog(self.one_dm_ps2, self.geminal, self.geminal, select)
        elif select == 'pqpq':
            cached_two_dm = self.init_two_dm(select)
            cached_two_dm.compute_ps2_two_dm_ap1rog(self.one_dm_ps2, self.geminal, self.geminal, select)
        elif select == 'ps2':
            cached_two_dm1 = self.init_two_dm('ppqq')
            cached_two_dm1.compute_ps2_two_dm_ap1rog(self.one_dm_ps2, self.geminal, self.geminal, 'ppqq')
            cached_two_dm2 = self.init_two_dm('pqpq')
            cached_two_dm2.compute_ps2_two_dm_ap1rog(self.one_dm_ps2, self.geminal, self.geminal, 'pqpq')
            return cached_two_dm1, cached_two_dm2
        else:
            raise NotImplementedError
        return cached_two_dm

    def update_three_dm(self, select, three_dm=None):
        '''Update 3-RDM

           **Optional arguments:**

           three_dm
                When provided, this 2-RDM is stored.
        '''
        if three_dm is not None:
            raise NotImplementedError
        if select == 'uuu':
            cached_three_dm = self.init_three_dm('uuu')
            cached_three_dm.compute_response_three_dm_ap1rog(self.geminal, self.lagrange, 'uuu')
        elif select == 'uud':
            cached_three_dm = self.init_three_dm('uud')
            cached_three_dm.compute_response_three_dm_ap1rog(self.geminal, self.lagrange, 'uud')
        elif select == 'uudoff':
            cached_three_dm = self.init_three_dm('uudoff')
            cached_three_dm.compute_response_three_dm_ap1rog(self.geminal, self.lagrange, 'uudoff')
        elif select == 'udd':
            cached_three_dm = self.init_three_dm('udd')
            cached_three_dm.compute_respons_three_dm_ap1rog(self.geminal, self.lagrange, 'udd')
        elif select == 'uddoff':
            cached_three_dm = self.init_three_dm('uddoff')
            cached_three_dm.compute_response_three_dm_ap1rog(self.geminal, self.lagrange, 'uddoff')
        else:
            raise NotImplementedError
        return cached_three_dm

    def update_four_dm(self, select, four_dm=None):
        '''Update 4-RDM

           **Optional arguments:**

           four_dm
                When provided, this 4-RDM is stored.
        '''
        if four_dm is not None:
            raise NotImplementedError
        if select == 'udud':
            cached_four_dm = self.init_four_dm('udud')
            cached_four_dm.compute_response_four_dm_ap1rog(self.two_dm_rpqpq, self.geminal, self.lagrange, 'udud')
        elif select == 'all':
            cached_four_dm = self.init_four_dm('udud')
            cached_four_dm.compute_response_four_dm_ap1rog(self.two_dm_rpqpq, self.geminal, self.lagrange, 'udud')
            return cached_four_dm
        else:
            raise NotImplementedError
        return cached_four_dm

    def clear(self):
        '''Clear all wavefunction information'''
        self._cache.clear()

    def clear_matrix(self):
        '''Clear the matrices'''
        self._cache.clear(tags='m', dealloc=True)

    def solve_wfn(self, one, two, exps, **kwargs):
        '''Solve for wave function.

           **Arguments:**

           one, two
                One- and four-index integrals (some Hamiltonian matrix elements)
                expressed in the AO basis (primitives).

           exps
                An expansion instance. It contains the AO/MO coefficients.

           **Keywords:**
                guess: initial guess for wfn.
                guesslm: initial guess Lagrange multipliers.
                solver: wfn/Lagrange solver; dictionary.
                indextrans: 4-index Transformation (tensordot).
                maxiter: maximum number of iterations; dictionary.
                thresh: thresholds; dictionary
                orbitaloptimizer: orbital optimization method.
        '''
        guess = kwargs.get('guess', None)
        guesslm = kwargs.get('guesslm', None)
        solver = kwargs.get('solver', None)
        indextrans = kwargs.get('indextrans', 'tensordot')
        maxiter = kwargs.get('maxiter', None)
        thresh = kwargs.get('thresh', None)
        orbitaloptimizer = kwargs.get('orbitaloptimizer', 'variational')

        one_mo, two_mo = self.update_mo_integrals(one, two, indextrans, exps)

        # Generate all sorts of matrices
        self.clear_matrix()
        self.update_matrix('all-scf', two_mo, one_mo)
        del one_mo, two_mo

        # Optimize OAP1roG wavefunction amplitudes:
        coeff = self.solve_for_wavefunction(guess, solver, thresh['wfn'], maxiter['wfniter'])
        self.update_geminal(coeff, self.nocc, self.nvirt)

        # Optimize OAP1roG Lagrange multipliers (lambda equations):
        if orbitaloptimizer == 'variational':
            lcoeff = self.solve_for_lagrange(guesslm, solver, thresh['wfn'], maxiter['wfniter'])
            self.update_lagrange(lcoeff, self.nocc, self.nvirt)

    @timer.with_section('ProjectedSEq')
    def solve_for_wavefunction(self, guess, solver, wfnthreshold, wfnmaxiter):
       '''Solves for geminal matrix
       '''
       kmatrix = self.get_matrix('vpp')
       vmatrix = self.get_matrix('v')
       vmatrixa = self.get_matrix('va')
       one_mo = self.get_matrix('t')
       energy = self.get_reference_energy()
       sol = opt.root(self.vector_function_ap1rog,guess,args=(kmatrix,vmatrix,vmatrixa,one_mo,energy),
                      jac=self.jacobian_ap1rog,method=solver['wfn'],options={'xtol': wfnthreshold,'maxfev': wfnmaxiter},callback=None)
       return sol.x

    @timer.with_section('LagrangeMult')
    def solve_for_lagrange(self, guess, solver, wfnthreshold, wfnmaxiter):
       '''Solves for Lagrange multipliers
       '''
       kmatrix = self.get_matrix('vpp')
       vmatrix = self.get_matrix('v')
       vmatrixa = self.get_matrix('va')
       one_mo = self.get_matrix('t')
       energy = self.get_reference_energy()
       cicoeff = self.geminal.pass_array()
       sol = opt.root(self.vector_function_lambda, guess, args=(cicoeff,kmatrix,vmatrix,vmatrixa,one_mo,energy),
                      method=solver['lagrange'], jac=self.jacobian_lambda, callback=None, options={'xtol': wfnthreshold})
       return sol.x

    def init_matrix(self, select):
        '''Initialize matrices
        '''
        if select not in ['t', 'vpp', 'v', 'va', 'vpq', 'vpqrq', 'vpqrr']:
            raise ValueError('The select argument must be one of .')
        if select in ['vpqrq', 'vpqrr']:
            matrix, new = self._cache.load('matrix_%s' % select, alloc=(self._lf.create_three_index, self.nbasis), tags='m')
        else:
            matrix, new = self._cache.load('matrix_%s' % select, alloc=(self._lf.create_two_index, self.nbasis), tags='m')
        if not new:
            raise RuntimeError('The matrix matrix_%s already exists. Call wfn.clear prior to updating the wfn.' % select)
        return matrix

    @timer.with_section('IntegralSort')
    def update_matrix(self, select, two_mo, one_mo=None):
        '''Derive all matrices.
           Vpp:         <pp|qq>,
           V:           <pq|pq>,
           Va:          <pq|qp>,
           Vpq:         2<pq|pq>-<pq|qp>,
           Vpqrq:       2<pq|rq>-<pq|qr>,
           Vpqrr:       <pq|rr>

           **Arguments:**

           select
                'vpp', 'v', 'va', 'vpqrq', 'vpqrr', 'all',
                'all-scf'.

           two_mo
                two-electron integrals to be sorted.

           **Optional arguments:**

           one_mo
                one-electron integrals
        '''
        if select == 't':
            cached_matrix = self.init_matrix(select)
            cached_matrix.assign(one_mo[0])
        elif select == 'vpp':
            cached_matrix = self.init_matrix(select)
            cached_matrix.apply_sorting_aabb(two_mo[0])
        elif select == 'v':
            cached_matrix = self.init_matrix(select)
            cached_matrix.apply_sorting_abab(two_mo[0])
        elif select == 'va':
            cached_matrix = self.init_matrix(select)
            cached_matrix.apply_sorting_abba(two_mo[0])
        elif select == 'vpq':
            cached_matrix = self.init_matrix(select)
            cached_matrix.apply_sorting_pq(two_mo[0])
        elif select == 'vpqrq':
            cached_matrix = self.init_matrix(select)
            cached_matrix.apply_sorting_pqrq(two_mo[0])
        elif select == 'vpqrr':
            cached_matrix = self.init_matrix(select)
            cached_matrix.apply_sorting_ijkk(two_mo[0])
        elif select == 'all':
            cached_matrix = self.init_matrix('t')
            cached_matrix.assign(one_mo[0])
            cached_matrix1= self.init_matrix('vpp')
            cached_matrix1.apply_sorting_aabb(two_mo[0])
            cached_matrix2= self.init_matrix('v')
            cached_matrix2.apply_sorting_abab(two_mo[0])
            cached_matrix3= self.init_matrix('va')
            cached_matrix3.apply_sorting_abba(two_mo[0])
            return cached_matrix, cached_matrix1, cached_matrix2, cached_matrix3
        elif select == 'all-scf':
            cached_matrix = self.init_matrix('t')
            cached_matrix.assign(one_mo[0])
            cached_matrix1= self.init_matrix('vpp')
            cached_matrix1.apply_sorting_aabb(two_mo[0])
            cached_matrix2= self.init_matrix('v')
            cached_matrix2.apply_sorting_abab(two_mo[0])
            cached_matrix3= self.init_matrix('va')
            cached_matrix3.apply_sorting_abba(two_mo[0])
            cached_matrix4 = self.init_matrix('vpq')
            cached_matrix4.apply_sorting_pq(two_mo[0])
            cached_matrix5 = self.init_matrix('vpqrq')
            cached_matrix5.apply_sorting_pqrq(two_mo[0])
            cached_matrix6 = self.init_matrix('vpqrr')
            cached_matrix6.apply_sorting_ijkk(two_mo[0])
            return cached_matrix, cached_matrix1, cached_matrix2, cached_matrix3, cached_matrix4, cached_matrix5, cached_matrix6
        else:
            raise RuntimeError('The matrix matrix_%s is not supported.' % select)
        return cached_matrix

    def get_matrix(self, select):
        '''Get a matrix.

           **Arguments:**

           select
                't', 'vpp', 'v', 'va', 'vpp', 'vpq', 'vpqrq', or 'vpqrr'.
        '''
        if not 'matrix_%s' % select in self._cache:
            raise NotImplementedError
            self.update_matrix(select, two_mo)
        return self._cache.load('matrix_%s' % select)

    # Functions for energy evaluation:
    def get_correlation_energy(self, coeff=None):
        '''Get correlation energy of restricted AP1roG.
        '''
        if coeff is None:
            coeff = self.geminal
        if isinstance(coeff, TwoIndex):
            tmp = coeff.copy()
        else:
            tmp = self.lf.create_two_index(self.npairs, self.nvirt)
            tmp.assign_array(coeff, self.npairs, self.nvirt)
        kmat = self.get_matrix('vpp')
        kia = kmat.copyview(0, self.npairs, self.npairs, self.nbasis, 1.0)
        tmp.imultiply(kia)
        return tmp.contract('ab', 1.0)

    def get_reference_energy(self):
        '''Get energy of reference Slater Determinant for restricted AP1roG.
        '''
        one = self.get_matrix('t')
        two = self.get_matrix('v')
        twoa = self.get_matrix('va')
        energy = one.contract('aa', 2.0, 0, self.npairs)
        energy += two.contract('ab', 2.0, 0, self.npairs, 0, self.npairs)
        energy += twoa.contract('ab', -1.0, 0, self.npairs, 0, self.npairs)
        return energy

    def get_total_energy(self, coeff=None):
        '''Get total energy (reference+correlation) including nuclear-repulsion energy for restricted AP1roG.
        '''
        return (self.get_correlation_energy(coeff)+self.get_reference_energy()+self.ecore)

    def get_rotation_matrix(self, coeff):
        '''Determine orbital rotation matrix for (oo), (vo), and (vv) blocks.
        '''
        output = np.identity((self.nbasis))
        kappa = np.zeros((self.nbasis, self.nbasis))
        indl = np.tril_indices(self.nbasis, -1)
        kappa[indl] = coeff
        kappa -= kappa.T

        # Approximate Unitary rotation matrix U = exp(-K) by U = 1 - K + 1/2 K^2 + O(3)
        output = output - kappa + 0.5*np.dot(kappa,kappa) #- 0.25*np.dot(tmp,np.dot(tmp,tmp))
        # orthogonalization because approximate U matrix might not be Unitary/Orthogonal:
        q,r = np.linalg.qr(output, mode='full')
        return q

    def do_orbital_rotation_step(self, lshift=1e-8, pos=True,
                                 optimizer='variational'):
        '''Get orbital rotation step (Newton--Raphson)

           **Arguments:**


           **Optional arguments:**

           lshift
               Level shift for approximate Hessian added to small elements.
               Default: 1e-8.

           pos
               Make all elements of Hessian positive if set to True.

           optimizer
               Orbital otimization method.

        '''
        # Calculate orbital gradient and diagonal approximation to the Hessian
        if optimizer == 'variational':
            grad = self.get_orbital_gradient()
            ind = np.where(abs(grad) < 1e-8)
            grad[ind] = 0.0
            # We use a diagonal hessian
            hessian = self.get_orbital_hessian(lshift, pos)
            kappa = -(np.divide(grad,hessian))
        elif optimizer == 'ps2c':
            grad = self.get_orbital_gradient(True)
            ind = np.where(abs(grad) < 1e-8)
            grad[ind] = 0.0
            # We use a diagonal hessian
            hessian = self.get_orbital_hessian(lshift, pos, True)
            kappa = -(np.divide(grad,hessian))
        else:
            raise NotImplementedError
        return kappa, grad, hessian

    # Vector function for AP1roG:
    def vector_function_ap1rog(self, coeff, kmatrix, vmatrix, vmatrixa, one_mo, energy):
        '''Construct vector function for nonlinear optimization of coefficients
           for restricted AP1roG.
        '''
        gmat = self.lf.create_two_index(self.npairs, self.nvirt)
        gmat.assign_array(coeff, self.npairs, self.nvirt)

        en = energy+self.get_correlation_energy(gmat)

        kia = kmatrix.copyview(0, self.npairs, self.npairs, self.nbasis)
        kij = kmatrix.copyview(0, self.npairs, 0, self.npairs)
        kab = kmatrix.copyview(self.npairs, self.nbasis, self.npairs, self.nbasis)
        kca = self.lf.create_one_index(self.nvirt)
        kci = self.lf.create_one_index(self.npairs)
        kii = self.lf.create_one_index(self.npairs)
        kaa = self.lf.create_one_index(self.nvirt)
        oneaa = self.lf.create_one_index(self.nvirt)
        oneii = self.lf.create_one_index(self.npairs)
        vaa = self.lf.create_one_index(self.nvirt)
        vii = self.lf.create_one_index(self.npairs)
        va = self.lf.create_one_index(self.nvirt)
        vi = self.lf.create_one_index(self.npairs)
        va_a = self.lf.create_one_index(self.nvirt)
        vi_a = self.lf.create_one_index(self.npairs)
        kci.contract_2twoindex(kia, gmat, 2)
        kca.contract_2twoindex(kia, gmat, 1)
        kii.assign_diagonal(kij)
        kaa.assign_diagonal(kab)
        oneaa.assign_diagonal(one_mo, 1, self.npairs, self.nbasis)
        oneii.assign_diagonal(one_mo, 1, 0, self.npairs)
        vaa.assign_diagonal(vmatrix, 1, self.npairs, self.nbasis)
        vii.assign_diagonal(vmatrix, 1, 0, self.npairs)
        vai = vmatrix.copyview(0, self.npairs, self.npairs, self.nbasis)
        vaia = vmatrixa.copyview(0, self.npairs, self.npairs, self.nbasis)
        va.contract_twoindex(vmatrix, 1, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        vi.contract_twoindex(vmatrix, 2, 1.0, 0, self.npairs, 0, self.npairs)
        va_a.contract_twoindex(vmatrixa, 1, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        vi_a.contract_twoindex(vmatrixa, 2, 1.0, 0, self.npairs, 0, self.npairs)
        kgeminal = kia.contract_twoindex(gmat)

        tmp = gmat.copy()
        tmp.imultiply(gmat)
        tmp.imultiply(kia)
        tmp2 = self.lf.create_two_index(self.npairs, self.npairs)
        tmp2.iadddott(gmat, kia, 1.0)

        result = self.lf.create_two_index(self.npairs,self.nvirt)
        # Kai:
        result.iadd(kia, 1.0)
        result.iadddot(kij, gmat, 1.0)
        result.iadddot(gmat, kab, 1.0)
        result.iadd(gmat, (kgeminal-en+energy))
        result.iadd(tmp, 2.0)
        result.multiplyt(gmat, kca, -2.0)
        result.multiply(gmat, kci, -2.0)
        result.multiply(gmat, kii, -1.0)
        result.multiplyt(gmat, kaa, -1.0)
        result.iadddot(tmp2, gmat)
        result.multiplyt(gmat, oneaa, 2.0)
        result.multiply(gmat, oneii, -2.0)
        result.multiplyt(gmat, vaa, 2.0)
        result.multiplyt(gmat, va, 4.0)
        result.multiply(gmat, vi, -4.0)
        result.multiply(gmat, vii, 2.0)
        result.multiply(gmat, vai, -4.0)
        result.multiplyt(gmat, va_a, -2.0)
        result.multiplyt(gmat, vaa, -1.0)
        result.multiply(gmat, vi_a, 2.0)
        result.multiply(gmat, vii, -1.0)
        result.multiply(gmat, vaia, 2.0)

        return result._array.ravel(order='C')

    # Jacobian for AP1roG:
    def jacobian_ap1rog(self, coeff, kmatrix, vmatrix, vmatrixa, one_mo, energy):
        '''Construct Jacobian for nonlinear optimization of coefficients for restricted AP1roG.
        '''
        jacobian = self.lf.create_two_index((self.npairs*self.nvirt),(self.npairs*self.nvirt))

        geminal = self.lf.create_two_index(self.npairs, self.nvirt)
        geminal.assign_array(coeff, self.npairs, self.nvirt)
        en = energy+self.get_correlation_energy(geminal)
        kmat = kmatrix.copyview(0, self.npairs, self.npairs, self.nbasis)
        kij = kmatrix.copyview(0, self.npairs, 0, self.npairs)
        kab = kmatrix.copyview(self.npairs, self.nbasis, self.npairs, self.nbasis)
        gmat = self.lf.create_two_index(self.npairs, self.nvirt)
        gmat.set_value(1.0)
        oneaa = self.lf.create_one_index(self.nvirt)
        oneii = self.lf.create_one_index(self.npairs)
        vaa = self.lf.create_one_index(self.nvirt)
        vii = self.lf.create_one_index(self.npairs)
        va = self.lf.create_one_index(self.nvirt)
        vi = self.lf.create_one_index(self.npairs)
        va_a = self.lf.create_one_index(self.nvirt)
        vi_a = self.lf.create_one_index(self.npairs)
        oneaa.assign_diagonal(one_mo, 1, self.npairs, self.nbasis)
        oneii.assign_diagonal(one_mo, 1, 0, self.npairs)
        vaa.assign_diagonal(vmatrix, 1, self.npairs, self.nbasis)
        vii.assign_diagonal(vmatrix, 1, 0, self.npairs)
        vai = vmatrix.copyview(0, self.npairs, self.npairs, self.nbasis)
        vaia = vmatrixa.copyview(0, self.npairs, self.npairs, self.nbasis)
        va.contract_twoindex(vmatrix, 1, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        vi.contract_twoindex(vmatrix, 2, 1.0, 0, self.npairs, 0, self.npairs)
        va_a.contract_twoindex(vmatrixa, 1, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        vi_a.contract_twoindex(vmatrixa, 2, 1.0, 0, self.npairs, 0, self.npairs)
        kcab = self.lf.create_two_index(self.nvirt, self.nvirt)
        kcij = self.lf.create_two_index(self.npairs, self.npairs)
        kcib = self.lf.create_one_index(self.npairs)
        kcja = self.lf.create_one_index(self.nvirt)
        kcib.contract_2twoindex(kmat, geminal, 2)
        kcja.contract_2twoindex(kmat, geminal, 1)
        kcab.iaddtdot(geminal, kmat)
        kcij.iadddott(geminal, kmat)
        eyeocc = self.lf.create_two_index(self.npairs, self.npairs)
        eyevir = self.lf.create_two_index(self.nvirt, self.nvirt)
        onesocc = self.lf.create_two_index(self.npairs, self.npairs)
        onesvir = self.lf.create_two_index(self.nvirt, self.nvirt)
        kronvir = self.lf.create_two_index((self.npairs*self.nvirt),(self.npairs*self.nvirt))
        kronocc = self.lf.create_two_index((self.npairs*self.nvirt),(self.npairs*self.nvirt))
        eyeocc.set_diagonal(1.0)
        eyevir.set_diagonal(1.0)
        onesocc.set_value(1.0)
        onesvir.set_value(1.0)
        kronocc.iaddkron(eyeocc, onesvir)
        kronvir.iaddkron(onesocc, eyevir)
        tmp = self.lf.create_two_index((self.npairs*self.nvirt),(self.npairs*self.nvirt))
        tmp.iaddouter(geminal, kmat, -2.0)

        result = self.lf.create_two_index(self.npairs,self.nvirt)
        factor = kmat.contract_twoindex(geminal)
        result.multiplyt(gmat, oneaa, 2.0)
        result.multiply(gmat, oneii, -2.0)
        result.multiplyt(gmat, vaa, 2.0)
        result.multiplyt(gmat, va, 4.0)
        result.multiply(gmat, vi, -4.0)
        result.multiply(gmat, vii, 2.0)
        result.multiply(gmat, vai, -4.0)
        result.multiplyt(gmat, va_a, -2.0)
        result.multiplyt(gmat, vaa, -1.0)
        result.multiply(gmat, vi_a, 2.0)
        result.multiply(gmat, vii, -1.0)
        result.multiply(gmat, vaia, 2.0)
        result.iaddnumber((energy))
        result.iaddnumber(factor-en)
        result.iadd(kcib, -1.0)
        result.iaddt(kcja, -1.0)

        jacobian.multiply(tmp, kronocc)
        jacobian.multiply(tmp, kronvir)
        jacobian.iaddkron(eyeocc, kcab)
        jacobian.iaddkron(kcij, eyevir)
        jacobian.iaddkron(eyeocc, kab)
        jacobian.iaddkron(kij, eyevir)
        jacobian.set_diagonal(result._array.ravel())

        return jacobian._array

    # Vector function for Lagrange multipliers of OAP1roG:
    def vector_function_lambda(self, lambdacoeff, coeff, kmatrix, vmatrix, vmatrixa, one_mo, energy):
        '''Construct vector function for optimization of Lagrange multipliers
           for restricted OAP1roG.
        '''
        gmat = self.lf.create_two_index(self.npairs, self.nvirt)
        gmat.assign_array(coeff, self.npairs, self.nvirt)
        en = energy+self.get_correlation_energy(coeff)

        lmat = self.lf.create_two_index(self.npairs, self.nvirt)
        lmat.assign_array(lambdacoeff, self.npairs, self.nvirt)

        kia = kmatrix.copyview(0, self.npairs, self.npairs, self.nbasis)
        kij = kmatrix.copyview(0, self.npairs, 0, self.npairs)
        kab = kmatrix.copyview(self.npairs, self.nbasis, self.npairs, self.nbasis)
        kca = self.lf.create_one_index(self.nvirt)
        kci = self.lf.create_one_index(self.npairs)
        kii = self.lf.create_one_index(self.npairs)
        kaa = self.lf.create_one_index(self.nvirt)
        lca = self.lf.create_one_index(self.nvirt)
        lci = self.lf.create_one_index(self.npairs)
        oneaa = self.lf.create_one_index(self.nvirt)
        oneii = self.lf.create_one_index(self.npairs)
        vaa = self.lf.create_one_index(self.nvirt)
        vii = self.lf.create_one_index(self.npairs)
        va = self.lf.create_one_index(self.nvirt)
        vi = self.lf.create_one_index(self.npairs)
        va_a = self.lf.create_one_index(self.nvirt)
        vi_a = self.lf.create_one_index(self.npairs)
        kci.contract_2twoindex(kia, gmat, 2)
        kca.contract_2twoindex(kia, gmat, 1)
        kii.assign_diagonal(kij)
        kaa.assign_diagonal(kab)
        oneaa.assign_diagonal(one_mo, 1, self.npairs, self.nbasis)
        oneii.assign_diagonal(one_mo, 1, 0, self.npairs)
        vaa.assign_diagonal(vmatrix, 1, self.npairs, self.nbasis)
        vii.assign_diagonal(vmatrix, 1, 0, self.npairs)
        vai = vmatrix.copyview(0, self.npairs, self.npairs, self.nbasis)
        vaia = vmatrixa.copyview(0, self.npairs, self.npairs, self.nbasis)
        va.contract_twoindex(vmatrix, 1, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        vi.contract_twoindex(vmatrix, 2, 1.0, 0, self.npairs, 0, self.npairs)
        va_a.contract_twoindex(vmatrixa, 1, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        vi_a.contract_twoindex(vmatrixa, 2, 1.0, 0, self.npairs, 0, self.npairs)
        lci.contract_2twoindex(gmat, lmat, 2, 1.0)
        lca.contract_2twoindex(gmat, lmat, 1, 1.0)
        kgeminal = kia.contract_twoindex(gmat)

        tmp = lmat.copy()
        tmp.imultiply(gmat)
        tmp.imultiply(kia)
        tmp2 = self.lf.create_two_index(self.npairs, self.npairs)
        tmp2.iadddott(kia, gmat, 1.0)
        tmp3 = self.lf.create_two_index(self.npairs, self.npairs)
        tmp3.iadddott(lmat, gmat, 1.0)

        result = self.lf.create_two_index(self.npairs,self.nvirt)
        # Kai:
        result.iadd(kia, 1.0)
        result.iadddot(kij, lmat, 1.0)
        result.iadddot(lmat, kab, 1.0)
        result.iadd(lmat, (kgeminal-en))
        result.iadd(tmp, 4.0)
        result.multiplyt(kia, lca, -2.0)
        result.multiply(kia, lci, -2.0)
        result.multiplyt(lmat, kca, -2.0)
        result.multiply(lmat, kci, -2.0)
        result.multiply(lmat, kii, -1.0)
        result.multiplyt(lmat, kaa, -1.0)
        result.iadddot(tmp2, lmat)
        result.iadddot(tmp3, kia)

        result.iadd(lmat, (energy))
        result.multiplyt(lmat, oneaa, 2.0)
        result.multiply(lmat, oneii, -2.0)
        result.multiplyt(lmat, vaa, 2.0)
        result.multiplyt(lmat, va, 4.0)
        result.multiply(lmat, vi, -4.0)
        result.multiply(lmat, vii, 2.0)
        result.multiply(lmat, vai, -4.0)
        result.multiplyt(lmat, va_a, -2.0)
        result.multiplyt(lmat, vaa, -1.0)
        result.multiply(lmat, vi_a, 2.0)
        result.multiply(lmat, vii, -1.0)
        result.multiply(lmat, vaia, 2.0)

        lambdacoeff = lambdacoeff.ravel(order='C')
        return result._array.ravel(order='C')

    # Jacobian for Lagrange multipliers of OAP1roG:
    def jacobian_lambda(self, lambdacoeff, cicoeff, kmatrix, vmatrix, vmatrixa, one_mo, energy):
        '''Construct Jacobian for nonlinear optimization of Lagrange multipliers for restricted OAP1roG.
        '''
        jacobian = self.lf.create_two_index((self.npairs*self.nvirt),(self.npairs*self.nvirt))

        geminal = self.lf.create_two_index(self.npairs, self.nvirt)
        geminal.assign_array(cicoeff, self.npairs, self.nvirt)
        en = energy+self.get_correlation_energy()
        kmat = kmatrix.copyview(0, self.npairs, self.npairs, self.nbasis)
        kij = kmatrix.copyview(0, self.npairs, 0, self.npairs)
        kab = kmatrix.copyview(self.npairs, self.nbasis, self.npairs, self.nbasis)
        gmat = self.lf.create_two_index(self.npairs, self.nvirt)
        gmat.set_value(1.0)
        oneaa = self.lf.create_one_index(self.nvirt)
        oneii = self.lf.create_one_index(self.npairs)
        vaa = self.lf.create_one_index(self.nvirt)
        vii = self.lf.create_one_index(self.npairs)
        va = self.lf.create_one_index(self.nvirt)
        vi = self.lf.create_one_index(self.npairs)
        va_a = self.lf.create_one_index(self.nvirt)
        vi_a = self.lf.create_one_index(self.npairs)
        oneaa.assign_diagonal(one_mo, 1, self.npairs, self.nbasis)
        oneii.assign_diagonal(one_mo, 1, 0, self.npairs)
        vaa.assign_diagonal(vmatrix, 1, self.npairs, self.nbasis)
        vii.assign_diagonal(vmatrix, 1, 0, self.npairs)
        vai = vmatrix.copyview(0, self.npairs, self.npairs, self.nbasis)
        vaia = vmatrixa.copyview(0, self.npairs, self.npairs, self.nbasis)
        va.contract_twoindex(vmatrix, 1, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        vi.contract_twoindex(vmatrix, 2, 1.0, 0, self.npairs, 0, self.npairs)
        va_a.contract_twoindex(vmatrixa, 1, 1.0, 0, self.npairs, self.npairs, self.nbasis)
        vi_a.contract_twoindex(vmatrixa, 2, 1.0, 0, self.npairs, 0, self.npairs)
        kcab = self.lf.create_two_index(self.nvirt, self.nvirt)
        kcij = self.lf.create_two_index(self.npairs, self.npairs)
        kcib = self.lf.create_one_index(self.npairs)
        kcja = self.lf.create_one_index(self.nvirt)
        kcib.contract_2twoindex(kmat, geminal, 2)
        kcja.contract_2twoindex(kmat, geminal, 1)
        kcab.iaddtdot(kmat, geminal)
        kcij.iadddott(kmat, geminal)
        eyeocc = self.lf.create_two_index(self.npairs, self.npairs)
        eyevir = self.lf.create_two_index(self.nvirt, self.nvirt)
        onesocc = self.lf.create_two_index(self.npairs, self.npairs)
        onesvir = self.lf.create_two_index(self.nvirt, self.nvirt)
        kronvir = self.lf.create_two_index((self.npairs*self.nvirt),(self.npairs*self.nvirt))
        kronocc = self.lf.create_two_index((self.npairs*self.nvirt),(self.npairs*self.nvirt))
        eyeocc.set_diagonal(1.0)
        eyevir.set_diagonal(1.0)
        onesocc.set_value(1.0)
        onesvir.set_value(1.0)
        kronocc.iaddkron(eyeocc, onesvir)
        kronvir.iaddkron(onesocc, eyevir)
        tmp2 = self.lf.create_two_index((self.npairs*self.nvirt),(self.npairs*self.nvirt))
        tmp2.iaddouter(kmat, geminal, -2.0)

        result = self.lf.create_two_index(self.npairs,self.nvirt)
        factor = kmat.contract_twoindex(geminal)
        result.multiplyt(gmat, oneaa, 2.0)
        result.multiply(gmat, oneii, -2.0)
        result.multiplyt(gmat, vaa, 2.0)
        result.multiplyt(gmat, va, 4.0)
        result.multiply(gmat, vi, -4.0)
        result.multiply(gmat, vii, 2.0)
        result.multiply(gmat, vai, -4.0)
        result.multiplyt(gmat, va_a, -2.0)
        result.multiplyt(gmat, vaa, -1.0)
        result.multiply(gmat, vi_a, 2.0)
        result.multiply(gmat, vii, -1.0)
        result.multiply(gmat, vaia, 2.0)
        result.iaddnumber((energy))
        result.iaddnumber(factor-en)
        result.iadd(kcib, -1.0)
        result.iaddt(kcja, -1.0)

        jacobian.multiply(tmp2, kronocc)
        jacobian.multiply(tmp2, kronvir)
        jacobian.iaddkron(eyeocc, kcab)
        jacobian.iaddkron(kcij, eyevir)
        jacobian.iaddkron(eyeocc, kab)
        jacobian.iaddkron(kij, eyevir)
        jacobian.set_diagonal(result._array.ravel())

        return jacobian._array

    # Function to calculate the orbital gradient for OAP1roG (oo,vo,vv):
    @timer.with_section('OOAP1roG grad')
    def get_orbital_gradient(self, ps2c=False):
        '''Determine orbital gradient for all non-reduntant orbital rotations
           (oo,vo,vv).
        '''
        one_mo = self.get_matrix('t')
        vpqrq = self.get_matrix('vpqrq')
        vpqrr = self.get_matrix('vpqrr')

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

        gradient = self.lf.create_two_index(self.nbasis, self.nbasis)

        two_dm_av = twodmppqq.symmetrize()
        gradient.contract_three2one(vpqrq, twodmpqpq, '13', 4)
        gradient.contract_three2one(vpqrr, two_dm_av, '21', 4)
        gradient.contract_three2one(vpqrq, twodmpqpq, '31', -4)
        gradient.contract_three2one(vpqrr, two_dm_av, '12', -4)
        gradient.multiply(one_mo, onedm, 4)
        gradient.multiplyt(one_mo, onedm, -4)

        self.clear_dm()

        ind = np.tril_indices(self.nbasis, -1)

        return gradient.pass_array(ind)

    # Approximate diagonal Hessian for OAP1roG:
    def get_orbital_hessian(self, lshift=1e-8, pos=True, ps2c=False):
        '''Construct diagonal approximation to Hessian for orbital optimization
           of restricted OAP1roG.
        '''
        one_mo = self.get_matrix('t')
        vpq = self.get_matrix('vpq')
        vpp = self.get_matrix('vpp')
        vmat = self.get_matrix('v')
        vmata = self.get_matrix('va')

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

        two_dm_av = twodmppqq.symmetrize()
        two_dm_av.set_diagonal(0.0)
        twodmpqpq.set_diagonal(0.0)
        one = onedm.copy()

        one_diag = self.lf.create_one_index(self.nbasis)
        two_c = self.lf.create_one_index(self.nbasis)
        two_ca = self.lf.create_one_index(self.nbasis)
        vmatdiag = self.lf.create_one_index(self.nbasis)
        vmatdiag.assign_diagonal(vmat)
        one_diag.assign_diagonal(one_mo)
        two_c.contract_2twoindex(vpq, twodmpqpq, 2, 1.0)
        two_ca.contract_2twoindex(vpp, two_dm_av, 2, 1.0)

        hessian = self.lf.create_two_index(self.nbasis, self.nbasis)
        # <qt> G_pt
        hessian.iadddot(vpq, twodmpqpq, 4.0)
        # <pt> G_qt
        hessian.iadddot(twodmpqpq, vpq, 4.0)
        # <qt> G_qt
        hessian.iadd(two_c, -4.0)
        # <pt> G_pt
        hessian.iaddt(two_c, -4.0)
        # <qt> G_pt
        hessian.iadddot(two_dm_av, vpp, 4.0)
        # <pt> G_qt
        hessian.iadddot(vpp, two_dm_av, 4.0)
        # <qt> G_qt
        hessian.iadd(two_ca, -4.0)
        # <pt> G_pt
        hessian.iaddt(two_ca, -4.0)
        # <pq> G_pq
        hessian.multiply(vmat, twodmpqpq, 8.0)
        hessian.multiply(vmata, twodmpqpq, -8.0)
        hessian.multiply(vpp, twodmpqpq, -16.0)
        hessian.multiply(vmat, two_dm_av, -8.0)
        hessian.multiplyt(vpp, two_dm_av, -8.0)
        # <ppqq> G_pp
        hessian.multiply(vpp, one, 8.0)
        hessian.multiplyt(vpp, one, 8.0)
        # <qq> g_pp
        hessian.multiplyt(one_diag, onedm, 4.0)
        # <pp> g_qq
        hessian.multiplyt(onedm, one_diag, 4.0)
        # <qq> g_qq
        hessian.multiply(one_diag, onedm, -4.0)
        # <pp> g_pp
        hessian.multiplytt(onedm, one_diag, -4.0)
        # missing terms due to zeros in DMs
        hessian.multiply(vmat, one, 4.0)
        hessian.multiplyt(vmat, one, 4.0)
        hessian.multiply(vmatdiag, one, -4.0)
        hessian.multiplytt(vmatdiag, one, -4.0)

        self.clear_dm()

        # Make everything positive Add levelshift:
        if pos:
            hessian.iabsolute()
        if lshift:
            hessian.iaddshift(lshift)
        ind = np.tril_indices(self.nbasis, -1)

        return hessian.pass_array(ind)

    @timer.with_section('exact Hessian')
    def get_exact_hessian(self, twomo):
        '''Construct exact Hessian for orbital optimization of restricted OAP1roG.
        '''
        one_mo = self.get_matrix('t')
        twomoa = twomo.copy()
        twomoa.add_exchange_part()
        vpqr = self.lf.create_three_index()
        vpqra = self.lf.create_three_index()
        vppr = self.lf.create_three_index()
        # Prepare 2el integrals
        vpqr.apply_sorting_pqrq_woexchange(twomo)
        vpqra.apply_sorting_pqrq_exchange(twomo)
        vppr.apply_sorting_ijkk(twomo)

        self.clear_dm()
        self.update_one_dm('response')
        self.update_two_dm('response')
        onedm = self.one_dm_response
        twodmpqpq = self.two_dm_rpqpq
        twodmppqq = self.two_dm_rppqq
        # Symmetrize 2DM
        two_dm_av = twodmppqq.symmetrize()
        twodmpqpqaa = twodmpqpq.copy()
        # Reset diagonal elements of DMs
        twodmpqpqaa.set_diagonal(0)
        twodmpqpq.set_diagonal(onedm)
        two_dm_av.set_diagonal(0)

        hessian = self.lf.create_four_index()
        tmp2i = self.lf.create_two_index()
        tmp3i = self.lf.create_three_index()
        tmp3i2 = self.lf.create_three_index()
        tmp3i3 = self.lf.create_three_index()
        tmp3i4 = self.lf.create_three_index()

        # (1)
        # aa
        hessian.contract_twoindex(twomoa, twodmpqpqaa, 'cd-acbd', 4)
        hessian.contract_twoindex(twomoa, twodmpqpqaa, 'ab-acdb', -4)
        hessian.contract_twoindex(twomoa, twodmpqpqaa, 'cd-acdb', -4)
        hessian.contract_twoindex(twomoa, twodmpqpqaa, 'ab-acbd', 4)
        # ab
        hessian.contract_twoindex(twomo, twodmpqpq, 'cd-acbd', 4)
        hessian.contract_twoindex(twomo, twodmpqpq, 'ab-acdb', -4)
        hessian.contract_twoindex(twomo, twodmpqpq, 'cd-acdb', -4)
        hessian.contract_twoindex(twomo, twodmpqpq, 'ab-acbd', 4)
        # (2)
        # aa
        hessian.contract_twoindex(twomoa, twodmpqpqaa, 'cb-acdb', 4)
        hessian.contract_twoindex(twomoa, twodmpqpqaa, 'ad-acbd', -4)
        hessian.contract_twoindex(twomoa, twodmpqpqaa, 'cb-acbd', -4)
        hessian.contract_twoindex(twomoa, twodmpqpqaa, 'ad-acdb', 4)
        # ab
        hessian.contract_twoindex(twomo, twodmpqpq, 'cb-acdb', 4)
        hessian.contract_twoindex(twomo, twodmpqpq, 'ad-acbd', -4)
        hessian.contract_twoindex(twomo, twodmpqpq, 'cb-acbd', -4)
        hessian.contract_twoindex(twomo, twodmpqpq, 'ad-acdb', 4)
        # (3)
        hessian.contract_twoindex(twomo, two_dm_av, 'bd-abcd', 4)
        hessian.contract_twoindex(twomo, two_dm_av, 'ad-abcd', -4)
        hessian.contract_twoindex(twomo, two_dm_av, 'bd-abdc', -4)
        hessian.contract_twoindex(twomo, two_dm_av, 'ac-abcd', 4)
        hessian.contract_twoindex(twomo, two_dm_av, 'bc-abdc', 4)
        hessian.contract_twoindex(twomo, two_dm_av, 'ac-abdc', -4)
        hessian.contract_twoindex(twomo, two_dm_av, 'bc-abcd', -4)
        hessian.contract_twoindex(twomo, two_dm_av, 'ad-abdc', 4)
        # Apq,qw (pw) (qv) (-qw) (-pv)
        tmp2i.contract_oneindex(one_mo, onedm, 'b-ab', 2)
        tmp2i.contract_oneindex(one_mo, onedm, 'a-ab', 2)
        # aa
        tmp2i.contract_three2one(vpqra, twodmpqpqaa, '132', 2)
        tmp2i.contract_three2one(vpqra, twodmpqpqaa, '13', 2)
        # ab
        tmp2i.contract_three2one(vpqr, twodmpqpq, '132', 2)
        tmp2i.contract_three2one(vpqr, twodmpqpq, '13', 2)
        tmp2i.contract_three2one(vppr, two_dm_av, '12', 2)
        tmp2i.contract_three2one(vppr, two_dm_av, '21', 2)
        # Apq,qw (pq,qw) (-pq,vq)
        tmp3i.expand_tothreeindex(one_mo, onedm, 'cab', -4)
        # aa
        tmp3i3.contract_twoindex(vpqra, twodmpqpqaa, 'db-adc', -4)
        # ab
        tmp3i3.contract_twoindex(vpqr, twodmpqpq, 'db-adc', -4)
        tmp3i3.contract_twoindex(vppr, two_dm_av, 'dc-adb', -4)
        # Apq,qw (pq,vp) (-pq,pw)
        tmp3i2.expand_tothreeindex(one_mo, onedm, 'acb', -4)
        # aa
        tmp3i4.contract_twoindex(vpqra, twodmpqpqaa, 'db-dac', -4)
        # ab
        tmp3i4.contract_twoindex(vpqr, twodmpqpq, 'db-dac', -4)
        tmp3i4.contract_twoindex(vppr, two_dm_av, 'dc-dab', -4)

        hessian.assign_hessian(tmp2i, tmp3i, tmp3i2, tmp3i3, tmp3i4)
        dim = (self.nbasis*(self.nbasis-1))/2
        output = self.lf.create_two_index(dim, dim)
        output.assign_tril_fourindex(hessian, self.nbasis)

        return output.pass_array()

    def prod(self, lst):
        return reduce(mul, lst, 1)

    def perm(self, a):
        n = len(a)
        r = range(n)
        s = permutations(r)
        return fsum(self.prod(a[i][sigma[i]] for i in r) for sigma in s)

    def print_solution(self, cithresh=0.01, excitationlevel=2,
                       amplitudestofile=False, filename="./ap1rog_amplitudes"):
        '''Print coefficients {ci} and Slater Determinant if |ci| > threshold.
           Prints up to heptuple excited pairs.
        '''
        it = []
        coeff = self.geminal._array.copy()
        matrix = np.identity(self._nocc)
        for i in range(excitationlevel):
            it.append(np.nditer(coeff,flags=['multi_index']))
        if amplitudestofile is True:
            filea = open(filename, 'w')
            filea.write('{0:30} {1:20}'.format('Determinant','Amplitude\n'))
            filea.write('{0:1} {1:20.16f}'.format(self.get_slater_determinant((-1,-1)), 1.0))
            filea.close()
            filea = open(filename, 'a')
            filea.write("\n")
        else:
            print '{0:10} {1:20.16f}'.format(self.get_slater_determinant((-1,-1)), 1.0)
        for i in range(excitationlevel):
            if i==0:
                for ci in it[0]:
                    if math.fabs(ci) >= (cithresh):
                        if amplitudestofile is True:
                            filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index), float(ci)))
                            filea.write("\n")
                        else:
                            print '{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index), float(ci))
            if self._nocc > 1:
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
                                                filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index), amplitude))
                                                filea.write("\n")
                                            else:
                                                print '{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index), amplitude)
                                        matrix = np.identity(self._nocc)
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
                                                            filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                           it[2].multi_index),
                                                                                               amplitude))
                                                            filea.write("\n")
                                                        else:
                                                            print '{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                      it[2].multi_index),
                                                                                          amplitude)
                                                    matrix = np.identity(self._nocc)
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
                                                                        filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                                       it[2].multi_index,it[3].multi_index),
                                                                                                           amplitude))
                                                                        filea.write("\n")
                                                                    else:
                                                                        print '{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                                  it[2].multi_index,it[3].multi_index),
                                                                                                      amplitude)
                                                                matrix = np.identity(self._nocc)
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
                                                                    filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                                   it[2].multi_index,it[3].multi_index,
                                                                                                                                   it[4].multi_index),
                                                                                                       amplitude))
                                                                    filea.write("\n")
                                                                else:
                                                                    print '{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                              it[2].multi_index,it[3].multi_index,
                                                                                                                              it[4].multi_index),
                                                                                                  amplitude)
                                                            matrix = np.identity(self._nocc)
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
                                                                            filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                                           it[2].multi_index,it[3].multi_index,
                                                                                                                                           it[4].multi_index,it[5].multi_index),
                                                                                                               amplitude))
                                                                            filea.write("\n")
                                                                        else:
                                                                            print '{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                                      it[2].multi_index,it[3].multi_index,
                                                                                                                                      it[4].multi_index,it[5].multi_index),
                                                                                                          amplitude)
                                                                    matrix = np.identity(self._nocc)
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
                                                                                    filea.write('{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                                                      it[2].multi_index,it[3].multi_index,
                                                                                                                                                      it[4].multi_index,it[5].multi_index,
                                                                                                                                                      it[6].multi_index),
                                                                                                                          amplitude))
                                                                                    filea.write("\n")
                                                                                else:
                                                                                    print '{0:10} {1:20.16f}'.format(self.get_slater_determinant(it[0].multi_index,it[1].multi_index,
                                                                                                                                              it[2].multi_index,it[3].multi_index,
                                                                                                                                              it[4].multi_index,it[5].multi_index,
                                                                                                                                              it[6].multi_index),
                                                                                                                  amplitude)
                                                                            matrix = np.identity(self._nocc)
            if (amplitudestofile is True):
                filea.close()

    def get_slater_determinant(self, indices, indices2=None, indices3=None, indices4=None, indices5=None, indices6=None, indices7=None):
        '''Return excited Slater Determinant.
        '''
        ref_orb = indices[0]
        excited = self._nocc+indices[1]
        if indices2 is None:
            ref_orb2 = -1
            excited2 = -1
        else:
            ref_orb2 = indices2[0]
            excited2 = int(self._nocc+indices2[1])

        if indices3 is None:
            ref_orb3 = -1
            excited3 = -1
        else:
            ref_orb3 = indices3[0]
            excited3 = int(self._nocc+indices3[1])

        if indices4 is None:
            ref_orb4 = -1
            excited4 = -1
        else:
            ref_orb4 = indices4[0]
            excited4 = int(self._nocc+indices4[1])

        if indices5 is None:
            ref_orb5 = -1
            excited5 = -1
        else:
            ref_orb5 = indices5[0]
            excited5 = int(self._nocc+indices5[1])

        if indices6 is None:
            ref_orb6 = -1
            excited6 = -1
        else:
            ref_orb6 = indices6[0]
            excited6 = int(self._nocc+indices6[1])

        if indices7 is None:
            ref_orb7 = -1
            excited7 = -1
        else:
            ref_orb7 = indices7[0]
            excited7 = int(self._nocc+indices7[1])

        sd = []
        for i in range(self.nbasis):
            if i < self._nocc:
                if (i != ref_orb) and (i != ref_orb2) and (i != ref_orb3) and (i != ref_orb4) and (i != ref_orb5) and (i != ref_orb6) and (i != ref_orb7):
                    sd.append(2)
                else:
                    sd.append(0)
            else:
                if (i != excited) and (i != excited2) and (i != excited3) and (i != excited4) and (i != excited5) and (i != excited6) and (i != excited7):
                    sd.append(0)
                else:
                    sd.append(2)

        return str(sd).translate(None, ", ")

    def do_checkpoint(self, exps, olp):
        '''Dump orbitals for restart
        '''
        exps.coeffs.tofile('./orbitals.dat')
        olp._array.tofile('./overlap.dat')

    def print_final(self, s='Final'):
        '''Print energies
        '''
        log.hline('-')
        log('%s reference Energy:     %16.12f' %(s,self.get_reference_energy()+self.ecore))
        log('%s correlation Energy:   %16.12f' %(s,self.get_correlation_energy()))
        log.hline('=')
        log('%s total Energy:         %16.12f' %(s, (self.get_total_energy())))
        log.hline('=')

    def dump_final(self, exps, olp, printoptions, dumpci, checkpoint):
        '''Dump final solution (orbitals, expansion)
        '''
        if checkpoint > 0:
            if log.do_medium:
                log('Writing orbitals to file')
            olp._array.tofile('./overlap.dat')
            exps.coeffs.tofile('./orbitals.dat')
            if log.do_medium:
                log.hline('-')
                log('Final solution for coefficients: (only printed if |c_i| > %f)' %printoptions['ci'])
        self.print_solution(printoptions['ci'],
                           printoptions['excitationlevel'],
                           amplitudestofile=dumpci['amplitudestofile'],
                           filename=dumpci['amplitudesfilename'])
        if printoptions['geminal'] == True:
            log.hline('-')
            log('Geminal matrix: (rows=occupied, columns=virtual)')
            self.print_geminal_coefficients()
        if log.do_medium:
            log.hline('=')

    def print_geminal_coefficients(self):
        '''Print geminal coefficients
        '''
        nocc = self.geminal._array.shape[0]
        nvirt = self.geminal._array.shape[1]
        count = 1
        print "      ",
        for nvirti in range(nvirt):
            print "%9i" %(nvirti+1+nocc),
        print
        for line in self.geminal._array:
            print "%5i %s" %(count,":  "),
            for element in line:
                print "%9.6f" %element,
            print
            count = count+1

    def print_options(self, guess, solver, maxiter, thresh,
                      printoptions, indextrans):
        '''Print optimization options.
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
        if (guess['guess']=="random"):
            log('Initial guess:')
            log('  type:                      %s' %guess['guess'])
            log('  factor:                    %2.3f' %guess['factor'])
        elif (guess['guess']=="const"):
            log('Initial guess:')
            log('  type:                      %s' %guess['guess'])
            log('  factor:                    %2.3f' %guess['factor'])
        else:
            log('Initial guess:')
            log('  type:                      %s' %guess['guess'])
        log('Solvers:')
        log('  wavefunction amplitudes:   %s' %solver['wfn'])
        log('Optimization thresholds:')
        log('  wavefunction:              %1.2e' %thresh['wfn'])
        log('Printing options:')
        log('  c_i:                       %1.2e' %printoptions['ci'])
        log('  excitation level:          %i' %printoptions['excitationlevel'])

    def print_options_scf(self, guess, solver, maxiter, lshift, stepsearch,
                      thresh, printoptions, checkpoint, indextrans,
                      orbitaloptimizer):
        '''Print optimization options.
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
            log('  type:                        %s' %guess['guess'])
            log('  factor:                      %2.3f' %guess['factor'])
            log('Solvers:')
            log('  wavefunction amplitudes:     %s' %solver['wfn'])
            log('  Lagrange multipliers:        %s' %solver['lagrange'])
            log('  orbital optimization:        %s' %orbitaloptimizer)
            log('Number of iterations:          %i' %maxiter['orbiter'])
            log('Checkpointing:                 %i' %checkpoint)
            log('4-index transformation:        %s' %indextrans)
            log('Level shift:                   %3.3e' %lshift)
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
            elif stepsearch['method']=='wolfe':
                log('Apply line search:')
                log('  line search method:          %s' %stepsearch['method'])
                log('  initial step size:           %1.3f' %stepsearch['stepa'])
                log('  maximum step size:           %2.3f' %stepsearch['maxstep'])
                log('  contraction factor:          %1.3f' %stepsearch['downscale'])
                log('  c1 factor:                   %1.3e' %stepsearch['c1'])
                log('  c2 factor:                   %1.3e' %stepsearch['c2'])
            elif stepsearch['method']=='backtracking':
                log('Apply line search:')
                log('  line search method:          %s' %stepsearch['method'])
                log('  initial step size:           %1.3f' %stepsearch['stepa'])
                log('  contraction factor:          %1.3f' %stepsearch['downscale'])
                log('  c1 factor:                   %1.3e' %stepsearch['c1'])
                log('  minimum step size:           %1.3e' %stepsearch['minstep'])
            else:
                log('No linesearch selected:')
                log('  step size:                   %1.3f' %stepsearch['stepa'])
            log('Optimization thresholds:')
            log('  wavefunction:                %1.2e' %thresh['wfn'])
            log('  energy:                      %1.2e' %thresh['energy'])
            log('  gradient:                    %1.2e' %thresh['gradientmax'])
            log('  gradient norm:               %1.2e' %thresh['gradientnorm'])
            log('Printing options:')
            log('  c_i:                         %1.2e' %printoptions['ci'])
            log('  excitation level:            %i' %printoptions['excitationlevel'])

    def start_up(self, exps, swapa, givensrot=None):
        '''Update orbitals prior to optimization

           **Arguments:**

           swapa
               Orbital swapping according to list swapa.

           **Optional arguments:**

           givensrot
               Givens rotation of orbital pairs according
               to list givensrot containing orbital indices and angle
               (index1, index2, angle in deg).

        '''
        if swapa:
            if log.do_medium:
                log('Swap orbitals:')
            exps.apply_swapping(swapa)
        if givensrot:
            if log.do_medium:
                 log('Apply Givens rotation:')
            exps.apply_givens_rotation(givensrot)
