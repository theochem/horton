# -*- coding: utf-8 -*-
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
'''Linesearch Methods'''


import numpy as np
from scipy import optimize as opt
from scipy import linalg as linalg
import math as math
from copy import deepcopy
from random import random

from horton.log import log, timer
from horton.correlatedwfn.trustregionopt import Dogleg, TruncatedCG, LevelShiftedNewton

__all__ = [
    'StepSearch', 'RStepSearch',
]


class StepSearch(object):
    def __init__(self, **kw):
        """
           **Arguments:**


           **Keywords:**
                :method: type of search method
                :optimizer: optimization method for line search
                :maxstep: maximum step size
                :maxeta: upper bound for trust region
                :mineta: lower bound for trust region
                :upscale: scale factor to increase trust radius
                :downscale: scale factor to decrease trust radius
                :maxiterouter: maxiter for outer iterations
                :maxiterinner: maxiter for inner iterations
                :stepa:
                :stepb:
                :minstep: minimum step size
                :c1:
                :c2:
                :trustradius: initial trust radius
                :maxtrustradius: maximum trust radius
                :threshold:
        """
        names = ["method", "maxstep", "maxeta", "mineta",
                 "upscale", "downscale", "maxiterouter", "maxiterinner",
                 "stepa", "stepb", "c1", "c2", "minstep", "trustradius",
                 "threshold", "optimizer", "maxtrustradius"]
        for name, value in kw.items():
            if name not in names:
                raise ValueError("Unknown keyword argument %s" % name)
            if value is not None:
                setattr(self, name, kw[name])

        def _get_method(self):
            return kw.pop('method')

        method = property(_get_method)

        def _get_maxstep(self):
            return kw.pop('maxstep')

        maxstep = property(_get_maxstep)

        def _get_maxeta(self):
            return kw.pop('maxeta')

        maxeta = property(_get_maxeta)

        def _get_mineta(self):
            return kw.pop('mineta')

        mineta = property(_get_mineta)

        def _get_upscale(self):
            return kw.pop('upscale')

        upscale = property(_get_upscale)

        def _get_downscale(self):
            return kw.pop('downscale')

        downscale = property(_get_downscale)

        def _get_maxiterouter(self):
            return kw.pop('maxiterouter')

        maxiterouter = property(_get_maxiterouter)

        def _get_maxiterinner(self):
            return kw.pop('maxiterinner')

        maxiterinner = property(_get_maxiterinner)

        def _get_stepa(self):
            return kw.pop('stepa')

        stepa = property(_get_stepa)

        def _get_minstep(self):
            return kw.pop('minstep')

        minstep = property(_get_minstep)

        def _get_c1(self):
            return kw.pop('c1')

        c1 = property(_get_c1)

        def _get_c2(self):
            return kw.pop('c2')

        c2 = property(_get_c2)

        def _get_maxtrustradius(self):
            return kw.pop('maxtrustradius')

        maxtrustradius = property(_get_maxtrustradius)

        def _get_trustradius(self):
            return kw.pop('trustradius')

        trustradius = property(_get_trustradius)

        def _get_threshold(self):
            return kw.pop('threshold')

        threshold = property(_get_threshold)

        def _get_optimizer(self):
            return kw.pop('optimizer')

        optimizer = property(_get_optimizer)

    def update_trustradius(self, new):
        '''Update current trastradius
        '''
        self.trustradius = new

    def update_stepa(self, new):
        '''Update current stepsize (for alpha-spin)
        '''
        self.stepa = new

    def __call__(self, obj, one, two, exps, **kwargs):
        raise NotImplementedError


class RStepSearch(StepSearch):
    @timer.with_section('Linesearch')
    def __call__(self, obj, one, two, exps, **kwargs):
        '''Optimize Newton-Rapshon step.

           **Arguments:**


           **Keywords:**

           kappa
                Initial stepsize.

           gradient, gradientb, hessian, hessianb
                Gradient and Hessian to recalculate step for alpha and beta electrons

        '''
        if self.method == 'None':
            self.do_no_linesearch(obj, one, two, exps, **kwargs)
        elif self.method == 'backtracking':
            self.do_backtracking(obj, one, two, exps, **kwargs)
        elif self.method == 'trust-region':
            self.do_trust_region(obj, one, two, exps, **kwargs)
        else:
            raise NotImplementedError

    def do_no_linesearch(self, obj, one, two, exps, **kwargs):
        '''Scale step size with factor self.stepa
        '''
        kappa = kwargs.get('kappa')

        rotation = obj.get_rotation_matrix((kappa*self.stepa))
        # Transform MOs:
        obj.rotate_orbitals(exps, rotation)

        obj.solve_wfn(one, two, exps, **kwargs)

    def do_backtracking(self, obj, one, two, exps, **kwargs):
        '''Backtracking line search.

           **Arguments:**

           obj
                A wfn instance.

           one/two
                Two and Four-index Hamiltonian.

           exps
                A AO/MO expansion.

        '''
        kappa = kwargs.get('kappa')
        gradient = kwargs.get('gradient', None)

        self.update_stepa(self.maxstep)
        totale_old = obj.get_total_energy()
        # Form transformation matrix for AOs:
        orb = exps.copy()
        rotation = obj.get_rotation_matrix((kappa*self.stepa))
        # Transform MOs:
        obj.rotate_orbitals(orb, rotation)

        obj.solve_wfn(one, two, orb, **kwargs)

        # Get energies
        totale = obj.get_total_energy()

        while self.check_line_search_condition(totale, totale_old, kappa, gradient):
            self.update_stepa(self.stepa*self.downscale)
            # Form transformation matrix for AOs:
            rotation = obj.get_rotation_matrix((kappa*self.stepa))
            # Transform MOs:
            orb = exps.copy()
            obj.rotate_orbitals(orb, rotation)

            obj.solve_wfn(one, two, orb, **kwargs)

            # Get energies
            totale = obj.get_total_energy()

            if self.stepa < self.minstep:
                exps.assign(orb)
                break
        exps.assign(orb)

    def do_trust_region(self, obj, one, two, exps, **kwargs):
        '''Do trust-region optimization.

           **Arguments:**

           obj
                A wfn instance.

           one/two
                Two and Four-index Hamiltonian.

           exps
                A AO/MO expansion.

        '''
        kappa = kwargs.get('kappa')
        gradient = kwargs.get('gradient', None)
        hessian = kwargs.get('hessian', None)

        tmp_orb = exps.copy()
        stepn = deepcopy(kappa)
        iteri = 1
        totale_ref = obj.get_total_energy()
        De = 0.0

        while True:
            norm = math.sqrt(np.dot(stepn, stepn))
            # If ||kappa||_2 is outside the trust region, find a new Newton
            # step inside the trust region:
            if norm > self.trustradius:
                # Preconditioned conjugate gradient
                if self.optimizer in ['pcg']:
                    optimizer = TruncatedCG(gradient, hessian, self.trustradius)
                    optimizer(**{'niter': self.maxiterinner, 'abstol': self.threshold})
                    stepn = optimizer.step
                # Powell's dogleg optimization:
                elif self.optimizer in ['dogleg']:
                    optimizer = Dogleg(kappa, gradient, hessian, self.trustradius)
                    optimizer()
                    stepn = optimizer.step
                # Levelshifted Newton equations:
                elif self.optimizer in ['lsn']:
                    optimizer = LevelShiftedNewton(gradient, hessian, self.trustradius)
                    optimizer()
                    stepn = optimizer.step
                else:
                    raise NotImplementedError
            # Form transformation matrix for AOs:
            rotation = obj.get_rotation_matrix(stepn)
            # Transform MOs:
            tmp_orb = exps.copy()
            obj.rotate_orbitals(tmp_orb, rotation)

            obj.solve_wfn(one, two, tmp_orb, **kwargs)

            # Get energies
            totale = obj.get_total_energy()

            # Determine ratio for a given step:
            De = totale-totale_ref
            Destimate = (np.dot(gradient, stepn)+0.5*np.dot(stepn, (hessian*stepn)))
            rho = De/Destimate
            small = False
            if abs(De) < 1e-8 and abs(Destimate) < 1e-8:
                small = True
            else:
                small = False

            # For a given ratio, change trust radius accordingly and accept or reject Newton step:
            if small and De <= 0.0:
                new = min(self.upscale*self.trustradius,self.maxtrustradius)
                self.update_trustradius(new)
                exps.assign(tmp_orb)
                break
            elif rho > self.maxeta and De <= 0.0:
                # Enlarge trust radius:
                new = min(self.upscale*self.trustradius,self.maxtrustradius)
                self.update_trustradius(new)
                exps.assign(tmp_orb)
                break
            elif rho >= self.mineta and rho <= self.maxeta and De <= 0.0:
                # Do nothing with trust radius:
                exps.assign(tmp_orb)
                break
            elif rho > 0 and rho < self.mineta and De <= 0.0:
                # Decrease trust radius:
                self.update_trustradius(self.downscale*self.trustradius)
                exps.assign(tmp_orb)
                break
            else:
                # Bad step! Reject Newton step and repeat with smaller trust region
                self.update_trustradius(self.downscale*self.trustradius)

            iteri = iteri+1
            if iteri > self.maxiterouter:
                log('Warning: Trust region search not converged after %i iterations. Trust region search aborted.' %self.maxiterouter)
                self.update_stepa(math.sqrt(np.dot(stepn, stepn)))
                exps.assign(tmp_orb)
                break
        self.update_stepa(math.sqrt(np.dot(stepn, stepn)))

    def check_line_search_condition(self, cetot, petot, kappa, grad):
        '''
        '''
        term = self.c1*self.stepa*np.dot(kappa, grad)
        if ((cetot - petot - term) <= 0):
            return False
        return True
