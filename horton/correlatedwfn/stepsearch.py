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
'''Linesearch Methods

   Optimize step length for orbital optimization
'''


import numpy as np

from horton.log import log, timer
from horton.orbital_utils import rotate_orbitals
from horton.correlatedwfn.trustregionopt import Dogleg, DoubleDogleg, TruncatedCG

__all__ = [
    'StepSearch',
    'RStepSearch',
]


class StepSearch(object):
    def __init__(self, lf, **kw):
        """
           **Arguments:**

           lf
                A LinalgFactory instance

           **Keywords:**
                :method: linesearch method used, one of ``trust-region``
                         (default), ``None``, ``backtracking``
                :stepa: scaling factor for step, used in ``backtracking`` and
                        ``None`` method (default 0.75)
                :c1: parameter used in ``backtracking`` (default 1e-4)
                :maxstep: maximum step length used in ``backtracking``
                          (default 0.75)
                :minstep: minimum step length used in ``backracking``
                          (default 1e-6)
                :maxiterouter: maximum number of line search steps (default 10)
                :maxiterinner: maximum number of optimization steps in each
                               line search step (used only in pcg, default 500)
                :maxeta: upper bound for estimated vs actual change in
                         ``trust-region`` (default 0.75)
                :mineta: lower bound for estimated vs actual change in
                         ``trust-region`` (default 0.25)
                :upscale: scaling factor to increase trustradius in
                          ``trust-region`` (default 2.0)
                :downscale: scaling factor to decrease trustradius in
                            ``trust-region`` (default 0.25)
                :trustradius: initial trustradius (default 0.75)
                :maxtrustradius: maximum trustradius (default 0.75)
                :threshold: trust-region optimization threshold, only used in
                            ``pcg`` method of ``trust-region``
                :optimizer: optimizes step to boundary of trustradius. One of
                            ``pcg``, ``dogleg``, ``ddl`` (default ddl)
        """
        self.lf = lf
        names = ["method", "maxstep", "maxeta", "mineta",
                 "upscale", "downscale", "maxiterouter", "maxiterinner",
                 "stepa", "stepb", "c1", "minstep", "trustradius",
                 "threshold", "optimizer", "maxtrustradius"]
        for name, value in kw.items():
            if name not in names:
                raise ValueError("Unknown keyword argument %s" % name)
            if value is not None:
                setattr(self, name, kw[name])

        def _get_lf(self):
            return self.lf

        lf = property(_get_lf)

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

           obj
                A class instance containing the wfn.

           one, two
                One- and Two-electron integrals

           exps
                An Expansion instance. Contains the AO/MO coefficients

           **Keywords:**

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

           **Arguments:**

           obj
                A class instance containing the wfn.

           one, two
                One- and Two-electron integrals

           exps
                An Expansion instance. Contains the AO/MO coefficients

           **Keywords:**
                :kappa: Initial step size

        '''
        kappa = kwargs.get('kappa')

        #
        # Scale current optimization step, rotate orbitals
        #
        kappa.iscale(self.stepa)
        rotation = obj.compute_rotation_matrix(kappa)
        rotate_orbitals(exps, rotation)

        #
        # Solve for wfn
        #
        obj.solve_wfn(one, two, exps, **kwargs)

    def do_backtracking(self, obj, one, two, exps, **kwargs):
        '''Backtracking line search.

           **Arguments:**

           obj
                A class instance containing the wfn.

           one, two
                One- and Two-electron integrals

           exps
                An Expansion instance. Contains the AO/MO coefficients

           **Keywords:**
                :kappa: Initial step size
                :gradient: Orbital gradient
        '''
        kappa = kwargs.get('kappa')
        gradient = kwargs.get('gradient', None)

        self.update_stepa(self.maxstep)
        ofun_ref = obj.compute_objective_function()
        #
        # Copy current orbitals
        #
        orb = exps.copy()
        #
        # Initial rotation
        #
        kappa.iscale(self.stepa)
        rotation = obj.compute_rotation_matrix(kappa)
        rotate_orbitals(orb, rotation)

        #
        # Solve for wfn
        #
        obj.solve_wfn(one, two, orb, **kwargs)

        ofun = obj.compute_objective_function()

        #
        # reduce step size until line search condition is satisfied
        #
        while self.check_line_search_condition(ofun, ofun_ref, kappa, gradient):
            self.update_stepa(self.stepa*self.downscale)
            #
            # New rotation
            #
            kappa.iscale(self.stepa)
            rotation = obj.compute_rotation_matrix(kappa)
            orb = exps.copy()
            rotate_orbitals(orb, rotation)

            #
            # Solve for wfn
            #
            obj.solve_wfn(one, two, orb, **kwargs)

            ofun = obj.compute_objective_function()

            #
            # Abort if scaling factor falls below threshold
            #
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
                One and Two-body Hamiltonian.

           exps
                A AO/MO expansion.


           **Keywords:**
                :kappa: Initial step size
                :gradient: Orbital gradient
                :hessian: Orbital Hessian
        '''
        kappa = kwargs.get('kappa')
        gradient = kwargs.get('gradient', None)
        hessian = kwargs.get('hessian', None)

        iteri = 1
        ofun_ref = obj.compute_objective_function()
        De = 0.0

        stepn = kappa.copy()
        while True:
            norm = stepn.norm()
            #
            # If ||kappa||_2 is outside the trust region, find a new Newton
            # step inside the trust region:
            #
            if norm > self.trustradius:
                #
                # Preconditioned conjugate gradient
                #
                if self.optimizer in ['pcg']:
                    optimizer = TruncatedCG(self.lf, gradient, hessian, self.trustradius)
                    optimizer(**{'niter': self.maxiterinner, 'abstol': self.threshold})
                    stepn = optimizer.step
                #
                # Powell's dogleg optimization:
                #
                elif self.optimizer in ['dogleg']:
                    optimizer = Dogleg(kappa, gradient, hessian, self.trustradius)
                    optimizer()
                    stepn = optimizer.step
                #
                # Powell's double dogleg optimization:
                #
                elif self.optimizer in ['ddl']:
                    optimizer = DoubleDogleg(kappa, gradient, hessian, self.trustradius)
                    optimizer()
                    stepn = optimizer.step
                else:
                    raise NotImplementedError
            #
            # New rotation
            #
            rotation = obj.compute_rotation_matrix(stepn)
            orb = exps.copy()
            rotate_orbitals(orb, rotation)

            #
            # Solve for wfn
            #
            obj.solve_wfn(one, two, orb, **kwargs)

            ofun = obj.compute_objective_function()

            #
            # Determine ratio for a given step:
            #
            hstep = hessian.new()
            hessian.mult(stepn, hstep)
            Destimate = gradient.dot(stepn)+0.5*stepn.dot(hstep)
            De = ofun-ofun_ref
            rho = De/Destimate

            #
            # For a given ratio, adjust trust radius and accept or reject
            # Newton step:
            #
            if rho > self.maxeta and De <= 0.0:
                #
                # Enlarge trust radius:
                #
                new = min(self.upscale*self.trustradius,self.maxtrustradius)
                self.update_trustradius(new)
                exps.assign(orb)
                break
            elif rho >= self.mineta and rho <= self.maxeta and De <= 0.0:
                #
                # Do nothing with trust radius:
                #
                exps.assign(orb)
                break
            elif rho > 0 and rho < self.mineta and De <= 0.0:
                #
                # Decrease trust radius:
                #
                self.update_trustradius(self.downscale*self.trustradius)
                exps.assign(orb)
                break
            else:
                #
                # Bad step! Reject Newton step and repeat with smaller trust
                # radius
                #
                self.update_trustradius(self.downscale*self.trustradius)

            iteri = iteri+1
            if iteri > self.maxiterouter:
                log('Warning: Trust region search not converged after %i iterations. Trust region search aborted.' %self.maxiterouter)
                self.update_stepa(stepn.norm())
                exps.assign(orb)
                break
        self.update_stepa(stepn.norm())

    def check_line_search_condition(self, cetot, petot, kappa, grad):
        '''
        '''
        term = self.c1*self.stepa*kappa.dot(grad)
        if ((cetot - petot - term) <= 0):
            return False
        return True
