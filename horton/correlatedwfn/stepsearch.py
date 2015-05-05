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
'''Step search methods for orbital rotations

   Optimize step length for orbital optimization (``kappa``) and rotate orbitals.

   Works only with diagonal approximation to the Hessian which is stored as a
   OneIndex instance.
'''


import numpy as np

from horton.log import log, timer
from horton.orbital_utils import rotate_orbitals
from horton.correlatedwfn.trustregionopt import Dogleg, DoubleDogleg, TruncatedCG
from horton.utils import check_type, check_options

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
                :method: step search method used, one of ``trust-region``
                         (default), ``None``, ``backtracking``
                :alpha: scaling factor for step, used in ``backtracking`` and
                        ``None`` method (default 0.75)
                :c1: parameter used in ``backtracking`` (default 1e-4)
                :minalpha: minimum step length used in ``backracking``
                           (default 1e-6)
                :maxiterouter: maximum number of step search steps (default 10)
                :maxiterinner: maximum number of optimization steps in each
                               step search step (used only in pcg, default 500)
                :maxeta: upper bound for estimated vs actual change in
                         ``trust-region`` (default 0.75)
                :mineta: lower bound for estimated vs actual change in
                         ``trust-region`` (default 0.25)
                :upscale: scaling factor to increase trust radius in
                          ``trust-region`` (default 2.0)
                :downscale: scaling factor to decrease trust radius in
                            ``trust-region`` and step length in ``backtracking``
                            (float) (default 0.25)
                :trustradius: initial trust radius (default 0.75)
                :maxtrustradius: maximum trust radius (default 0.75)
                :threshold: trust-region optimization threshold, only used in
                            ``pcg`` method of ``trust-region``
                :optimizer: optimizes step to boundary of trust radius. One of
                            ``pcg``, ``dogleg``, ``ddl`` (default ddl)
        """
        self.lf = lf
        #
        # Check keywords and set default arguments, types and options are also
        # checked
        #
        names = []
        def _helper(x,y):
            names.append(x)
            kw.setdefault(x,y)
        _helper('method', 'trust-region')
        _helper('alpha', 1.0)
        _helper('c1', 0.0001)
        _helper('minalpha', 1e-6)
        _helper('maxiterouter', 10)
        _helper('maxiterinner', 500)
        _helper('maxeta', 0.75)
        _helper('mineta', 0.25)
        _helper('upscale', 2.0)
        _helper('downscale', 0.25)
        _helper('trustradius', 0.75)
        _helper('maxtrustradius', 0.75)
        _helper('threshold', 1e-8)
        _helper('optimizer', 'ddl')

        for name, value in kw.items():
            if name not in names:
                raise ValueError("Unknown keyword argument %s" % name)
            if value is not None:
                if value < 0:
                    raise ValueError('Cannot set attribute %s because of illegal value %s' %(name, value))
                setattr(self, name, kw[name])

        check_options('method', self.method, 'None', 'backtracking', 'trust-region')
        check_options('optimizer', self.optimizer, 'pcg', 'dogleg', 'ddl')
        check_type('alpha', self.alpha, int, float)
        check_type('c1', self.c1, int, float)
        check_type('minalpha', self.minalpha, int, float)
        check_type('maxiterouter', self.maxiterouter, int)
        check_type('maxiterinner', self.maxiterinner, int)
        check_type('maxeta', self.maxeta, int, float)
        check_type('mineta', self.mineta, int, float)
        check_type('upscale', self.upscale, float)
        check_type('downscale', self.downscale, float)
        check_type('trustradius', self.trustradius, float)
        check_type('maxtrustradius', self.maxtrustradius, float)
        check_type('threshold', self.threshold, float)

        self.alpha0 = self.alpha

        def _get_lf(self):
            return self.lf

        lf = property(_get_lf)

        def _get_method(self):
            return kw.pop('method')

        method = property(_get_method)

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

        def _get_alpha(self):
            return kw.pop('alpha')

        alpha = property(_get_alpha)

        def _get_minalpha(self):
            return kw.pop('minalpha')

        minalpha = property(_get_minalpha)

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

    def update_alpha(self, new):
        '''Update current scaling factor
        '''
        self.alpha = new

    def __call__(self, obj, one, two, orb, **kwargs):
        raise NotImplementedError


class RStepSearch(StepSearch):
    @timer.with_section('Linesearch')
    def __call__(self, obj, one, two, orb, **kwargs):
        '''Optimize Newton-Rapshon step.

           **Arguments:**

           obj
                A class instance containing the objective function.

           one, two
                One- and Two-electron integrals (TwoIndex and FourIndex/Cholesky
                instances)

           orb
                An Expansion instance. Contains the MO coefficients

           **Keywords:**

        '''
        if self.method == 'None':
            self.do_no_stepsearch(obj, one, two, orb, **kwargs)
        elif self.method == 'backtracking':
            self.do_backtracking(obj, one, two, orb, **kwargs)
        elif self.method == 'trust-region':
            self.do_trust_region(obj, one, two, orb, **kwargs)

    def do_no_stepsearch(self, obj, one, two, orb, **kwargs):
        '''Scale step size with factor self.alpha

           **Arguments:**

           obj
                A class instance containing the objective function.

           one, two
                One- and Two-electron integrals (TwoIndex and FourIndex/Cholesky
                instances)

           orb
                An Expansion instance. Contains the MO coefficients

           **Keywords:**
                :kappa: Initial step size (OneIndex instance)
        '''
        kappa = kwargs.get('kappa')

        #
        # Scale current optimization step, rotate orbitals
        #
        kappa.iscale(self.alpha)
        rotation = obj.compute_rotation_matrix(kappa)
        rotate_orbitals(orb, rotation)

        #
        # Solve for wfn/population model
        #
        try:
            obj.solve_model(one, two, orb, **kwargs)
        except:
            pass

    def do_backtracking(self, obj, one, two, orb, **kwargs):
        '''Backtracking line search.

           **Arguments:**

           obj
                A class instance containing the objctive function.

           one, two
                One- and Two-electron integrals (TwoIndex and FourIndex/Cholesky
                instances)

           orb
                An Expansion instance. Contains the MO coefficients

           **Keywords:**
                :kappa: Initial step size (OneIndex instance)
                :gradient: Orbital gradient (OneIndex instance)
        '''
        kappa = kwargs.get('kappa')
        gradient = kwargs.get('gradient')

        #
        # Update scaling factor to initial value
        #
        self.update_alpha(self.alpha0)

        #
        # Calculate objective function
        #
        ofun_ref = obj.compute_objective_function()

        #
        # Copy current orbitals
        #
        orb_ = orb.copy()
        #
        # Initial rotation
        #
        kappa.iscale(self.alpha)
        rotation = obj.compute_rotation_matrix(kappa)
        rotate_orbitals(orb_, rotation)

        #
        # Solve for wfn/population model if required
        #
        try:
            obj.solve_model(one, two, orb_, **kwargs)
        except:
            pass

        #
        # Calculate objective function for rotated orbitals
        #
        ofun = obj.compute_objective_function()

        #
        # reduce step size until line search condition is satisfied
        #
        while self.check_line_search_condition(ofun, ofun_ref, kappa, gradient):
            self.update_alpha(self.alpha*self.downscale)
            #
            # New rotation
            #
            kappa.iscale(self.alpha)
            rotation = obj.compute_rotation_matrix(kappa)
            orb_ = orb.copy()
            rotate_orbitals(orb_, rotation)
            #
            # Solve for wfn/population model if required
            #
            try:
                obj.solve_model(one, two, orb_, **kwargs)
            except:
                pass
            #
            # Calculate objective function for scaled rotation
            #
            ofun = obj.compute_objective_function()
            #
            # Abort if scaling factor falls below threshold
            #
            if self.alpha < self.minalpha:
                break
        orb.assign(orb_)

    def do_trust_region(self, obj, one, two, orb, **kwargs):
        '''Do trust-region optimization.

           **Arguments:**

           obj
                A class instance containing the objective function.

           one/two
                One and Two-body Hamiltonian (TwoIndex and FourIndex/Cholesky
                instances)

           orb
                (Expansion instance) An MO expansion

           **Keywords:**
                :kappa: Initial step size (OneIndex instance)
                :gradient: Orbital gradient (OneIndex instance)
                :hessian: Orbital Hessian (OneIndex instance)
        '''
        kappa = kwargs.get('kappa')
        gradient = kwargs.get('gradient')
        hessian = kwargs.get('hessian')

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
                if self.optimizer == 'pcg':
                    optimizer = TruncatedCG(self.lf, gradient, hessian, self.trustradius)
                    optimizer(**{'niter': self.maxiterinner, 'abstol': self.threshold})
                    stepn = optimizer.step
                #
                # Powell's dogleg optimization:
                #
                elif self.optimizer == 'dogleg':
                    optimizer = Dogleg(kappa, gradient, hessian, self.trustradius)
                    optimizer()
                    stepn = optimizer.step
                #
                # Powell's double dogleg optimization:
                #
                elif self.optimizer == 'ddl':
                    optimizer = DoubleDogleg(kappa, gradient, hessian, self.trustradius)
                    optimizer()
                    stepn = optimizer.step
            #
            # New rotation
            #
            rotation = obj.compute_rotation_matrix(stepn)
            orb_ = orb.copy()
            rotate_orbitals(orb_, rotation)

            #
            # Solve for wfn/population model if required
            #
            try:
                obj.solve_model(one, two, orb_, **kwargs)
            except:
                pass

            #
            # Calculate objective function after rotation
            #
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
                orb.assign(orb_)
                break
            elif rho >= self.mineta and rho <= self.maxeta and De <= 0.0:
                #
                # Do nothing with trust radius:
                #
                orb.assign(orb_)
                break
            elif rho > 0 and rho < self.mineta and De <= 0.0:
                #
                # Decrease trust radius:
                #
                self.update_trustradius(self.downscale*self.trustradius)
                orb.assign(orb_)
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
                orb.assign(orb_)
                break

    def check_line_search_condition(self, cetot, petot, kappa, grad):
        '''Check if Armijo condition is satisfied
           (e_tot(0+alpha*kappa)-e_tot(0) <= c1*alpha*(kappa,grad))
        '''
        term = self.c1*self.alpha*kappa.dot(grad)
        return (cetot - petot - term) > 0
