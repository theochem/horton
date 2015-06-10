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
'''Optimization of Newton step

   Module contains different trust radius optimization algorithms
'''


import numpy as np

__all__ = [
    'Dogleg',
    'DoubleDogleg',
    'TruncatedCG',
]


class Dogleg():
    def __init__(self, pn, g, H, radius):
        '''
        Powell’s single dogleg steps approximate the trust region step over
        the two dimensional subspace spanned by the Cauchy step and any Newton
        step. Depending on the size of current step length, the step is
        adjusted yielding one of the following steps:

          * if p(newton) < radius:   Do full Newton step (not calculated here)
          * elif p(cauchy) > radius: Do Cauchy step
          * else:                    Do dogleg step


        ** Arguements **

        pn
            The full Newton step

        g
            The gradient, a OneIndex instance

        H
            The Hessian, a OneIndex instance

        radius
            The current trustradius

        ** Returns **

          :step:       final step, a OneIndex instance
          :stepNorm:   Euclidian norm of the step

        '''
        self.pn = pn
        self.H = H
        self.g = g
        self.Delta = radius

        self.step = None
        self.stepNorm = 0.0

    def __call__(self):
        """The dogleg algorithm."""

        #
        # Calculate Cauchy step pc.
        #
        tmp = self.H.mult(self.g)
        curv = self.g.dot(tmp)
        factor = -self.g.dot(self.g)/curv
        pc = self.g.copy()
        pc.iscale(factor)
        dotpc = pc.dot(pc)
        normpc = pc.norm()

        #
        # Test if the Cauchy step pc exits the trust region.
        #
        if normpc >= self.Delta:
            self.step = pc
            self.step.iscale(self.Delta/normpc)
            self.stepNorm = np.sqrt(np.dot(self.step, self.step))
            return

        #
        # Otherwise do a dogleg step.
        # Find the solution to the scalar quadratic equation.
        #
        d_pn_pc = self.pn.copy()
        d_pn_pc.iadd(pc, -1.0)
        dot_d_pn_pc = d_pn_pc.dot(d_pn_pc)
        dot_pc_d_pn_pc = self.pn.dot(d_pn_pc)
        term = dot_pc_d_pn_pc*dot_pc_d_pn_pc-dot_d_pn_pc*(dotpc-self.Delta*self.Delta)
        tau = (-dot_pc_d_pn_pc+np.sqrt(term))/dot_d_pn_pc

        #
        # Return Dogleg step
        #
        self.step = pc
        self.step.iadd(d_pn_pc, tau)
        self.stepNorm = self.step.norm()
        return


class DoubleDogleg():
    def __init__(self, pn, g, H, radius):
        '''
        Double dogleg steps approximate the trust region step

          * if ||p(g)|| >= radius:      Do step in g direction
          * elif c*||p(n)|| >= radius:  Do ddl step
          * elif ||p(n)|| >= radius:    Do scaled Newton step
          * else:                       Do Newton step


        ** Arguements **

        pn
            The full Newton step

        g
            The gradient, a OneIndex instance

        H
            The Hessian, a OneIndex instance

        radius
            The current trustradius

        ** Returns **

          :step:       final step, a OneIndex instance
          :stepNorm:   Euclidian norm of the step

        '''
        self.pn = pn
        self.H = H
        self.g = g
        self.Delta = radius

        self.step = None
        self.stepNorm = 0.0

    def __call__(self):
        """The double dogleg algorithm."""

        #
        # Calculate gradient step pg.
        #
        tmp = self.H.mult(self.g)
        curv = self.g.dot(tmp)
        factor = -self.g.dot(self.g)/curv
        pg = self.g.copy()
        pg.iscale(factor)

        #
        # Calculate norms.
        #
        normpg = pg.norm()
        normpn = self.pn.norm()

        #
        # Calculate gamma.
        #
        gamma = pg.dot(pg)/pg.dot(self.pn)

        if self.Delta <= normpg:
            #
            # p = Delta/norm(g)*g
            #
            pg.iscale(self.Delta/normpg)
            self.step = pg
            self.stepNorm = self.step.norm()
        elif self.Delta <= (gamma*normpn):
            #
            # p = a*p(g) + (1-a)*gamma*p(n)
            #
            a = gamma*gamma*normpn*normpn-normpg*normpg
            b = self.Delta*self.Delta-normpg*normpg
            alpha = (a-b)/(a+np.sqrt(a*b))

            pg.iscale(alpha)
            self.pn.iscale((1-alpha))
            self.step = pg.copy()
            self.step.iadd(self.pn)
            self.stepNorm = self.step.norm()
        elif self.Delta <= normpn:
            #
            # p = Delta*p(g)/norm(p(n))
            #
            self.pn.iscale(self.Delta/normpn)
            self.step = self.pn.copy()
            self.stepNorm = self.step.norm()
        else:
            #
            # p = p(n)
            #
            self.step = self.pn.copy()
            self.stepNorm = self.step.norm()

        return


class TruncatedCG:
    def __init__(self, lf, g, H, radius, **kwargs):
        '''
        Solve the quadratic trust-region subproblem

          minimize    < g, s > + 1/2 < s, Hs >
          subject to  < s, s >  <=  radius

        by means of the truncated conjugate gradient algorithm as described in
        Trond Steihaug, SIAM Journal on Numerical Analysis (1983), 20, 626-637.

        The algorithm stops as soon as the preconditioned norm of the gradient
        falls under

            max( abstol, reltol * g0 )

        where g0 is the preconditioned norm of the initial gradient (or the
        Euclidian norm if no preconditioner is given), or as soon as the
        iterates cross the boundary of the trust region.

        ** Arguements **

        lf
            A linalg factory instance

        g
            The gradient, a OneIndex instance

        H
            The Hessian, a OneIndex instance

        ** Returns **

          :step:       final step
          :niter:      number of iterations
          :stepNorm:   Euclidian norm of the step

        '''
        self.lf = lf
        self.H = H
        self.g = g
        self.n = g.shape[0]
        self.Delta = radius

        self.step = None
        self.stepNorm = 0.0
        self.niter = 0

    def __call__(self, **kwargs):
        '''
        Solve the trust-region subproblem.

        ** Keywords **

          :s0:         initial guess (default: [0,0,...,0]),
          :abstol:     absolute stopping tolerance (default: 1.0e-8),
          :reltol:     relative stopping tolerance (default: 1.0e-6),
          :maxiter:    maximum number of iterations (default: 2n),
          :prec:       a user-defined preconditioner.
        '''
        abstol  = kwargs.get('absol', 1.0e-8)
        reltol  = kwargs.get('reltol', 1.0e-6)
        maxiter = kwargs.get('maxiter', 2*self.n)
        prec    = kwargs.get('prec', lambda v: v)

        #
        # Initialization (Step 1)
        #
        r = self.g.copy()
        if 's0' in kwargs:
            s = kwargs['s0']
            snorm2 = s.norm()
            tmp = self.H.mul(s)
            r.iadd(tmp)
        else:
            s = self.lf.create_one_index(self.g.shape[0])
            snorm2 = 0.0

        y = prec(r).copy()
        ry = r.dot(y)
        sqrtry = np.sqrt(ry)

        #
        # Initialize r as a copy of g not to alter the original g
        #
        p = y.copy()                       # p = - preconditioned residual
        p.iscale(-1)
        i = 0

        tolerance = max(abstol, reltol*sqrtry)
        maxThresh = True
        maxIter = True
        pastBoundary = True

        #
        # while sqrtry > stopTol and k < maxiter
        #
        while (maxThresh or maxIter) and pastBoundary:

            i += 1
            #
            # Step 2:
            #
            Hp = self.H.mult(p)
            dot_pHp = p.dot(Hp)

            #
            # Compute steplength to the boundary.
            #
            sigma = self.to_boundary(s, p, self.Delta, ss=snorm2)

            #
            # Step 3:
            # Compute CG steplength.
            #
            alpha = ry/dot_pHp

            #
            # (p, Hp) < 0, go to boundary
            #
            if (dot_pHp <= 0 or alpha > sigma):
                #
                # p leads past the trust-region boundary. Move to the boundary.
                #
                s.iadd(p, sigma)
                snorm2 = s.norm()
                pastBoundary = False
                continue

            #
            # Step 4:
            # Move to next iteration step.
            #
            s.iadd(p, alpha)

            #
            # Step 5:
            #
            r.iadd(Hp, alpha)
            y = prec(r).copy()
            ry_next = r.dot(y)
            beta = ry_next/ry
            p.iscale(beta)
            p.iadd(y, -1.0)
            ry = ry_next
            sqrtry = np.sqrt(ry)
            snorm2 = s.dot(s)

            if i >= maxiter:
                maxIter = False
            if sqrtry <= tolerance:
                maxThresh = False

        #
        # Output info about the last iteration.
        #
        self.step = s
        self.niter = i
        self.stepNorm = np.sqrt(snorm2)
        return

    def to_boundary(self, s, p, radius, ss=None):
        '''
        Given vectors `s` and `p` and a trust-region radius `radius` > 0,
        return the positive scalar `sigma` such that

           || s + sigma * p || = radius

        in Euclidian norm. If known, supply optional argument `ss` whose value
        should be the squared Euclidian norm of `s`.
        '''
        sp = s.dot(p)
        pp = p.dot(p)
        if not ss:
            ss = s.dot(s)
        sigma = (-sp + np.sqrt(abs(sp*sp + pp * (radius*radius - ss))))
        sigma = sigma/pp
        return sigma
