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
'''Optimization of the Newton step'''


import numpy as np
from math import sqrt

__all__ = [
    'Dogleg',
    'TruncatedCG',
    'LevelShiftedNewton',
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
            The gradient

        H
            The Hessian. The Hessian is assumed to be a vector!

        radius
            The current trustradius

        ** Returns **

          :step:       final step,
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

        # Calculate Cauchy step pc.
        curv = np.dot(self.g, self.H*self.g)
        pc = -np.dot(self.g, self.g)/curv*self.g
        dotpc = np.dot(pc, pc)
        normpc = sqrt(dotpc)

        # Test if the Cauchy step pc exits the trust region.
        if normpc >= self.Delta:
            self.step = self.Delta*pc/normpc
            self.stepNorm = sqrt(np.dot(self.step, self.step))
            return

        # Otherwise do a dogleg step.
        # Find the solution to the scalar quadratic equation.
        d_pn_pc = self.pn-pc
        dot_d_pn_pc = np.dot(d_pn_pc, d_pn_pc)
        dot_pc_d_pn_pc = np.dot(self.pn, d_pn_pc)
        term = dot_pc_d_pn_pc*dot_pc_d_pn_pc-dot_d_pn_pc*(dotpc-self.Delta*self.Delta)
        tau = (-dot_pc_d_pn_pc+sqrt(term))/dot_d_pn_pc

        # Return Dogleg step
        self.step = pc+tau*d_pn_pc
        self.stepNorm = sqrt(np.dot(self.step, self.step))
        return


class TruncatedCG:
    def __init__(self, g, H, radius, **kwargs):
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

        g
            The gradient

        H
            The Hessian. The Hessian is assumed to be a vector!

        radius
            The current trustradius

        ** Returns **

          :step:       final step
          :niter:      number of iterations
          :stepNorm:   Euclidian norm of the step

        '''
        self.H = H
        self.g = g
        self.n = len(g)
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

        # Initialization (Step 1)
        r = self.g.copy()
        if 's0' in kwargs:
            s = kwargs['s0']
            snorm2 = np.linalg.norm(s)
            r += self.H*s
        else:
            s = np.zeros(self.n)
            snorm2 = 0.0

        y = prec(r)
        ry = np.dot(r, y)
        sqrtry = sqrt(ry)

        # Initialize r as a copy of g not to alter the original g
        p = -y                       # p = - preconditioned residual
        i = 0

        tolerance = max(abstol, reltol*sqrtry)
        maxThresh = True
        maxIter = True
        pastBoundary = True

        #while sqrtry > stopTol and k < maxiter
        while (maxThresh or maxIter) and pastBoundary:

            i += 1
            # Step 2:
            Hp = self.H*p
            dot_pHp = np.dot(p, Hp)

            # Compute steplength to the boundary.
            sigma = self.to_boundary(s, p, self.Delta, ss=snorm2)

            # Step 3:
            # Compute CG steplength.
            alpha = ry/dot_pHp

            # (p, Hp) < 0, go to boundary
            if (dot_pHp <= 0 or alpha > sigma):
                # p leads past the trust-region boundary. Move to the boundary.
                s += sigma*p
                snorm2 = sqrt(np.dot(s, s))
                pastBoundary = False
                continue

            # Step 4:
            # Move to next iteration step.
            s += alpha*p

            # Step 5:
            r += alpha*Hp
            y = prec(r)
            ry_next = np.dot(r, y)
            beta = ry_next/ry
            p = -y + beta * p
            ry = ry_next
            sqrtry = sqrt(ry)
            snorm2 = np.dot(s,s)

            if i >= maxiter:
                maxIter = False
            if sqrtry <= tolerance:
                maxThresh = False

        # Output info about the last iteration.
        self.step = s
        self.niter = i
        self.stepNorm = sqrt(snorm2)
        return

    def to_boundary(self, s, p, radius, ss=None):
        '''
        Given vectors `s` and `p` and a trust-region radius `radius` > 0,
        return the positive scalar `sigma` such that

           || s + sigma * p || = radius

        in Euclidian norm. If known, supply optional argument `ss` whose value
        should be the squared Euclidian norm of `s`.
        '''
        sp = np.dot(s,p)
        pp = np.dot(p,p)
        if not ss:
            ss = np.dot(s,s)
        sigma = (-sp + sqrt(abs(sp*sp + pp * (radius*radius - ss))))
        sigma = sigma/pp
        return sigma


class LevelShiftedNewton:
    def __init__(self, g, H, radius, **kwargs):
        '''
        Solve level shifted Newton equations

        The algorithm stops as soon as the residual norm is smaller than the
        residual norm of R1 by a factor of k

        ** Arguements **

        g
            The gradient

        H
            The Hessian. The Hessian is assumed to be a vector!

        radius
            The current trustradius

        ** Returns **

          :step:       final step
          :niter:      number of iterations

        '''
        self.H = H
        self.g = g
        self.n = len(g)
        self.Delta = radius

        self.step = None
        self.niter = 0

    def __call__(self, **kwargs):
        '''
        Solve the trust-region subproblem.

        ** Keywords **

          :s0:         initial guess (default: [0,0,...,0]),
          :tol:        absolute stopping tolerance (default: 5.0e-3),
          :maxiter:    maximum number of iterations (default: 20),
          :prec:       a user-defined preconditioner.
        '''
        tol  = kwargs.get('tol', 5.0e-2)
        maxiter = kwargs.get('maxiter', 30)
        prec    = kwargs.get('prec', lambda v, w: np.divide(v, (self.H-w)))
        mu = np.amin(self.H)
        b = []

        # Set up first trial vectors
        b0 = np.zeros(self.n)
        normg = np.linalg.norm(self.g)
        gnorm = self.g/normg
        b.append(gnorm)
        b.append(self.residual(0.0, gnorm))
        tmp = np.zeros(self.n)
        tmp[np.argmin(self.H)] = np.amin(self.H)
        b.append(tmp/np.sqrt(2))
        b = self.gramschmidt(*b)

        maxThresh = True
        maxIter = True
        alpha = normg
        r1 = self.residual(0.0, gnorm)
        if self.Delta < tol:
            tol = 5e-3

        it = 0
        #while sqrtry > stopTol and k < maxiter
        while maxIter and maxThresh:

            # Construct augmented Hessian
            augH = np.zeros((4+it,4+it))
            augH[1, 0] = alpha*normg
            augH[0, 1] = augH[1, 0]
            for i in range(1, 4+it):
                for j in range(1, 4+it):
                    augH[i,j] = self.augH(b[i-1], b[j-1])

            eig, evec = np.linalg.eigh(augH)
            minind = np.argmin(eig)
            mu = eig[minind]
            s = b0
            for i in range(len(b)):
                s += evec[i+1, minind]*b[i]
            alpha = self.get_alpha(s)

            # New trial vector
            rn = self.residual(mu, s, alpha)
            prn = prec(rn, mu)
            b.append(prn)

            # Perform Gram-Schmidt
            b = self.gramschmidt(*b)

            it += 1
            maxIter = (it<maxiter)
            maxThresh = np.linalg.norm(rn)>(np.linalg.norm(r1)*tol)
        alpha = self.get_alpha(s)
        # Output info about the last iteration.
        self.step = s*alpha
        self.niter = it
        return

    def residual(self, mu, step, alpha=1.0):
        '''
        '''
        return (-self.g-(self.H-mu)*step*alpha)

    def get_alpha(self, vec):
        return (self.Delta/np.linalg.norm(vec))

    def gs_cofficient(self, v1, v2):
        return np.dot(v2, v1)/np.dot(v1, v1)

    def multiply(self, cofficient, v):
        return map((lambda x : x*cofficient), v)

    def proj(self, v1, v2):
        return self.multiply(self.gs_cofficient(v1, v2), v1)

    def normalize(self, *vecs):
        result = []
        for vec in vecs:
             result.append(vec/np.linalg.norm(vec))
        return result

    def gramschmidt(self, *vecs):
        result = []
        for vec in vecs:
            temp_vec = vec
            if vec.any():
                for j in result:
                    proj_vec = self.proj(j, vec)
                    temp_vec = map(lambda x, y: x-y, temp_vec, proj_vec)
                result.append(temp_vec)
        return self.normalize(*result)

    def augH(self, v1, v2):
        return np.dot(v1, self.H*v2)
