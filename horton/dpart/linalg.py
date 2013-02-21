# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import numpy as np


__all__ = ['solve_positive', 'quadratic_solver', 'constrained_solver']


def solve_positive(A, B, bindings_eq=[], upper_limits=None, rcond=0):
    npar = A.shape[1]
    bindings_ineq = []
    for i in xrange(npar):
        row = np.zeros(npar, float)
        row[i] = 1
        bindings_ineq.append((row, 0))
    if upper_limits is not None:
        for i in xrange(npar):
            row = np.zeros(npar, float)
            row[i] = -1
            bindings_ineq.append((row, -upper_limits[i]))

    return quadratic_solver(A, B, bindings_eq, bindings_ineq, rcond)


def quadratic_solver(A, b, binding_eq, binding_ineq, rcond=0):
    """A quadratic program solver.

       **Arguments:**

       A
            The coefficients of the equations

       b
            The right-hand side

       binding_eq
            Linear equality constraint equations. These equations are of the
            form. [(row_i, val_i), ...]. A solution x must satisfy all equations
            dot(row_i, x) = val_i. Can be None.

       binding_ineq
            Linear inequality constraint equations.  These equations are of the
            form. [(row_i, val_i), ...]. A solution x must satisfy all equations
            dot(row_i, x) >= val_i.

       rcond
            When different from zero, rcond**2 is added to the diagonal of the
            normal equations to regularize the fit.

       The implementation is based on the normal equations in order to keep
       things simple. The constraint equations are not inverted separately,
       which avoids a lot of trouble.
    """
    if False:
        print A
        print b
        print binding_eq
        print binding_ineq
        # When this routine fails for a certain system. This printout can be
        # used to add a test_quadratic_solver_hard* to the unit tests.
        from cPickle import dumps
        print (dumps(A),)
        print (dumps(b),)
        print (dumps(binding_eq),)
        print (dumps(binding_ineq),)
        print rcond

    # The matrix A must be positive definite
    evals = np.linalg.eigvalsh(A)
    if evals.min() < 0:
        raise ValueError('The matrix A must be positive definite.')

    # Thresholds, should probably become arguments
    eps_x = 1e-10 # the allowed error on a constraint.
    eps_l = 1e-10 # the allowed error on a Lagrange_multiplier.
    eps_g = 1e-10 # the allowed error on the projected gradient.
    maxiter = A.shape[1]*4

    # Useful dimensions
    npar = len(A) # number of unknowns
    ne = len(binding_eq) # number of equality constraints
    ni = len(binding_ineq) # number of inequality constraints

    # Set initial status of the inequality constraints.
    active = np.zeros(ni, dtype=bool)
    # Keep track of the worst condition number.
    cn_max = 0

    for counter in xrange(maxiter):
        #print 'iter', counter, active.astype(int)
        # Setup the constrained problem and run the constrained solver.
        R = []
        s = []
        for i in xrange(ni): # inequality constraints come first
            if active[i]:
                R.append(binding_ineq[i][0])
                s.append(binding_ineq[i][1])
        for i in xrange(ne): # then equality constraints
            R.append(binding_eq[i][0])
            s.append(binding_eq[i][1])
        R = np.array(R)
        if len(R.shape) == 1:
            R = R.reshape(1,-1)
        s = np.array(s)
        x, l, cn = constrained_solver(A, b, R, s, rcond)
        cn_max = max(cn, cn_max) # condition number

        if ni == 0:
            return x

        # A) Activate the inequality constraints that is most violated.

        # Compute all the inequality constraints:
        dev_ineq = np.zeros(ni, float)
        for i in xrange(ni):
            dev_ineq[i] = np.dot(binding_ineq[i][0], x) - binding_ineq[i][1]
        #print 'x         ', x
        #print 'active    ', active.astype(int)
        #print 'deviations', dev_ineq
        #print 'lagrange  ', l
        # Select a constraint for activation.
        need_deact = False
        iactivate = dev_ineq.argmin()
        if dev_ineq[iactivate] < -eps_x:
            #print '  activate', i
            if active.sum() + len(binding_eq) >= npar:
                # If the number of constraints exceeds the number of parameters,
                # we need to deactivate a constraint in order to avoid an over-
                # constrained system.
                need_deact = True
        else:
            iactivate = None
            need_deact = True

        # B) If no new inequality constraints were violated or if the number
        # of constraints has become too large, find the constraint with the
        # largest positive Lagrange multiplier and deactivate it.
        lmax = None
        ideactivate = None
        if need_deact:
            ia = 0
            for i in xrange(ni):
                if active[i]:
                    li = l[ia]
                    if li > eps_l and (lmax is None or li > lmax):
                        lmax = li
                        ideactivate = i
                    ia += 1
            if ideactivate is not None:
                active[ideactivate] = False
                #print '  deactivate', imax

        # A) Finish activation
        if iactivate is not None:
            active[iactivate] = True

        if iactivate is not None or ideactivate is not None:
            # If something happened to the constraints, redo the constrained fit.
            continue

        # C) If no constraints were activated or deactivated, we have a solution.
        # Double check whether solution satisfies all conditions.

        # Compute all equality constraints.
        dev_eq = np.zeros(ne, float)
        for i in xrange(ne):
            dev_eq[i] = np.dot(binding_eq[i][0], x) - binding_eq[i][1]

        # Raise assertion errors where needed.
        if (dev_ineq < -eps_x).any():
            raise AssertionError('Violation of inequality constraint in final result.')
        if (dev_eq < -eps_x).any():
            raise AssertionError('Violation of equality constraint in final result.')

        # Check for projected gradient.
        grad = np.dot(A, x) - b + rcond**2*x
        if R.size > 0:
            U, S, Vt = np.linalg.svd(R, full_matrices=True)
            ortho = Vt[:len(R)]
            grad -= np.dot(ortho.T, np.dot(ortho, grad))
        if (abs(grad) > eps_g).any():
            raise AssertionError('Projected gradient is not close enough to zero.')

        return x

    raise AssertionError('Qaudratic solver did not converge in %i iterations. Worst condition number was %.1e' % (maxiter, cn))


def constrained_solver(A, b, R, s, rcond=0):
    """Minimizes 1/2 x^t A x - b^t x, subject to R x = s.

       The implementation is based on an eigen decomposition, which is also used
       to check the condition number and optionally regularize the equations.
       This is not the most efficient approach, neither does it optimally make
       use of numerical precision. However, the approach is robust and simple.
    """
    ncon = len(R)
    npar = len(A)
    nbig = npar + ncon
    if R.size == 0:
        big_A = A + rcond**2*np.identity(npar)
        big_b = b
    else:
        big_A = np.zeros((nbig, nbig), float)
        big_A[:npar,:npar] = A + rcond**2*np.identity(npar)
        big_A[npar:,:npar] = R
        big_A[:npar,npar:] = R.T
        big_b = np.concatenate([b, s])
    evals, evecs = np.linalg.eigh(big_A)
    big_x = np.dot(evecs, np.dot(evecs.T, big_b)/evals)
    return big_x[:npar], big_x[npar:], abs(evals).max()/abs(evals).min()
