# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
# --
r'''A light-weight quadratic programming solver

   Problems of the following type can be solved:

   .. math::
        min_{x} \frac{1}{2} x^T A x - b^T x \\
        R x = s
        x \ge 0

   This is not the most general type of quadratic programming but it is
   sufficient for the needs in HORTON.

   The equality constraints are optional. When A is positive definite, a
   polynomial-scaling interior point algorithm can be used:
   ``QPSolver.find_local``. If A has negative eigenvalues, only the brute force
   solver is reliable but it's cost scales exponentially with problem size.

   No a priori feasibility tests are carried out and it is assumed that in
   ``QPSolver.find_local`` one is able to provide a feasible initial guess. The
   routine ``QPSolver.find_brute`` will complain a posteriori if the problem
   is infeasible.

   When A is negative definite, the solution may be unbounded. In the
   ``QPSolver.find_local`` algorithm, this will lead to a convergence failure
   in a limited number of steps. In `QPSolver.find_brute``, this situation is
   properly detected.

   In summary, it is assumed that one only works with well-behaved problems.
   Issues with feasibilty and boundedness are not always handled smoothly.

   Notation
   ........

   The following naming conventions are used in the quadprog code:

   free
        Refers to the components of the solution that is not constrained to
        zero. The variable free is always an array with booleans indicating
        which components are free (True) or frozen (False).

   frozen
        Opposite of free.

   x
        A solution vector

   nx
        The dimension of the solution vector

   nl
        The number of constraints
'''


# How to implement the primal feasibility test
# --------------------------------------------
# Checking the linear feasibility can be done by solving the following linear
# programming problem. Rewrite the equality constraints such that all components
# of s are positive. Then solve the following problem:
#
# min_y W  with W = d^T.y
# subject to  r.x + y = s  and  y >= 0  and  x >= 0.
#
# The initial feasible point for this problem is y=s and x=0. If after
# optimization, W is zero, the problem is feasible and the corresponding x is
# a feasible point. If W is positive, the linear problem is infeasible.


import numpy as np


__all__ = ['QPSolver']


class FeasibilityError(Exception):
    '''Error raised when the problem appears to be infeasible'''
    pass


class BoundedError(Exception):
    '''Error raised when the problem appears to be unbounded'''
    pass


class ConvergenceError(Exception):
    '''Error raised when the maximum number of iterations is exceeded'''
    pass


def _counter_to_free(counter, free):
    '''Convert a positive integer to a decomposition in bits

       **Arguments:**

       counter
            The integer

       free
            An output array with booleans to store the bits.
    '''
    power = 1
    for j in xrange(len(free)):
        free[j] = (counter & power) != 0
        power *= 2


def find_1d_root(fn, (x0, y0), (x2, y2), eps):
    '''Find the root of a 1D function

       **Arguments:**

       fn
            The function to be zeroed.

       (x0, y0), (x2, y2)
            Argument and function value pairs for the initial bracket.

       eps
            The allowed error on the function.
    '''
    # we want y0 < 0 and y2 > 0
    if y2 < 0:
        x0, y0, x2, y2 = x2, y2, x0, y0
    assert y0 < 0
    assert y2 > 0
    # root finder loop
    while True:
        # When t would be allowed to be close to 0 or 1, slow convergence may
        # occur with exponential-like (very non-linear) functions.
        t = np.clip(y0/(y0-y2), 0.1, 0.9)
        # compute the new point
        x1 = x0 + t*(x2-x0)
        y1 = fn(x1)
        # decide on convergence or which edge of the bracket to replace
        if abs(y1) < eps:
            return x1, y1
        elif y1 > 0:
            x2, y2 = x1, y1
        else:
            x0, y0 = x1, y1


def solve_safe(a, b):
    '''Try to solve with numpy.linalg.solve. Use SVD as fallback.

       **Arguments:**

       a, b
            Arrays for the matrix equation a x = b.
    '''
    try:
        # The default is to use the standard solver as it is faster and more
        # accurate than the SVD code below.
        return np.linalg.solve(a, b)
    except np.linalg.LinAlgError:
        # If the standard procedure breaks, fall back to an SVD-based solver
        # that selects the least-norm solution in case of trouble.
        u, s, vt = np.linalg.svd(a, full_matrices=False)
        rank = (s != 0).sum()
        u = u[:,:rank]
        s = s[:rank]
        vt = vt[:rank]
        assert s.min() > 0
        return np.dot(vt.T, np.dot(u.T, b)/s)


def check_constrained_problem(a, b, r=None, s=None):
    '''Check the matrices of a quadratic problem with linear constraints

       **Arguments:**

       a
            Quadratic coefficients. Array with shape (nx, nx).

       b
            Linear coefficients. Array with shape (nx,).

       **Optional arguments:**

       r
            Coefficients of equality constraints. Array with shape (nl, nx).

       s
            Targets of equality constraints. Array with shape (nl,).

    '''
    if r is not None and r.size == 0:
        r = None
    if r is None:
        s = None
    nx = len(a)
    if nx == 0:
        raise TypeError('The array a can not be empty.')
    if a.shape != (nx, nx):
        raise TypeError('The array a must be square.')
    if not (a == a.T).all():
        raise ValueError('The matrix a must be symmetric.')
    if b.shape != (nx,):
        raise TypeError('The size fo a and b do not match.')
    if r is not None:
        nl = len(r)
        if nl > nx:
            raise TypeError('There must not be more constraints than unknowns.')
        if r.shape != (nl, nx):
            raise TypeError('The shape of the equality constraints coefficients array is wrong.')
        if s.shape != (nl,):
            raise TypeError('The target vector of the inequality constraints has the wrong shape.')
    return a, b, r, s


def diagonal_form(a, b, r=None, s=None):
    '''Transform the constrained quadratic problem to an unconstraind diagonal form

       **Arguments:**

       a
            Quadratic coefficients. Array with shape (nx, nx).

       b
            Linear coefficients. Array with shape (nx,).

       **Optional arguments:**

       r
            Coefficients of equality constraints. Array with shape (nl, nx).

       s
            Targets of equality constraints. Array with shape (nl,).

       **Returns:**

       x0
            A particular solution of the constraints. Array with shape (nx,).

       basis
            A basis for the nullspace of the constraint equations (nb, nx),
            where nb >= nx - nl.

       a_diag
            The matrix a in the diagonal basis. Array with shape (nb,)

       b_diag
            The vector b transformed to the diagonal. Array with shape (nb,)

       evecs
            The eigen vectors of the diagonalization. Array with shape (nb, nb)

       A vector in diagonal basis can be transformed back as follows::

           x = x0 + np.dot(basis.T, np.dot(evecs, x_diag))

       The solution can be found as follows::

           x_diag = b_diag/evals
    '''
    a, b, r, s = check_constrained_problem(a, b, r, s)

    if r is not None:
        # parameterize the constraints
        u, singvals, vt = np.linalg.svd(r, full_matrices=True)
        rank = (singvals != 0).sum()
        if rank == 0:
            raise FeasibilityError('Constraints could not be satisfied in diagonal_form.')
        singvals = singvals[:rank]
        u = u[:rank]
        basis = vt[rank:] # rows are basis of the orthogonal complement
        vt = vt[:rank]
        x0 = np.dot(vt.T, np.dot(u.T, s)/singvals) # particular solution

        # it may be that there is not orthogonal complement. handle this nicely
        if basis.size == 0:
            return x0, None, None, None, None

        # transform a and b to unconstrained basis
        a_unconstrained = np.dot(np.dot(basis, a), basis.T)
        b_unconstrained = np.dot(basis, b - np.dot(a, x0))
    else:
        x0 = None
        basis = None
        a_unconstrained = a
        b_unconstrained = b

    # diagonalize
    evals, evecs = np.linalg.eigh(a_unconstrained)
    b_diag = np.dot(evecs.T, b_unconstrained)

    return x0, basis, evals, evecs, b_diag


def solve_constrained(a, b, r=None, s=None, eps=1e-10):
    '''Solve the quadratic problem, with linear constraints.

       The equality constraints are honored.

       **Arguments:**

       a
            Quadratic coefficients. Array with shape (nx, nx).

       b
            Linear coefficients. Array with shape (nx,).

       **Optional arguments:**

       r
            Coefficients of equality constraints. Array with shape (nl, nx).

       s
            Targets of equality constraints. Array with shape (nl,).
    '''
    x0, basis, evals, evecs, b_diag = diagonal_form(a, b, r, s)
    if evals is None:
        return x0
    x_diag = np.zeros(b_diag.shape)
    for i in xrange(len(x_diag)):
        if evals[i] != 0:
            x_diag[i] = b_diag[i]/evals[i]
    if x0 is None:
        return np.dot(evecs, x_diag)
    else:
        return x0 + np.dot(basis.T, np.dot(evecs, x_diag))

    # This was the old code. It is not as robust when a is ill-conditioned.
    ############################################################################
    # a, b, r, s = check_constrained_problem(a, b, r, s)
    # nx = len(a)
    # nl = 0 if r is None else len(r)
    # # construct the augmented equations
    # a2 = np.zeros((nx+nl, nx+nl))
    # b2 = np.zeros(nx+nl)
    # a2[:nx, :nx] = a
    # b2[:nx] = b
    # if r is not None:
    #     # fill in the constraints if present
    #     a2[nx:, :nx] = r
    #     a2[:nx, nx:] = r.T
    #     b2[nx:] = s
    # # solve the augmented problem
    # x2 = solve_safe(a2, b2)
    # x = x2[:nx] # [:nx] = drop the lagrange multipliers from the result
    # # check if the constraints are satisfied
    # if nl > 0 and abs(np.dot(r, x) - s).max() > eps:
    #     raise FeasibilityError('Constraints could not be satisfied in solve_constrained.')
    # # done
    # return x
    ############################################################################


def solve_radius(a, b, center, radius, r=None, s=None):
    '''Solve the quadratic problem, with linear constraints and maximum radius.

       The equality constraints are honored.

       **Arguments:**

       a
            Quadratic coefficients. Array with shape (nx, nx).

       b
            Linear coefficients. Array with shape (nx,).

       center
            The center around which the solution is to be found.

       radius
            The maximum Euclidean distance of the returned solution x from
            the given center.

       **Optional arguments:**

       r
            Coefficients of equality constraints. Array with shape (nl, nx).

       s
            Targets of equality constraints. Array with shape (nl,).
    '''
    # A) convert the problem into an unconstrained diagonal form
    x0, basis, evals, evecs, b_diag = diagonal_form(a, b, r, s)
    if evals is None:
        # only one solution is allowed
        min_radius = np.linalg.norm(x0 - center)
        if min_radius > radius:
            raise FeasibilityError('The constraints to not allow a solution closer to the center than %f > radius.' % min_radius)
        return x0
    if basis is None:
        b_shift = np.dot(center, evecs)
    else:
        # check if the distance from the constraints to the center is smaller than
        # the radius.
        yc = solve_safe(np.dot(basis, basis.T), np.dot(center - x0, basis.T))
        xc = x0 + np.dot(basis.T, yc)
        min_radius = np.linalg.norm(xc - center)
        if min_radius > radius:
            raise FeasibilityError('The constraints to not allow a solution closer to the center than %f > radius.' % min_radius)
        b_shift = np.dot(np.dot(center - x0, basis.T), evecs)

    # B) define an error function
    def solve(ridge):
        '''solution in diagional basis'''
        if ridge == 0.0:
            return b_diag/evals
        else:
            return b_diag/(evals + ridge) + b_shift/(evals/ridge+1)

    def to_orig(x_diag):
        '''transformation from diagonal to original basis'''
        x_free = np.dot(evecs, x_diag)
        if basis is not None:
            return x0 + np.dot(basis.T, x_free)
        else:
            return x_free

    def compute_error(ridge):
        '''compute the error on the radius. negative means small than trust radius'''
        x_diag = solve(ridge)
        x = to_orig(x_diag)
        return radius - np.linalg.norm(x - center)

    # Find a good starting point for the ridge parameter(s)
    evals_min = evals.min()
    if evals_min > 0:
        ridge0 = 0.0
        ridge1 = evals_min
    else:
        ridge0 = -1.1*evals_min
        ridge1 = -2.2*evals_min

    error0 = compute_error(ridge0)
    if error0 < 0:
        # In this case, we have to find the ridge parameter that sets the
        # error to zero

        # the error on this root finding problem may be relatively large
        eps1d = 1e-4*radius

        # fix right end of the bracket
        error1 = compute_error(ridge1)
        while error1 < 0:
            ridge1 *= 2
            assert ridge1 < 1e10
            error1 = compute_error(ridge1)

        # 1d root finder, it guarantees that the last call to error is done
        # with the returned ridge parameter
        ridge, error = find_1d_root(compute_error, (ridge0, error0), (ridge1, error1), eps1d)
    else:
        ridge = ridge0

    return to_orig(solve(ridge))



class QPSolver(object):
    '''A Quadratic Programming Solver'''
    def __init__(self, a, b, r=None, s=None, eps=1e-10):
        r'''The problem is defined as follows;

        .. math::
            \min_x \frac{1}{2} x^T A x - b^T x \\
            R x = s \\
            x \ge 0

        Parameters
        ----------

        a : np.ndarray
            The symmetric matrix :math:`A \in \mathbb{R}^{n \times n}`
        b : np.ndarray
            The column vector :math:`b \in \mathbb{R}^n`.
        r : np.ndarray or None
            The matrix with constraint coefficients, :math:`R \in \mathbb{R}^{l \times n}`.
        s : np.ndarray or None
            The matrix with constraint targets, :math:`s \in \mathbb{R}^l`.
        eps : float
              A general threshold used for several purposes, e.g. the validity of a
              solution.
        '''
        a, b, r, s = check_constrained_problem(a, b, r, s)

        # assign attributes
        self.a = a
        self.b = b
        self.r = r
        self.s = s
        self.eps = eps
        # intialize internals
        self._free = None

    def _get_nx(self):
        '''The number of unkowns'''
        return self.a.shape[0]

    nx = property(_get_nx)

    def _get_nl(self):
        '''The number of constraints'''
        return 0 if self.r is None else self.r.shape[0]

    nl = property(_get_nl)

    def compute_cost(self, x):
        '''Compute the function to be minimized for the given x'''
        return np.dot(x, 0.5*np.dot(self.a, x) - self.b)

    def check_feasible(self, x, free=None):
        '''Check if a solution, x, matches the constraints

           **Arguments:**

           x
                The solution to test

           **Optional arguments:**

           free
                When given, it is assumed that the fixed variables must be zero
                and that the rest can be anything (also negative). When not
                given, all components of x must be zero or positive.
        '''
        if self.nl > 0:
            if abs(np.dot(self.r, x) - self.s).max() > self.eps:
                raise FeasibilityError('Point does not satisfy equality constraints.')
        if free is None:
            if (x < 0).any():
                raise FeasibilityError('Solution has negative components.')
        else:
            (x[~free] == 0).all()

    def check_solution(self, x):
        '''Check if a solution, x, is a valid local minimum'''
        self.check_feasible(x)
        gradient, rmsd_free, rmsd_frozen, rmsd_neg = self.get_rmsds(x)
        if rmsd_free > self.eps:
            raise ValueError('The error on the gradient of the free components is too large.')
        if rmsd_frozen > self.eps:
            raise ValueError('The error on the gradient of the frozen components is too large.')
        if rmsd_neg > self.eps:
            raise ValueError('Some components are too negative.')

    def get_free_problem(self, free):
        '''Return the matrix a, b, r and s without the frozen columns/rows

           **Arguments:**

           free
                Boolean array with the free components

           **Returns:**

           a2
                The coefficients of the free system, including rows and colums
                for the lagrange multipliers.

           b2
                The right-hand side of the equations, including the targets of
                the equations.

           nfree
                The number of equations due to free components of x.
        '''
        assert free.dtype == bool
        # When the free argument has changed, we have to reconstruct the sliced
        # matrices.
        if self._free is None or (free != self._free).any():
            self._free = free.copy()
            self._nfree = free.sum()
            nfree = self._nfree
            assert nfree >= self.nl # more constraints than free params does not make sense
            # Slice the matrices an keep them as attributes for later
            self._a_free = self.a[free,:][:,free]
            self._b_free = self.b[free]
            self._r_free = self.r[:,free] if self.nl > 0 else None
            self._s_free = self.s if self.nl > 0 else None
        return self._a_free, self._b_free, self._r_free, self._s_free

    def get_rmsds(self, x):
        '''Quantify how far x deviates from a local solution

           **Returns:**

           gradient
                The gradient projected onto the equality constraints, if any.

           rmsd_free
                The gradient RMSD for the non-zero components

           rmsd_frozen
                The gradient RMSD for the zero components. Only contributions
                from negative components are included.

           rmsd_neg
                The rmsd of the negative components of the solution.
        '''
        rmsd_neg = np.sqrt((x*x*(x<0)).mean())
        gradient = np.dot(self.a, x) - self.b
        free = x > 0
        if self.nl > 0:
            # subtract the components of the gradient orthogonal to the constraints
            a_free, b_free, r_free, s_free = self.get_free_problem(free)
            g_free = np.dot(a_free, x[free]) - b_free
            lagrange = np.linalg.lstsq(r_free.T, -g_free)[0]
            gradient += np.dot(self.r.T, lagrange)
        rmsd_free = np.sqrt((gradient**2*free).mean())
        rmsd_frozen = np.sqrt((gradient**2*(~free)*(gradient<0)).mean())
        return gradient, rmsd_free, rmsd_frozen, rmsd_neg

    def compute_cn(self, free):
        '''Compute the condition number of the problem.

           **Arguments:**

           free
                Boolean array with the free components.
        '''
        a_free, b_free, r_free, s_free = self.get_free_problem(free)
        x0, basis, evals, evecs, b_diag = diagonal_form(a_free, b_free, r_free, s_free)
        if evals is None:
            return 0.0 # this happens when the solution is fully constrained.
        abs_evals = abs(evals)
        if abs_evals.min() == 0:
            return np.inf
        else:
            return abs_evals.max()/abs_evals.min()

    def _free_to_full(self, x_free):
        '''Convert a solution of get_free_system to the original space\

           **Arguments:**

           x_free
               A solution of the equations returned by get_free_system, obtained
               with solve_constrained or solve_radius

           **Returns:** a solution vector x with length nx.
        '''
        x = np.zeros(self.nx)
        x[self._free] = x_free
        return x

    def solve(self, free):
        '''Solve the problem, keeping certain components fixed at zero

           The equality constraints are honored. However, free components may
           become negative. Fixed components are internally excluded from the
           equations and their result is manually set to zero.

           **Arguments:**

           free
                Boolean array with the free components
        '''
        a_free, b_free, r_free, s_free = self.get_free_problem(free)
        x_free = solve_constrained(a_free, b_free, r_free, s_free, self.eps)
        return self._free_to_full(x_free)

    def solve_radius(self, free, center, radius):
        '''Solve the equations, keeping certain components fixed at zero, with maximum radius.

           The equality constraints are honored. However, free components may
           become negative. Fixed components are internally excluded from the
           equations and their result is manually set to zero.

           **Arguments:**

           free
                Boolean array with the free components

           center
                The center of the restraint, a vector with length nx. (:math:`x_0`)

           radius
                The maximum Euclidean distance of the returned solution x from
                the given center
        '''
        a_free, b_free, r_free, s_free = self.get_free_problem(free)
        x_free = solve_radius(a_free, b_free, center[free], radius, r_free, s_free)
        return self._free_to_full(x_free)

    def find_brute(self):
        '''Brute force solution of the quadratic programming problem

           **Returns:**

           cost
                The value of the cost function at the solution

           x
                The solution vector
        '''
        best = None
        free = np.zeros(self.nx, dtype=bool)
        for counter in xrange(1, 2**self.nx):
            _counter_to_free(counter, free)
            if free.sum() < self.nl:
                # skip overdetermined constraints
                continue
            x = self.solve(free)
            if x.min() >= 0.0:
                cost = self.compute_cost(x)
                if best is None or best[0] > cost:
                    best = cost, x, free.copy()

        # Is there any solution at all?
        if best is None:
            raise FeasibilityError('No solutions found in brute force search.')

        cost, x, free = best

        # Check if the solution is valid
        gradient, rmsd_free, rmsd_frozen, rmsd_neg = self.get_rmsds(x)
        assert rmsd_neg < self.eps
        if rmsd_frozen > self.eps or rmsd_free > self.eps:
            raise BoundedError('No bounded solution could be found.')

        # Check if the solution is not a saddle point
        a_free, b_free, r_free, s_free = self.get_free_problem(free)
        x0, basis, evals, evecs, b_diag = diagonal_form(a_free, b_free, r_free, s_free)
        if evals is not None and evals.min() < 0:
            raise BoundedError('No bounded solution could be found.')

        self.check_solution(x)
        return cost, x

    def find_local(self, x, trust_radius, maxiter=None):
        '''A local solver for the quadratic programming problem

           **Arguments:**

           x
                The initial guess for the solution. This must be a feasible
                point.

           trust_radius
                The maximum step size. This can be fairly large.

           **Optional arguments:**

           maxiter
                The maximum number of iterations, this is 10*nx by default.

           **Returns:**

           cost
                The value of the cost function at the solution

           x
                The solution vector


           Note that this local search will always lead to the global optimum if
           the matrix :math:`A` is positive definite.
        '''
        # keep a copy of the initial guess for debugging.
        guess = x.copy()

        # test if the initial state is feasible
        self.check_feasible(x)

        # set maxiter
        if maxiter is None:
            maxiter = self.nx*10

        # get all results associated with x
        free = x != 0.0
        cost = self.compute_cost(x)

        # iterative part
        for counter in xrange(maxiter):
            # check for convergence
            gradient, rmsd_free, rmsd_frozen, rmsd_neg = self.get_rmsds(x)
            if (rmsd_free < self.eps) and (rmsd_frozen < self.eps) and (rmsd_neg < self.eps):
                self.check_solution(x)
                return cost, x

            # deactivate the constraint with the most negative gradient, if any
            if rmsd_free < self.eps:
                ilow = (gradient*(~free)).argmin()
                if gradient[ilow] < 0:
                    free[ilow] = True

            # take a step with the current constraints, taking into account the
            # trust radius
            new_x = self.solve_radius(free, x, trust_radius)

            # find the first inequality constraint on the line segment x--new_x
            # that is violated, if any. only consider violations in this test
            tmin = 1.0
            imin = None
            for i in xrange(self.nx):
                if free[i] and new_x[i] < 0:
                    t = x[i]/(x[i] - new_x[i])
                    if tmin > t:
                        tmin = t
                        imin = i
            if tmin < 1.0:
                # in case of violation, step back the last feasible point on
                # the line segment.
                assert tmin > 0.0
                assert imin is not None
                new_x = (new_x - x)*tmin + x
                # fix the corresponding violation
                free[imin] = False
                # make the corresponding component a hard zero
                assert abs(new_x[imin]) < 1e-10
                new_x[imin] = 0.0

            # perform some checks
            new_cost = self.compute_cost(new_x)
            assert new_cost - cost < self.eps


            # let the new x become the current x and update the corresponding
            # quantities
            x = new_x
            cost = new_cost

        self.log(guess)
        raise ConvergenceError('Local search failed to converge')

    def log(self, x=None):
        '''Print out the qps for debugging purposes

           **Optional arguments:**

           x
                An solution vector (that causes problems).
        '''
        def _print_array(ar):
            print '    np.array(['
            for row in ar:
                print '        [%s],' % (', '.join(repr(v) for v in row))
            print '    ]),'

        def _print_vector(ve):
            print '    np.array([%s]),' % (', '.join(repr(v) for v in ve))

        print '#'*80
        print 'qps = QPSolver('
        _print_array(self.a)
        _print_vector(self.b)
        if self.nl > 0:
            _print_array(self.r)
            _print_vector(self.s)
        else:
            print 'None, None,'
        print ')'
        if x is not None:
            print 'x = np.array([%s])' % (', '.join(repr(v) for v in x))
        print '#'*80
