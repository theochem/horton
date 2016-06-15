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
'''Optimal Damping SCF algorithm'''


import numpy as np

from horton.log import log, timer
from horton.exceptions import NoSCFConvergence
from horton.meanfield.utils import check_dm, compute_commutator
from horton.meanfield.convergence import convergence_error_commutator


__all__ = ['ODASCFSolver', 'check_cubic']


def find_min_cubic(f0, f1, g0, g1):
    '''Find the minimum of a cubic polynomial in the range [0,1]

       **Arguments:**

       f0
            The function at argument 0
       f1
            The function at argument 1
       g0
            The derivative at argument 0
       g1
            The derivative at argument 1
    '''
    # coefficients of the polynomial a*x**3 + b*x**2 + c*x + d
    d = f0
    c = g0
    a = g1 - 2*f1 + c + 2*d
    b = f1 - a - c - d

    # find the roots of the derivative
    discriminant = b**2 - 3*a*c # simplified expression, not a mistake!
    if discriminant >= 0:
        if b*b < abs(3*a*c)*1e5:
            if log.do_high:
                log('               cubic')
            xa = (-b + np.sqrt(discriminant))/(3*a)
            xb = (-b - np.sqrt(discriminant))/(3*a)
            # test the solutions
            for x in xa, xb:
                if x >= 0 and x <= 1:
                    # compute the curvature at the solution
                    curv = 6*a*x+2*b
                    if curv > 0:
                        # Only one of two solutions has the right curvature, no
                        # need to compare the energy of both solutions
                        return x
        elif b > 0: # only b > 0 because b is also the curvature
            if log.do_high:
                log('               quadratic')
            x = -0.5*c/b
            if x >= 0 and x <= 1:
                return x

    # If we get here, no solution was found in the interval. One of the
    # boundaries is then the minimum
    if log.do_high:
        log('               edge')
    if f0 < f1:
        return 0.0
    else:
        return 1.0


def find_min_quadratic(g0, g1):
    '''Find the minimum of a quadratic polynomial in the range [0,1]

       **Arguments:**

       g0
            The derivative at argument 0
       g1
            The derivative at argument 1
    '''
    if g0 > 0:
        if g1 < -g0:
            return 1.0
        else:
            return 0.0
    else:
        if g1 > 0:
            # the regular case:
            return g0 / (g0 - g1)
        else:
            return 1.0


class ODASCFSolver(object):
    '''Optimal damping SCF algorithm (with cubic interpolation)'''
    kind = 'dm' # input/output variable is the density matrix

    def __init__(self, threshold=1e-8, maxiter=128, skip_energy=False, debug=False):
        '''
           **Optional arguments:**

           maxiter
                The maximum number of iterations. When set to None, the SCF loop
                will go one until convergence is reached.

           threshold
                The convergence threshold for the wavefunction

           skip_energy
                When set to True, the final energy is not computed.

           debug
                When set to True, for each mixing step, a plot is made to double
                check the cubic interpolation.
        '''
        self.maxiter = maxiter
        self.threshold = threshold
        self.skip_energy = skip_energy
        self.debug = debug

    @timer.with_section('SCF')
    def __call__(self, ham, lf, overlap, occ_model, *dm0s):
        '''Find a self-consistent set of density matrices.

           **Arguments:**

           ham
                An effective Hamiltonian.

           lf
                The linalg factory to be used.

           overlap
                The overlap operator.

           occ_model
                Model for the orbital occupations.

           dm1, dm2, ...
                The initial density matrices. The number of dms must match
                ham.ndm.
        '''
        # Some type checking
        if ham.ndm != len(dm0s):
            raise TypeError('The number of initial density matrices does not match the Hamiltonian.')

        # Check input density matrices.
        for i in xrange(ham.ndm):
            check_dm(dm0s[i], overlap, lf)
        occ_model.check_dms(overlap, *dm0s)

        if log.do_medium:
            log('Starting SCF with optimal damping. ham.ndm=%i' % ham.ndm)
            log.hline()
            log(' Iter               Energy         Error      Mixing')
            log.hline()

        fock0s = [lf.create_two_index() for i in xrange(ham.ndm)]
        fock1s = [lf.create_two_index() for i in xrange(ham.ndm)]
        dm1s = [lf.create_two_index() for i in xrange(ham.ndm)]
        exps = [lf.create_expansion() for i in xrange(ham.ndm)]
        work = dm0s[0].new()
        commutator = dm0s[0].new()
        converged = False
        counter = 0
        mixing = None
        error = None
        while self.maxiter is None or counter < self.maxiter:
            # feed the latest density matrices in the hamiltonian
            ham.reset(*dm0s)
            # Construct the Fock operators in point 0
            ham.compute_fock(*fock0s)
            # Compute the energy in point 0
            energy0 = ham.compute_energy()

            if log.do_medium:
                if mixing is None:
                    log('%5i %20.13f' % (counter, energy0))
                else:
                    log('%5i %20.13f  %12.5e  %10.5f' % (counter, energy0, error, mixing))

            # go to point 1 by diagonalizing the fock matrices
            for i in xrange(ham.ndm):
                exps[i].from_fock(fock0s[i], overlap)
            # Assign new occupation numbers.
            occ_model.assign(*exps)
            # Construct the density matrices
            for i in xrange(ham.ndm):
                exps[i].to_dm(dm1s[i])

            # feed the latest density matrices in the hamiltonian
            ham.reset(*dm1s)
            # Compute the fock matrices in point 1
            ham.compute_fock(*fock1s)
            # Compute the energy in point 1
            energy1 = ham.compute_energy()

            # Compute the derivatives of the energy at point 0 and 1 for a
            # linear interpolation of the density matrices
            deriv0 = 0.0
            deriv1 = 0.0
            for i in xrange(ham.ndm):
                deriv0 += fock0s[i].contract_two('ab,ab', dm1s[i])
                deriv0 -= fock0s[i].contract_two('ab,ab', dm0s[i])
                deriv1 += fock1s[i].contract_two('ab,ab', dm1s[i])
                deriv1 -= fock1s[i].contract_two('ab,ab', dm0s[i])
            deriv0 *= ham.deriv_scale
            deriv1 *= ham.deriv_scale

            # find the lambda that minimizes the cubic polynomial in the range [0,1]
            if log.do_high:
                log('           E0: % 10.5e     D0: % 10.5e' % (energy0, deriv0))
                log('        E1-E0: % 10.5e     D1: % 10.5e' % (energy1-energy0, deriv1))
            mixing = find_min_cubic(energy0, energy1, deriv0, deriv1)

            if self.debug:
                check_cubic(ham, dm0s, dm1s, energy0, energy1, deriv0, deriv1)

            # compute the mixed density and fock matrices (in-place in dm0s and fock0s)
            for i in xrange(ham.ndm):
                dm0s[i].iscale(1-mixing)
                dm0s[i].iadd(dm1s[i], mixing)
                fock0s[i].iscale(1-mixing)
                fock0s[i].iadd(fock1s[i], mixing)

            # Compute the convergence criterion.
            errorsq = 0.0
            for i in xrange(ham.ndm):
                compute_commutator(dm0s[i], fock0s[i], overlap, work, commutator)
                errorsq += commutator.contract_two('ab,ab', commutator)
            error = errorsq**0.5
            if error < self.threshold:
                converged = True
                break
            elif mixing == 0.0:
                raise NoSCFConvergence('The ODA algorithm made a zero step without reaching convergence.')

            # counter
            counter += 1

        if log.do_medium:
            ham.log()

        if not converged:
            raise NoSCFConvergence

        return counter

    def error(self, ham, lf, overlap, *dms):
        '''See :py:func:`horton.meanfield.convergence.convergence_error_commutator`.'''
        return convergence_error_commutator(ham, lf, overlap, *dms)


def check_cubic(ham, dm0s, dm1s, e0, e1, g0, g1, do_plot=True):
    '''Method to test the correctness of the cubic interpolation

       **Arguments:**

       ham
            An effective Hamiltonian.

       dm0s
            A list of density matrices for the initial state.

       dm1s
            A list of density matrices for the final state.

       e0
            The energy of the initial state.

       e1
            The energy of the final state.

       g0
            The derivative of the energy for the initial state.

       g1
            The derivative of the energy for the initial state.

       This method interpolates between two state with a parameter lambda going
       from 0 to 1. The arguments g0 and g1 are derivatives toward this
       parameter lambda.

       **Optional arguments:**

       do_plot
            When set to True, a plot is made to visualize the interpolation.
            Otherwise, an AssertionError is made if the error on the
            interpolation. is too large.

       This function is mainly useful as a tool to double check the
       implementation of Fock matrices.
    '''
    # coefficients of the polynomial a*x**3 + b*x**2 + c*x + d
    d = e0
    c = g0
    a = g1 - 2*e1 + c + 2*d
    b = e1 - a - c - d

    ndm = len(dm0s)
    assert ndm == len(dm1s)
    dm2s = [dm1.new() for dm1 in dm1s]
    xs = np.array([0.001, 0.002, 0.003, 0.004, 0.005, 0.995, 0.996, 0.997, 0.998, 0.999])
    energies = []
    for x in xs:
        for i in xrange(ndm):
            dm2s[i].assign(dm0s[i])
            dm2s[i].iscale(1-x)
            dm2s[i].iadd(dm1s[i], x)
        ham.reset(*dm2s)
        e2 = ham.compute_energy()
        energies.append(e2)
    energies = np.array(energies)

    if do_plot:
        # make a nice figure
        xxs = np.concatenate([np.linspace(0, 0.006, 60), np.linspace(0.994, 1.0, 60)])
        poly = a*xxs**3+b*xxs**2+c*xxs+d
        import matplotlib.pyplot as pt
        pt.clf()
        pt.subplot(121)
        pt.plot(xxs, poly, 'k-', label='cubic')
        pt.plot(xs, energies, 'ro', label='ref')
        pt.plot([0,1],[e0,e1], 'b--', label='linear')
        pt.xlim(0, 0.006)
        pt.ylim(min(energies[:5].min(), poly[:60].min()), max(energies[:5].max(), poly[:60].max()))
        pt.legend(loc=0)
        pt.subplot(122)
        pt.plot(xxs, poly, 'k-', label='cubic')
        pt.plot(xs, energies, 'ro', label='ref')
        pt.plot([0,1],[e0,e1], 'b--', label='linear')
        pt.xlim(0.994, 1.0)
        pt.ylim(min(energies[-5:].min(), poly[-60:].min()), max(energies[-5:].max(), poly[-60:].max()))
        pt.legend(loc=0)
        pt.savefig('check_cubic_%+024.17f.png' % e0)
    else:
        # if not plotting, check that the errors are not too large
        poly = a*xs**3+b*xs**2+c*xs+d
        relative_errors = abs(poly[:5]-energies[:5])/(energies[:5] - e0)
        assert relative_errors.max() < 0.03
        relative_errors = abs(poly[-5:]-energies[-5:])/(energies[-5:] - e1)
        assert relative_errors.max() < 0.03
