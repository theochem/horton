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
'''Optimal Damping Self-Consistent Field algorithm'''


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
    kind = 'dm' # basic variable is the density matrix

    def __init__(self, threshold=1e-8, maxiter=128, skip_energy=False, debug=False):
        self.maxiter = maxiter
        self.threshold = threshold
        self.skip_energy = skip_energy
        self.debug = debug

    @timer.with_section('SCF')
    def __call__(self, ham, lf, overlap, occ_model, *dm0s):
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

        fock0s = [lf.create_one_body() for i in xrange(ham.ndm)]
        fock1s = [lf.create_one_body() for i in xrange(ham.ndm)]
        dm1s = [lf.create_one_body() for i in xrange(ham.ndm)]
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
            for fock0 in fock0s:
                fock0.clear()
            ham.compute_fock(*fock0s)
            # Compute the energy in point 0
            energy0 = ham.compute()

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
            for fock1 in fock1s:
                fock1.clear()
            ham.compute_fock(*fock1s)
            # Compute the energy in point 1
            energy1 = ham.compute()

            # Compute the derivatives of the energy at point 0 and 1 for a
            # linear interpolation of the density matrices
            deriv0 = 0.0
            deriv1 = 0.0
            for i in xrange(ham.ndm):
                deriv0 += fock0s[i].expectation_value(dm1s[i])
                deriv0 -= fock0s[i].expectation_value(dm0s[i])
                deriv1 += fock1s[i].expectation_value(dm1s[i])
                deriv1 -= fock1s[i].expectation_value(dm0s[i])
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
                errorsq += commutator.expectation_value(commutator)
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
        return convergence_error_commutator(ham, lf, overlap, *dms)


def check_cubic(ham, dm0s, dm1s, e0, e1, g0, g1, do_plot=True):
    # coefficients of the polynomial a*x**3 + b*x**2 + c*x + d
    d = e0
    c = g0
    a = g1 - 2*e1 + c + 2*d
    b = e1 - a - c - d

    ndm = len(dm0s)
    assert ndm == len(dm1s)
    dm2s = [dm1.new() for dm1 in dm1s]
    xs = np.arange(0, 1.001, 0.1)
    energies = []
    for x in xs:
        for i in xrange(ndm):
            dm2s[i].assign(dm0s[i])
            dm2s[i].iscale(1-x)
            dm2s[i].iadd(dm1s[i], x)
        ham.reset(*dm2s)
        e2 = ham.compute()
        energies.append(e2)
    energies = np.array(energies)

    if do_plot:
        # make a nice figure
        xxs = np.arange(0, 1.00001, 0.001)
        poly = a*xxs**3+b*xxs**2+c*xxs+d
        import matplotlib.pyplot as pt
        pt.clf()
        pt.plot(xxs, poly, 'k-', label='cubic')
        pt.plot(xs, energies, 'ro', label='ref')
        pt.plot([0,1],[e0,e1], 'b--', label='linear')
        pt.ylim(min(energies.min(), poly.min()), max(energies.max(), poly.max()))
        pt.legend(loc=0)
        pt.savefig('check_cubic_cs_%+024.17f.png' % e0)
    else:
        # if not plotting, check that the errors are not too large
        poly = a*xs**3+b*xs**2+c*xs+d
        error = abs(poly-energies).max()
        oom = energies.max() - energies.min()
        assert error < 0.01*oom
