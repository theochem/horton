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
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN, check_dm


__all__ = ['check_cubic_cs', 'check_cubic_os', 'converge_scf_oda']


@timer.with_section('SCF')
def converge_scf_oda(ham, wfn, lf, overlap, maxiter=128, threshold=1e-6, debug=False):
    '''Minimize the energy of the wavefunction with optimal-damping SCF

       **Arguments:**

       ham
            A Hamiltonian instance.

       wfn
            The wavefunction object to be optimized.

       lf
            The linalg factory to be used.

       overlap
            The overlap operator.

       **Optional arguments:**

       maxiter
            The maximum number of iterations. When set to None, the SCF loop
            will go one until convergence is reached.

       threshold
            The convergence threshold for the wavefunction

       debug
            Make debug plots with matplotlib of the cubic interpolation

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.
    '''
    log.cite('cances2001', 'using the optimal damping algorithm (ODA) SCF')
    if isinstance(wfn, RestrictedWFN):
        return converge_scf_oda_cs(ham, wfn, lf, overlap, maxiter, threshold, debug)
    elif isinstance(wfn, UnrestrictedWFN):
        return converge_scf_oda_os(ham, wfn, lf, overlap, maxiter, threshold, debug)
    else:
        raise NotImplementedError


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
    '''Find the minimum of a cubic polynomial in the range [0,1]

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


def check_cubic_cs(ham, wfn, dm0, dm1, e0, e1, g0, g1, do_plot=True):
    dm1 = dm1.copy()

    # coefficients of the polynomial a*x**3 + b*x**2 + c*x + d
    d = e0
    c = g0
    a = g1 - 2*e1 + c + 2*d
    b = e1 - a - c - d

    dm2 = dm0.new()
    xs = np.arange(0, 1.001, 0.1)
    energies = []
    for x in xs:
        dm2.assign(dm0)
        dm2.iscale(1-x)
        dm2.iadd(dm1, x)

        wfn.clear()
        wfn.update_dm('alpha', dm2)
        ham.clear()
        e2 = ham.compute()
        energies.append(e2)
    energies = np.array(energies)

    if do_plot:
        # make a nice figure
        ys = np.arange(0, 1.00001, 0.001)
        poly = a*ys**3+b*ys**2+c*ys+d
        import matplotlib.pyplot as pt
        pt.clf()
        pt.plot(ys, poly, 'k-', label='cubic')
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


def converge_scf_oda_cs(ham, wfn, lf, overlap, maxiter=128, threshold=1e-6, debug=False):
    '''Minimize the energy of the closed-shell wavefunction with optimal damping

       **Arguments:**

       ham
            A Hamiltonian instance.

       wfn
            The wavefunction object to be optimized.

       lf
            The linalg factory to be used.

       overlap
            The overlap operator.

       **Optional arguments:**

       maxiter
            The maximum number of iterations. When set to None, the SCF loop
            will go one until convergence is reached.

       threshold
            The convergence threshold for the wavefunction.

       debug
            Make debug plots with matplotlib of the cubic interpolation

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.
    '''

    if log.do_medium:
        log('Starting restricted closed-shell SCF with optimal damping')
        log.hline()
        log(' Iter               Energy  Error(alpha)      Mixing')
        log.hline()

    # initialization of variables and datastructures
    # suffixes
    #    0 = current or initial state
    #    1 = state after conventional SCF step
    #    2 = state after optimal damping
    fock0 = lf.create_one_body()
    fock1 = lf.create_one_body()
    dm0 = lf.create_one_body()
    dm2 = lf.create_one_body()
    converged = False
    mixing = None
    error = None

    # Get rid of outdated stuff
    ham.clear()

    # If an input density matrix is present, check if it sensible. This avoids
    # redundant testing of a density matrix that is derived here from the
    # orbitals.
    if 'dm_alpha' in wfn._cache:
        check_dm(wfn.dm_alpha, overlap, lf, 'alpha')

    counter = 0
    while maxiter is None or counter < maxiter:
        # A) Construct Fock operator, compute energy and keep dm at current/initial point
        fock0.clear()
        ham.compute_fock(fock0, None)
        energy0 = ham.compute()
        dm0.assign(wfn.dm_alpha)

        if log.do_medium:
            if mixing is None:
                log('%5i %20.13f' % (counter, energy0))
            else:
                log('%5i %20.13f  %12.5e  %10.5f' % (counter, energy0, error, mixing))

        # B) Diagonalize fock operator and go to the next point
        wfn.clear()
        wfn.update_exp(fock0, overlap)
        # Let the hamiltonian know that the wavefunction has changed.
        ham.clear()

        # C) Compute Fock matrix at new point (lambda=1)
        fock1.clear()
        ham.compute_fock(fock1, None)
        # Compute energy at new point
        energy1 = ham.compute()
        # take the density matrix
        dm1 = wfn.dm_alpha

        # D) Find the optimal damping
        # Compute the derivatives of the energy towards lambda at edges 0 and 1
        ev_00 = fock0.expectation_value(dm0)
        ev_10 = fock1.expectation_value(dm0)
        ev_01 = fock0.expectation_value(dm1)
        ev_11 = fock1.expectation_value(dm1)
        energy0_deriv = 2*(ev_01-ev_00)
        energy1_deriv = 2*(ev_11-ev_10)

        # find the lambda that minimizes the cubic polynomial in the range [0,1]
        if log.do_high:
            log('           E0: % 10.5e     D0: % 10.5e' % (energy0, energy0_deriv))
            log('        E1-E0: % 10.5e     D1: % 10.5e' % (energy1-energy0, energy1_deriv))
        mixing = find_min_cubic(energy0, energy1, energy0_deriv, energy1_deriv)

        if debug:
            check_cubic_cs(ham, wfn, dm0, dm1, energy0, energy1, energy0_deriv, energy1_deriv)

        # E) Construct the new dm
        # Put the mixed dm in dm_old, which is local in this routine.
        dm2.clear()
        dm2.iadd(dm1, factor=mixing)
        dm2.iadd(dm0, factor=1-mixing)

        # Wipe the caches and use the interpolated density matrix
        wfn.clear()
        wfn.update_dm('alpha', dm2)
        ham.clear()

        error = dm2.distance(dm0)
        if error < threshold:
            #if abs(energy0_deriv) > threshold:
            #    raise RuntimeError('The ODA algorithm stopped a point with non-zero gradient.')
            converged = True
            break

        # counter
        counter += 1

    # compute orbitals, energies and occupation numbers
    # Note: suffix 0 is used for final state here
    dm0.assign(wfn.dm_alpha)
    fock0.clear()
    ham.compute_fock(fock0, None)
    wfn.clear()
    wfn.update_exp(fock0, overlap, dm0)
    energy0 = ham.compute()
    if log.do_medium:
        log('%5i %20.13f  %12.5e  %10.5f' % (counter+1, energy0, error, mixing))
        log.blank()

    if log.do_medium:
        ham.log_energy()

    if not converged:
        raise NoSCFConvergence

    return counter


def check_cubic_os(ham, wfn, dm0a, dm0b, dm1a, dm1b, e0, e1, g0, g1, do_plot=True):
    dm1a = dm1a.copy()
    dm1b = dm1b.copy()

    # coefficients of the polynomial a*x**3 + b*x**2 + c*x + d
    d = e0
    c = g0
    a = g1 - 2*e1 + c + 2*d
    b = e1 - a - c - d

    dm2a = dm0a.new()
    dm2b = dm0b.new()
    xs = np.arange(0, 1.001, 0.1)
    energies = []
    for x in xs:
        dm2a.assign(dm0a)
        dm2a.iscale(1-x)
        dm2a.iadd(dm1a, x)
        dm2b.assign(dm0b)
        dm2b.iscale(1-x)
        dm2b.iadd(dm1b, x)

        wfn.clear()
        wfn.update_dm('alpha', dm2a)
        wfn.update_dm('beta', dm2b)
        ham.clear()
        e2 = ham.compute()
        energies.append(e2)
    energies = np.array(energies)

    if do_plot:
        # make a nice figure
        ys = np.arange(0, 1.00001, 0.001)
        poly = a*ys**3+b*ys**2+c*ys+d
        import matplotlib.pyplot as pt
        pt.clf()
        pt.plot(ys, poly, 'k-', label='cubic')
        pt.plot(xs, energies, 'ro', label='ref')
        pt.plot([0,1],[e0,e1], 'b--', label='linear')
        pt.ylim(min(energies.min(), poly.min()), max(energies.max(), poly.max()))
        pt.legend(loc=0)
        pt.savefig('check_cubic_os_%+024.17f.png' % e0)
    else:
        # if not plotting, check that the errors are not too large
        poly = a*xs**3+b*xs**2+c*xs+d
        error = abs(poly-energies).max()
        oom = energies.max() - energies.min()
        assert error < 0.01*oom


def converge_scf_oda_os(ham, wfn, lf, overlap, maxiter=128, threshold=1e-6, debug=False):
    '''Minimize the energy of the open-shell wavefunction with optimal damping

       **Arguments:**

       ham
            A Hamiltonian instance.

       wfn
            The wavefunction object to be optimized.

       lf
            The linalg factory to be used.

       overlap
            The overlap operator.

       **Optional arguments:**

       maxiter
            The maximum number of iterations. When set to None, the SCF loop
            will go one until convergence is reached.

       threshold
            The convergence threshold for the wavefunction.

       debug
            Make debug plots with matplotlib of the cubic interpolation

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.

       **Returns:** the number of iterations
    '''

    if log.do_medium:
        log('Starting restricted closed-shell SCF with optimal damping')
        log.hline()
        log(' Iter               Energy  Error(alpha)  Error(beta)       Mixing')
        log.hline()

    # initialization of variables and datastructures
    # suffixes
    #    0 = current or initial state
    #    1 = state after conventional SCF step
    #    2 = state after optimal damping
    #    a = alpha
    #    b = beta
    fock0a = lf.create_one_body()
    fock1a = lf.create_one_body()
    fock0b = lf.create_one_body()
    fock1b = lf.create_one_body()
    dm0a = lf.create_one_body()
    dm2a = lf.create_one_body()
    dm0b = lf.create_one_body()
    dm2b = lf.create_one_body()
    converged = False
    mixing = None
    errora = None
    errorb = None

    # Get rid of outdated stuff
    ham.clear()

    # If an input density matrix is present, check if it sensible. This avoids
    # redundant testing of a density matrix that is derived here from the
    # orbitals.
    if 'dm_alpha' in wfn._cache:
        check_dm(wfn.dm_alpha, overlap, lf, 'alpha')
    if 'dm_beta' in wfn._cache:
        check_dm(wfn.dm_beta, overlap, lf, 'beta')

    counter = 0
    while maxiter is None or counter < maxiter:
        # A) Construct Fock operator, compute energy and keep dm at current/initial point
        fock0a.clear()
        fock0b.clear()
        ham.compute_fock(fock0a, fock0b)
        energy0 = ham.compute()
        dm0a.assign(wfn.dm_alpha)
        dm0b.assign(wfn.dm_beta)

        if log.do_medium:
            if mixing is None:
                log('%5i %20.13f' % (counter, energy0))
            else:
                log('%5i %20.13f  %12.5e  %12.5e  %10.5f' % (counter, energy0, errora, errorb, mixing))

        # B) Diagonalize fock operator and go to state 1
        wfn.clear()
        wfn.update_exp(fock0a, fock0b, overlap)
        # Let the hamiltonian know that the wavefunction has changed.
        ham.clear()

        # C) Compute Fock matrix at state 1
        fock1a.clear()
        fock1b.clear()
        ham.compute_fock(fock1a, fock1b)
        # Compute energy at new point
        energy1 = ham.compute()
        # take the density matrix
        dm1a = wfn.dm_alpha
        dm1b = wfn.dm_beta

        # D) Find the optimal damping
        # Compute the derivatives of the energy towards lambda at edges 0 and 1
        ev_00 = fock0a.expectation_value(dm0a) + fock0b.expectation_value(dm0b)
        ev_01 = fock0a.expectation_value(dm1a) + fock0b.expectation_value(dm1b)
        ev_10 = fock1a.expectation_value(dm0a) + fock1b.expectation_value(dm0b)
        ev_11 = fock1a.expectation_value(dm1a) + fock1b.expectation_value(dm1b)
        energy0_deriv = (ev_01-ev_00)
        energy1_deriv = (ev_11-ev_10)

        # find the lambda that minimizes the cubic polynomial in the range [0,1]
        if log.do_high:
            log('           E0: % 10.5e     D0: % 10.5e' % (energy0, energy0_deriv))
            log('        E1-E0: % 10.5e     D1: % 10.5e' % (energy1-energy0, energy1_deriv))
        mixing = find_min_cubic(energy0, energy1, energy0_deriv, energy1_deriv)
        if log.do_high:
            log('       mixing: % 10.5f' % mixing)

        if debug:
            check_cubic_os(ham, wfn, dm0a, dm0b, dm1a, dm1b, energy0, energy1, energy0_deriv, energy1_deriv)

        # E) Construct the new dm
        # Put the mixed dm in dm_old, which is local in this routine.
        dm2a.clear()
        dm2a.iadd(dm1a, factor=mixing)
        dm2a.iadd(dm0a, factor=1-mixing)
        dm2b.clear()
        dm2b.iadd(dm1b, factor=mixing)
        dm2b.iadd(dm0b, factor=1-mixing)

        # Wipe the caches and use the interpolated density matrix
        wfn.clear()
        wfn.update_dm('alpha', dm2a)
        wfn.update_dm('beta', dm2b)
        ham.clear()

        errora = dm2a.distance(dm0a)
        errorb = dm2b.distance(dm0b)
        if errora < threshold and errorb < threshold:
            if abs(energy0_deriv) > threshold:
                raise RuntimeError('The ODA algorithm stopped a point with non-zero gradient.')
            converged = True
            break

        # counter
        counter += 1

    # compute orbitals, energies and occupation numbers
    # Note: suffix 0 is used for final state here
    dm0a.assign(wfn.dm_alpha)
    dm0b.assign(wfn.dm_beta)
    fock0a.clear()
    fock0b.clear()
    ham.compute_fock(fock0a, fock0b)
    wfn.clear()
    wfn.update_exp(fock0a, fock0b, overlap, dm0a, dm0b)
    energy0 = ham.compute()
    if log.do_medium:
        log('%5i %20.13f  %12.5e  %12.5e  %10.5f' % (counter+1, energy0, errora, errorb, mixing))
        log.blank()

    if log.do_medium:
        ham.log_energy()

    if not converged:
        raise NoSCFConvergence

    return counter
