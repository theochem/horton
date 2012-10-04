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

from horton.log import log, timer
from horton.wfn import ClosedShellWFN, OpenShellWFN


__all__ = ['converge_scf', 'converge_scf_oda_cs', 'convergence_error']


def converge_scf(ham, max_iter=128, threshold=1e-8):
    '''Minimize the energy of the wavefunction with basic SCF

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       max_iter
            The maximum number of iterations

       threshold
            The convergence threshold for the wavefunction

       **Returns:**

       converged
            True of converged, False if not
    '''
    with timer.section('SCF'):
        if isinstance(ham.system.wfn, ClosedShellWFN):
            return converge_scf_cs(ham, max_iter, threshold)
        elif isinstance(ham.system.wfn, OpenShellWFN):
            return converge_scf_os(ham, max_iter, threshold)
        else:
            raise NotImplementedError


def converge_scf_cs(ham, max_iter=128, threshold=1e-8):
    '''Minimize the energy of the wavefunction with basic closed-shell SCF

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       max_iter
            The maximum number of iterations.

       threshold
            The convergence threshold for the wavefunction.

       **Returns:**

       converged
            True of converged, False if not.
    '''
    if log.do_medium:
        log('Starting restricted closed-shell SCF')
        log.hline()
        log('Iter  Error(alpha)')
        log.hline()

    lf = ham.system.lf
    wfn = ham.system.wfn
    fock = lf.create_one_body()
    converged = False
    for i in xrange(max_iter):
        # Construct the Fock operator
        fock.reset()
        ham.compute_fock(fock, None)
        # Check for convergence
        error = lf.error_eigen(fock, ham.overlap, wfn.exp_alpha)

        if log.do_medium:
            log('%4i  %12.5e' % (i, error))

        if error < threshold:
            converged = True
            break
        # Diagonalize the fock operator
        wfn.invalidate() # discard previous wfn state
        wfn.update_exp(fock, ham.overlap)
        # Let the hamiltonian know that the wavefunction has changed.
        ham.invalidate()
        # Write intermediate results to checkpoint
        ham.system.update_chk('wfn')
        # Take a backup copy of the Fock matrix

    if log.do_medium:
        log.hline()
        log('SCF converged: %s' % converged)

    return converged


def converge_scf_os(ham, max_iter=128, threshold=1e-8):
    '''Minimize the energy of the wavefunction with basic open-shell SCF

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       max_iter
            The maximum number of iterations.

       threshold
            The convergence threshold for the wavefunction.

       **Returns:**

       converged
            True of converged, False if not.
    '''
    if log.do_medium:
        log('Starting unrestricted open-shell SCF')
        log.hline()
        log('Iter  Error(alpha)  Error(beta)')
        log.hline()

    lf = ham.system.lf
    wfn = ham.system.wfn
    fock_alpha = lf.create_one_body()
    fock_beta = lf.create_one_body()
    converged = False
    for i in xrange(max_iter):
        # Construct the Fock operators
        fock_alpha.reset()
        fock_beta.reset()
        ham.compute_fock(fock_alpha, fock_beta)
        # Check for convergence
        error_alpha = lf.error_eigen(fock_alpha, ham.overlap, wfn.exp_alpha)
        error_beta = lf.error_eigen(fock_beta, ham.overlap, wfn.exp_beta)

        if log.do_medium:
            log('%4i  %12.5e  %12.5e' % (i, error_alpha, error_beta))

        if error_alpha < threshold and error_beta < threshold:
            converged = True
            break
        # Diagonalize the fock operators
        wfn.invalidate()
        wfn.update_exp(fock_alpha, fock_beta, ham.overlap)
        # Let the hamiltonian know that the wavefunction has changed.
        ham.invalidate()
        # Write intermediate results to checkpoint
        ham.system.update_chk('wfn')
        # Take backup copies of the Fock matrices.

    if log.do_medium:
        log.hline()
        log('SCF converged: %s' % converged)

    return converged


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
        if abs(a) > 0:
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


def converge_scf_oda_cs(ham, max_iter=128, threshold=1e-8):
    '''Minimize the energy of the wavefunction with basic closed-shell SCF

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       max_iter
            The maximum number of iterations.

       threshold
            The convergence threshold for the wavefunction.

       **Returns:**

       converged
            True of converged, False if not.
    '''
    log.cite('cances2001', 'using the optimal damping algorithm (ODA)')

    if log.do_medium:
        log('Starting restricted closed-shell SCF with optimal damping')
        log.hline()
        log('Iter     Energy      Error(alpha)   Mixing')
        log.hline()

    # initialization of variables and datastructures
    lf = ham.system.lf
    wfn = ham.system.wfn
    fock = lf.create_one_body()
    fock_old = lf.create_one_body()
    dm_old = lf.create_one_body()
    dm_new = lf.create_one_body()
    converged = False
    mixing = None
    error = None

    for i in xrange(max_iter):
        # A) Construct Fock operator, compute energy and keep dm at old point
        fock_old.reset()
        ham.compute_fock(fock_old, None)
        energy_old = ham.compute_energy()
        dm_old.assign(wfn.dm_alpha)

        if log.do_medium:
            if mixing is None:
                log('%5i %15.10f' % (i, energy_old, ))
            else:
                log('%5i %15.10f %10.5e %10.5f' % (i, energy_old, error, mixing))

        # B) Diagonalize fock operator and go to the next point
        wfn.invalidate()
        wfn.update_exp(fock_old, ham.overlap)
        # Let the hamiltonian know that the wavefunction has changed.
        ham.invalidate()

        # C) Compute Fock matrix at new point (lambda=1)
        fock.reset()
        ham.compute_fock(fock, None)
        # Compute energy at new point
        energy = ham.compute_energy()
        # take the density matrix
        dm = wfn.dm_alpha

        # D) Find the optimal damping
        # Compute the derivatives of the energy towards lambda at edges 0 and 1
        ev_11 = fock.expectation_value(dm)
        ev_01 = fock_old.expectation_value(dm)
        ev_10 = fock.expectation_value(dm_old)
        ev_00 = fock_old.expectation_value(dm_old)
        energy_deriv = 2*(ev_11-ev_10)
        energy_old_deriv = 2*(ev_01-ev_00)

        # find the lambda that minimizes the cubic polynomial in the range [0,1]
        if log.do_high:
            log('           E0: % 10.5e     D0: % 10.5e' % (energy_old, energy_old_deriv))
            log('        E1-E0: % 10.5e     D1: % 10.5e' % (energy-energy_old, energy_deriv))
        mixing = find_min_cubic(energy_old, energy, energy_old_deriv, energy_deriv)

        # E) Construct the new dm
        # Put the mixed dm in dm_old, which is local in this routine.
        dm_new.reset()
        dm_new.iadd(dm, factor=mixing)
        dm_new.iadd(dm_old, factor=1-mixing)

        # Wipe the caches and use the interpolated density matrix
        wfn.invalidate()
        wfn.update_dm('alpha', dm_new)
        ham.invalidate()

        error = dm_new.distance(dm_old)
        if error < threshold:
            converged = True
            break

        # Write intermediate wfn to checkpoint.
        ham.system.update_chk('wfn')

    if converged:
        # compute orbitals, energies and occupation numbers
        fock_old.reset()
        ham.compute_fock(fock_old, None)
        wfn.invalidate()
        wfn.update_exp(fock_old, ham.overlap, dm_new)
        energy = ham.compute_energy()
        # Write final wfn to checkpoint.
        ham.system.update_chk('wfn')
        if log.do_medium:
            log('%5i %15.10f %10.5e %10.5f' % (i, energy, error, mixing))

    if log.do_medium:
        log.hline()
        log('SCF converged: %s' % converged)

    return converged


def convergence_error(ham):
    '''Compute the self-consistency error

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Returns:**

       error
            The SCF error. This measure (not this function) is also used
            in the SCF algorithm to check for convergence.
    '''
    if isinstance(ham.system.wfn, ClosedShellWFN):
        return convergence_error_cs(ham)
    elif isinstance(ham.system.wfn, OpenShellWFN):
        return convergence_error_os(ham)
    else:
        raise NotImplementedError


def convergence_error_cs(ham):
    '''Compute the self-consistency error for a close-shell wavefunction

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Returns:**

       error
            The SCF error. This measure (not this function) is also used
            in the SCF algorithm to check for convergence.
    '''
    lf = ham.system.lf
    wfn = ham.system.wfn
    fock = lf.create_one_body()
    # Construct the Fock operator
    fock.reset()
    ham.compute_fock(fock, None)
    # Compute error
    return lf.error_eigen(fock, ham.overlap, wfn.exp_alpha)


def convergence_error_os(ham):
    '''Compute the self-consistency error for an open-shell wavefunction

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Returns:**

       error
            The maximum of the alpha and beta SCF errors. This measure
            (not this function) is also used in the SCF algorithm to check
            for convergence.
    '''
    lf = ham.system.lf
    wfn = ham.system.wfn
    fock_alpha = lf.create_one_body()
    fock_beta = lf.create_one_body()
    # Construct the Fock operators
    fock_alpha.reset()
    fock_beta.reset()
    ham.compute_fock(fock_alpha, fock_beta)
    # Compute errors
    error_alpha = lf.error_eigen(fock_alpha, ham.overlap, wfn.exp_alpha)
    error_beta = lf.error_eigen(fock_beta, ham.overlap, wfn.exp_beta)
    return max(error_alpha, error_beta)
