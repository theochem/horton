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
'''Basic Self-Consistent Field algorithm'''


from horton.log import log, timer
from horton.exceptions import NoSCFConvergence
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN
from horton.meanfield.hamiltonian import RestrictedEffectiveHamiltonian, \
    UnrestrictedEffectiveHamiltonian


__all__ = ['converge_scf']


@timer.with_section('SCF')
def converge_scf(ham, wfn, lf, overlap, maxiter=128, threshold=1e-8, skip_energy=False):
    '''Minimize the energy of the wavefunction with basic SCF

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

       skip_energy
            When set to True, the final energy is not computed.

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.

       **Returns:** the number of iterations
    '''
    if isinstance(wfn, RestrictedWFN):
        assert isinstance(ham, RestrictedEffectiveHamiltonian)
        return converge_scf_cs(ham, wfn, lf, overlap, maxiter, threshold)
    elif isinstance(wfn, UnrestrictedWFN):
        assert isinstance(ham, UnrestrictedEffectiveHamiltonian)
        return converge_scf_os(ham, wfn, lf, overlap, maxiter, threshold)
    else:
        raise NotImplementedError


def converge_scf_cs(ham, wfn, lf, overlap, maxiter=128, threshold=1e-8, skip_energy=False):
    '''Minimize the energy of the wavefunction with basic closed-shell SCF

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

       skip_energy
            When set to True, the final energy is not computed.

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.

       **Returns:** the number of iterations
    '''
    if log.do_medium:
        log('Starting restricted closed-shell SCF')
        log.hline()
        log('Iter  Error(alpha)')
        log.hline()

    fock = lf.create_one_body()
    converged = False
    counter = 0
    while maxiter is None or counter < maxiter:
        # feed the latest density matrices in the hamiltonian
        ham.reset(wfn.dm_alpha)
        # Construct the Fock operator
        fock.clear()
        ham.compute_fock(fock)
        # Check for convergence
        error = lf.error_eigen(fock, overlap, wfn.exp_alpha)

        if log.do_medium:
            log('%4i  %12.5e' % (counter, error))

        if error < threshold:
            converged = True
            break
        # Diagonalize the fock operator
        wfn.clear() # discard previous wfn state
        wfn.update_exp(fock, overlap)
        # counter
        counter += 1

    if log.do_medium:
        log.blank()

    if not skip_energy:
        ham.compute()
        if log.do_medium:
            ham.log()

    if not converged:
        raise NoSCFConvergence

    return counter


def converge_scf_os(ham, wfn, lf, overlap, maxiter=128, threshold=1e-8, skip_energy=False):
    '''Minimize the energy of the wavefunction with basic open-shell SCF

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

       skip_energy
            When set to True, the final energy is not computed.

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.
    '''
    if log.do_medium:
        log('Starting unrestricted open-shell SCF')
        log.hline()
        log('Iter  Error(alpha)  Error(beta)')
        log.hline()
    fock_alpha = lf.create_one_body()
    fock_beta = lf.create_one_body()
    converged = False
    counter = 0
    while maxiter is None or counter < maxiter:
        # feed the latest density matrices in the hamiltonian
        ham.reset(wfn.dm_alpha, wfn.dm_beta)
        # Construct the Fock operators
        fock_alpha.clear()
        fock_beta.clear()
        ham.compute_fock(fock_alpha, fock_beta)
        # Check for convergence
        error_alpha = lf.error_eigen(fock_alpha, overlap, wfn.exp_alpha)
        error_beta = lf.error_eigen(fock_beta, overlap, wfn.exp_beta)

        if log.do_medium:
            log('%4i  %12.5e  %12.5e' % (counter, error_alpha, error_beta))

        if error_alpha < threshold and error_beta < threshold:
            converged = True
            break
        # Diagonalize the fock operators
        wfn.clear()
        wfn.update_exp(fock_alpha, fock_beta, overlap)
        # counter
        counter += 1

    if log.do_medium:
        log.blank()

    if not skip_energy:
        ham.compute()
        if log.do_medium:
            ham.log()

    if not converged:
        raise NoSCFConvergence

    return counter
