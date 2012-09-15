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


from horton.wfn import ClosedShellWFN, OpenShellWFN


__all__ = ['converge_scf', 'convergence_error']


def converge_scf(ham, max_iter=128, threshold=1e-8, mixing=0.0):
    '''Minimize the energy of the wavefunction with basic SCF

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       max_iter
            The maximum number of iterations

       threshold
            The convergence threshold for the wavefunction

       mixing
            The amount of mixing between the old and the new Fock operator(s).

       **Returns:**

       converged
            True of converged, False if not
    '''
    ham.system.wfn.check_normalization(ham.system.get_overlap())
    if isinstance(ham.system.wfn, ClosedShellWFN):
        return converge_scf_cs(ham, max_iter, threshold)
    elif isinstance(ham.system.wfn, OpenShellWFN):
        return converge_scf_os(ham, max_iter, threshold)
    else:
        raise NotImplementedError


def converge_scf_cs(ham, max_iter=128, threshold=1e-8, mixing=0.0):
    '''Minimize the energy of the wavefunction with basic closed-shell SCF

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       max_iter
            The maximum number of iterations.

       threshold
            The convergence threshold for the wavefunction.

       mixing
            The amount of mixing between the old and the new Fock operator.

       **Returns:**

       converged
            True of converged, False if not.
    '''
    if mixing >= 1 or mixing < 0:
        raise ValueError('The mixing argument should be in the range [0,1[.')
    lf = ham.system.lf
    wfn = ham.system.wfn
    nbasis = ham.system.obasis.nbasis
    fock = lf.create_one_body(nbasis)
    old_fock = None
    converged = False
    for i in xrange(max_iter):
        # Construct the Fock operator
        fock.reset()
        ham.compute_fock(fock, None)
        if mixing > 0.0 and old_fock is not None:
            old_fock.iscale(mixing)
            fock.iscale(1-mixing)
            fock.iadd(old_fock)
        # Check for convergence
        error = lf.error_eigen(fock, ham.overlap, wfn.expansion)
        if error < threshold:
            converged = True
            break
        # Diagonalize the fock operator
        lf.diagonalize(fock, ham.overlap, wfn.expansion)
        # Let the hamiltonian know that the wavefunction has changed.
        ham.invalidate()
        # Write intermediate results to checkpoint
        ham.system.update_chk('wfn')
        # Take a backup copy of the Fock matrix
        if mixing > 0:
            if old_fock is None:
                old_fock = lf.create_one_body()
            old_fock.reset()
            old_fock.iadd(fock)
    return converged


def converge_scf_os(ham, max_iter=128, threshold=1e-8, mixing=0.0):
    '''Minimize the energy of the wavefunction with basic open-shell SCF

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       max_iter
            The maximum number of iterations.

       threshold
            The convergence threshold for the wavefunction.

       mixing
            The amount of mixing between the old and the new Fock operators.

       **Returns:**

       converged
            True of converged, False if not.
    '''
    if mixing >= 1 or mixing < 0:
        raise ValueError('The mixing argument should be in the range [0,1[.')
    lf = ham.system.lf
    wfn = ham.system.wfn
    nbasis = ham.system.obasis.nbasis
    fock_alpha = lf.create_one_body(nbasis)
    fock_beta = lf.create_one_body(nbasis)
    old_fock_alpha = None
    old_fock_beta = None
    converged = False
    for i in xrange(max_iter):
        # Construct the Fock operators
        fock_alpha.reset()
        fock_beta.reset()
        ham.compute_fock(fock_alpha, fock_beta)
        if mixing > 0.0 and old_fock_alpha is not None:
            old_fock_alpha.iscale(mixing)
            fock_alpha.iscale(1-mixing)
            fock_alpha.iadd(old_fock_alpha)
            old_fock_beta.iscale(mixing)
            fock_beta.iscale(1-mixing)
            fock_beta.iadd(old_fock_beta)
        # Check for convergence
        error_alpha = lf.error_eigen(fock_alpha, ham.overlap, wfn.alpha_expansion)
        error_beta = lf.error_eigen(fock_beta, ham.overlap, wfn.beta_expansion)
        if error_alpha < threshold and error_beta < threshold:
            converged = True
            break
        # Diagonalize the fock operators
        lf.diagonalize(fock_alpha, ham.overlap, wfn.alpha_expansion)
        lf.diagonalize(fock_beta, ham.overlap, wfn.beta_expansion)
        # Let the hamiltonian know that the wavefunction has changed.
        ham.invalidate()
        # Write intermediate results to checkpoint
        ham.system.update_chk('wfn')
        # Take backup copies of the Fock matrices.
        if mixing > 0:
            if old_fock is None:
                old_fock_alpha = lf.create_one_body()
                old_fock_beta = lf.create_one_body()
            old_fock_alpha.reset()
            old_fock_alpha.iadd(fock)
            old_fock_beta.reset()
            old_fock_beta.iadd(fock)
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
    nbasis = ham.system.obasis.nbasis
    fock = lf.create_one_body(nbasis)
    # Construct the Fock operator
    fock.reset()
    ham.compute_fock(fock, None)
    # Compute error
    return lf.error_eigen(fock, ham.overlap, wfn.expansion)


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
    nbasis = ham.system.obasis.nbasis
    fock_alpha = lf.create_one_body(nbasis)
    fock_beta = lf.create_one_body(nbasis)
    # Construct the Fock operators
    fock_alpha.reset()
    fock_beta.reset()
    ham.compute_fock(fock_alpha, fock_beta)
    # Compute errors
    error_alpha = lf.error_eigen(fock_alpha, ham.overlap, wfn.alpha_expansion)
    error_beta = lf.error_eigen(fock_beta, ham.overlap, wfn.beta_expansion)
    return max(error_alpha, error_beta)
