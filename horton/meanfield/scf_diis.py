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
'''Abstract DIIS code used by the different DIIS implementations'''


import numpy as np

from horton.log import log
from horton.exceptions import NoSCFConvergence
from horton.meanfield.convergence import compute_commutator
from horton.meanfield.scf_oda import find_min_cubic, find_min_quadratic, check_cubic_cs


__all__ = [
    'CSSCFStep', 'CSCubicODAStep', 'CSQuadraticODAStep'
]



class CSSCFStep(object):
    def __init__(self, ham):
        self.ham = ham

    def __call__(self, energy0, fock, fock_interpolated=False, debug=False):
        '''
           Note that ``fock`` is also used as output argument.
        '''
        wfn = self.ham.system.wfn

        # TODO, avoid reallocation by putting this in a class
        overlap = self.ham.system.get_overlap()

        #if fock_interpolated:
        #    # Construct the Fock operator for the new DM
        #    fock.clear()
        #    self.ham.compute_fock(fock, None)

        # Construct the new DM (regular SCF step)
        wfn.clear()
        wfn.update_exp(fock, overlap)
        self.ham.clear()

        # Construct the Fock operator for the new DM
        fock.clear()
        self.ham.compute_fock(fock, None)

        # the mixing coefficient
        return 1.0



class CSCubicODAStep(CSSCFStep):
    def __call__(self, energy0, fock, fock_interpolated=False, debug=False):
        wfn = self.ham.system.wfn

        # TODO, avoid reallocation by putting this in a class
        overlap = self.ham.system.get_overlap()

        # keep copies of current state:
        dm0 = wfn.dm_alpha.copy()
        fock0 = fock.copy()
        if fock_interpolated:
            # if the fock matrix was interpolated, recompute it from the
            # interpolated density matrix.
            fock0.clear()
            self.ham.compute_fock(fock0, None)
        if energy0 is None:
            energy0 = self.ham.compute()

        wfn.clear()
        wfn.update_exp(fock0, overlap)
        # Let the hamiltonian know that the wavefunction has changed.
        self.ham.clear()

        # second point
        fock1 = fock0.copy()
        fock1.clear()
        self.ham.compute_fock(fock1, None)
        # Compute energy at new point
        energy1 = self.ham.compute()
        # take the density matrix
        dm1 = wfn.dm_alpha.copy()

        # Compute the derivatives of the energy towards lambda at edges 0 and 1
        ev_00 = fock0.expectation_value(dm0)
        ev_10 = fock1.expectation_value(dm0)
        ev_01 = fock0.expectation_value(dm1)
        ev_11 = fock1.expectation_value(dm1)
        energy0_deriv = 2*(ev_01-ev_00)
        energy1_deriv = 2*(ev_11-ev_10)

        # find the lambda that minimizes the cubic polynomial in the range [0,1]
        if log.do_high:
            log('           E0: % 10.3e     D0: % 10.3e' % (energy0, energy0_deriv))
            log('        E1-E0: % 10.3e     D1: % 10.3e' % (energy1-energy0, energy1_deriv))
        mixing = find_min_cubic(energy0, energy1, energy0_deriv, energy1_deriv)
        if mixing == 0:
            mixing = 0.1
        if log.do_high:
            log('       mixing: % 10.7f' % mixing)

        if debug:
            check_cubic_cs(self.ham, dm0, dm1, energy0, energy1, energy0_deriv, energy1_deriv)

        # E) Construct the new dm
        # Put the mixed dm in dm_old, which is local in this routine.
        dm1.iscale(mixing)
        dm1.iadd(dm0, factor=1-mixing)

        wfn.clear()
        wfn.update_dm('alpha', dm1)
        self.ham.clear()

        fock.clear()
        self.ham.compute_fock(fock, None)
        wfn.update_exp(fock, overlap, dm1)

        # the mixing coefficient
        return mixing


class CSQuadraticODAStep(CSSCFStep):
    def __call__(self, energy0, fock, fock_interpolated=False, debug=False):
        wfn = self.ham.system.wfn

        # TODO, avoid reallocation by putting this in a class
        overlap = self.ham.system.get_overlap()

        # keep copies of current state:
        dm0 = wfn.dm_alpha.copy()
        fock0 = fock.copy()
        if fock_interpolated:
            # if the fock matrix was interpolated, recompute it from the
            # interpolated density matrix.
            fock0.clear()
            self.ham.compute_fock(fock0, None)

        wfn.clear()
        wfn.update_exp(fock0, overlap)
        # Let the hamiltonian know that the wavefunction has changed.
        self.ham.clear()

        # second point
        fock1 = fock0.copy()
        fock1.clear()
        self.ham.compute_fock(fock1, None)
        # take the density matrix
        dm1 = wfn.dm_alpha.copy()

        # Compute the derivatives of the energy towards lambda at edges 0 and 1
        ev_00 = fock0.expectation_value(dm0)
        ev_10 = fock1.expectation_value(dm0)
        ev_01 = fock0.expectation_value(dm1)
        ev_11 = fock1.expectation_value(dm1)
        energy0_deriv = 2*(ev_01-ev_00)
        energy1_deriv = 2*(ev_11-ev_10)

        # find the lambda that minimizes the cubic polynomial in the range [0,1]
        if log.do_high:
            log('           D0: % 10.3e' % (energy0_deriv))
            log('           D1: % 10.3e' % (energy1_deriv))
        mixing = find_min_quadratic(energy0_deriv, energy1_deriv)
        if mixing == 0:
            mixing = 0.1
        if log.do_high:
            log('       mixing: % 10.7f' % mixing)

        # E) Construct the new dm
        # Put the mixed dm in dm_old, which is local in this routine.
        dm1.iscale(mixing)
        dm1.iadd(dm0, factor=1-mixing)

        wfn.clear()
        wfn.update_dm('alpha', dm1)
        self.ham.clear()

        fock.clear()
        self.ham.compute_fock(fock, None)
        wfn.update_exp(fock, overlap, dm1)

        # the mixing coefficient
        return mixing


def converge_scf_diis_cs(ham, DIISHistoryClass, maxiter=128, threshold=1e-6, nvector=6, prune_old_states=False, skip_energy=False, scf_step='regular'):
    '''Minimize the energy of the closed-shell wavefunction with EDIIS

       **Arguments:**

       ham
            A Hamiltonian instance.

       DIISHistoryClass
            A DIIS history class.

       **Optional arguments:**

       maxiter
            The maximum number of iterations. When set to None, the SCF loop
            will go one until convergence is reached.

       threshold
            The convergence threshold for the wavefunction

       prune_old_states
            When set to True, old states are pruned from the history when their
            coefficient is zero. Pruning starts at the oldest state and stops
            as soon as a state is encountered with a non-zero coefficient. Even
            if some newer states have a zero coefficient.

       skip_energy
            When set to True, the final energy is not computed. Note that some
            DIIS variants need to compute the energy anyway. for these methods
            this option is irrelevant.

       scf_step
            The type of SCF step to take after the interpolated states was
            create from the DIIS history. This can be 'regular', 'oda2' or
            'oda3'.

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.
    '''
    # allocated and define some one body operators
    lf = ham.system.lf
    wfn = ham.system.wfn
    overlap = ham.system.get_overlap()
    history = DIISHistoryClass(lf, nvector, overlap)
    fock = lf.create_one_body()
    dm = lf.create_one_body()

    # The scf step
    if scf_step == 'regular':
        cs_scf_step = CSSCFStep(ham)
    elif scf_step == 'oda2':
        cs_scf_step = CSQuadraticODAStep(ham)
    elif scf_step == 'oda3':
        cs_scf_step = CSCubicODAStep(ham)
    else:
        raise ValueError('scf_step argument not recognized: %s' % scf_step)

    # Get rid of outdated stuff
    ham.clear()

    if log.do_medium:
        log('Starting restricted closed-shell %s-SCF' % history.name)
        log.hline()
        log('Iter  Error(alpha) CN(alpha)  Last(alpha) nv Method          Energy       Change')
        log.hline()

    converged = False
    counter = 0
    while maxiter is None or counter < maxiter:
        # Construct the Fock operator from scratch:
        if history.nused == 0:
            # Update the Fock matrix
            fock.clear()
            ham.compute_fock(fock, None)
            # Put this state also in the history
            energy = ham.compute() if history.need_energy else None
            # Add the current fock+dm pair to the history
            error = history.add(energy, wfn.dm_alpha, fock)
            if log.do_high:
                log('          DIIS add')
            if error < threshold:
                converged = True
                break
            if log.do_high:
                log.blank()
            if log.do_medium:
                energy_str = ' '*20 if energy is None else '% 20.13f' % energy
                log('%4i %12.5e                         %2i   %20s' % (
                    counter, error, history.nused, energy_str
                ))
            if log.do_high:
                log.blank()
            fock_interpolated = False
        else:
            energy = None
            fock_interpolated = True

        # Construct a new density matrix and fock matrix (based on the current
        # density matrix or fock matrix). The fock matrix is modified in-place
        # while the density matrix is stored in ham.system.wfn
        mixing = cs_scf_step(energy, fock, fock_interpolated)
        if mixing == 0.0:
            converged = True
            break

        energy = ham.compute() if history.need_energy else None

        # Write intermediate results to checkpoint
        ham.system.update_chk('wfn')

        # Add the current (dm, fock) pair to the history
        if log.do_high:
            log('          DIIS add')
        error = history.add(energy, wfn.dm_alpha, fock)

        # break when converged
        if error < threshold:
            converged = True
            break

        if log.do_high:
            log.blank()
        if log.do_medium:
            energy_str = ' '*20 if energy is None else '% 20.13f' % energy
            log('%4i %12.5e                         %2i   %20s' % (
                counter, error, history.nused, energy_str
            ))
        if log.do_high:
            log.blank()

        # get extra/intra-polated Fock matrix
        while True:
            energy_approx, coeffs, cn, method = history.solve(dm, fock)
            #if coeffs[coeffs<0].sum() < -1:
            #    if log.do_high:
            #        log('          DIIS (coeffs too negative) -> drop %i and retry' % history.stack[0].identity)
            #    history.shrink()
            if history.nused <= 2:
                break
            if coeffs[-1] == 0.0:
                if log.do_high:
                    log('          DIIS (last coeff zero) -> drop %i and retry' % history.stack[0].identity)
                history.shrink()
            else:
                break
        
        if False and len(coeffs) == 2:
            tmp = dm.copy()
            import matplotlib.pyplot as pt
            xs = np.linspace(0.0, 1.0, 25)
            a, b = history._setup_equations()
            energies1 = []
            energies2 = []
            for x in xs:
                x_coeffs = np.array([1-x, x])
                energies1.append(np.dot(x_coeffs, 0.5*np.dot(a, x_coeffs) - b))
                history._build_combinations(x_coeffs, tmp, None)
                wfn.clear()
                wfn.update_dm('alpha', tmp)
                ham.clear()
                energies2.append(ham.compute())
                print x, energies1[-1], energies2[-1]
            pt.clf()
            pt.plot(xs, energies1, label='est')
            pt.plot(xs, energies2, label='ref')
            pt.axvline(coeffs[1], color='k')
            pt.legend(loc=0)
            pt.savefig('ediis2_test_%05i.png' % counter)

        wfn.clear()
        wfn.update_dm('alpha', dm)
        ham.clear()

        if energy_approx is not None:
            energy_change = energy_approx - min(state.energy for state in history.stack)
        else:
            energy_change = None

        # log
        if log.do_high:
            history.log(coeffs)

        if log.do_medium:
            change_str = ' '*10 if energy_change is None else '% 12.7f' % energy_change
            log('%4i              %10.3e %12.7f %2i %s                      %12s' % (
                counter, cn, coeffs[-1], history.nused, method,
                change_str
            ))

        if log.do_high:
            log.blank()

        if prune_old_states:
            # get rid of old states with zero coeff
            for i in xrange(history.nused):
                if coeffs[i] == 0.0:
                    if log.do_high:
                        log('          DIIS insignificant -> drop %i' % history.stack[0].identity)
                    history.shrink()
                else:
                    break

        # counter
        counter += 1

    if log.do_medium:
        if converged:
            log('%4i %12.5e (converged)' % (counter, error))
        log.blank()

    if not skip_energy or history.need_energy:
        if not history.need_energy:
            ham.compute()
        if log.do_medium:
            ham.log_energy()

    if not converged:
        raise NoSCFConvergence


class DIISState(object):
    '''A single record (vector) in a DIIS history object.'''
    def __init__(self, lf, work, overlap):
        '''
           **Arguments:**

           lf
                The LinalgFactor used to create the one-body operators.

           work
                A one body operator to be used as a temporary variable. This
                object is allocated by the history object.

           overlap
                The overlap matrix.
        '''
        # Not all of these need to be used.
        self.work = work
        self.overlap = overlap
        self.energy = np.nan
        self.norm = np.nan
        self.dm = lf.create_one_body()
        self.fock = lf.create_one_body()
        self.commutator = lf.create_one_body()
        self.identity = None # every state has a different id.

    def clear(self):
        '''Reset this record.'''
        self.energy = np.nan
        self.norm = np.nan
        self.dm.clear()
        self.fock.clear()
        self.commutator.clear()

    def assign(self, identity, energy, dm, fock):
        '''Assign a new state.

           **Arguments:**

           identity
                A unique id for the new state.

           energy
                The energy of the new state.

           dm
                The density matrix of the new state.

           fock
                The Fock matrix of the new state.
        '''
        self.identity = identity
        self.energy = energy
        self.dm.assign(dm)
        self.fock.assign(fock)
        compute_commutator(dm, fock, self.overlap, self.work, self.commutator)
        self.norm = self.commutator.expectation_value(self.commutator)


class DIISHistory(object):
    '''A base class of DIIS histories'''
    name = None
    need_energy = None

    def __init__(self, lf, nvector, overlap, dots_matrices):
        '''
           **Arguments:**

           lf
                The LinalgFactor used to create the one-body operators.

           nvector
                The maximum size of the history.

           overlap
                The overlap matrix of the system.

           dots_matrices
                Matrices in which dot products will be stored

           **Useful attributes:**

           used
                The actual number of vectors in the history.
        '''
        self.work = lf.create_one_body()
        self.stack = [DIISState(lf, self.work, overlap) for i in xrange(nvector)]
        self.overlap = overlap
        self.dots_matrices = dots_matrices
        self.nused = 0
        self.idcounter = 0
        self.commutator = lf.create_one_body()

    def _get_nvector(self):
        '''The maximum size of the history'''
        return len(self.stack)

    nvector = property(_get_nvector)

    def log(self, coeffs):
        eref = min(state.energy for state in self.stack[:self.nused])
        if eref is None:
            log('          DIIS history          norm         coeff         id')
            for i in xrange(self.nused):
                state = self.stack[i]
                log('          DIIS history  %12.5e  %12.7f   %8i' % (state.norm, coeffs[i], state.identity))
        else:
            log('          DIIS history          norm        energy         coeff         id')
            for i in xrange(self.nused):
                state = self.stack[i]
                log('          DIIS history  %12.5e  %12.5e  %12.7f   %8i' % (state.norm, state.energy-eref, coeffs[i], state.identity))
        log.blank()

    def shrink(self):
        '''Remove the oldest item from the history'''
        self.nused -= 1
        state = self.stack.pop(0)
        state.clear()
        self.stack.append(state)
        for dots in self.dots_matrices:
            dots[:-1] = dots[1:]
            dots[:,:-1] = dots[:,1:]
            dots[-1] = np.nan
            dots[:,-1] = np.nan

    def add(self, energy, dm, fock):
        '''Add new state to the history.

           **Arguments:**

           energy
                The energy of the new state.

           dm
                The density matrix of the new state.

           fock
                The Fock matrix of the new state.
        '''
        # There must be a free spot
        if self.nused == self.nvector:
            self.shrink()

        # assign dm and fock
        state = self.stack[self.nused]
        state.assign(self.idcounter, energy, dm, fock)
        self.idcounter += 1

        # prepare for next iteration
        self.nused += 1
        return np.sqrt(state.norm)

    def _build_combinations(self, coeffs, dm_output, fock_output):
        '''Construct a linear combination of density/fock matrices'''
        if dm_output is not None:
            self._linear_combination(coeffs, [self.stack[i].dm for i in xrange(self.nused)], dm_output)
        if fock_output is not None:
            self._linear_combination(coeffs, [self.stack[i].fock for i in xrange(self.nused)], fock_output)

    def _linear_combination(self, coeffs, ops, output):
        '''Make a linear combination of one-body objects'''
        output.clear()
        for i in xrange(self.nused):
            output.iadd(ops[i], factor=coeffs[i])
