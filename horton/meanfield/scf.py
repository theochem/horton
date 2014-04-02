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
from horton.meanfield.convergence import convergence_error_eigen


__all__ = ['PlainSCFSolver']


class PlainSCFSolver(object):
    kind = 'exp' # basic variable is the wfn expansion

    def __init__(self, threshold=1e-8, maxiter=128, skip_energy=False):
        self.maxiter = maxiter,
        self.threshold = threshold
        self.skip_energy = skip_energy

    @timer.with_section('SCF')
    def __call__(self, ham, lf, overlap, occ_model, *exps):
        # Some type checking
        if ham.ndm != len(exps):
            raise TypeError('The number of initial orbital expansions does not match the Hamiltonian.')
        # Impose the requested occupation numbers
        occ_model.assign(*exps)
        # Check the orthogonality of the orbitals
        for exp in exps:
            exp.check_normalization(overlap)

        if log.do_medium:
            log('Starting plain SCF solver. ndm=%i' % ham.ndm)
            log.hline()
            log('Iter         Error')
            log.hline()

        focks = [lf.create_one_body() for i in xrange(ham.ndm)]
        dms = [lf.create_one_body() for i in xrange(ham.ndm)]
        converged = False
        counter = 0
        while self.maxiter is None or counter < self.maxiter:
            # convert the orbital expansions to density matrices
            for i in xrange(ham.ndm):
                exps[i].to_dm(dms[i])
            # feed the latest density matrices in the hamiltonian
            ham.reset(*dms)
            # Construct the Fock operator
            for fock in focks:
                fock.clear()
            ham.compute_fock(*focks)
            # Check for convergence
            error = 0.0
            for i in xrange(ham.ndm):
                error += lf.error_eigen(focks[i], overlap, exps[i])
            if log.do_medium:
                log('%4i  %12.5e' % (counter, error))
            if error < self.threshold:
                converged = True
                break
            # Diagonalize the fock operators to obtain new orbitals and
            for i in xrange(ham.ndm):
                exps[i].from_fock(focks[i], overlap)
            # Assign new occupation numbers.
            occ_model.assign(*exps)
            # counter
            counter += 1

        if log.do_medium:
            log.blank()

        if not self.skip_energy:
            ham.compute()
            if log.do_medium:
                ham.log()

        if not converged:
            raise NoSCFConvergence

        return counter

    def error(self, ham, lf, overlap, *exps):
        return convergence_error_eigen(ham, lf, overlap, *exps)
