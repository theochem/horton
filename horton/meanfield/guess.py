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
'''Initial guesses for wavefunctions'''


from horton.log import log, timer
from horton.matrix.base import TwoIndex, Expansion


__all__ = ['guess_core_hamiltonian']


@timer.with_section('Initial Guess')
def guess_core_hamiltonian(overlap, *args, **kwargs):
    '''Guess the orbitals by diagonalizing a core Hamiltonian

       **Arguments:**

       overlap
            The overlap operator.

       core1, core2, ...
            A number of operators that add up to the core Hamiltonian. Any set
            of operators whose sum resembles a Fock operator is fine. Usually,
            one passes the kinetic energy and nuclear attraction integrals.

       exp1, exp2, ...
            A list of wavefunction expansion objects (output arguments)

       This method only modifies the expansion coefficients and the orbital
       energies.
    '''
    if len(kwargs) != 0:
        raise TypeError('Unknown keyword arguments: %s' % kwargs.keys())

    if log.do_medium:
        log('Performing a core Hamiltonian guess.')
        log.blank()

    core = []
    exps = []
    for arg in args:
        if isinstance(arg, TwoIndex):
            core.append(arg)
        elif isinstance(arg, Expansion):
            exps.append(arg)
        else:
            raise TypeError('argument of unsupported type: %s' % arg)

    if len(core) == 0:
        raise TypeError('At least one term is needed for the core Hamiltonian.')
    if len(exps) == 0:
        raise TypeError('At least one wavefunction expansion is needed.')

    # Take sum of operators for core hamiltonian
    hamcore = core[0].copy()
    for term in core[1:]:
        hamcore.iadd(term)

    # Compute orbitals.
    exps[0].from_fock(hamcore, overlap)
    # Copy to other expansions.
    for i in xrange(1, len(exps)):
        exps[i].coeffs[:] = exps[0].coeffs
        exps[i].energies[:] = exps[0].energies
