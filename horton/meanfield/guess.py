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
'''Initial guesses for wavefunctions'''


from horton.log import log, timer
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = ['guess_core_hamiltonian']


@timer.with_section('Initial Guess')
def guess_core_hamiltonian(wfn, overlap, *core):
    '''Guess the wavefunction from a core hamiltonian

       **Arguments:**

       wfn
            The wavefunction in which the orbitals are stored.

       overlap
            The overlap operator.

       core1, core2, ...
            A number of operators that add up to the core Hamiltonian. Any sum
            of operators that resembles a Fock operator is fine.
    '''
    if log.do_medium:
        log('Performing a hamiltonian core guess.')
        log.blank()
    if len(core) == 0:
        raise TypeError('At least one term is needed for the core Hamiltonian.')

    # Take sum
    hamcore = core[0].copy()
    for term in core[1:]:
        hamcore.iadd(term)

    # Compute orbitals
    wfn.clear()
    if isinstance(wfn, RestrictedWFN):
        wfn.update_exp(hamcore, overlap)
    elif isinstance(wfn, UnrestrictedWFN):
        wfn.update_exp(hamcore, hamcore, overlap)
    else:
        raise NotImplementedError
