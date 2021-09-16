# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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


__all__ = ['guess_core_hamiltonian']


@timer.with_section('Initial Guess')
def guess_core_hamiltonian(overlap, core, *orbs):
    '''Guess the orbitals by diagonalizing a core Hamiltonian.

    Parameters
    ----------
    overlap : np.ndarray, shape=(nbasis, nbasis), dtype=float
        The overlap operator.
    core : np.ndarray, shape=(nbasis, nbasis), dtype=float
        The core Hamiltonian. operator that resembles a Fock operator is fine. Usually,
        one adds the kinetic energy and nuclear attraction integrals.
    orb1, orb2, ... : Orbitals
        A list of Orbitals objects (output arguments)

    This method only modifies the expansion coefficients and the orbital energies.
    '''
    if log.do_medium:
        log('Performing a core Hamiltonian guess.')
        log.blank()

    if len(orbs) == 0:
        raise TypeError('At least one set of orbitals.')

    # Compute orbitals.
    orbs[0].from_fock(core, overlap)
    # Copy to other Orbitals objects.
    for i in range(1, len(orbs)):
        orbs[i].coeffs[:] = orbs[0].coeffs
        orbs[i].energies[:] = orbs[0].energies
