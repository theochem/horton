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


__all__ = ['guess_hamiltonian_core']


@timer.with_section('Initial Guess')
def guess_hamiltonian_core(system):
    if log.do_medium:
        log('Performing a hamiltonian core guess.')
        log.blank()
    if isinstance(system.wfn, RestrictedWFN):
        guess_hamiltonian_core_cs(system)
    elif isinstance(system.wfn, UnrestrictedWFN):
        guess_hamiltonian_core_os(system)
    else:
        raise NotImplementedError


def guess_hamiltonian_core_cs(system):
    overlap = system.get_overlap()
    hamcore = system.lf.create_one_body()
    hamcore.iadd(system.get_kinetic())
    hamcore.iadd(system.get_nuclear_attraction())
    system.wfn.clear()
    system.wfn.update_exp(hamcore, overlap)


def guess_hamiltonian_core_os(system):
    overlap = system.get_overlap()
    hamcore = system.lf.create_one_body()
    hamcore.iadd(system.get_kinetic())
    hamcore.iadd(system.get_nuclear_attraction())
    system.wfn.clear()
    system.wfn.update_exp(hamcore, hamcore, overlap)
