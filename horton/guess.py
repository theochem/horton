# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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
"""Initial guesses for wavefunctions"""


from horton.log import log, timer
from horton.wfn import ClosedShellWFN, OpenShellWFN


__all__ = ['guess_hamiltonian_core']


def guess_hamiltonian_core(system):
    if log.do_medium:
        log('Performing a hamiltonian core guess.')
    with timer.section('Initial Guess'):
        if isinstance(system.wfn, ClosedShellWFN):
            guess_hamiltonian_core_cs(system)
        elif isinstance(system.wfn, OpenShellWFN):
            guess_hamiltonian_core_os(system)
        else:
            raise NotImplementedError


def guess_hamiltonian_core_cs(system):
    overlap = system.get_overlap()
    hamcore = system.lf.create_one_body()
    hamcore.iadd(system.get_kinetic(), 1)
    hamcore.iadd(system.get_nuclear_attraction(), -1)
    system.wfn.invalidate()
    system.wfn.update_exp(hamcore, overlap)
    system.update_chk('wfn')


def guess_hamiltonian_core_os(system):
    overlap = system.get_overlap()
    hamcore = system.lf.create_one_body()
    hamcore.iadd(system.get_kinetic(), 1)
    hamcore.iadd(system.get_nuclear_attraction(), -1)
    system.wfn.invalidate()
    system.wfn.update_exp(hamcore, hamcore, overlap)
    system.update_chk('wfn')
