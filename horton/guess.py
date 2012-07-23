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
"""Initial guesses for wavefunctions"""


from horton.wfn import ClosedShellWFN, OpenShellWFN


__all__ = ['guess_hamiltionian_core']


def guess_hamiltionian_core(system):
    if isinstance(system.wfn, ClosedShellWFN):
        guess_hamiltonian_core_cs(system)
    elif isinstance(system.wfn, OpenShellWFN):
        guess_hamiltonian_core_os(system)
    else:
        raise NotImplementedError


def guess_hamiltonian_core_cs(system):
    overlap = system.get_kinetic()
    hamcore = system.lf.create_one_body(overlap.nbasis)
    hamcore.iadd(system.get_kinetic(), 1)
    hamcore.iadd(system.get_nuclear_attraction(), -1)
    system.lf.diagonalize(hamcore, overlap, system.wfn.expansion)


def guess_hamiltonian_core_os(system):
    overlap = system.get_kinetic()
    hamcore = system.lf.create_one_body(overlap.nbasis)
    hamcore.iadd(system.get_kinetic(), 1)
    hamcore.iadd(system.get_nuclear_attraction(), -1)
    system.lf.diagonalize(hamcore, overlap, system.wfn.alpha_expansion)
    # TODO: just copy the alpha guess to the beta
    system.lf.diagonalize(hamcore, overlap, system.wfn.beta_expansion)
