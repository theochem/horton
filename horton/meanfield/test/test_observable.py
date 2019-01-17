# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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
"""Unit tests for horton/meanfield/observable.py."""


from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_dot_hessian, \
    check_dot_hessian_polynomial, check_dot_hessian_cache


def setup_rhf_case():
    """Prepare datastructures for R-HF calculation on Water."""
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    mol = IOData.from_file(fn_fchk)
    mol.dm_alpha = mol.exp_alpha.to_dm()

    # RHF Effective Hamiltonian
    olp = mol.obasis.compute_overlap(mol.lf)
    core = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    core.iadd(na)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RTwoIndexTerm(core, 'core'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
    ]
    ham = REffHam(terms)

    return mol, olp, core, ham


def test_dot_hessian_rhf_polynomial():
    mol, olp, core, ham = setup_rhf_case()
    check_dot_hessian_polynomial(olp, core, ham, [mol.exp_alpha])


def test_dot_hessian_rhf_fd():
    mol, _olp, _core, ham = setup_rhf_case()
    check_dot_hessian(ham, mol.lf, mol.dm_alpha)


def test_cache_dot_hessian_rhf():
    mol, _olp, _core, ham = setup_rhf_case()
    check_dot_hessian_cache(ham, mol.dm_alpha)


def setup_uhf_case():
    """Prepare datastructures for UHF calculation."""
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    mol.dm_alpha = mol.exp_alpha.to_dm()
    mol.dm_beta = mol.exp_beta.to_dm()

    # UHF Effective Hamiltonian
    olp = mol.obasis.compute_overlap(mol.lf)
    core = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    core.iadd(na)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        UTwoIndexTerm(core, 'core'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
    ]
    ham = UEffHam(terms)

    return mol, olp, core, ham


def test_dot_hessian_uhf_polynomial():
    mol, olp, core, ham = setup_uhf_case()
    check_dot_hessian_polynomial(olp, core, ham, [mol.exp_alpha, mol.exp_beta])


def test_dot_hessian_uhf_fd():
    mol, _olp, _core, ham = setup_uhf_case()
    check_dot_hessian(ham, mol.lf, mol.dm_alpha, mol.dm_beta)


def test_cache_dot_hessian_uhf():
    mol, _olp, _core, ham = setup_uhf_case()
    check_dot_hessian_cache(ham, mol.dm_alpha, mol.dm_beta)
