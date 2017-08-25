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
"""Unit tests for horton/meanfield/observable.py."""

from .common import check_dot_hessian, \
    check_dot_hessian_polynomial, check_dot_hessian_cache, load_orbsa_dms, load_orbsb_dms, load_olp, \
    load_kin, load_na, load_er, load_er_chol, load_orbs_alpha, load_orbs_beta

from .. import RTwoIndexTerm, RDirectTerm, RExchangeTerm, REffHam, UTwoIndexTerm, UDirectTerm, \
    UExchangeTerm, UEffHam


def setup_rhf_case(cholesky=False):
    """Prepare data structures for R-HF calculation on Water."""
    fname = 'water_sto3g_hf_g03_fchk'
    dm_alpha = load_orbsa_dms(fname)
    orb_alpha = load_orbs_alpha(fname)

    # RHF Effective Hamiltonian
    olp = load_olp(fname)
    core = load_kin(fname)
    core += load_na(fname)
    if cholesky:
        er = load_er_chol(fname)
    else:
        er = load_er(fname)
    terms = [
        RTwoIndexTerm(core, 'core'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
    ]
    ham = REffHam(terms)

    return dm_alpha, olp, core, ham, orb_alpha


def test_dot_hessian_rhf_polynomial():
    dma, olp, core, ham, orb_alpha = setup_rhf_case()
    check_dot_hessian_polynomial(olp, core, ham, [orb_alpha])


def test_dot_hessian_rhf_fd():
    dma, olp, core, ham, orb_alpha = setup_rhf_case()
    check_dot_hessian(ham, dma)


def test_cache_dot_hessian_rhf():
    dma, olp, core, ham, orb_alpha = setup_rhf_case()
    check_dot_hessian_cache(ham, dma)


def test_dot_hessian_rhf_polynomial_cholesky():
    dma, olp, core, ham, orb_alpha = setup_rhf_case(True)
    check_dot_hessian_polynomial(olp, core, ham, [orb_alpha])


def test_dot_hessian_rhf_fd_cholesky():
    dma, olp, core, ham, orb_alpha = setup_rhf_case(True)
    check_dot_hessian(ham, dma)


def test_cache_dot_hessian_rhf_cholesky():
    dma, olp, core, ham, orb_alpha = setup_rhf_case(True)
    check_dot_hessian_cache(ham, dma)


def setup_uhf_case(cholesky=False):
    """Prepare data structures for UHF calculation."""
    fname = 'h3_hfs_321g_fchk'
    dma = load_orbsa_dms(fname)
    dmb = load_orbsb_dms(fname)

    orb_alpha = load_orbs_alpha(fname)
    orb_beta = load_orbs_beta(fname)

    # UHF Effective Hamiltonian
    olp = load_olp(fname)
    core = load_kin(fname)
    core += load_na(fname)
    er = load_er(fname)
    terms = [
        UTwoIndexTerm(core, 'core'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
    ]
    ham = UEffHam(terms)

    return dma, dmb, olp, core, ham, orb_alpha, orb_beta


def test_dot_hessian_uhf_polynomial():
    dma, dmb, olp, core, ham, orb_alpha, orb_beta = setup_uhf_case()
    check_dot_hessian_polynomial(olp, core, ham, [orb_alpha, orb_beta])


def test_dot_hessian_uhf_fd():
    dma, dmb, olp, core, ham, orb_alpha, orb_beta = setup_uhf_case()
    check_dot_hessian(ham, dma, dmb)


def test_cache_dot_hessian_uhf():
    dma, dmb, olp, core, ham, orb_alpha, orb_beta = setup_uhf_case()
    check_dot_hessian_cache(ham, dma, dmb)


def test_dot_hessian_uhf_polynomial_cholesky():
    dma, dmb, olp, core, ham, orb_alpha, orb_beta = setup_uhf_case(True)
    check_dot_hessian_polynomial(olp, core, ham, [orb_alpha, orb_beta])


def test_dot_hessian_uhf_fd_cholesky():
    dma, dmb, olp, core, ham, orb_alpha, orb_beta = setup_uhf_case(True)
    check_dot_hessian(ham, dma, dmb)


def test_cache_dot_hessian_uhf_cholesky():
    dma, dmb, olp, core, ham, orb_alpha, orb_beta = setup_uhf_case(True)
    check_dot_hessian_cache(ham, dma, dmb)
