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


import numpy as np
from nose.plugins.attrib import attr

from .common import check_hf_cs_hf, check_lih_os_hf, \
    check_water_cs_hfs, check_n2_cs_hfs, check_h3_os_hfs, check_h3_os_pbe, \
    check_co_cs_pbe, check_water_cs_m05, \
    check_methyl_os_tpss, load_mdata, load_olp, load_kin, load_na, load_er, load_nn, \
    load_orbs_alpha, load_orbs_beta
from .. import EDIISSCFSolver, NoSCFConvergence, RTwoIndexTerm, RDirectTerm, RExchangeTerm, REffHam, \
    AufbauOccModel, UTwoIndexTerm, UDirectTerm, UExchangeTerm, UEffHam, guess_core_hamiltonian


def test_hf_cs_hf():
    check_hf_cs_hf(EDIISSCFSolver(threshold=1e-7))


def test_lih_os_hf():
    check_lih_os_hf(EDIISSCFSolver(threshold=1e-7))


def test_water_cs_hfs():
    check_water_cs_hfs(EDIISSCFSolver(threshold=1e-6))


@attr('slow')
def test_n2_cs_hfs():
    check_n2_cs_hfs(EDIISSCFSolver(threshold=1e-6))


@attr('slow')
def test_h3_os_hfs():
    check_h3_os_hfs(EDIISSCFSolver(threshold=1e-6))


@attr('slow')
def test_co_cs_pbe():
    check_co_cs_pbe(EDIISSCFSolver(threshold=1e-5))


@attr('slow')
def test_h3_os_pbe():
    check_h3_os_pbe(EDIISSCFSolver(threshold=1e-6))


# TODO: Move to higher level test
# def test_vanadium_sc_hf():
#     with assert_raises(NoSCFConvergence):
#         check_vanadium_sc_hf(EDIISSCFSolver(threshold=1e-10, maxiter=10))


def test_water_cs_m05():
    check_water_cs_m05(EDIISSCFSolver(threshold=1e-6))


def test_methyl_os_tpss():
    check_methyl_os_tpss(EDIISSCFSolver(threshold=1e-3))


def test_interpol_hf_cs_hf():
    fname = 'hf_sto3g_fchk'
    mdata = load_mdata(fname)

    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    external = {'nn': load_nn(fname)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    occ_model = AufbauOccModel(5)

    orbs = [load_orbs_alpha(fname)]
    check_interpol_hf(ham, orbs, olp, kin, na, occ_model)


def test_interpol_lih_os_hf():
    fname = 'li_h_3_21G_hf_g09_fchk'
    mdata = load_mdata(fname)

    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    external = {'nn': load_nn(fname)}
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham = UEffHam(terms, external)
    occ_model = AufbauOccModel(2, 1)

    orbs = [load_orbs_alpha(fname), load_orbs_beta(fname)]
    check_interpol_hf(ham, orbs, olp, kin, na, occ_model)


def check_interpol_hf(ham, orbs, olp, kin, na, occ_model):
    guess_core_hamiltonian(olp, kin + na, *orbs)
    dms = [exp.to_dm() for exp in orbs]
    scf_solver = EDIISSCFSolver(maxiter=4)
    try:
        scf_solver(ham, olp, occ_model, *dms)
    except NoSCFConvergence:
        pass

    # test harmonic approximation of the energy. This requires access to the
    # internals of the ediis solver.
    b, e = scf_solver._history._setup_equations()
    x = np.zeros(len(e))
    alphas = np.arange(0.0, 1.00001, 0.01)
    npt = len(alphas)
    energies_approx = np.zeros(npt)
    energies_hf = np.zeros(npt)
    for ipt in xrange(npt):
        x[0] = 1 - alphas[ipt]
        x[1] = alphas[ipt]
        energies_approx[ipt] = np.dot(x, 0.5 * np.dot(b, x) - e)
        # compute the hf energy
        scf_solver._history._build_combinations(x, dms, None)
        ham.reset(*dms)
        energies_hf[ipt] = ham.compute_energy()
    if False:
        import matplotlib.pyplot as pt
        pt.clf()
        pt.plot(alphas, energies_approx, 'k-', label='approx')
        pt.plot(alphas, energies_hf, 'r-', label='hf')
        pt.legend(loc=0)
        pt.savefig('foo.png')
    assert abs(energies_approx - energies_hf).max() < 1e-6
