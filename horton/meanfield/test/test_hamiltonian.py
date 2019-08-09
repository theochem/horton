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


from horton.grid import BeckeMolGrid
from .common import check_interpolation, helper_compute, load_mdata, load_kin, load_na, load_er, \
    load_nn, load_orbs_alpha, load_orbs_beta, get_obasis, load_olp
from .. import UTwoIndexTerm, UDirectTerm, UExchangeTerm, UEffHam, RTwoIndexTerm, RDirectTerm, \
    RGridGroup, RDiracExchange, REffHam, RExchangeTerm, PlainSCFSolver, AufbauOccModel, \
    convergence_error_eigen


def test_energy_hydrogen():
    fname = 'h_sto3g_fchk'
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
        UTwoIndexTerm(na, 'ne'),
    ]
    external = {'nn': load_nn(fname)}
    ham = UEffHam(terms, external)
    helper_compute(ham, load_orbs_alpha(fname), load_orbs_beta(fname))
    print ham.cache['energy'] - -4.665818503844346E-01
    assert abs(ham.cache['energy'] - -4.665818503844346E-01) < 1e-8


def test_cubic_interpolation_hfs_cs():
    fname = 'water_hfs_321g_fchk'
    mdata = load_mdata(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RGridGroup(get_obasis(fname), grid, [
            RDiracExchange(),
        ]),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)

    check_interpolation(ham, olp, kin, na, [load_orbs_alpha(fname)])


def test_perturbation():
    fname = 'n2_hfs_sto3g_fchk'
    mdata = load_mdata(fname)
    scf_solver = PlainSCFSolver(maxiter=1024)

    # Without perturbation
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)
    occ_model = AufbauOccModel(7)

    orb_alpha = load_orbs_alpha(fname)

    assert convergence_error_eigen(ham, olp, orb_alpha) > 1e-8
    scf_solver(ham, olp, occ_model, orb_alpha)
    assert convergence_error_eigen(ham, olp, orb_alpha) < 1e-8
    energy0 = ham.compute_energy()

    # Construct a perturbation based on the Mulliken AIM operator
    assert get_obasis(fname).nbasis % 2 == 0
    nfirst = get_obasis(fname).nbasis / 2
    operator = load_olp(fname).copy()
    operator[:nfirst, nfirst:] *= 0.5
    operator[nfirst:, :nfirst] *= 0.5
    operator[nfirst:, nfirst:] = 0.0

    # Apply the perturbation with opposite signs and check that, because of
    # symmetry, the energy of the perturbed wavefunction is the same in both
    # cases, and higher than the unperturbed.
    energy1_old = None
    for scale in 0.1, -0.1:
        # Perturbation
        tmp = scale * operator
        perturbation = RTwoIndexTerm(tmp, 'pert')
        # Hamiltonian
        terms = [
            RTwoIndexTerm(kin, 'kin'),
            RDirectTerm(er, 'hartree'),
            RExchangeTerm(er, 'x_hf'),
            RTwoIndexTerm(na, 'ne'),
            perturbation,
        ]
        ham = REffHam(terms)
        orb_alpha = load_orbs_alpha(fname)
        assert convergence_error_eigen(ham, olp, orb_alpha) > 1e-8
        scf_solver(ham, olp, occ_model, orb_alpha)
        assert convergence_error_eigen(ham, olp, orb_alpha) < 1e-8
        energy1 = ham.compute_energy()
        energy1 -= ham.cache['energy_pert']

        assert energy1 > energy0
        if energy1_old is None:
            energy1_old = energy1
        else:
            assert abs(energy1 - energy1_old) < 1e-7


def test_ghost_hf():
    fname = 'water_dimer_ghost_fchk'
    mdata = load_mdata(fname)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)
    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham, olp, load_orbs_alpha(fname)) < 1e-5
