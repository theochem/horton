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


import numpy as np

from horton.cext import compute_nucnuc
from horton.context import context
from horton.gbasis.gobasis import get_gobasis
from horton.grid.molgrid import BeckeMolGrid
from horton.io.iodata import IOData
from horton.log import log
from horton.matrix.dense import DenseLinalgFactory
from horton.meanfield.builtin import RDiracExchange, UDiracExchange
from horton.meanfield.convergence import convergence_error_eigen
from horton.meanfield.gridgroup import RGridGroup, UGridGroup
from horton.meanfield.guess import guess_core_hamiltonian
from horton.meanfield.hamiltonian import REffHam, UEffHam
from horton.meanfield.libxc import RLibXCLDA, ULibXCLDA, RLibXCGGA, ULibXCGGA
from horton.meanfield.observable import RTwoIndexTerm, RDirectTerm, RExchangeTerm
from horton.meanfield.observable import UTwoIndexTerm, UDirectTerm, UExchangeTerm
from horton.meanfield.occ import AufbauOccModel, FixedOccModel
from horton.meanfield.scf_oda import check_cubic


__all__ = [
    'check_cubic_wrapper', 'check_interpolation', 'check_solve', 'helper_compute',
    'check_hf_cs_hf', 'check_lih_os_hf', 'check_water_cs_hfs',
    'check_n2_cs_hfs', 'check_h3_os_hfs', 'check_h3_os_pbe', 'check_co_cs_pbe',
    'check_vanadium_sc_hf',
]


def check_cubic_wrapper(ham, dm0s, dm1s, do_plot=False):
    focks = [dm0.new() for dm0 in dm0s]

    # evaluate stuff at dm0
    ham.reset(*dm0s)
    e0 = ham.compute_energy()
    ham.compute_fock(*focks)
    g0 = 0.0
    for i in xrange(ham.ndm):
        g0 += focks[i].contract_two('ab,ba', dm1s[i])
        g0 -= focks[i].contract_two('ab,ba', dm0s[i])
    g0 *= ham.deriv_scale

    # evaluate stuff at dm1
    ham.reset(*dm1s)
    e1 = ham.compute_energy()
    ham.compute_fock(*focks)
    g1 = 0.0
    for i in xrange(ham.ndm):
        g1 += focks[i].contract_two('ab,ba', dm1s[i])
        g1 -= focks[i].contract_two('ab,ba', dm0s[i])
    g1 *= ham.deriv_scale

    check_cubic(ham, dm0s, dm1s, e0, e1, g0, g1, do_plot)


def check_interpolation(ham, lf, olp, kin, na, exps, do_plot=False):
    dm0s = [exp.to_dm() for exp in exps]
    guess_core_hamiltonian(olp, kin, na, *exps)
    dm1s = [exp.to_dm() for exp in exps]
    check_cubic_wrapper(ham, dm0s, dm1s, do_plot)


def check_solve(ham, scf_solver, occ_model, lf, olp, kin, na, *exps):
    guess_core_hamiltonian(olp, kin, na, *exps)
    if scf_solver.kind == 'exp':
        occ_model.assign(*exps)
        assert scf_solver.error(ham, lf, olp, *exps) > scf_solver.threshold
        scf_solver(ham, lf, olp, occ_model, *exps)
        assert scf_solver.error(ham, lf, olp, *exps) < scf_solver.threshold
    else:
        occ_model.assign(*exps)
        dms = [exp.to_dm() for exp in exps]
        assert scf_solver.error(ham, lf, olp, *dms) > scf_solver.threshold
        scf_solver(ham, lf, olp, occ_model, *dms)
        assert scf_solver.error(ham, lf, olp, *dms) < scf_solver.threshold
        focks = [lf.create_two_index() for i in xrange(ham.ndm)]
        ham.compute_fock(*focks)
        for i in xrange(ham.ndm):
            exps[i].from_fock(focks[i], olp)
        occ_model.assign(*exps)


def helper_compute(ham, lf, *exps):
    # Test energy before scf
    dms = [exp.to_dm() for exp in exps]
    ham.reset(*dms)
    ham.compute_energy()
    focks = [lf.create_two_index() for exp in exps]
    ham.compute_fock(*focks)
    return ham.cache['energy'], focks


@log.with_level(log.high)
def check_hf_cs_hf(scf_solver):
    fn_fchk = context.get_fn('test/hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)

    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    occ_model = AufbauOccModel(5)

    check_solve(ham, scf_solver, occ_model, mol.lf, olp, kin, na, mol.exp_alpha)

    # test orbital energies
    expected_energies = np.array([
        -2.59083334E+01, -1.44689996E+00, -5.57467136E-01, -4.62288194E-01,
        -4.62288194E-01, 5.39578910E-01,
    ])
    assert abs(mol.exp_alpha.energies - expected_energies).max() < 1e-5

    ham.compute_energy()
    # compare with g09
    assert abs(ham.cache['energy'] - -9.856961609951867E+01) < 1e-8
    assert abs(ham.cache['energy_kin'] - 9.766140786239E+01) < 2e-7
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_x_hf'] - 4.561984106482E+01) < 1e-6
    assert abs(ham.cache['energy_ne'] - -2.465756615329E+02) < 1e-6
    assert abs(ham.cache['energy_nn'] - 4.7247965053) < 1e-8


@log.with_level(log.high)
def check_lih_os_hf(scf_solver):
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = IOData.from_file(fn_fchk)

    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham = UEffHam(terms, external)
    occ_model = AufbauOccModel(2, 1)

    check_solve(ham, scf_solver, occ_model, mol.lf, olp, kin, na, mol.exp_alpha, mol.exp_beta)

    expected_alpha_energies = np.array([
        -2.76116635E+00, -7.24564188E-01, -1.79148636E-01, -1.28235698E-01,
        -1.28235698E-01, -7.59817520E-02, -1.13855167E-02, 6.52484445E-03,
        6.52484445E-03, 7.52201895E-03, 9.70893294E-01,
    ])
    expected_beta_energies = np.array([
        -2.76031162E+00, -2.08814026E-01, -1.53071066E-01, -1.25264964E-01,
        -1.25264964E-01, -1.24605870E-02, 5.12761388E-03, 7.70499854E-03,
        7.70499854E-03, 2.85176080E-02, 1.13197479E+00,
    ])
    assert abs(mol.exp_alpha.energies - expected_alpha_energies).max() < 1e-5
    assert abs(mol.exp_beta.energies - expected_beta_energies).max() < 1e-5

    ham.compute_energy()
    # compare with g09
    assert abs(ham.cache['energy'] - -7.687331212191962E+00) < 1e-8
    assert abs(ham.cache['energy_kin'] - 7.640603924034E+00) < 2e-7
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_x_hf'] - 2.114420907894E+00) < 1e-7
    assert abs(ham.cache['energy_ne'] - -1.811548789281E+01) < 2e-7
    assert abs(ham.cache['energy_nn'] - 0.6731318487) < 1e-8


@log.with_level(log.high)
def check_water_cs_hfs(scf_solver):
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RGridGroup(mol.obasis, grid, [
            RDiracExchange(),
        ]),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file and different integration grids:
    assert convergence_error_eigen(ham, mol.lf, olp, mol.exp_alpha) < 3e-5

    # Recompute the orbitals and orbital energies. This should be reasonably OK.
    dm_alpha = mol.exp_alpha.to_dm()
    ham.reset(dm_alpha)
    ham.compute_energy()
    fock_alpha = mol.lf.create_two_index()
    ham.compute_fock(fock_alpha)
    mol.exp_alpha.from_fock(fock_alpha, olp)

    expected_energies = np.array([
        -1.83691041E+01, -8.29412411E-01, -4.04495188E-01, -1.91740814E-01,
        -1.32190590E-01, 1.16030419E-01, 2.08119657E-01, 9.69825207E-01,
        9.99248500E-01, 1.41697384E+00, 1.47918828E+00, 1.61926596E+00,
        2.71995350E+00
    ])

    assert abs(mol.exp_alpha.energies - expected_energies).max() < 2e-4
    assert abs(ham.cache['energy_ne'] - -1.977921986200E+02) < 1e-7
    assert abs(ham.cache['energy_kin'] - 7.525067610865E+01) < 1e-9
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_x_dirac'] - 3.864299848058E+01) < 2e-4
    assert abs(ham.cache['energy'] - -7.474134898935590E+01) < 2e-4
    assert abs(ham.cache['energy_nn'] - 9.1571750414) < 2e-8

    # Converge from scratch and check energies
    occ_model = AufbauOccModel(5)
    check_solve(ham, scf_solver, occ_model, mol.lf, olp, kin, na, mol.exp_alpha)

    ham.compute_energy()
    assert abs(ham.cache['energy_ne'] - -1.977921986200E+02) < 1e-4
    assert abs(ham.cache['energy_kin'] - 7.525067610865E+01) < 1e-4
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_x_dirac'] - 3.864299848058E+01) < 2e-4
    assert abs(ham.cache['energy'] - -7.474134898935590E+01) < 2e-4


@log.with_level(log.high)
def check_n2_cs_hfs(scf_solver):
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'veryfine', random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}

    libxc_term = RLibXCLDA('x')
    terms1 = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RGridGroup(mol.obasis, grid, [libxc_term]),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham1 = REffHam(terms1, external)

    builtin_term = RDiracExchange()
    terms2 = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RGridGroup(mol.obasis, grid, [builtin_term]),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham2 = REffHam(terms2, external)

    # Compare the potential computed by libxc with the builtin implementation
    energy1, focks1 = helper_compute(ham1, mol.lf, mol.exp_alpha)
    energy2, focks2 = helper_compute(ham2, mol.lf, mol.exp_alpha)
    libxc_pot = ham1.cache.load('pot_libxc_lda_x_alpha')
    builtin_pot = ham2.cache.load('pot_x_dirac_alpha')
    # Libxc apparently approximates values of the potential below 1e-4 with zero.
    assert abs(libxc_pot - builtin_pot).max() < 1e-4
    # Check of the libxc energy matches our implementation
    assert abs(energy1 - energy2) < 1e-10
    ex1 = ham1.cache['energy_libxc_lda_x']
    ex2 = ham2.cache['energy_x_dirac']
    assert abs(ex1 - ex2) < 1e-10

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham1, mol.lf, olp, mol.exp_alpha) < 1e-5
    assert convergence_error_eigen(ham2, mol.lf, olp, mol.exp_alpha) < 1e-5

    occ_model = AufbauOccModel(7)
    for ham in ham1, ham2:
        # Converge from scratch
        check_solve(ham, scf_solver, occ_model, mol.lf, olp, kin, na, mol.exp_alpha)

        # test orbital energies
        expected_energies = np.array([
            -1.37107053E+01, -1.37098006E+01, -9.60673085E-01, -3.57928483E-01,
            -3.16017655E-01, -3.16017655E-01, -2.12998316E-01, 6.84030479E-02,
            6.84030479E-02, 7.50192517E-01,
        ])
        assert abs(mol.exp_alpha.energies - expected_energies).max() < 3e-5

        ham.compute_energy()
        assert abs(ham.cache['energy_ne'] - -2.981579553570E+02) < 1e-5
        assert abs(ham.cache['energy_kin'] - 1.061620887711E+02) < 1e-5
        assert abs(ham.cache['energy'] - -106.205213597) < 1e-4
        assert abs(ham.cache['energy_nn'] - 23.3180604505) < 1e-8
    assert abs(ham1.cache['energy_hartree'] + ham1.cache['energy_libxc_lda_x'] - 6.247259253877E+01) < 1e-4
    assert abs(ham2.cache['energy_hartree'] + ham2.cache['energy_x_dirac'] - 6.247259253877E+01) < 1e-4


@log.with_level(log.high)
def check_h3_os_hfs(scf_solver):
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'veryfine', random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}

    libxc_term = ULibXCLDA('x')
    terms1 = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UGridGroup(mol.obasis, grid, [libxc_term]),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham1 = UEffHam(terms1, external)

    builtin_term = UDiracExchange()
    terms2 = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UGridGroup(mol.obasis, grid, [builtin_term]),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham2 = UEffHam(terms2, external)

    # Compare the potential computed by libxc with the builtin implementation
    energy1, focks1 = helper_compute(ham1, mol.lf, mol.exp_alpha, mol.exp_beta)
    energy2, focks2 = helper_compute(ham2, mol.lf, mol.exp_alpha, mol.exp_beta)
    libxc_pot = ham1.cache.load('pot_libxc_lda_x_both')[:,0]
    builtin_pot = ham2.cache.load('pot_x_dirac_alpha')
    # Libxc apparently approximates values of the potential below 1e-4 with zero.
    assert abs(libxc_pot - builtin_pot).max() < 1e-4
    # Check of the libxc energy matches our implementation
    assert abs(energy1 - energy2) < 1e-10
    ex1 = ham1.cache['energy_libxc_lda_x']
    ex2 = ham2.cache['energy_x_dirac']
    assert abs(ex1 - ex2) < 1e-10

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham1, mol.lf, olp, mol.exp_alpha, mol.exp_beta) < 1e-5
    assert convergence_error_eigen(ham2, mol.lf, olp, mol.exp_alpha, mol.exp_beta) < 1e-5

    occ_model = AufbauOccModel(2, 1)
    for ham in ham1, ham2:
        # Converge from scratch
        check_solve(ham, scf_solver, occ_model, mol.lf, olp, kin, na, mol.exp_alpha, mol.exp_beta)

        # test orbital energies
        expected_energies = np.array([
            -4.93959157E-01, -1.13961330E-01, 2.38730924E-01, 7.44216538E-01,
            8.30143356E-01, 1.46613581E+00
        ])
        assert abs(mol.exp_alpha.energies - expected_energies).max() < 1e-5
        expected_energies = np.array([
            -4.34824166E-01, 1.84114514E-04, 3.24300545E-01, 7.87622756E-01,
            9.42415831E-01, 1.55175481E+00
        ])
        assert abs(mol.exp_beta.energies - expected_energies).max() < 1e-5

        ham.compute_energy()
        # compare with g09
        assert abs(ham.cache['energy_ne'] - -6.832069993374E+00) < 1e-5
        assert abs(ham.cache['energy_kin'] - 1.870784279014E+00) < 1e-5
        assert abs(ham.cache['energy'] - -1.412556114057104E+00) < 1e-5
        assert abs(ham.cache['energy_nn'] - 1.8899186021) < 1e-8

    assert abs(ham1.cache['energy_hartree'] + ham1.cache['energy_libxc_lda_x'] - 1.658810998195E+00) < 1e-6
    assert abs(ham2.cache['energy_hartree'] + ham2.cache['energy_x_dirac'] - 1.658810998195E+00) < 1e-6


@log.with_level(log.high)
def check_co_cs_pbe(scf_solver):
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'fine', random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RGridGroup(mol.obasis, grid, [
            RLibXCGGA('x_pbe'),
            RLibXCGGA('c_pbe'),
        ]),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)

    # Test energy before scf
    energy, focks = helper_compute(ham, mol.lf, mol.exp_alpha)
    assert abs(energy - -1.116465967841901E+02) < 1e-4

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham, mol.lf, olp, mol.exp_alpha) < 1e-5

    # Converge from scratch
    occ_model = AufbauOccModel(7)
    check_solve(ham, scf_solver, occ_model, mol.lf, olp, kin, na, mol.exp_alpha)

    # test orbital energies
    expected_energies = np.array([
         -1.86831122E+01, -9.73586915E+00, -1.03946082E+00, -4.09331776E-01,
         -3.48686522E-01, -3.48686522E-01, -2.06049056E-01, 5.23730418E-02,
         5.23730418E-02, 6.61093726E-01
    ])
    assert abs(mol.exp_alpha.energies - expected_energies).max() < 1e-2

    ham.compute_energy()
    # compare with g09
    assert abs(ham.cache['energy_ne'] - -3.072370116827E+02) < 1e-2
    assert abs(ham.cache['energy_kin'] - 1.103410779827E+02) < 1e-2
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_libxc_gga_x_pbe'] + ham.cache['energy_libxc_gga_c_pbe'] - 6.273115782683E+01) < 1e-2
    assert abs(ham.cache['energy'] - -1.116465967841901E+02) < 1e-4
    assert abs(ham.cache['energy_nn'] - 22.5181790889) < 1e-7


@log.with_level(log.high)
def check_h3_os_pbe(scf_solver):
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'veryfine', random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UGridGroup(mol.obasis, grid, [
            ULibXCGGA('x_pbe'),
            ULibXCGGA('c_pbe'),
        ]),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham = UEffHam(terms, external)

    # compute the energy before converging
    dm_alpha = mol.exp_alpha.to_dm()
    dm_beta = mol.exp_beta.to_dm()
    ham.reset(dm_alpha, dm_beta)
    ham.compute_energy()
    assert abs(ham.cache['energy'] - -1.593208400939354E+00) < 1e-5

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham, mol.lf, olp, mol.exp_alpha, mol.exp_beta) < 2e-6

    # Converge from scratch
    occ_model = AufbauOccModel(2, 1)
    check_solve(ham, scf_solver, occ_model, mol.lf, olp, kin, na, mol.exp_alpha, mol.exp_beta)

    # test orbital energies
    expected_energies = np.array([
        -5.41141676E-01, -1.56826691E-01, 2.13089637E-01, 7.13565167E-01,
        7.86810564E-01, 1.40663544E+00
    ])
    assert abs(mol.exp_alpha.energies - expected_energies).max() < 2e-5
    expected_energies = np.array([
        -4.96730336E-01, -5.81411249E-02, 2.73586652E-01, 7.41987185E-01,
        8.76161160E-01, 1.47488421E+00
    ])
    assert abs(mol.exp_beta.energies - expected_energies).max() < 2e-5

    ham.compute_energy()
    # compare with g09
    assert abs(ham.cache['energy_ne'] - -6.934705182067E+00) < 1e-5
    assert abs(ham.cache['energy_kin'] - 1.948808793424E+00) < 1e-5
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_libxc_gga_x_pbe'] + ham.cache['energy_libxc_gga_c_pbe'] - 1.502769385597E+00) < 1e-5
    assert abs(ham.cache['energy'] - -1.593208400939354E+00) < 1e-5
    assert abs(ham.cache['energy_nn'] - 1.8899186021) < 1e-8


@log.with_level(log.high)
def check_vanadium_sc_hf(scf_solver):
    """Try to converge the SCF for the neutral vanadium atom with fixe fractional occupations.

    Parameters
    ----------
    scf_solver : one of the SCFSolver types in HORTON
                 A configured SCF solver that must be tested.
    """
    # vanadium atoms
    numbers = np.array([23])
    pseudo_numbers = numbers.astype(float)
    coordinates = np.zeros((1, 3), float)

    # Simple basis set
    obasis = get_gobasis(coordinates, numbers, 'def2-tzvpd')

    # Dense matrices
    lf = DenseLinalgFactory(obasis.nbasis)

    # Compute integrals
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(coordinates, pseudo_numbers, lf)
    er = obasis.compute_electron_repulsion(lf)

    # Setup of restricted HF Hamiltonian
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)

    # Define fractional occupations of interest. (Spin-compensated case)
    occ_model = FixedOccModel(np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                        1.0, 1.0, 0.5]))

    # Allocate orbitals and make the initial guess
    exp_alpha = lf.create_expansion(obasis.nbasis)
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    # SCF test
    check_solve(ham, scf_solver, occ_model, lf, olp, kin, na, exp_alpha)
