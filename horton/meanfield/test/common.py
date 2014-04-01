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
#pylint: skip-file


from horton import *


__all__ = [
    'check_cubic_cs_wrapper', 'check_cubic_os_wrapper',
    'check_scf_hf_cs_hf', 'check_scf_water_cs_hfs',
]


def check_cubic_cs_wrapper(ham, wfn, dm0, dm1, do_plot=False):
    fock = dm0.copy()
    fock.clear()

    # evaluate stuff at dm0
    ham.clear()
    wfn.clear()
    wfn.update_dm('alpha', dm0)
    e0 = ham.compute()
    fock.clear()
    ham.compute_fock(fock, None)
    ev_00 = fock.expectation_value(dm0)
    ev_01 = fock.expectation_value(dm1)
    g0 = 2*(ev_01 - ev_00)

    # evaluate stuff at dm1
    ham.clear()
    wfn.clear()
    wfn.update_dm('alpha', dm1)
    e1 = ham.compute()
    fock.clear()
    ham.compute_fock(fock, None)
    ev_10 = fock.expectation_value(dm0)
    ev_11 = fock.expectation_value(dm1)
    g1 = 2*(ev_11 - ev_10)

    check_cubic_cs(ham, wfn, dm0, dm1, e0, e1, g0, g1, do_plot)


def check_cubic_os_wrapper(ham, wfn, dma0, dmb0, dma1, dmb1, do_plot=False):
    focka = dma0.copy()
    focka.clear()
    fockb = dmb0.copy()
    fockb.clear()

    # evaluate stuff at 0
    ham.clear()
    wfn.clear()
    wfn.update_dm('alpha', dma0)
    wfn.update_dm('beta', dmb0)
    e0 = ham.compute()
    focka.clear()
    fockb.clear()
    ham.compute_fock(focka, fockb)
    ev_00 = focka.expectation_value(dma0) + fockb.expectation_value(dmb0)
    ev_01 = focka.expectation_value(dma1) + fockb.expectation_value(dmb1)
    g0 = (ev_01 - ev_00)

    # evaluate stuff at 1
    ham.clear()
    wfn.clear()
    wfn.update_dm('alpha', dma1)
    wfn.update_dm('beta', dmb1)
    e1 = ham.compute()
    focka.clear()
    fockb.clear()
    ham.compute_fock(focka, fockb)
    ev_10 = focka.expectation_value(dma0) + fockb.expectation_value(dmb0)
    ev_11 = focka.expectation_value(dma1) + fockb.expectation_value(dmb1)
    g1 = (ev_11 - ev_10)

    check_cubic_os(ham, wfn, dma0, dmb0, dma1, dmb1, e0, e1, g0, g1, do_plot)


@log.with_level(log.high)
def check_scf_hf_cs_hf(scf_wrapper):
    fn_fchk = context.get_fn('test/hf_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)

    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        OneBodyTerm(kin, mol.wfn, 'kin'),
        DirectTerm(er, mol.wfn, 'hartree'),
        ExchangeTerm(er, mol.wfn, 'x_hf'),
        OneBodyTerm(na, mol.wfn, 'ne'),
    ]
    ham = Hamiltonian(terms, external)

    guess_core_hamiltonian(mol.wfn, olp, kin, na)
    assert scf_wrapper.convergence_error(ham, mol.wfn, mol.lf, olp) > scf_wrapper.kwargs['threshold']
    scf_wrapper(ham, mol.wfn, mol.lf, olp)
    assert scf_wrapper.convergence_error(ham, mol.wfn, mol.lf, olp) < scf_wrapper.kwargs['threshold']

    # test orbital energies
    expected_energies = np.array([
        -2.59083334E+01, -1.44689996E+00, -5.57467136E-01, -4.62288194E-01,
        -4.62288194E-01, 5.39578910E-01,
    ])
    assert abs(mol.wfn.exp_alpha.energies - expected_energies).max() < 1e-5

    ham.compute()
    # compare with g09
    assert abs(ham.cache['energy'] - -9.856961609951867E+01) < 1e-8
    assert abs(ham.cache['energy_kin'] - 9.766140786239E+01) < 2e-7
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_x_hf'] - 4.561984106482E+01) < 1e-6
    assert abs(ham.cache['energy_ne'] - -2.465756615329E+02) < 1e-6
    assert abs(ham.cache['energy_nn'] - 4.7247965053) < 1e-8


@log.with_level(log.high)
def check_scf_water_cs_hfs(scf_wrapper):
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = Molecule.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        OneBodyTerm(kin, mol.wfn, 'kin'),
        DirectTerm(er, mol.wfn, 'hartree'),
        GridGroup(mol.obasis, grid, mol.wfn, [
            DiracExchange(mol.wfn),
        ]),
        OneBodyTerm(na, mol.wfn, 'ne'),
    ]
    ham = Hamiltonian(terms, external)

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file and different integration grids:
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 3e-5

    # The energies should also be in reasonable agreement. Repeated to check for
    # stupid bugs
    for i in xrange(2):
        ham.clear()
        ham.compute()
        expected_energies = np.array([
            -1.83691041E+01, -8.29412411E-01, -4.04495188E-01, -1.91740814E-01,
            -1.32190590E-01, 1.16030419E-01, 2.08119657E-01, 9.69825207E-01,
            9.99248500E-01, 1.41697384E+00, 1.47918828E+00, 1.61926596E+00,
            2.71995350E+00
        ])

        assert abs(mol.wfn.exp_alpha.energies - expected_energies).max() < 2e-4
        assert abs(ham.cache['energy_ne'] - -1.977921986200E+02) < 1e-7
        assert abs(ham.cache['energy_kin'] - 7.525067610865E+01) < 1e-9
        assert abs(ham.cache['energy_hartree'] + ham.cache['energy_exchange_dirac'] - 3.864299848058E+01) < 2e-4
        assert abs(ham.cache['energy'] - -7.474134898935590E+01) < 2e-4
        assert abs(ham.cache['energy_nn'] - 9.1571750414) < 2e-8

    # Converge from scratch
    guess_core_hamiltonian(mol.wfn, olp, kin, na)
    assert scf_wrapper.convergence_error(ham, mol.wfn, mol.lf, olp) > scf_wrapper.kwargs['threshold']
    scf_wrapper(ham, mol.wfn, mol.lf, olp)
    assert scf_wrapper.convergence_error(ham, mol.wfn, mol.lf, olp) < scf_wrapper.kwargs['threshold']

    assert abs(ham.cache['energy_ne'] - -1.977921986200E+02) < 1e-4
    assert abs(ham.cache['energy_kin'] - 7.525067610865E+01) < 1e-4
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_exchange_dirac'] - 3.864299848058E+01) < 2e-4
    assert abs(ham.cache['energy'] - -7.474134898935590E+01) < 2e-4
