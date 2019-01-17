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


from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_interpolation, helper_compute


def test_energy_hydrogen():
    fn_fchk = context.get_fn('test/h_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
        UTwoIndexTerm(na, 'ne'),
    ]
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    ham = UEffHam(terms, external)
    helper_compute(ham, mol.lf, mol.exp_alpha, mol.exp_beta)
    assert abs(ham.cache['energy'] - -4.665818503844346E-01) < 1e-8


def test_cubic_interpolation_hfs_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RGridGroup(mol.obasis, grid, [
            RDiracExchange(),
        ]),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)

    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha])


def test_perturbation():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    scf_solver = PlainSCFSolver(maxiter=1024)

    # Without perturbation
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)
    occ_model = AufbauOccModel(7)

    assert convergence_error_eigen(ham, mol.lf, olp, mol.exp_alpha) > 1e-8
    scf_solver(ham, mol.lf, olp, occ_model, mol.exp_alpha)
    assert convergence_error_eigen(ham, mol.lf, olp, mol.exp_alpha) < 1e-8
    energy0 = ham.compute_energy()

    # Construct a perturbation based on the Mulliken AIM operator
    assert mol.obasis.nbasis % 2 == 0
    nfirst = mol.obasis.nbasis / 2
    operator = mol.obasis.compute_overlap(mol.lf).copy()
    operator._array[:nfirst,nfirst:] *= 0.5
    operator._array[nfirst:,:nfirst] *= 0.5
    operator._array[nfirst:,nfirst:] = 0.0

    # Apply the perturbation with oposite signs and check that, because of
    # symmetry, the energy of the perturbed wavefunction is the same in both
    # cases, and higher than the unperturbed.
    energy1_old = None
    for scale in 0.1, -0.1:
        # Perturbation
        tmp = operator.copy()
        tmp.iscale(scale)
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
        assert convergence_error_eigen(ham, mol.lf, olp, mol.exp_alpha) > 1e-8
        scf_solver(ham, mol.lf, olp, occ_model, mol.exp_alpha)
        assert convergence_error_eigen(ham, mol.lf, olp, mol.exp_alpha) < 1e-8
        energy1 = ham.compute_energy()
        energy1 -= ham.cache['energy_pert']

        assert energy1 > energy0
        if energy1_old is None:
            energy1_old = energy1
        else:
            assert abs(energy1 - energy1_old) < 1e-7


def test_ghost_hf():
    fn_fchk = context.get_fn('test/water_dimer_ghost.fchk')
    mol = IOData.from_file(fn_fchk)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)
    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham, mol.lf, olp, mol.exp_alpha) < 1e-5
