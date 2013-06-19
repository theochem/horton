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


from horton import *


def test_dm_water_sto3g_hf():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    dm = sys.wfn.dm_full
    assert abs(dm.get_element(0, 0) - 2.10503807) < 1e-7
    assert abs(dm.get_element(0, 1) - -0.439115917) < 1e-7
    assert abs(dm.get_element(1, 1) - 1.93312061) < 1e-7


def test_conversion_dm_exp():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    oes = sys.wfn.exp_alpha.energies.copy()
    dm = sys.lf.create_one_body()
    dm.assign(sys.wfn.dm_alpha)

    ham = Hamiltonian(sys, [HartreeFockExchange()])
    fock = sys.lf.create_one_body()
    ham.compute_fock(fock, None)
    sys.wfn.invalidate()
    sys.wfn.update_exp(fock, sys.get_overlap(), dm)

    assert abs(sys.wfn.exp_alpha.occupations.sum() - 5) < 1e-8
    assert abs(oes - sys.wfn.exp_alpha.energies).max() < 3e-8
    assert (dm._array == sys.wfn.dm_alpha._array).all()
    # ugly hack
    sys.wfn._cache.clear_item('dm_alpha')
    assert abs(dm._array - sys.wfn.dm_alpha._array).max() < 1e-9


def test_dm_lih_sto3g_hf():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)

    dm = sys.wfn.dm_full
    assert abs(dm.get_element(0, 0) - 1.96589709) < 1e-7
    assert abs(dm.get_element(0, 1) - 0.122114249) < 1e-7
    assert abs(dm.get_element(1, 1) - 0.0133112081) < 1e-7
    assert abs(dm.get_element(10, 10) - 4.23924688E-01) < 1e-7

    dm = sys.wfn.dm_spin
    assert abs(dm.get_element(0, 0) - 1.40210760E-03) < 1e-9
    assert abs(dm.get_element(0, 1) - -2.65370873E-03) < 1e-9
    assert abs(dm.get_element(1, 1) - 5.38701212E-03) < 1e-9
    assert abs(dm.get_element(10, 10) - 4.23889148E-01) < 1e-7


def test_spin_li_h():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    sz, ssq = sys.wfn.get_spin(sys.get_overlap())
    assert sz == 0.5
    assert abs(ssq - 0.75) < 1e-7
    # swap the spins and test again
    wfn = UnrestrictedWFN(sys.wfn.occ_model, sys.lf, sys.wfn.nbasis, sys.wfn.norb)
    exp_alpha = wfn.init_exp('alpha')
    exp_beta = wfn.init_exp('beta')
    exp_alpha.assign(sys.wfn.exp_beta)
    exp_beta.assign(sys.wfn.exp_alpha)
    sz, ssq = wfn.get_spin(sys.get_overlap())
    assert sz == -0.5
    assert abs(ssq - 0.75) < 1e-7


def test_spin_h3_hfs():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)
    sz, ssq = sys.wfn.get_spin(sys.get_overlap())
    assert sz == 0.5
    assert abs(ssq - 0.7530) < 1e-4


def test_spin_h3_pbe():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    sys = System.from_file(fn_fchk)
    sz, ssq = sys.wfn.get_spin(sys.get_overlap())
    assert sz == 0.5
    assert abs(ssq - 0.7530) < 1e-4


def test_spin_ch3_hf():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    sz, ssq = sys.wfn.get_spin(sys.get_overlap())
    assert sz == 0.5
    assert abs(ssq - 0.7632) < 1e-4


def test_homo_lumo_water_hf():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    assert abs(sys.wfn.homo_energy - -3.87671783E-01) < 1e-8
    assert abs(sys.wfn.lumo_energy - 6.03082408E-01) < 1e-8


def test_homo_lumo_ch3_hf():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    assert abs(sys.wfn.homo_energy - -3.63936540E-01) < 1e-8
    assert abs(sys.wfn.lumo_energy - 3.28562907E-01) < 1e-8


def test_homo_lumo_h():
    fn_fchk = context.get_fn('test/atom_001_001_hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    assert abs(sys.wfn.homo_energy - -4.66581850E-01) < 1e-8
    assert abs(sys.wfn.lumo_energy - 3.08024094E-01) < 1e-8


def test_homo_lumo_he():
    fn_fchk = context.get_fn('test/helium_hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    assert abs(sys.wfn.homo_energy - -8.76035508E-01) < 1e-8
    assert sys.wfn.lumo_energy is None # All orbitals are occupied


def test_setup_wfn_cs():
    sys = System(np.zeros((1,3), float), np.array([6]), obasis='3-21g')
    setup_mean_field_wfn(sys, 0, 1)
    assert isinstance(sys.wfn, RestrictedWFN)
    assert sys.wfn.occ_model.nalpha == 3
    assert sys.wfn.occ_model.nbeta == 3

    sys = System(np.zeros((1,3), float), np.array([6]), obasis='3-21g')
    setup_mean_field_wfn(sys)
    assert isinstance(sys.wfn, RestrictedWFN)
    assert sys.wfn.occ_model.nalpha == 3
    assert sys.wfn.occ_model.nbeta == 3

    sys = System(np.zeros((1,3), float), np.array([6]), obasis='3-21g')
    setup_mean_field_wfn(sys, 2)
    assert isinstance(sys.wfn, RestrictedWFN)
    assert sys.wfn.occ_model.nalpha == 2
    assert sys.wfn.occ_model.nbeta == 2

    try:
        sys = System(np.zeros((1,3), float), np.array([6]), obasis='3-21g')
        setup_mean_field_wfn(sys, 0, 2)
        assert False
    except ValueError:
        pass


def test_setup_wfn_os():
    sys = System(np.zeros((1,3), float), np.array([7]), obasis='3-21g')
    setup_mean_field_wfn(sys, 0, 2)
    assert isinstance(sys.wfn, UnrestrictedWFN)
    assert sys.wfn.occ_model.nalpha == 4
    assert sys.wfn.occ_model.nbeta == 3

    sys = System(np.zeros((1,3), float), np.array([8]), obasis='3-21g')
    setup_mean_field_wfn(sys, 1)
    assert isinstance(sys.wfn, UnrestrictedWFN)
    assert sys.wfn.occ_model.nalpha == 4
    assert sys.wfn.occ_model.nbeta == 3
    setup_mean_field_wfn(sys, 1)

    try:
        sys = System(np.zeros((1,3), float), np.array([7]), obasis='3-21g')
        setup_mean_field_wfn(sys, 0, 1)
        assert False
    except ValueError:
        pass


def test_setup_wfn_cs_fractional():
    sys = System(np.zeros((1,3), float), np.array([6]), obasis='3-21g')
    setup_mean_field_wfn(sys, 0.2, 1.0)
    assert isinstance(sys.wfn, RestrictedWFN)
    assert abs(sys.wfn.occ_model.nalpha - 2.9) < 1e-10
    assert abs(sys.wfn.occ_model.nbeta - 2.9) < 1e-10
    sys.wfn.init_exp('alpha')
    sys.wfn.occ_model.assign(sys.wfn.exp_alpha)
    assert abs(sys.wfn.exp_alpha.occupations[:4] - [1.0, 1.0, 0.9, 0.0]).max() < 1e-10
    setup_mean_field_wfn(sys, 0.4, 1.0)
    assert isinstance(sys.wfn, RestrictedWFN)
    assert abs(sys.wfn.occ_model.nalpha - 2.8) < 1e-10
    assert abs(sys.wfn.occ_model.nbeta - 2.8) < 1e-10
    assert abs(sys.wfn.exp_alpha.occupations[:4] - [1.0, 1.0, 0.8, 0.0]).max() < 1e-10


def test_setup_wfn_os_fractional():
    sys = System(np.zeros((1,3), float), np.array([7]), obasis='3-21g')
    setup_mean_field_wfn(sys, 0.4, 2.5)
    assert isinstance(sys.wfn, UnrestrictedWFN)
    assert abs(sys.wfn.occ_model.nalpha - (7.0-0.4+2.5-1)/2) < 1e-10
    assert abs(sys.wfn.occ_model.nbeta - (7.0-0.4-2.5+1)/2) < 1e-10
    sys.wfn.init_exp('alpha')
    sys.wfn.init_exp('beta')
    sys.wfn.occ_model.assign(sys.wfn.exp_alpha, sys.wfn.exp_beta)
    assert abs(sys.wfn.exp_alpha.occupations[:6] - [1.0, 1.0, 1.0, 1.0, 0.05, 0.0]).max() < 1e-10
    assert abs(sys.wfn.exp_beta.occupations[:6] - [1.0, 1.0, 0.55, 0.0, 0.0, 0.0]).max() < 1e-10
    setup_mean_field_wfn(sys, 0.4, 2.0)
    assert isinstance(sys.wfn, UnrestrictedWFN)
    assert abs(sys.wfn.occ_model.nalpha - (7.0-0.4+2.0-1)/2) < 1e-10
    assert abs(sys.wfn.occ_model.nbeta - (7.0-0.4-2.0+1)/2) < 1e-10
    assert abs(sys.wfn.exp_alpha.occupations[:6] - [1.0, 1.0, 1.0, 0.8, 0.0, 0.0]).max() < 1e-10
    assert abs(sys.wfn.exp_beta.occupations[:6] - [1.0, 1.0, 0.8, 0.0, 0.0, 0.0]).max() < 1e-10
    setup_mean_field_wfn(sys, 0.4, 'free')
    assert isinstance(sys.wfn, UnrestrictedWFN)
    assert isinstance(sys.wfn.occ_model, AufbauSpinOccModel)
    assert abs(sys.wfn.occ_model.nel - (7.0-0.4)) < 1e-10
