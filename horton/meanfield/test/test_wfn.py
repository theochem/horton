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


from nose.tools import assert_raises
import h5py as h5
from horton import *


def test_dm_water_sto3g_hf():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    mol = Molecule.from_file(fn_fchk)
    dm = mol.wfn.dm_full
    assert abs(dm.get_element(0, 0) - 2.10503807) < 1e-7
    assert abs(dm.get_element(0, 1) - -0.439115917) < 1e-7
    assert abs(dm.get_element(1, 1) - 1.93312061) < 1e-7


def test_clear_exp():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    mol = Molecule.from_file(fn_fchk)
    assert 'exp_alpha' in mol.wfn._cache
    mol.wfn.dm_alpha
    mol.wfn.clear_exp()
    assert 'exp_alpha' not in mol.wfn._cache
    assert 'dm_alpha' in mol.wfn._cache


def test_clear_dm():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    mol = Molecule.from_file(fn_fchk)
    assert 'exp_alpha' in mol.wfn._cache
    mol.wfn.dm_alpha
    mol.wfn.clear_dm()
    assert 'exp_alpha' in mol.wfn._cache
    assert 'dm_alpha' not in mol.wfn._cache


def test_conversion_dm_exp():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    oes = sys.wfn.exp_alpha.energies.copy()
    dm = sys.lf.create_one_body()
    dm.assign(sys.wfn.dm_alpha)

    olp = sys.get_overlap()
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
        DirectTerm(er, sys.lf, sys.wfn),
        ExchangeTerm(er, sys.lf, sys.wfn),
        OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
    ]
    ham = Hamiltonian(terms)
    fock = sys.lf.create_one_body()
    ham.compute_fock(fock, None)
    sys.wfn.clear()
    sys.wfn.update_exp(fock, olp, dm)

    assert abs(sys.wfn.exp_alpha.occupations.sum() - 5) < 1e-8
    assert abs(oes - sys.wfn.exp_alpha.energies).max() < 3e-8
    assert (dm._array == sys.wfn.dm_alpha._array).all()
    # ugly hack
    sys.wfn.clear_dm()
    assert abs(dm._array - sys.wfn.dm_alpha._array).max() < 1e-9


def test_dm_lih_sto3g_hf():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = Molecule.from_file(fn_fchk)

    dm = mol.wfn.dm_full
    assert abs(dm.get_element(0, 0) - 1.96589709) < 1e-7
    assert abs(dm.get_element(0, 1) - 0.122114249) < 1e-7
    assert abs(dm.get_element(1, 1) - 0.0133112081) < 1e-7
    assert abs(dm.get_element(10, 10) - 4.23924688E-01) < 1e-7

    dm = mol.wfn.dm_spin
    assert abs(dm.get_element(0, 0) - 1.40210760E-03) < 1e-9
    assert abs(dm.get_element(0, 1) - -2.65370873E-03) < 1e-9
    assert abs(dm.get_element(1, 1) - 5.38701212E-03) < 1e-9
    assert abs(dm.get_element(10, 10) - 4.23889148E-01) < 1e-7


def check_spin(fn_fchk, sz0, ssq0, eps):
    path_fchk = context.get_fn('test/%s' % fn_fchk)
    mol = Molecule.from_file(path_fchk)
    olp = mol.lf.create_one_body()
    mol.obasis.compute_overlap(olp)
    sz, ssq = mol.wfn.get_spin(olp)
    assert abs(sz - sz0) < eps
    assert abs(ssq - ssq0) < eps
    if isinstance(mol.wfn, RestrictedWFN):
        sz, ssq = mol.wfn.get_spin()
        assert sz == 0.0
        assert ssq == 0.0
    else:
        # swap the spins and test again
        wfn = UnrestrictedWFN(mol.lf, mol.wfn.nbasis, mol.wfn.occ_model, mol.wfn.norb)
        exp_alpha = wfn.init_exp('alpha')
        exp_beta = wfn.init_exp('beta')
        exp_alpha.assign(mol.wfn.exp_beta)
        exp_beta.assign(mol.wfn.exp_alpha)
        sz, ssq = wfn.get_spin(olp)
        assert abs(sz + sz0) < eps
        assert abs(ssq - ssq) < eps

def test_spin_li_h():
    check_spin('li_h_3-21G_hf_g09.fchk', 0.5, 0.75, 1e-7)


def test_spin_h3_hfs():
    check_spin('h3_hfs_321g.fchk', 0.5, 0.7530, 1e-4)


def test_spin_h3_pbe():
    check_spin('h3_pbe_321g.fchk', 0.5, 0.7530, 1e-4)


def test_spin_ch3_hf():
    check_spin('ch3_hf_sto3g.fchk', 0.5, 0.7632, 1e-4)


def test_spin_water_hf():
    check_spin('water_sto3g_hf_g03.fchk', 0.0, 0.0, 1e-10)


def test_homo_lumo_water_hf():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    mol = Molecule.from_file(fn_fchk)
    assert abs(mol.wfn.homo_energy - -3.87671783E-01) < 1e-8
    assert abs(mol.wfn.lumo_energy - 6.03082408E-01) < 1e-8


def test_homo_lumo_ch3_hf():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    assert abs(mol.wfn.homo_energy - -3.63936540E-01) < 1e-8
    assert abs(mol.wfn.lumo_energy - 3.28562907E-01) < 1e-8


def test_homo_lumo_h():
    fn_fchk = context.get_fn('test/atom_001_001_hf_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    assert abs(mol.wfn.homo_energy - -4.66581850E-01) < 1e-8
    assert abs(mol.wfn.lumo_energy - 3.08024094E-01) < 1e-8


def test_homo_lumo_he():
    fn_fchk = context.get_fn('test/helium_hf_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    assert abs(mol.wfn.homo_energy - -8.76035508E-01) < 1e-8
    assert mol.wfn.lumo_energy is None # All orbitals are occupied


def test_setup_wfn_cs():
    numbers = np.array([6])
    nbasis = 10
    lf = DenseLinalgFactory

    wfn = setup_mean_field_wfn(nbasis, numbers, lf, 0, 1)
    assert isinstance(wfn, RestrictedWFN)
    assert wfn.nbasis == nbasis
    assert wfn.nel == 6
    assert wfn.occ_model.nalpha == 3
    assert wfn.occ_model.nbeta == 3

    wfn = setup_mean_field_wfn(nbasis, numbers, lf)
    assert isinstance(wfn, RestrictedWFN)
    assert wfn.nbasis == nbasis
    assert wfn.nel == 6
    assert wfn.occ_model.nalpha == 3
    assert wfn.occ_model.nbeta == 3

    wfn = setup_mean_field_wfn(nbasis, numbers, lf, 2)
    assert isinstance(wfn, RestrictedWFN)
    assert wfn.occ_model.nalpha == 2
    assert wfn.occ_model.nbeta == 2

    with assert_raises(ValueError):
        wfn = setup_mean_field_wfn(nbasis, numbers, lf, 0, 2)


def test_setup_wfn_os():
    nbasis = 10
    lf = DenseLinalgFactory

    wfn = setup_mean_field_wfn(nbasis, np.array([7]), lf, 0, 2)
    assert isinstance(wfn, UnrestrictedWFN)
    assert wfn.nbasis == nbasis
    assert wfn.nalpha == 4
    assert wfn.nbeta == 3
    assert wfn.occ_model.nalpha == 4
    assert wfn.occ_model.nbeta == 3

    wfn = setup_mean_field_wfn(nbasis, np.array([8]), lf, 1)
    assert isinstance(wfn, UnrestrictedWFN)
    assert wfn.nbasis == nbasis
    assert wfn.nalpha == 4
    assert wfn.nbeta == 3
    assert wfn.occ_model.nalpha == 4
    assert wfn.occ_model.nbeta == 3

    with assert_raises(ValueError):
        wfn = setup_mean_field_wfn(nbasis, np.array([7]), lf, 0, 1)


def test_setup_wfn_cs_fractional():
    nbasis = 10
    lf = DenseLinalgFactory()
    wfn = setup_mean_field_wfn(nbasis, np.array([6]), lf, 0.2, 1.0)
    assert isinstance(wfn, RestrictedWFN)
    assert abs(wfn.occ_model.nalpha - 2.9) < 1e-10
    assert abs(wfn.occ_model.nbeta - 2.9) < 1e-10
    wfn.init_exp('alpha')
    wfn.occ_model.assign(wfn.exp_alpha)
    assert abs(wfn.exp_alpha.occupations[:4] - [1.0, 1.0, 0.9, 0.0]).max() < 1e-10
    wfn = setup_mean_field_wfn(nbasis, np.array([6]), lf, 0.4, 1.0)
    assert isinstance(wfn, RestrictedWFN)
    assert abs(wfn.occ_model.nalpha - 2.8) < 1e-10
    assert abs(wfn.occ_model.nbeta - 2.8) < 1e-10
    wfn.init_exp('alpha')
    wfn.occ_model.assign(wfn.exp_alpha)
    assert abs(wfn.exp_alpha.occupations[:4] - [1.0, 1.0, 0.8, 0.0]).max() < 1e-10


def test_setup_wfn_os_fractional():
    nbasis = 10
    lf = DenseLinalgFactory()
    wfn = setup_mean_field_wfn(nbasis, np.array([7]), lf, 0.4, 2.5)
    assert isinstance(wfn, UnrestrictedWFN)
    assert abs(wfn.occ_model.nalpha - (7.0-0.4+2.5-1)/2) < 1e-10
    assert abs(wfn.occ_model.nbeta - (7.0-0.4-2.5+1)/2) < 1e-10
    wfn.init_exp('alpha')
    wfn.init_exp('beta')
    wfn.occ_model.assign(wfn.exp_alpha, wfn.exp_beta)
    assert abs(wfn.exp_alpha.occupations[:6] - [1.0, 1.0, 1.0, 1.0, 0.05, 0.0]).max() < 1e-10
    assert abs(wfn.exp_beta.occupations[:6] - [1.0, 1.0, 0.55, 0.0, 0.0, 0.0]).max() < 1e-10
    wfn = setup_mean_field_wfn(nbasis, np.array([7]), lf, 0.4, 2.0)
    assert isinstance(wfn, UnrestrictedWFN)
    assert abs(wfn.occ_model.nalpha - (7.0-0.4+2.0-1)/2) < 1e-10
    assert abs(wfn.occ_model.nbeta - (7.0-0.4-2.0+1)/2) < 1e-10
    wfn.init_exp('alpha')
    wfn.init_exp('beta')
    wfn.occ_model.assign(wfn.exp_alpha, wfn.exp_beta)
    assert abs(wfn.exp_alpha.occupations[:6] - [1.0, 1.0, 1.0, 0.8, 0.0, 0.0]).max() < 1e-10
    assert abs(wfn.exp_beta.occupations[:6] - [1.0, 1.0, 0.8, 0.0, 0.0, 0.0]).max() < 1e-10
    wfn = setup_mean_field_wfn(nbasis, np.array([7]), lf, 0.4, 'free')
    assert isinstance(wfn, UnrestrictedWFN)
    assert isinstance(wfn.occ_model, AufbauSpinOccModel)
    assert abs(wfn.occ_model.nel - (7.0-0.4)) < 1e-10


def test_fermi_occ_model_cs_helium():
    fn_fchk = context.get_fn('test/helium_hf_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    mol.wfn.occ_model = FermiOccModel(1.0)
    mol.wfn.occ_model.assign(mol.wfn.exp_alpha)
    assert (mol.wfn.exp_alpha.occupations == [1.0]).all()


def test_fermi_occ_model_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = Molecule.from_file(fn_fchk)
    for temperature in 300, 3000, 10000, 30000:
        mol.wfn.occ_model = FermiOccModel(5.0, temperature=temperature)
        mol.wfn.occ_model.assign(mol.wfn.exp_alpha)
        occ = mol.wfn.exp_alpha.occupations
        assert abs(occ.sum() - 5.0) < 1e-8
        assert (occ[1:] <= occ[:-1]).all()


def test_fermi_occ_model_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = Molecule.from_file(fn_fchk)
    for temperature in 300, 3000, 10000, 30000:
        mol.wfn.occ_model = FermiOccModel(2.0, 1.0, temperature=temperature)
        mol.wfn.occ_model.assign(mol.wfn.exp_alpha)
        occ_a = mol.wfn.exp_alpha.occupations
        assert abs(occ_a.sum() - 2) < 1e-8
        assert (occ_a[1:] <= occ_a[:-1]).all()
        occ_b = mol.wfn.exp_beta.occupations
        assert abs(occ_b.sum() - 1) < 1e-8
        assert (occ_b[1:] <= occ_b[:-1]).all()


def test_fermi_occ_model_hdf5():
    with h5.File('horton.meanfield.test.test_wfn.test_fermi_occ_model_hdf5', driver='core', backing_store=False) as f:
        om1 = FermiOccModel(2.0, 1.0, temperature=1025, eps=1e-10)
        om1.to_hdf5(f)
        om2 = FermiOccModel.from_hdf5(f, None)
        assert om1.nalpha == om2.nalpha
        assert om1.nbeta == om2.nbeta
        assert om1.temperature == om2.temperature
        assert om1.eps == om2.eps


def test_level_shift():
    fn_fchk = context.get_fn('test/helium_hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    overlap = sys.get_overlap()
    dm_alpha1 = sys.wfn.dm_alpha.copy()
    ls_alpha = sys.wfn.get_level_shift('alpha', overlap)
    sys.wfn.clear()
    sys.wfn.update_exp(ls_alpha, overlap)
    dm_alpha2 = sys.wfn.dm_alpha.copy()
    assert abs(dm_alpha1._array - dm_alpha2._array) < 1e-5
