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


import tempfile, os, h5py as h5

from horton import *


def compare_expansions(wfn1, wfn2, spin):
    if wfn1._cache.has('exp_%s' % spin):
        assert wfn2._cache.has('exp_%s' % spin)
        e1 = wfn1.get_exp(spin)
        e2 = wfn2.get_exp(spin)
        assert e1.nbasis == e2.nbasis
        assert e1.nfn == e2.nfn
        assert (e1.coeffs == e2.coeffs).all()
        assert (e1.energies == e2.energies).all()
        assert (e1.occupations == e2.occupations).all()
    else:
        assert not wfn2._cache.has('exp_%s' % spin)


def compare_all_expansions(wfn1, wfn2):
    compare_expansions(wfn1, wfn2, 'alpha')
    compare_expansions(wfn1, wfn2, 'beta')


def compare_dms(wfn1, wfn2, select):
    if wfn1._cache.has('dm_%s' % select):
        assert wfn2._cache.has('dm_%s' % select)
        dm1 = wfn1.get_dm(select)
        dm2 = wfn2.get_dm(select)
        assert dm1.nbasis == dm2.nbasis
        assert (dm1._array == dm2._array).all()
    else:
        assert not wfn2._cache.has('dm_%s' % select)


def compare_all_dms(wfn1, wfn2):
    compare_dms(wfn1, wfn2, 'alpha')
    compare_dms(wfn1, wfn2, 'beta')
    compare_dms(wfn1, wfn2, 'full')
    compare_dms(wfn1, wfn2, 'spin')


def compare_one_body(sys1, sys2, key):
    if key in sys1.operators:
        assert key in sys2.operators
        op1 = sys1.operators[key]
        op2 = sys2.operators[key]
        if isinstance(op1, DenseOneBody):
            assert isinstance(op2, DenseOneBody)
            assert op1.nbasis == op2.nbasis
            assert (op1._array == op2._array).all()
        else:
            raise NotImplementedError
    else:
        assert key not in sys2.operators


def compare_two_body(sys1, sys2, key):
    if key in sys1.operators:
        assert key in sys2.operators
        op1 = sys1.operators[key]
        op2 = sys2.operators[key]
        if isinstance(op1, DenseTwoBody):
            assert isinstance(op2, DenseTwoBody)
            assert op1.nbasis == op2.nbasis
            assert (op1._array == op2._array).all()
        else:
            raise NotImplementedError
    else:
        assert key not in sys2.operators


def compare_occ_model(occ_model1, occ_model2):
    assert occ_model1.__class__ == occ_model2.__class__
    if isinstance(occ_model1, AufbauOccModel):
        assert occ_model1.nalpha == occ_model2.nalpha
        assert occ_model1.nbeta == occ_model2.nbeta
    else:
        raise NotImplementedError


def compare_wfns(wfn1, wfn2):
    if isinstance(wfn1, ClosedShellWFN):
        assert isinstance(wfn2, ClosedShellWFN)
        assert wfn1.nbasis == wfn2.nbasis
        assert wfn1.norb == wfn2.norb
        compare_all_expansions(wfn1, wfn2)
        compare_all_dms(wfn1, wfn2)
        compare_occ_model(wfn1.occ_model, wfn2.occ_model)
    elif isinstance(wfn1, OpenShellWFN):
        assert isinstance(wfn2, OpenShellWFN)
        assert wfn1.nbasis == wfn2.nbasis
        assert wfn1.norb == wfn2.norb
        compare_all_expansions(wfn1, wfn2)
        compare_all_dms(wfn1, wfn2)
        compare_occ_model(wfn1.occ_model, wfn2.occ_model)
    elif wfn1 is None:
        assert wfn2 is None
    else:
        raise NotImplementedError


def compare_systems(sys1, sys2):
    assert (sys1.numbers == sys2.numbers).all()
    assert (sys1.coordinates == sys2.coordinates).all()
    # orbital basis
    if sys1.obasis is not None:
        assert (sys1.obasis.centers == sys2.obasis.centers).all()
        assert (sys1.obasis.shell_map == sys2.obasis.shell_map).all()
        assert (sys1.obasis.nprims == sys2.obasis.nprims).all()
        assert (sys1.obasis.shell_types == sys2.obasis.shell_types).all()
        assert (sys1.obasis.alphas == sys2.obasis.alphas).all()
        assert (sys1.obasis.con_coeffs == sys2.obasis.con_coeffs).all()
    else:
        assert sys2.obasis is None
    # wfn
    compare_wfns(sys1.wfn, sys2.wfn)
    # one-body operators
    compare_one_body(sys1, sys2, 'olp')
    compare_one_body(sys1, sys2, 'kin')
    compare_one_body(sys1, sys2, 'na')
    # two-body operators
    compare_two_body(sys1, sys2, 'er')


def test_chk_initialization_filename_cs():
    tmpdir = tempfile.mkdtemp('horton.test.test_checkpoint.test_chk_initialization_filename_cs')
    try:
        fn_chk = '%s/chk.h5' % tmpdir
        fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
        fn_log = context.get_fn('test/water_sto3g_hf_g03.log')
        sys1 = System.from_file(fn_fchk, fn_log, chk=fn_chk)
        del sys1
        sys1 = System.from_file(fn_fchk, fn_log)
        sys2 = System.from_file(fn_chk)
        compare_systems(sys1, sys2)
    finally:
        if os.path.isfile(fn_chk):
            os.remove(fn_chk)
        os.rmdir(tmpdir)


def test_chk_initialization_filename_os():
    tmpdir = tempfile.mkdtemp('horton.test.test_checkpoint.test_chk_initialization_filename_os')
    try:
        fn_chk = '%s/chk.h5' % tmpdir
        sys1 = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'), chk=fn_chk)
        del sys1
        sys1 = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
        sys2 = System.from_file(fn_chk)
        compare_systems(sys1, sys2)
    finally:
        if os.path.isfile(fn_chk):
            os.remove(fn_chk)
        os.rmdir(tmpdir)


def test_chk_initialization_file():
    chk = h5.File('horton.test.test_checkpoint.test_chk_initialization_file', driver='core', backing_store=False)
    sys1 = System.from_file(context.get_fn('test/hf_sto3g.fchk'), chk=chk)
    del sys1
    sys1 = System.from_file(context.get_fn('test/hf_sto3g.fchk'))
    sys2 = System.from_file(chk)
    compare_systems(sys1, sys2)
    chk.close()


def test_chk_initialization_override():
    tmpdir = tempfile.mkdtemp('horton.test.test_checkpoint.test_chk_override')
    try:
        fn_chk1 = '%s/chk1.h5' % tmpdir
        fn_chk2 = '%s/chk2.h5' % tmpdir
        sys1 = System.from_file(context.get_fn('test/hf_sto3g.fchk'), chk=fn_chk1)
        del sys1
        sys1 = System.from_file(context.get_fn('test/hf_sto3g.fchk'))
        sys2 = System.from_file(fn_chk1, chk=fn_chk2)

        compare_systems(sys1, sys2)
        assert os.path.isfile(fn_chk2)

        sys3 = System.from_file(fn_chk2, chk=None)
        compare_systems(sys1, sys3)
        sys3.numbers[:] = 0
        sys3.update_chk('numbers')

        sys4 = System.from_file(fn_chk2, chk=None)
        compare_systems(sys1, sys4)
    finally:
        if os.path.isfile(fn_chk1):
            os.remove(fn_chk1)
        if os.path.isfile(fn_chk2):
            os.remove(fn_chk2)
        os.rmdir(tmpdir)


def test_chk_update1():
    chk = h5.File('horton.test.test_checkpoint.test_chk_update1', driver='core', backing_store=False)
    sys1 = System.from_file(context.get_fn('test/hf_sto3g.fchk'), chk=chk)
    sys1.numbers[:] = [3, 2]
    sys1.coordinates[0,2] = 0.25
    sys1.update_chk()
    del sys1
    sys1 = System.from_file(chk)
    assert (sys1.numbers == [3, 2]).all()
    assert sys1.coordinates[0,2] == 0.25
    chk.close()


def test_chk_update2():
    chk = h5.File('horton.test.test_checkpoint.test_chk_update2', driver='core', backing_store=False)
    sys1 = System.from_file(context.get_fn('test/hf_sto3g.fchk'), chk=chk)
    sys1.numbers[:] = [3, 2]
    sys1.coordinates[0,2] = 0.25
    sys1.update_chk('coordinates')
    del sys1
    sys1 = System.from_file(chk)
    assert (sys1.numbers != [3, 2]).all()
    assert sys1.coordinates[0,2] == 0.25
    chk.close()


def test_chk_update3():
    chk = h5.File('horton.test.test_checkpoint.test_chk_update2', driver='core', backing_store=False)
    sys1 = System.from_file(context.get_fn('test/hf_sto3g.fchk'), chk=chk)
    sys1.numbers[:] = [3, 2]
    sys1.update_chk()
    sys1.coordinates[0,2] = 0.25
    del sys1
    sys1 = System.from_file(chk)
    assert (sys1.numbers == [3, 2]).all()
    assert sys1.coordinates[0,2] != 0.25
    chk.close()


def test_chk_operators():
    chk = h5.File('horton.test.test_checkpoint.test_chk_operators', driver='core', backing_store=False)
    sys1 = System.from_file(context.get_fn('test/hf_sto3g.fchk'), chk=chk)
    ar_olp = sys1.get_overlap()._array
    ar_kin = sys1.get_kinetic()._array
    ar_na = sys1.get_nuclear_attraction()._array
    ar_er = sys1.get_electron_repulsion()._array
    # manually put er integrals in checkpoint, not done by default
    sys1.update_chk('operators.er')
    del sys1
    sys1 = System.from_file(chk)
    assert 'olp' in sys1.operators
    assert 'kin' in sys1.operators
    assert 'na' in sys1.operators
    assert 'er' in sys1.operators
    assert (ar_olp == sys1.get_overlap()._array).all()
    assert (ar_kin == sys1.get_kinetic()._array).all()
    assert (ar_na == sys1.get_nuclear_attraction()._array).all()
    assert (ar_er == sys1.get_electron_repulsion()._array).all()


def test_chk_guess_scf_cs():
    chk = h5.File('horton.test.test_checkpoint.test_chk_guess_scf_cs', driver='core', backing_store=False)
    fn_fchk = context.get_fn('test/hf_sto3g.fchk')
    sys = System.from_file(fn_fchk, chk=chk)

    guess_hamiltonian_core(sys)
    c = sys.wfn.exp_alpha._coeffs
    e = sys.wfn.exp_alpha._energies
    dma = sys.wfn.dm_alpha._array
    del sys
    sys = System.from_file(chk)
    assert (sys.wfn.exp_alpha._coeffs == c).all()
    assert (sys.wfn.exp_alpha._energies == e).all()
    assert (sys.wfn.dm_alpha._array == dma).all()

    ham = Hamiltonian(sys, [HartreeFock()])
    converge_scf(ham, 5)
    c = sys.wfn.exp_alpha._coeffs
    e = sys.wfn.exp_alpha._energies
    dma = sys.wfn.dm_alpha._array
    ham.compute_energy()
    energy = sys.props['energy']
    energy_kin = sys.props['energy_kin']
    energy_hartree = sys.props['energy_hartree']
    energy_exchange_fock = sys.props['energy_exchange_fock']
    energy_ne = sys.props['energy_ne']
    energy_nn = sys.props['energy_nn']
    del sys
    del ham
    sys = System.from_file(chk)
    assert (sys.wfn.exp_alpha._coeffs == c).all()
    assert (sys.wfn.exp_alpha._energies == e).all()
    assert (sys.wfn.dm_alpha._array == dma).all()
    assert sys.props['energy'] == energy
    assert sys.props['energy_kin'] == energy_kin
    assert sys.props['energy_hartree'] == energy_hartree
    assert sys.props['energy_exchange_fock'] == energy_exchange_fock
    assert sys.props['energy_ne'] == energy_ne
    assert sys.props['energy_nn'] == energy_nn


def test_chk_guess_scf_os():
    chk = h5.File('horton.test.test_checkpoint.test_chk_guess_scf_os', driver='core', backing_store=False)
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk, chk=chk)

    guess_hamiltonian_core(sys)
    ac = sys.wfn.exp_alpha._coeffs
    bc = sys.wfn.exp_beta._coeffs
    ae = sys.wfn.exp_alpha._energies
    be = sys.wfn.exp_beta._energies
    dma = sys.wfn.dm_alpha._array
    dmb = sys.wfn.dm_beta._array
    del sys
    sys = System.from_file(chk)
    assert (sys.wfn.exp_alpha._coeffs == ac).all()
    assert (sys.wfn.exp_beta._coeffs == bc).all()
    assert (sys.wfn.exp_alpha._energies == ae).all()
    assert (sys.wfn.exp_beta._energies == be).all()
    assert (sys.wfn.dm_alpha._array == dma).all()
    assert (sys.wfn.dm_beta._array == dmb).all()

    ham = Hamiltonian(sys, [HartreeFock()])
    converge_scf(ham, 5)
    ac = sys.wfn.exp_alpha._coeffs
    bc = sys.wfn.exp_beta._coeffs
    ae = sys.wfn.exp_alpha._energies
    be = sys.wfn.exp_beta._energies
    dma = sys.wfn.dm_alpha._array
    dmb = sys.wfn.dm_beta._array
    ham.compute_energy()
    energy = sys.props['energy']
    energy_kin = sys.props['energy_kin']
    energy_hartree = sys.props['energy_hartree']
    energy_exchange_fock = sys.props['energy_exchange_fock']
    energy_ne = sys.props['energy_ne']
    energy_nn = sys.props['energy_nn']
    del sys
    del ham
    sys = System.from_file(chk)
    assert (sys.wfn.exp_alpha._coeffs == ac).all()
    assert (sys.wfn.exp_beta._coeffs == bc).all()
    assert (sys.wfn.exp_alpha._energies == ae).all()
    assert (sys.wfn.exp_beta._energies == be).all()
    assert (sys.wfn.dm_alpha._array == dma).all()
    assert (sys.wfn.dm_beta._array == dmb).all()
    assert sys.props['energy'] == energy
    assert sys.props['energy_kin'] == energy_kin
    assert sys.props['energy_hartree'] == energy_hartree
    assert sys.props['energy_exchange_fock'] == energy_exchange_fock
    assert sys.props['energy_ne'] == energy_ne
    assert sys.props['energy_nn'] == energy_nn
