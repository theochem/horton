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


def compare_expansions(e1, e2):
    assert (e1.coeffs == e2.coeffs).all()
    if e1.energies is None:
        assert e2.energies is None
    else:
        assert (e1.energies == e2.energies).all()


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


def compare_systems(sys1, sys2):
    assert (sys1.numbers == sys2.numbers).all()
    assert (sys1.coordinates == sys2.coordinates).all()
    # basis
    if sys1.basis is not None:
        assert (sys1.basis.centers == sys2.basis.centers).all()
        assert (sys1.basis.shell_map == sys2.basis.shell_map).all()
        assert (sys1.basis.nprims == sys2.basis.nprims).all()
        assert (sys1.basis.shell_types == sys2.basis.shell_types).all()
        assert (sys1.basis.alphas == sys2.basis.alphas).all()
        assert (sys1.basis.con_coeffs == sys2.basis.con_coeffs).all()
    else:
        assert sys2.basis is None
    # wfn
    if isinstance(sys1.wfn, ClosedShellWFN):
        assert isinstance(sys2.wfn, ClosedShellWFN)
        assert sys1.wfn.nep == sys2.wfn.nep
        assert sys1.wfn.nbasis == sys2.wfn.nbasis
        assert sys1.wfn.norb == sys2.wfn.norb
        assert (sys1.wfn.expansion.coeffs == sys2.wfn.expansion.coeffs).all()
        compare_expansions(sys1.wfn.expansion, sys2.wfn.expansion)
    elif isinstance(sys1.wfn, OpenShellWFN):
        assert isinstance(sys2.wfn, OpenShellWFN)
        assert sys1.wfn.nalpha == sys2.wfn.nalpha
        assert sys1.wfn.nbeta == sys2.wfn.nbeta
        assert sys1.wfn.nbasis == sys2.wfn.nbasis
        assert sys1.wfn.norb == sys2.wfn.norb
        compare_expansions(sys1.wfn.alpha_expansion, sys2.wfn.alpha_expansion)
        compare_expansions(sys1.wfn.beta_expansion, sys2.wfn.beta_expansion)
    elif sys1.wfn is None:
        assert sys2.wfn is None
    else:
        raise NotImplementedError
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
        assert sys2.wfn is not None
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
        assert sys2.wfn is not None
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
