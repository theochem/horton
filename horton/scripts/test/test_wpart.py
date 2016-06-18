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


import os, h5py as h5

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.test.common import check_script, tmpdir
from horton.part.test.common import get_proatomdb_hf_sto3g
from horton.scripts.test.common import copy_files, check_files
from horton.scripts.wpart import wpart_schemes


def test_wpart_schemes():
    assert 'b' in wpart_schemes
    assert 'h' in wpart_schemes
    assert 'hi' in wpart_schemes
    assert wpart_schemes['hi'] is HirshfeldIWPart
    assert wpart_schemes['hi'].options == ['lmax', 'threshold', 'maxiter']
    assert not wpart_schemes['hi'].linear
    assert wpart_schemes['h'].linear
    assert wpart_schemes['b'].linear

    for WPartClass in wpart_schemes.itervalues():
        assert hasattr(WPartClass, 'options')


def write_atomdb_sto3g(dn, do_deriv=True):
    padb = get_proatomdb_hf_sto3g()
    if not do_deriv:
        # remove the derivatives of the pro-atoms (ugly hack to test backward compatibility).
        for record in padb._records:
            record._deriv = None
    padb.to_file(os.path.join(dn, 'atoms.h5'))


def check_script_water_sto3g(scheme, do_deriv=True):
    with tmpdir('horton.scripts.test.test_wpart.test_script_water_sto3g_%s' % scheme) as dn:
        fn_fchk = 'water_sto3g_hf_g03.fchk'
        copy_files(dn, [fn_fchk])
        if scheme == 'b':
            check_script('horton-wpart.py %s water_sto3g_hf_g03_wpart.h5:wpart %s --debug' % (fn_fchk, scheme), dn)
        else:
            write_atomdb_sto3g(dn, do_deriv)
            check_script('horton-wpart.py %s water_sto3g_hf_g03_wpart.h5:wpart %s atoms.h5 --slow' % (fn_fchk, scheme), dn)
        fn_h5 = 'water_sto3g_hf_g03_wpart.h5'
        check_files(dn, [fn_h5])
        with h5.File(os.path.join(dn, fn_h5)) as f:
            assert 'wpart' in f
            assert abs(f['wpart/charges'][:].sum()) < 1e-2
            assert 'wpart/spin_charges' not in f
            assert 'atom_00000' in f['wpart']
            for s in 'density_decomposition', 'hartree_decomposition':
                assert s in f['wpart']['atom_00000']
                assert 'spline_00000' in f['wpart']['atom_00000'][s]
                assert 'rtransform' in f['wpart']['atom_00000'][s]['spline_00000'].attrs
                assert 'extrapolation' in f['wpart']['atom_00000'][s]['spline_00000'].attrs
                assert 'y' in f['wpart']['atom_00000'][s]['spline_00000']
                assert 'd' in f['wpart']['atom_00000'][s]['spline_00000']
            if scheme != 'b':
                assert 'spline_prodensity' in f['wpart']['atom_00000']
                assert 'spline_prohartree' in f['wpart']['atom_00000']


def test_script_water_sto3g_b():
    check_script_water_sto3g('b')


def test_script_water_sto3g_h():
    check_script_water_sto3g('h')


def test_script_water_sto3g_h_noderiv():
    check_script_water_sto3g('h', do_deriv=False)


def test_script_water_sto3g_hi():
    check_script_water_sto3g('hi')


def test_script_water_sto3g_hi_noderiv():
    check_script_water_sto3g('hi', do_deriv=False)


def check_script_ch3_rohf_sto3g(scheme, do_deriv=True):
    with tmpdir('horton.scripts.test.test_wpart.test_script_ch3_rohf_sto3g_%s' % scheme) as dn:
        fn_fchk = 'ch3_rohf_sto3g_g03.fchk'
        copy_files(dn, [fn_fchk])
        if scheme == 'b':
            check_script('horton-wpart.py %s foo.h5:wpart %s --debug' % (fn_fchk, scheme), dn)
        else:
            write_atomdb_sto3g(dn, do_deriv)
            check_script('horton-wpart.py %s foo.h5:wpart %s atoms.h5' % (fn_fchk, scheme), dn)
        fn_h5 = 'foo.h5'
        check_files(dn, [fn_h5])
        with h5.File(os.path.join(dn, fn_h5)) as f:
            assert 'wpart' in f
            assert abs(f['wpart/charges'][:].sum()) < 1e-2
            assert abs(f['wpart/spin_charges'][:].sum() - 1) < 1e-2


def test_script_ch3_rohf_sto3g_b():
    check_script_ch3_rohf_sto3g('b')


def test_script_ch3_rohf_sto3g_h():
    check_script_ch3_rohf_sto3g('h')


def test_script_ch3_rohf_sto3g_h_noderiv():
    check_script_ch3_rohf_sto3g('h', do_deriv=False)


def test_script_ch3_rohf_sto3g_hi():
    check_script_ch3_rohf_sto3g('hi')


def test_script_ch3_rohf_sto3g_hi_noderiv():
    check_script_ch3_rohf_sto3g('hi', do_deriv=False)
