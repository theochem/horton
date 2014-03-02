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


import os, h5py as h5

from horton.test.common import check_script, tmpdir
from horton.scripts.test.common import copy_files, check_files
from horton.part.test.common import get_proatomdb_hf_sto3g


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
            check_script('horton-wpart.py %s %s' % (fn_fchk, scheme), dn)
        else:
            write_atomdb_sto3g(dn, do_deriv)
            check_script('horton-wpart.py %s %s atoms.h5' % (fn_fchk, scheme), dn)
        fn_h5 = 'water_sto3g_hf_g03_wpart.h5'
        check_files(dn, [fn_h5])
        with h5.File(os.path.join(dn, fn_h5)) as f:
            assert 'wpart' in f
            assert scheme in f['wpart']


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


def test_script_water_sto3g_he():
    check_script_water_sto3g('he')


def test_script_water_sto3g_he_noderiv():
    check_script_water_sto3g('he', do_deriv=False)
