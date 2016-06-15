# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2015 The HORTON Development Team
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
from nose.plugins.attrib import attr

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.test.common import check_script, tmpdir
from horton.scripts.test.common import copy_files, check_files, write_random_lta_cube
from horton.scripts.cpart import cpart_schemes


def test_cpart_schemes():
    assert 'h' in cpart_schemes
    assert 'hi' in cpart_schemes
    assert 'he' in cpart_schemes
    assert cpart_schemes['hi'] is HirshfeldICPart

    for CPartClass in cpart_schemes.itervalues():
        assert hasattr(CPartClass, 'options')


def write_atomdb_refatoms(dn):
    padb = ProAtomDB.from_refatoms(numbers=[8,14], max_kation=3, max_anion=3)
    padb.to_file(os.path.join(dn, 'atoms.h5'))


def check_script_jbw_coarse(scheme):
    with tmpdir('horton.scripts.test.test_cpart.test_script_jbw_coarse_%s' % scheme) as dn:
        fn_cube = 'jbw_coarse_aedens.cube'
        copy_files(dn, [fn_cube])
        write_atomdb_refatoms(dn)
        fn_h5 = 'foobar.h5'
        check_script('horton-cpart.py %s foobar.h5:cpart/%s_r1 %s atoms.h5' % (fn_cube, scheme, scheme), dn)
        check_files(dn, [fn_h5])
        with h5.File(os.path.join(dn, fn_h5)) as f:
            assert 'cpart' in f
            assert scheme + '_r1' in f['cpart']


@attr('slow')
def test_script_jbw_coarse_h():
    # other schemes don't work because the cube file is too crappy.
    check_script_jbw_coarse('h')


def check_script_lta(fn_sym, suffix, do_spin=False):
    with tmpdir('horton.scripts.test.test_cpart.test_script_lta_coarse_h_%s' % suffix) as dn:
        # prepare files
        if fn_sym is not None:
            copy_files(dn, [fn_sym])
        write_atomdb_refatoms(dn)

        # write a random cube file
        fn_cube = 'dens.cube'
        mol = write_random_lta_cube(dn, fn_cube)

        # if needed, write a random spin cube file
        if do_spin:
            fn_spin = 'spin.cube'
            molspin = write_random_lta_cube(dn, fn_spin)

        # run the script
        fn_h5 = '%s_cpart.h5' % fn_cube[:-5]
        opts = ''
        if not (fn_sym is None):
            opts += ' --symmetry=%s' % fn_sym
        if do_spin:
            opts += ' --spindens=%s' % fn_spin
        check_script('horton-cpart.py %s %s:cpart/h_r1 h atoms.h5 %s' % (fn_cube, fn_h5, opts), dn)

        # check the output
        check_files(dn, [fn_h5])
        with h5.File(os.path.join(dn, fn_h5)) as f:
            assert 'cpart' in f
            assert 'h_r1' in f['cpart']
            assert 'charges' in f['cpart/h_r1']
            if do_spin:
                assert 'spin_charges' in f['cpart/h_r1']
            if fn_sym is not None:
                assert 'symmetry' in f['cpart/h_r1']
                assert 'charges' in f['cpart/h_r1/symmetry']
                if do_spin:
                    assert 'spin_charges' in f['cpart/h_r1/symmetry']
                assert 'cartesian_multipoles' in f['cpart/h_r1/symmetry']
                for name, ds in f['cpart/h_r1/symmetry'].iteritems():
                    assert ds.shape[0] == mol.symmetry.natom
                    assert ds.shape[1] == 2


@attr('slow')
def test_script_lta():
    check_script_lta(None, 'nosym')


@attr('slow')
def test_script_lta_spin():
    check_script_lta(None, 'nosym_spin', True)


@attr('slow')
def test_script_lta_sym():
    check_script_lta('lta_gulp.cif', 'sym')


@attr('slow')
def test_script_lta_sym_spin():
    check_script_lta('lta_gulp.cif', 'sym_spin', True)
