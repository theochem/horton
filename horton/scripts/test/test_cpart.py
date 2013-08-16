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


import os, h5py as h5

from horton import *
from horton.test.common import check_script, tmpdir
from horton.scripts.test.common import copy_files, check_files, write_random_lta_cube


def write_atomdb_refatoms(dn):
    padb = ProAtomDB.from_refatoms(numbers=[8,14], max_kation=3, max_anion=3)
    padb.to_file(os.path.join(dn, 'atoms.h5'))


def check_script_jbw_coarse(scheme):
    with tmpdir('horton.scripts.test.test_cpart.test_script_jbw_coarse_%s' % scheme) as dn:
        fn_cube = 'jbw_coarse_aedens.cube'
        copy_files(dn, [fn_cube])
        write_atomdb_refatoms(dn)
        check_script('horton-cpart.py %s %s atoms.h5' % (fn_cube, scheme), dn)
        fn_h5 = '%s_cpart.h5' % fn_cube[:-5]
        check_files(dn, [fn_h5])
        with h5.File(os.path.join(dn, fn_h5)) as f:
            assert 'cpart' in f
            assert scheme + '_r1' in f['cpart']


def test_script_jbw_coarse_h():
    # other schemes don't work because the cube file is too crappy.
    check_script_jbw_coarse('h')


def check_script_lta(fn_sym, suffix):
    with tmpdir('horton.scripts.test.test_cpart.test_script_lta_coarse_h_%s' % suffix) as dn:
        # prepare files
        if fn_sym is not None:
            copy_files(dn, [fn_sym])
        write_atomdb_refatoms(dn)

        # write a random cube file
        fn_cube = 'dens.cube'
        sys = write_random_lta_cube(dn, fn_cube)

        # run the script
        if fn_sym is None:
            check_script('horton-cpart.py %s h atoms.h5' % (fn_cube), dn)
        else:
            check_script('horton-cpart.py %s h atoms.h5 --symmetry=%s' % (fn_cube, fn_sym), dn)

        # check the output
        fn_h5 = '%s_cpart.h5' % fn_cube[:-5]
        check_files(dn, [fn_h5])
        with h5.File(os.path.join(dn, fn_h5)) as f:
            assert 'cpart' in f
            assert 'h_r1' in f['cpart']
            if fn_sym is not None:
                assert 'symmetry' in f['system/extra']
                assert 'symmetry' in f['cpart/h_r1']
                assert 'charges' in f['cpart/h_r1/symmetry']
                assert 'cartesian_multipoles' in f['cpart/h_r1/symmetry']
                for name, ds in f['cpart/h_r1/symmetry'].iteritems():
                    assert ds.shape[0] == sys.extra['symmetry'].natom
                    assert ds.shape[1] == 2


def test_script_lta():
    check_script_lta(None, 'nosym')


def test_script_lta_sym():
    check_script_lta('lta_gulp.cif', 'sym')
