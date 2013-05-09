# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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


import tempfile, shutil, os, h5py as h5

from horton import *
from horton.test.common import check_script
from horton.scripts.test.common import copy_files, check_files, write_random_lta_cube


def write_atomdb_refatoms(tmpdir):
    rtf = ExpRTransform(1e-3, 1e1, 30)
    rgrid = RadialGrid(rtf)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rgrid, 26))
    padb = ProAtomDB.from_refatoms(atgrid, numbers=[8,14], max_kation=3, max_anion=3)
    padb.to_file(os.path.join(tmpdir, 'atoms.h5'))


def check_script_jbw_coarse(scheme):
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_cpart.test_script_jbw_coarse_%s' % scheme)
    try:
        fn_cube = 'jbw_coarse_aedens.cube'
        copy_files(tmpdir, [fn_cube])
        write_atomdb_refatoms(tmpdir)
        check_script('horton-cpart.py %s %s atoms.h5' % (fn_cube, scheme), tmpdir)
        fn_h5 = '%s.h5' % fn_cube
        check_files(tmpdir, [fn_h5])
        with h5.File(os.path.join(tmpdir, fn_h5)) as f:
            assert 'cpart' in f
            assert scheme + '_r1' in f['cpart']
    finally:
        shutil.rmtree(tmpdir)


def test_script_jbw_coarse_h():
    # other schemes don't work because the cube file is too crappy.
    check_script_jbw_coarse('h')


def check_script_lta(fn_sym, suffix):
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_cpart.test_script_lta_coarse_h_%s' % suffix)
    try:
        # prepare files
        if fn_sym is not None:
            copy_files(tmpdir, [fn_sym])
        write_atomdb_refatoms(tmpdir)

        # write a random cube file
        fn_cube = 'dens.cube'
        sys = write_random_lta_cube(tmpdir, fn_cube)

        # run the script
        if fn_sym is None:
            check_script('horton-cpart.py %s h atoms.h5' % (fn_cube), tmpdir)
        else:
            check_script('horton-cpart.py %s h atoms.h5 --symmetry=%s' % (fn_cube, fn_sym), tmpdir)

        # check the output
        fn_h5 = '%s.h5' % fn_cube
        check_files(tmpdir, [fn_h5])
        with h5.File(os.path.join(tmpdir, fn_h5)) as f:
            assert 'cpart' in f
            assert 'h_r1' in f['cpart']
            if fn_sym is not None:
                assert 'symmetry' in f['system/props']
                assert 'symmetry' in f['cpart/h_r1']
                assert 'charges' in f['cpart/h_r1/symmetry']
                assert 'cartesian_moments' in f['cpart/h_r1/symmetry']
                for name, ds in f['cpart/h_r1/symmetry'].iteritems():
                    assert ds.shape[0] == sys.props['symmetry'].natom
                    assert ds.shape[1] == 2
    finally:
        shutil.rmtree(tmpdir)


def test_script_lta():
    check_script_lta(None, 'nosym')


def test_script_lta_sym():
    check_script_lta('lta_gulp.cif', 'sym')
