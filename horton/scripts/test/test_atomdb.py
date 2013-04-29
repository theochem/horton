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


import tempfile, shutil, os
from horton.context import context
from horton.periodic import periodic
from horton.part.proatomdb import ProAtomDB
from horton.test.common import check_script
from horton.scripts.test.common import copy_files, check_files

from horton.scripts.atomdb import *


def test_iter_mults():
    assert list(iter_mults(4, True)) == [1]
    assert list(iter_mults(4, False)) == [1, 3]
    assert list(iter_mults(42, True)) == [7]
    assert list(iter_mults(42, False)) == [7, 5, 3, 1]
    assert list(iter_mults(47, True)) == [2]
    assert list(iter_mults(47, False)) == [2]


def test_iter_states():
    l = list(iter_states('1-3', 2, 2, True))
    assert l == [(1, -2, 2), (1, -1, 1), (1, 0, 2), (2, -2, 1), (2, -1, 2),
                 (2, 0, 1), (2, 1, 2), (3, -2, 2), (3, -1, 1), (3, 0, 2),
                 (3, 1, 1), (3, 2, 2)]
    l = list(iter_states('1, 6-8', 3, 1, False))
    assert l == [(1, -1, 1), (1, -1, 3), (1, 0, 2), (6, -1, 4), (6, -1, 2),
                 (6, 0, 3), (6, 0, 5), (6, 0, 1), (6, 1, 2), (6, 1, 4),
                 (6, 2, 1), (6, 2, 3), (6, 3, 2), (6, 3, 4), (7, -1, 3),
                 (7, -1, 1), (7, 0, 4), (7, 0, 2), (7, 1, 3), (7, 1, 5),
                 (7, 1, 1), (7, 2, 2), (7, 2, 4), (7, 3, 1), (7, 3, 3),
                 (8, -1, 2), (8, 0, 3), (8, 0, 1), (8, 1, 4), (8, 1, 2),
                 (8, 2, 3), (8, 2, 5), (8, 2, 1), (8, 3, 2), (8, 3, 4)]


def test_script_input_adf():
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_atomdb.test_script_input_adf')
    try:
        fn_template = 'template_atomdb_adf.in'
        copy_files(tmpdir, [fn_template])
        check_script('horton-atomdb.py input adf H,13 %s -l 6' % fn_template, tmpdir)
        dirs = [
            '001__h_001_q+00/mult02', '001__h_002_q-01/mult01',
            '001__h_003_q-02/mult02', '013_al_011_q+02/mult02',
            '013_al_012_q+01/mult01', '013_al_013_q+00/mult02',
            '013_al_014_q-01/mult03', '013_al_015_q-02/mult04',
        ]
        fns = [
            'run_adf.sh', 'settings.txt'
        ] + ['%s/atom.in' % d for d in dirs] + ['%s/grid.in' % d for d in dirs]
        check_files(tmpdir, fns)
    finally:
        shutil.rmtree(tmpdir)


def test_script_input_cp2k():
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_atomdb.test_script_input_cp2k')
    try:
        fn_template = 'template_atomdb_cp2k.in'
        fn_valence = 'include_atomdb_cp2k_valence.inc'
        fn_ppot = 'include_atomdb_cp2k_ppot.inc'
        copy_files(tmpdir, [fn_template, fn_valence, fn_ppot])
        check_script('horton-atomdb.py input cp2k Ca,F %s' % fn_template, tmpdir)
        fns = [
            '020_ca_021_q-01/mult02/atom.in', '020_ca_020_q+00/mult01/atom.in',
            '020_ca_019_q+01/mult02/atom.in', '020_ca_018_q+02/mult01/atom.in',
            '020_ca_017_q+03/mult02/atom.in', '009__f_010_q-01/mult01/atom.in',
            '009__f_009_q+00/mult02/atom.in', '009__f_008_q+01/mult03/atom.in',
            '009__f_007_q+02/mult04/atom.in',
            'run_cp2k.sh', 'settings.txt'
        ]
        check_files(tmpdir, fns)
    finally:
        shutil.rmtree(tmpdir)


def check_script_input_gaussian(binary):
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_atomdb.test_script_input_%s' % binary)
    try:
        fn_template = 'template_atomdb_gaussian.in'
        fn_basis1 = 'include_atomdb_gaussian_basis.001_000_00'
        fn_basis8 = 'include_atomdb_gaussian_basis.008_000_00'
        copy_files(tmpdir, [fn_template, fn_basis1, fn_basis8])
        check_script('horton-atomdb.py input %s 1,O %s' % (binary, fn_template), tmpdir)
        fns = [
            '001__h_003_q-02/mult02/atom.in', '001__h_002_q-01/mult01/atom.in',
            '001__h_001_q+00/mult02/atom.in', '008__o_010_q-02/mult01/atom.in',
            '008__o_009_q-01/mult02/atom.in', '008__o_008_q+00/mult03/atom.in',
            '008__o_007_q+01/mult04/atom.in', '008__o_006_q+02/mult03/atom.in',
            '008__o_005_q+03/mult02/atom.in',
            'run_%s.sh' % binary, 'settings.txt'
        ]
        check_files(tmpdir, fns)
    finally:
        shutil.rmtree(tmpdir)


def test_script_input_g03():
    check_script_input_gaussian('g03')


def test_script_input_g09():
    check_script_input_gaussian('g09')


def test_script_input_orca():
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_atomdb.test_script_input_orca')
    try:
        fn_template = 'template_atomdb_orca.in'
        copy_files(tmpdir, [fn_template])
        check_script('horton-atomdb.py input orca H,13 %s --no-hund' % fn_template, tmpdir)
        fns = [
            '001__h_003_q-02/mult02/atom.in', '001__h_003_q-02/mult04/atom.in',
            '001__h_002_q-01/mult01/atom.in', '001__h_002_q-01/mult03/atom.in',
            '001__h_001_q+00/mult02/atom.in', '013_al_015_q-02/mult04/atom.in',
            '013_al_015_q-02/mult02/atom.in', '013_al_014_q-01/mult03/atom.in',
            '013_al_014_q-01/mult05/atom.in', '013_al_014_q-01/mult01/atom.in',
            '013_al_013_q+00/mult02/atom.in', '013_al_013_q+00/mult04/atom.in',
            '013_al_012_q+01/mult01/atom.in', '013_al_012_q+01/mult03/atom.in',
            '013_al_011_q+02/mult02/atom.in', '013_al_010_q+03/mult01/atom.in',
            'run_orca.sh', 'settings.txt'
        ]
        check_files(tmpdir, fns)
    finally:
        shutil.rmtree(tmpdir)


def copy_atom_output(fn, number, charge, mult, tmpdir, fn_out):
    pop = number - charge
    symbol = periodic[number].symbol.lower().rjust(2, '_')
    destination = os.path.join(
        tmpdir, '%03i_%s_%03i_q%+03i' % (number, symbol, pop, charge),
        'mult%02i' % mult
    )
    os.makedirs(destination)
    destination = os.path.join(destination, fn_out)
    shutil.copy(context.get_fn(os.path.join('test', fn)), destination)


def make_settings(srtf, nll, program, tmpdir):
    with open(os.path.join(tmpdir, 'settings.txt'), 'w') as f:
        print >> f, srtf
        print >> f, nll
        print >> f, program


def test_script_convert_cp2k():
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_atomdb.test_script_convert_cp2k')
    try:
        copy_atom_output('atom_op2.cp2k.out', 8, +2, 3, tmpdir, 'atom.cp2k.out')
        copy_atom_output('atom_op1.cp2k.out', 8, +1, 4, tmpdir, 'atom.cp2k.out')
        copy_atom_output('atom_o.cp2k.out',   8,  0, 2, tmpdir, 'atom.cp2k.out')
        copy_atom_output('atom_om1.cp2k.out', 8, -1, 1, tmpdir, 'atom.cp2k.out')
        copy_atom_output('atom_om2.cp2k.out', 8, -2, 0, tmpdir, 'atom.cp2k.out')
        make_settings('ExpRTransform 0.01 10.0 100', 26, 'cp2k', tmpdir)
        check_script('horton-atomdb.py convert', tmpdir)
        # check presence of files
        fns = ['atoms.h5', 'dens_008__o.png', 'rdens_008__o.png', 'fukui_008__o.png', 'rfukui_008__o.png']
        check_files(tmpdir, fns)
        # load proatomdb file and check some contents
        padb = ProAtomDB.from_file(os.path.join(tmpdir, 'atoms.h5'))
        assert padb.get_numbers() == [8]
        assert padb.get_charges(8) == [+2, +1, 0, -1, -2]
        assert not padb.get_record(8, -2).safe
    finally:
        shutil.rmtree(tmpdir)


def test_script_convert_g09():
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_atomdb.test_script_convert_g09')
    try:
        copy_atom_output('atom_014_013_hf_lan.fchk', 14, +1, 2, tmpdir, 'atom.fchk')
        make_settings('ExpRTransform 0.01 10.0 20', 26, 'g09', tmpdir)
        check_script('horton-atomdb.py convert', tmpdir)
        # check presence of files
        fns = ['atoms.h5', 'dens_014_si.png', 'rdens_014_si.png', 'fukui_014_si.png', 'rfukui_014_si.png']
        check_files(tmpdir, fns)
        # load proatomdb file and check some contents
        padb = ProAtomDB.from_file(os.path.join(tmpdir, 'atoms.h5'))
        assert padb.get_numbers() == [14]
        assert padb.get_charges(14) == [+1]
    finally:
        shutil.rmtree(tmpdir)


def test_script_convert_g03():
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_atomdb.test_script_convert_g03')
    try:
        copy_atom_output('atom_001_001_hf_sto3g.fchk', 1,  0, 2, tmpdir, 'atom.fchk')
        copy_atom_output('atom_001_002_hf_sto3g.fchk', 1, -1, 1, tmpdir, 'atom.fchk')
        copy_atom_output('atom_008_007_hf_sto3g.fchk', 8, +1, 4, tmpdir, 'atom.fchk')
        copy_atom_output('atom_008_008_hf_sto3g.fchk', 8,  0, 3, tmpdir, 'atom.fchk')
        copy_atom_output('atom_008_009_hf_sto3g.fchk', 8, -1, 2, tmpdir, 'atom.fchk')
        make_settings('ExpRTransform 0.01 10.0 20', 26, 'g03', tmpdir)
        check_script('horton-atomdb.py convert', tmpdir)
        # check presence of files
        fns = [
            'atoms.h5',
            'dens_001__h.png', 'rdens_001__h.png', 'fukui_001__h.png', 'rfukui_001__h.png',
            'dens_008__o.png', 'rdens_008__o.png', 'fukui_008__o.png', 'rfukui_008__o.png',
        ]
        check_files(tmpdir, fns)
        # load proatomdb file and check some contents
        padb = ProAtomDB.from_file(os.path.join(tmpdir, 'atoms.h5'))
        assert padb.get_numbers() == [1, 8]
        assert padb.get_charges(1) == [0, -1]
        assert padb.get_charges(8) == [+1, 0, -1]
    finally:
        shutil.rmtree(tmpdir)
