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


import os, shutil
from nose.plugins.attrib import attr

from horton.context import context
from horton.periodic import periodic
from horton.part.proatomdb import ProAtomDB
from horton.test.common import check_script, tmpdir
from horton.scripts.test.common import copy_files, check_files

from horton.scripts.atomdb import iter_mults, iter_states, plot_atoms


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


def test_script_input_cp2k():
    with tmpdir('horton.scripts.test.test_atomdb.test_script_input_cp2k') as dn:
        fn_template = 'template_atomdb_cp2k.in'
        fn_valence = 'include_atomdb_cp2k_valence.inc'
        fn_ppot = 'include_atomdb_cp2k_ppot.inc'
        copy_files(dn, [fn_template, fn_valence, fn_ppot])
        check_script('horton-atomdb.py input cp2k Ca,F %s' % fn_template, dn)
        fns = [
            '020_ca_021_q-01/mult02/atom.in', '020_ca_020_q+00/mult01/atom.in',
            '020_ca_019_q+01/mult02/atom.in', '020_ca_018_q+02/mult01/atom.in',
            '020_ca_017_q+03/mult02/atom.in', '009__f_010_q-01/mult01/atom.in',
            '009__f_009_q+00/mult02/atom.in', '009__f_008_q+01/mult03/atom.in',
            '009__f_007_q+02/mult04/atom.in',
            'run_cp2k.sh',
        ]
        check_files(dn, fns)


def check_script_input_gaussian(binary):
    with tmpdir('horton.scripts.test.test_atomdb.test_script_input_%s' % binary) as dn:
        fn_template = 'template_atomdb_gaussian.in'
        fn_basis1 = 'include_atomdb_gaussian_basis.001_000_00'
        fn_basis8 = 'include_atomdb_gaussian_basis.008_000_00'
        copy_files(dn, [fn_template, fn_basis1, fn_basis8])
        check_script('horton-atomdb.py input %s 1,O %s' % (binary, fn_template), dn)
        fns = [
            '001__h_003_q-02/mult02/atom.in', '001__h_002_q-01/mult01/atom.in',
            '001__h_001_q+00/mult02/atom.in', '008__o_010_q-02/mult01/atom.in',
            '008__o_009_q-01/mult02/atom.in', '008__o_008_q+00/mult03/atom.in',
            '008__o_007_q+01/mult04/atom.in', '008__o_006_q+02/mult03/atom.in',
            '008__o_005_q+03/mult02/atom.in',
            'run_%s.sh' % binary,
        ]
        check_files(dn, fns)


def test_script_input_g03():
    check_script_input_gaussian('g03')


def test_script_input_g09():
    check_script_input_gaussian('g09')


def test_script_input_orca():
    with tmpdir('horton.scripts.test.test_atomdb.test_script_input_orca') as dn:
        fn_template = 'template_atomdb_orca.in'
        copy_files(dn, [fn_template])
        check_script('horton-atomdb.py input orca H,13 %s --no-hund' % fn_template, dn)
        fns = [
            '001__h_003_q-02/mult02/atom.in', '001__h_003_q-02/mult04/atom.in',
            '001__h_002_q-01/mult01/atom.in', '001__h_002_q-01/mult03/atom.in',
            '001__h_001_q+00/mult02/atom.in', '013_al_015_q-02/mult04/atom.in',
            '013_al_015_q-02/mult02/atom.in', '013_al_014_q-01/mult03/atom.in',
            '013_al_014_q-01/mult05/atom.in', '013_al_014_q-01/mult01/atom.in',
            '013_al_013_q+00/mult02/atom.in', '013_al_013_q+00/mult04/atom.in',
            '013_al_012_q+01/mult01/atom.in', '013_al_012_q+01/mult03/atom.in',
            '013_al_011_q+02/mult02/atom.in', '013_al_010_q+03/mult01/atom.in',
            'run_orca.sh',
        ]
        check_files(dn, fns)


def copy_atom_output(fn, number, charge, mult, dn, fn_out):
    pop = number - charge
    symbol = periodic[number].symbol.lower().rjust(2, '_')
    destination = os.path.join(
        dn, '%03i_%s_%03i_q%+03i' % (number, symbol, pop, charge),
        'mult%02i' % mult
    )
    os.makedirs(destination)
    destination = os.path.join(destination, fn_out)
    shutil.copy(context.get_fn(os.path.join('test', fn)), destination)


def make_fake_run_script(program, dn):
    with open(os.path.join(dn, 'run_%s.sh' % program), 'w') as f:
        print >> f, '#!/bin/bash'
        print >> f, 'echo "Foo"'


@attr('slow')
def test_script_convert_cp2k():
    with tmpdir('horton.scripts.test.test_atomdb.test_script_convert_cp2k') as dn:
        copy_atom_output('atom_op2.cp2k.out', 8, +2, 3, dn, 'atom.cp2k.out')
        copy_atom_output('atom_op1.cp2k.out', 8, +1, 4, dn, 'atom.cp2k.out')
        copy_atom_output('atom_o.cp2k.out',   8,  0, 2, dn, 'atom.cp2k.out')
        copy_atom_output('atom_om1.cp2k.out', 8, -1, 1, dn, 'atom.cp2k.out')
        copy_atom_output('atom_om2.cp2k.out', 8, -2, 0, dn, 'atom.cp2k.out')
        make_fake_run_script('cp2k', dn)
        check_script('horton-atomdb.py convert', dn)
        # check presence of files
        fns = ['atoms.h5', 'dens_008__o.png', 'rdens_008__o.png', 'fukui_008__o.png', 'rfukui_008__o.png']
        check_files(dn, fns)
        # load proatomdb file and check some contents
        padb = ProAtomDB.from_file(os.path.join(dn, 'atoms.h5'))
        assert padb.get_numbers() == [8]
        assert padb.get_charges(8) == [+2, +1, 0, -1, -2]
        assert not padb.get_record(8, -2).safe
        assert padb.get_rgrid(8).size == 71


@attr('slow')
def test_script_convert_g09():
    with tmpdir('horton.scripts.test.test_atomdb.test_script_convert_g09') as dn:
        copy_atom_output('atom_014_013_hf_lan.fchk', 14, +1, 2, dn, 'atom.fchk')
        make_fake_run_script('g09', dn)
        check_script('horton-atomdb.py convert --grid medium', dn)
        # check presence of files
        fns = ['atoms.h5', 'dens_014_si.png', 'rdens_014_si.png', 'fukui_014_si.png', 'rfukui_014_si.png']
        check_files(dn, fns)
        # load proatomdb file and check some contents
        padb = ProAtomDB.from_file(os.path.join(dn, 'atoms.h5'))
        assert padb.get_numbers() == [14]
        assert padb.get_charges(14) == [+1]
        assert padb.get_rgrid(14).size == 49


@attr('slow')
def test_script_convert_g03():
    with tmpdir('horton.scripts.test.test_atomdb.test_script_convert_g03') as dn:
        copy_atom_output('atom_001_001_hf_sto3g.fchk', 1,  0, 2, dn, 'atom.fchk')
        copy_atom_output('atom_001_002_hf_sto3g.fchk', 1, -1, 1, dn, 'atom.fchk')
        copy_atom_output('atom_008_007_hf_sto3g.fchk', 8, +1, 4, dn, 'atom.fchk')
        copy_atom_output('atom_008_008_hf_sto3g.fchk', 8,  0, 3, dn, 'atom.fchk')
        copy_atom_output('atom_008_009_hf_sto3g.fchk', 8, -1, 2, dn, 'atom.fchk')
        make_fake_run_script('g03', dn)
        check_script('horton-atomdb.py convert', dn)
        # check presence of files
        fns = [
            'atoms.h5',
            'dens_001__h.png', 'rdens_001__h.png', 'fukui_001__h.png', 'rfukui_001__h.png',
            'dens_008__o.png', 'rdens_008__o.png', 'fukui_008__o.png', 'rfukui_008__o.png',
        ]
        check_files(dn, fns)
        # load proatomdb file and check some contents
        padb = ProAtomDB.from_file(os.path.join(dn, 'atoms.h5'))
        assert padb.get_numbers() == [1, 8]
        assert padb.get_charges(1) == [0, -1]
        assert padb.get_charges(8) == [+1, 0, -1]


def test_plot_atoms():
    padb = ProAtomDB.from_refatoms(numbers=[8, 1], max_cation=1, max_anion=1)
    with tmpdir('horton.scripts.test.test_atomdb.test_plot_atoms') as dn:
        plot_atoms(padb, dn)
        fns = [
            'dens_001__h.png', 'rdens_001__h.png', 'fukui_001__h.png', 'rfukui_001__h.png',
            'dens_008__o.png', 'rdens_008__o.png', 'fukui_008__o.png', 'rfukui_008__o.png',
        ]
        check_files(dn, fns)
