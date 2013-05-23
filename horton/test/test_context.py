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


import os, subprocess

from horton import context


def test_context():
    fn = context.get_fn('basis/sto-3g.nwchem')
    assert os.path.isfile(fn)
    fns = context.glob('basis/*.nwchem')
    assert fn in fns


def test_data_files():
    # Find files in data that were not checked in.
    # This test only makes sense if ran inside the source tree. The purpose is
    # to detect mistakes in the development process.
    if context.data_dir == os.path.abspath('data/'):
        lines = subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', 'data']).split('\n')
        for line in lines:
            line = line.strip()
            if len(line) != 0:
                raise ValueError('The following file is not checked in: %s' % line)
