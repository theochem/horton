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
#--
#pylint: skip-file


import os, subprocess, importlib
from glob import glob
from nose.tools import assert_raises

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
    if context.data_dir == os.path.abspath('data/') and os.path.isdir('.git'):
        lines = subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', 'data']).split('\n')
        for line in lines:
            line = line.strip()
            if len(line) != 0:
                raise ValueError('The following file is not checked in: %s' % line)


def test_do_not_use_iodata():
    # Check for the use of the IOData class outside horton.io, horton.scripts
    # or horton.part.proatomdb.

    # find packages
    packages = {'horton': []}
    for fn in glob('horton/*/__init__.py'):
        subpackage = fn.split('/')[1]
        if subpackage in ['test', 'io', 'scripts']:
            continue
        packages['horton.%s' % subpackage] = []
    # find modules
    for package, modules in packages.iteritems():
        stub = package.replace('.', '/')
        for fn in sorted(glob('%s/*.py' % stub) + glob('%s/*.so' % stub)):
            module = fn.split('/')[-1][:-3]
            if module in ['__init__', 'proatomdb']:
                continue
            modules.append(module)
    # try to import IOData
    for package, modules in packages.iteritems():
        for module in modules:
            m = importlib.import_module('%s.%s' % (package, module))
            with assert_raises(AttributeError):
                m.IOData
