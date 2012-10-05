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


import os, subprocess, sys


def run_example(dirname, fn_py):
    # fix python path
    env = dict(os.environ)
    python_path = env.get('PYTHONPATH')
    if python_path is None:
        python_path = os.getcwd()
    else:
        python_path += ':' + os.getcwd()
    env['PYTHONPATH'] = python_path
    env['HORTONDATA'] = os.path.join(os.getcwd(), 'data')

    # prepare Popen arguments
    root = os.path.join("examples", dirname)
    assert os.path.isdir(root)
    assert os.path.isfile(os.path.join(root, fn_py))

    # run example and pass through the output
    p = subprocess.Popen(['./%s' % fn_py], cwd=root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    p.wait()
    sys.stdout.write(p.stdout.read())
    sys.stderr.write(p.stderr.read())

    # final check
    assert p.returncode == 0


def test_example_001_hf_water():
    run_example('001_hf_water', 'run.py')

def test_example_002_hfs_water():
    run_example('002_hfs_water', 'run.py')
