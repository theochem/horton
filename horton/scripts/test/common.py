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


import subprocess, sys, os, shlex


__all__ = ['check_script']


def check_script(command, workdir):
    env = dict(os.environ)
    root_dir = os.getcwd()
    env['PYTHONPATH'] = root_dir + ':' + env['PYTHONPATH']
    env['HORTONDATA'] = os.path.join(root_dir, 'data')
    env['PATH'] = os.path.join(root_dir, 'scripts') + ':' + env['PATH']
    proc = subprocess.Popen(shlex.split(command), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=workdir, env=env)
    outdata, errdata = proc.communicate()
    print 'Standard output'
    print '+'*80
    print outdata
    print '+'*80
    print 'Standard error'
    print '+'*80
    print errdata
    print '+'*80
    assert proc.returncode == 0
