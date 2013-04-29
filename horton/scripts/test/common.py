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


import shutil, os
from horton import context


__all__ = ['copy_files', 'check_files']


def copy_files(tmpdir, fns):
    for fn in fns:
        shutil.copy(context.get_fn(os.path.join('test', fn)), os.path.join(tmpdir, fn))


def check_files(tmpdir, fns):
    for fn in fns:
        assert os.path.isfile(os.path.join(tmpdir, fn))
