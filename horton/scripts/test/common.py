# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2022 The HORTON Development Team
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


import os
import shutil

from horton import context


__all__ = ['copy_files', 'check_files']


def copy_files(dn, fns):
    for fn in fns:
        shutil.copy(context.get_fn(os.path.join('test', fn)), os.path.join(dn, fn))


def check_files(dn, fns):
    for fn in fns:
        assert os.path.isfile(os.path.join(dn, fn)), "Missing %s" % fn
