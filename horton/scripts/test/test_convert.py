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

from horton.test.common import check_script, tmpdir
from horton.scripts.test.common import copy_files, check_files

def test_scripts():
    with tmpdir('horton.scripts.test.test_convert.test_scripts') as dn:
        fn_fchk = 'water_sto3g_hf_g03.fchk'
        copy_files(dn, [fn_fchk])
        check_script('horton-convert.py %s test.xyz' % fn_fchk, dn)
