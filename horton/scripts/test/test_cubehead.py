# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2015 The Horton Development Team
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
#pylint: skip-file


from horton.test.common import check_script, tmpdir
from horton.scripts.test.common import check_files, copy_files


def test_script():
    with tmpdir('horton.scripts.test.test_cubehead.test_script') as dn:
        fn_fchk = 'water_sto3g_hf_g03.fchk'
        copy_files(dn, [fn_fchk])
        check_script('horton-cubehead.py %s cubehead.txt' % fn_fchk, dn)
        check_files(dn, ['cubehead.txt'])
