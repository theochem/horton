#!/usr/bin/env python
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


from horton import *


def test_becke_n2_hfs_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialIntGrid(rtf)
    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, keep_subgrids=True)
    bp = BeckeWPart(sys, grid)
    bp.do_populations()
    assert abs(bp['populations'] - 7).max() < 1e-4
    bp.do_charges()
    assert abs(bp['charges']).max() < 1e-4
    bp.invalidate()
    try:
        bp['charges']
        assert False
    except KeyError:
        pass
    bp.do_charges()
    assert abs(bp['populations'] - 7).max() < 1e-4
    assert abs(bp['charges']).max() < 1e-4


def test_becke_nonlocal_lih_hf_321g():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialIntGrid(rtf)

    grid1 = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, keep_subgrids=True)
    bp1 = BeckeWPart(sys, grid1)

    grid2 = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, keep_subgrids=False)
    bp2 = BeckeWPart(sys, grid2, local=False)

    bp1.do_charges()
    bp2.do_charges()
    assert abs(bp1['charges'] - bp2['charges']).max() < 5e-4
