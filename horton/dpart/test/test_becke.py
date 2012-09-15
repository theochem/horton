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
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 1e1, 100)
    grid = BeckeMolGrid(sys, (rtf, 110, 100), random_rotate=False, keep_subgrids=1)
    bdp = BeckeDPart(grid)
    bdp.do_populations()
    assert abs(bdp['populations'] - 7).max() < 1e-4
    bdp.do_charges()
    assert abs(bdp['charges']).max() < 1e-4
    bdp.invalidate()
    try:
        bdp['charges']
        assert False
    except KeyError:
        pass
    bdp.do_charges()
    assert abs(bdp['populations'] - 7).max() < 1e-4
    assert abs(bdp['charges']).max() < 1e-4

def test_becke_nonlocal_lih_hf_321g():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 1e1, 100)

    grid1 = BeckeMolGrid(sys, (rtf, 110, 100), random_rotate=False, keep_subgrids=1)
    bdp1 = BeckeDPart(grid1)

    grid2 = BeckeMolGrid(sys, (rtf, 110, 100), random_rotate=False, keep_subgrids=0)
    bdp2 = BeckeDPart(grid2, local=False)

    bdp1.do_charges()
    bdp2.do_charges()
    assert abs(bdp1['charges'] - bdp2['charges']).max() < 5e-4
