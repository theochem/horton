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


def test_base_exceptions():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    int1d = TrapezoidIntegrator1D()
    rtf = LogRTransform(1e-3, 1e1, 100)
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=0)

    try:
        bdp = BaseDPart(grid)
        assert False
    except ValueError:
        pass

    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=1)

    try:
        bdp = BaseDPart(grid)
        assert False
    except NotImplementedError:
        pass

    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=2)

    try:
        bdp = BaseDPart(grid)
        assert False
    except NotImplementedError:
        pass

    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=0)
    try:
        bdp = BaseDPart(grid, local=False)
        assert False
    except NotImplementedError:
        pass
