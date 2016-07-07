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
# --
"""Unit tests for horton.meanfield.gridgroup."""


from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_gridgroup_density_cutoff():
    # prepare some molecule
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)

    # make and populate a fake cache object (normally done by the effective
    # hamiltonian)
    cache = Cache()
    cache['dm_alpha'] = mol.exp_alpha.to_dm()

    # normal use case
    rgg = RGridGroup(mol.obasis, grid, [RLibXCMGGA('c_tpss')])
    alpha_basics = rgg._update_grid_basics(cache, 'alpha')
    mask1 = alpha_basics[:, 0] >= 1e-9
    mask2 = alpha_basics[:, 0] == 0.0
    assert (mask1 | mask2).all()
    assert (alpha_basics[mask2, :] == 0.0).all()

    # use all grid points
    rgg = RGridGroup(mol.obasis, grid, [RLibXCMGGA('c_tpss')], density_cutoff=0.0)
    alpha_basics = rgg._update_grid_basics(cache, 'alpha')
    assert (alpha_basics[:, 0] >= 0.0).all()
