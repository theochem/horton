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

import numpy as np
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def get_h2o_er():
    fn = context.get_fn('test/water.xyz')
    mol = IOData.from_file(fn)
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    lf = DenseLinalgFactory(obasis.nbasis)
    return obasis, obasis.compute_electron_repulsion(lf)._array


def test_select_2index():
    obasis, er = get_h2o_er()

    lookups = obasis.shell_lookup
    print lookups, obasis.basis_offsets
    for index0 in np.arange(obasis.nbasis):
        for index2 in np.arange(obasis.nbasis):
            pbegin0, pend0, pbegin2, pend2 = select_2index(obasis, index0, index2)

            shell0 = lookups[index0]
            shell2 = lookups[index2]

            assert pbegin0 <= index0 and index0 < pend0, (pbegin0, index0, pend0)
            assert pbegin2 <= index2 and index2 < pend2, (pbegin2, index2, pend2)

            assert shell0 == lookups[pbegin0] == lookups[pend0-1]
            assert shell2 == lookups[pbegin2] == lookups[pend2-1]

            assert pbegin0 == obasis.basis_offsets[shell0]
            if shell0+1 < obasis.nshell:
                assert pend0 == obasis.basis_offsets[shell0+1], pend0
                print pend0+1
            if shell2+1 < obasis.nshell:
                assert pend2 == obasis.basis_offsets[shell2+1], pend2

def test_compute_diagonal():
    obasis, er = get_h2o_er()

    ref_diag = np.einsum('iijj->ij', er)
    test_diag = np.zeros_like(ref_diag)
    compute_diagonal(obasis, test_diag)

    assert np.allclose(ref_diag, test_diag)

def test_get_2index_slice():
    obasis, er = get_h2o_er()

    for index0 in np.arange(obasis.nbasis):
        for index2 in np.arange(obasis.nbasis):
            ref_slice = er[index0,:,index2,:]
            test_slice = np.zeros_like(ref_slice)
            get_2index_slice(obasis, index0, index2, test_slice)
            assert np.allclose(ref_slice, test_slice), (index0, index2,
                    ref_slice,test_slice)
