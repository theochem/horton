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


def get_h2o_er(linalg_factory=DenseLinalgFactory):
    fn = context.get_fn('test/water.xyz')
    mol = IOData.from_file(fn)
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    lf = linalg_factory(obasis.nbasis)
    return obasis, obasis.compute_electron_repulsion(lf)._array

def pcholesky4(A, thresh=1e-8):
    for i in A.shape: #assumes square matrix
        assert i == A.shape[0]
    Ls = []
    d = np.inf
    counter=1
    while True:
        d = np.einsum('iijj->ij',A) - sum([i*i for i in Ls],np.zeros(A.shape[0:2]))
        idx_d = np.unravel_index(np.argmax(d),A.shape[0:2])
        print "Iteration " , counter, " selected d: ", d[idx_d]
        if d[idx_d] < thresh:
            print "Condition met. Exiting loop"
            break
        past_L = sum([i*i[idx_d] for i in Ls], np.zeros(A.shape[0:2]))
        Ls.append((d[idx_d]**-0.5)*(A[:, idx_d[0], :,idx_d[1]] - past_L)) #ugly

        counter += 1
        print ""
    return Ls

def test_cholesky_array():
    obasis, ref_er = get_h2o_er()

    vecs = compute_cholesky(obasis)
    test_er = np.einsum('kac,kbd->abcd', vecs, vecs)

    assert np.allclose(ref_er, test_er), abs(ref_er - test_er).max()

def test_cholesky_from_gbasis():
    obasis, ref_er = get_h2o_er()
    obasis2, vecs = get_h2o_er(CholeskyLinalgFactory)

    test_er = np.einsum('kac,kbd->abcd', vecs, vecs)

    assert np.allclose(ref_er, test_er), abs(ref_er - test_er).max()
