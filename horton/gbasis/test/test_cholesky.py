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


def get_h2o_obasis():
    fn = context.get_fn('test/water.xyz')
    mol = IOData.from_file(fn)
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    return obasis


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


def test_cholesky_coulomb_array():
    obasis = get_h2o_obasis()
    ref = obasis.compute_electron_repulsion(DenseLinalgFactory(obasis.nbasis))._array
    vecs = compute_cholesky(obasis, GB4ElectronRepulsionIntegralLibInt(obasis.max_shell_type))
    chol = np.einsum('kac,kbd->abcd', vecs, vecs)
    np.testing.assert_allclose(ref, chol, rtol=1e-5, atol=1e-8)


def test_cholesky_coulomb_from_gbasis():
    obasis = get_h2o_obasis()
    ref = obasis.compute_electron_repulsion(DenseLinalgFactory(obasis.nbasis))._array
    vecs = obasis.compute_electron_repulsion(CholeskyLinalgFactory(obasis.nbasis))._array
    chol = np.einsum('kac,kbd->abcd', vecs, vecs)
    np.testing.assert_allclose(ref, chol, rtol=1e-5, atol=1e-8)


def test_cholesky_erf_from_gbasis():
    obasis = get_h2o_obasis()
    mu = 1e4
    ref = obasis.compute_erf_repulsion(DenseLinalgFactory(obasis.nbasis), mu)._array
    vecs = obasis.compute_erf_repulsion(CholeskyLinalgFactory(obasis.nbasis), mu)._array
    chol = np.einsum('kac,kbd->abcd', vecs, vecs)
    np.testing.assert_allclose(ref, chol, rtol=1e-5, atol=1e-8)


def test_cholesky_gauss_from_gbasis():
    obasis = get_h2o_obasis()
    c = 1.2
    alpha = 0.5
    ref = obasis.compute_gauss_repulsion(DenseLinalgFactory(obasis.nbasis), c, alpha)._array
    vecs = obasis.compute_gauss_repulsion(CholeskyLinalgFactory(obasis.nbasis), c, alpha)._array
    chol = np.einsum('kac,kbd->abcd', vecs, vecs)
    np.testing.assert_allclose(ref, chol, rtol=1e-5, atol=1e-8)


def test_cholesky_ralpha_from_gbasis():
    obasis = get_h2o_obasis()
    alpha = -2.0
    ref = obasis.compute_ralpha_repulsion(DenseLinalgFactory(obasis.nbasis), alpha)._array
    vecs = obasis.compute_ralpha_repulsion(CholeskyLinalgFactory(obasis.nbasis), alpha)._array
    chol = np.einsum('kac,kbd->abcd', vecs, vecs)
    np.testing.assert_allclose(ref, chol, rtol=1e-5, atol=1e-8)
