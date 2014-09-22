# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import numpy as np
from horton import *


def get_four_cho_dens(nbasis=10, nvec=8):
    '''Create random cholesky vectors and matching dense four-index object'''
    vecs = []
    for ivec in xrange(nvec):
        vec = np.random.normal(0, 1, (nbasis, nbasis))
        vec = (vec+vec.T)/2
        vecs.append(vec)
    chob = CholeskyFourIndex(nbasis)
    chob._array = np.array(vecs)
    chob._array2 = chob._array
    chob._nvec = nvec
    erb = DenseFourIndex(nbasis)
    erb._array[:] = np.einsum('kac,kbd->abcd', chob._array, chob._array2)
    return chob, erb


def test_cholesky_get_slice():
    chob, erb = get_four_cho_dens()

    indices = "abab->ab"
    indices2 = "aabb-> ba"
    indices3 = "abba->ab"

    indices4 = 'abcc->bac'
    indices5 = 'abcc->abc'
    indices6 = 'abcb->abc'
    indices7 = 'abbc->abc'

    indices4 = 'abcc->bac'
    indices5 = 'abcc->abc'
    indices6 = 'abcb->abc'
    indices7 = 'abbc->abc'

    list_indices = [indices, indices2, indices3, indices4, indices5, indices6,
            indices7]

    for i in list_indices:
        assert np.allclose(erb.get_slice(i),chob.get_slice(i))

def test_cholesky_esum():
    chob, erb = get_four_cho_dens()

    assert np.allclose(erb.esum(), chob.esum())

def test_cholesky_apply_direct():
    chob, erb = get_four_cho_dens()

    A = np.random.random((erb.nbasis, erb.nbasis))
    dm = DenseTwoIndex(A.shape[0])
    dm._array = A

    out = DenseTwoIndex(A.shape[0])
    out2 = DenseTwoIndex(A.shape[0])

    erb.apply_direct(dm, out)
    chob.apply_direct(dm, out2)

    assert np.allclose(out._array, out2._array)

def test_cholesky_apply_exchange():
    chob, erb = get_four_cho_dens()

    A = np.random.random((erb.nbasis,erb.nbasis))
    dm = DenseTwoIndex(A.shape[0])
    dm._array = A

    out = DenseTwoIndex(A.shape[0])
    out2 = DenseTwoIndex(A.shape[0])

    erb.apply_exchange(dm, out)
    chob.apply_exchange(dm, out2)

    assert np.allclose(out._array, out2._array)

def test_cholesky_get_dense():
    chob, erb = get_four_cho_dens()

    assert np.allclose(erb._array, chob._get_dense())

def test_cholesky_four_index_transform_tensordot():
    chob, erb = get_four_cho_dens()

    A = np.random.random((erb.nbasis, erb.nbasis))
    de = DenseExpansion(A.shape[0])
    de._coeffs = A

    A2 = np.random.random((erb.nbasis, erb.nbasis))
    de2 = DenseExpansion(A2.shape[0])
    de2._coeffs = A2

    mo1 = DenseFourIndex(erb.nbasis)
    mo2 = CholeskyFourIndex(erb.nbasis)
    mo2._array = np.zeros_like(chob._array)
    mo2.reset_array2()
    mo2._nvec = chob._array.shape[2]

    mo1.apply_four_index_transform_tensordot(erb, de, de2, de, de2)
    mo2.apply_four_index_transform_tensordot(chob, de, de2, de, de2)

    assert np.allclose(mo1._array, mo2._get_dense())

def test_cholesky_four_index_transform_einsum():
    chob, erb = get_four_cho_dens()

    A = np.random.random((erb.nbasis,erb.nbasis))
    de = DenseExpansion(A.shape[0])
    de._coeffs = A

    A2 = np.random.random((erb.nbasis,erb.nbasis))
    de2 = DenseExpansion(A2.shape[0])
    de2._coeffs = A2

    mo1 = DenseFourIndex(erb.nbasis)
    mo2 = CholeskyFourIndex(erb.nbasis)
    mo2._array = np.zeros_like(chob._array)
    mo2.reset_array2()
    mo2._nvec = chob._array.shape[2]

    mo1.apply_four_index_transform_einsum(erb, de, de2, de, de2)
    mo2.apply_four_index_transform_einsum(chob, de, de2, de, de2)

    assert np.allclose(mo1._array, mo2._get_dense())
