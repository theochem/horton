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


import numpy as np, h5py as h5
from nose.tools import assert_raises

from horton import *


def test_linalg_factory_constructors():
    lf = CholeskyLinalgFactory(5)
    assert lf.default_nbasis == 5
    lf = CholeskyLinalgFactory()
    assert lf.default_nbasis is None
    lf.default_nbasis = 10

    # Four-index tests
    op4 = lf.create_four_index(nvec=8)
    assert isinstance(op4, CholeskyFourIndex)
    lf.create_four_index.__check_init_args__(lf, op4, nvec=8)
    assert op4.nbasis == 10
    assert op4.nvec == 8
    assert op4.shape == (10, 10, 10, 10)
    assert not op4.is_decoupled
    assert op4.hermitian_vecs

    op4 = lf.create_four_index(8, 4)
    lf.create_four_index.__check_init_args__(lf, op4, 8, 4)
    assert op4.nbasis == 8
    assert op4.nvec == 4
    assert not op4.is_decoupled
    assert op4.hermitian_vecs

    array = np.random.normal(0, 1, (5, 10, 10))
    op4 = lf.create_four_index(10, array=array, hermitian_vecs=False)
    lf.create_four_index.__check_init_args__(lf, op4, nvec=5, hermitian_vecs=False)
    assert op4._array is array
    assert op4._array2 is array
    assert op4.nbasis == 10
    assert op4.nvec == 5
    assert not op4.is_decoupled
    assert not op4.hermitian_vecs

    array2 = np.random.normal(0, 1, (5, 10, 10))
    op4 = lf.create_four_index(10, array=array, array2=array2, hermitian_vecs=False)
    lf.create_four_index.__check_init_args__(lf, op4, nvec=5, hermitian_vecs=False)
    assert op4._array is array
    assert op4._array2 is array2
    assert op4.nbasis == 10
    assert op4.nvec == 5
    assert op4.is_decoupled
    assert not op4.hermitian_vecs


def test_linalg_objects_del():
    lf = CholeskyLinalgFactory()
    with assert_raises(TypeError):
        lf.create_four_index()


#
# Tests on the CholeskyFourIndex stuff
#


def get_four_cho_dense(nbasis=10, nvec=8, sym=8):
    '''Create random 2-index Cholesky vectors and matching dense four-index object

       **Optional arguments:**

       nbasis
            The number of basis functions

       nvec
            The number of Cholesky vectors

       sym
            The amount of symmetries in the FourIndex object:

            * ``8``: all possible symmtries with respect to index permutations
                     apply.
            * ``4``: symmetric with respect to swapping two indexes of the same
                     electron.
            * ``2``: symmetric with respect to swapping two electrons.
            * ``1``: no symmetries.
    '''
    check_options('sym', sym, 1, 2, 4, 8)
    hermitian_vecs = sym>2
    cho = CholeskyFourIndex(nbasis, nvec, hermitian_vecs=sym>2)
    if sym == 1 or sym == 4:
        cho.decouple_array2()
    cho.randomize()
    dense = DenseFourIndex(nbasis, symmetry=cho.symmetry)
    dense._array[:] = np.einsum('kac,kbd->abcd', cho._array, cho._array2)
    return cho, dense


def test_four_index_hdf5():
    lf = CholeskyLinalgFactory(5)
    a = lf.create_four_index(5, 3)
    a.randomize()
    with h5.File('horton.matrix.test.test_cholesky.test_four_index_hdf5', driver='core', backing_store=False) as f:
        a.to_hdf5(f)
        b = CholeskyFourIndex.from_hdf5(f)
        assert a == b


def test_four_index_copy_new_randomize_clear_assign():
    lf = CholeskyLinalgFactory(5)
    for args in (None, 3), (4, 3):
        a = lf.create_four_index(*args)
        b = a.copy()
        b.randomize()
        assert a != b
        c = b.copy()
        c.new.__check_init_args__(c, b)
        assert b == c
        d = c.new()
        assert a == d
        b.assign(c)
        assert b == c


def test_four_index_iscale():
    lf = CholeskyLinalgFactory()
    op = lf.create_four_index(3, 2)
    op.randomize()
    tmp = op._array.copy()
    op.iscale(3.0)
    assert abs(op._array - 3**0.5*tmp).max() < 1e-10
    assert abs(op._array2 - 3**0.5*tmp).max() < 1e-10
    op.decouple_array2()
    op.randomize()
    tmp = op._array.copy()
    tmp2 = op._array2.copy()
    op.iscale(3.0)
    assert abs(op._array - 3**0.5*tmp).max() < 1e-10
    assert abs(op._array2 - 3**0.5*tmp2).max() < 1e-10


def test_four_index_get():
    cho, dense = get_four_cho_dense(nbasis=3, sym=1)
    for i0 in xrange(dense.shape[0]):
        for i1 in xrange(dense.shape[1]):
            for i2 in xrange(dense.shape[2]):
                for i3 in xrange(dense.shape[3]):
                    assert cho.get_element(i0, i1, i2, i3) == \
                           dense.get_element(i0, i1, i2, i3)


def test_four_index_check_symmetry():
    for sym in 1, 2, 4, 8:
        cho = get_four_cho_dense(sym=sym)[0]
        cho.check_symmetry()
        assert cho.symmetry == sym


def test_four_index_itranspose():
    for sym in 1, 2, 4, 8:
        cho = get_four_cho_dense(sym=sym)[0]
        cho.itranspose()
        cho.check_symmetry()


def check_four_sum(sym):
    cho, dense = get_four_cho_dense(sym=sym)
    assert np.allclose(dense.sum(), cho.sum())


def test_four_sum_1():
    check_four_sum(1)


def test_four_sum_2():
    check_four_sum(2)


def test_four_sum_4():
    check_four_sum(4)


def test_four_sum_8():
    check_four_sum(8)


def check_four_slice_to_X(sym):
    cho, dense = get_four_cho_dense(sym=sym)

    for subscripts in "abab->ab", "aabb->ab", "abba->ab":
        # Return value
        factor = np.random.uniform(1, 2)
        assert np.allclose(dense.slice_to_two(subscripts, factor=factor)._array,
                           cho.slice_to_two(subscripts, factor=factor)._array)
        # Output argument
        dense_out = DenseTwoIndex(dense.nbasis)
        cho_out = DenseTwoIndex(cho.nbasis)
        dense.slice_to_two(subscripts, out=dense_out, factor=factor)
        cho.slice_to_two(subscripts, out=cho_out, factor=factor)
        assert np.allclose(dense_out._array, cho_out._array)
        # Output argument without clear
        factor = np.random.uniform(1, 2)
        dense.slice_to_two(subscripts, out=dense_out, factor=factor, clear=False)
        cho.slice_to_two(subscripts, out=cho_out, factor=factor, clear=False)
        assert np.allclose(dense_out._array, cho_out._array)

    for subscripts in 'abcc->bac', 'abcc->abc', 'abcb->abc', 'abbc->abc':
        # Return value
        factor = np.random.uniform(1, 2)
        assert np.allclose(dense.slice_to_three(subscripts, factor=factor)._array,
                           cho.slice_to_three(subscripts, factor=factor)._array)
        # Output argument
        dense_out = DenseThreeIndex(dense.nbasis)
        cho_out = DenseThreeIndex(cho.nbasis)
        dense.slice_to_three(subscripts, out=dense_out, factor=factor)
        cho.slice_to_three(subscripts, out=cho_out, factor=factor)
        assert np.allclose(dense_out._array, cho_out._array)
        # Output argument without clear
        factor = np.random.uniform(1, 2)
        dense.slice_to_three(subscripts, out=dense_out, factor=factor, clear=False)
        cho.slice_to_three(subscripts, out=cho_out, factor=factor, clear=False)
        assert np.allclose(dense_out._array, cho_out._array)


def test_four_slice_to_X_1():
    check_four_slice_to_X(1)


def test_four_slice_to_X_2():
    check_four_slice_to_X(2)


def test_four_slice_to_X_4():
    check_four_slice_to_X(4)


def test_four_slice_to_X_8():
    check_four_slice_to_X(8)


def check_four_contract_two_to_two_direct(sym):
    cho, dense = get_four_cho_dense(sym=sym)
    dm = DenseTwoIndex(dense.nbasis)
    dm.randomize()
    factor = np.random.uniform(1, 2)
    # check return value
    assert np.allclose(dense.contract_two_to_two('abcd,bd->ac', dm, factor=factor)._array,
                       cho.contract_two_to_two('abcd,bd->ac', dm, factor=factor)._array)
    # check output argument
    out_dense = DenseTwoIndex(dense.nbasis)
    out_cho = DenseTwoIndex(cho.nbasis)
    dense.contract_two_to_two('abcd,bd->ac', dm, out_dense, factor)
    cho.contract_two_to_two('abcd,bd->ac', dm, out_cho, factor)
    assert np.allclose(out_dense._array, out_cho._array)
    # check output argument abd clear=False
    factor = np.random.uniform(1, 2)
    dense.contract_two_to_two('abcd,bd->ac', dm, out_dense, factor, clear=False)
    cho.contract_two_to_two('abcd,bd->ac', dm, out_cho, factor, clear=False)
    assert np.allclose(out_dense._array, out_cho._array)



def test_four_contract_two_to_two_direct_1():
    check_four_contract_two_to_two_direct(1)


def test_four_contract_two_to_two_direct_2():
    check_four_contract_two_to_two_direct(2)


def test_four_contract_two_to_two_direct_4():
    check_four_contract_two_to_two_direct(4)


def test_four_contract_two_to_two_direct_8():
    check_four_contract_two_to_two_direct(8)


def check_four_contract_two_to_two_exchange(sym):
    cho, dense = get_four_cho_dense(sym=sym)
    dm = DenseTwoIndex(dense.nbasis)
    dm.randomize()
    factor = np.random.uniform(1, 2)
    # check return value
    assert np.allclose(dense.contract_two_to_two('abcd,cb->ad', dm, factor=factor)._array,
                       cho.contract_two_to_two('abcd,cb->ad', dm, factor=factor)._array)
    # check output argument
    out_dense = DenseTwoIndex(dense.nbasis)
    out_cho = DenseTwoIndex(cho.nbasis)
    dense.contract_two_to_two('abcd,cb->ad', dm, out_dense, factor)
    cho.contract_two_to_two('abcd,cb->ad', dm, out_cho, factor)
    assert np.allclose(out_dense._array, out_cho._array)
    # check output argument abd clear=False
    factor = np.random.uniform(1, 2)
    dense.contract_two_to_two('abcd,cb->ad', dm, out_dense, factor, clear=False)
    cho.contract_two_to_two('abcd,cb->ad', dm, out_cho, factor, clear=False)
    assert np.allclose(out_dense._array, out_cho._array)


def test_four_contract_two_to_two_exchange_1():
    check_four_contract_two_to_two_exchange(1)


def test_four_contract_two_to_two_exchange_2():
    check_four_contract_two_to_two_exchange(2)


def test_four_contract_two_to_two_exchange_4():
    check_four_contract_two_to_two_exchange(4)


def test_four_contract_two_to_two_exchange_8():
    check_four_contract_two_to_two_exchange(8)


def check_four_index_transform(sym_in, sym_exp, method):
    cho, dense = get_four_cho_dense(sym=sym_in)
    dense_mo = dense.new()
    cho_mo = cho.new()
    if sym_exp == 8:
        exp0 = DenseExpansion(dense.nbasis)
        exp0.randomize()
        dense_mo.assign_four_index_transform(dense, exp0, method=method)
        cho_mo.assign_four_index_transform(cho, exp0, method=method)
        assert cho_mo.is_decoupled == cho.is_decoupled
        assert cho_mo.hermitian_vecs == cho.hermitian_vecs
    elif sym_exp == 4:
        exp0 = DenseExpansion(dense.nbasis)
        exp0.randomize()
        exp1 = DenseExpansion(dense.nbasis)
        exp1.randomize()
        dense_mo.assign_four_index_transform(dense, exp0, exp1, method=method)
        cho_mo.assign_four_index_transform(cho, exp0, exp1, method=method)
        assert cho_mo.is_decoupled
        assert cho_mo.hermitian_vecs == cho.hermitian_vecs
    elif sym_exp == 2:
        exp0 = DenseExpansion(dense.nbasis)
        exp0.randomize()
        exp2 = DenseExpansion(dense.nbasis)
        exp2.randomize()
        dense_mo.assign_four_index_transform(dense, exp0, exp2=exp2, method=method)
        cho_mo.assign_four_index_transform(cho, exp0, exp2=exp2, method=method)
        assert cho_mo.is_decoupled == cho.is_decoupled
        assert not cho_mo.hermitian_vecs
    elif sym_exp == 1:
        exp0 = DenseExpansion(dense.nbasis)
        exp0.randomize()
        exp1 = DenseExpansion(dense.nbasis)
        exp1.randomize()
        exp2 = DenseExpansion(dense.nbasis)
        exp2.randomize()
        exp3 = DenseExpansion(dense.nbasis)
        exp3.randomize()
        dense_mo.assign_four_index_transform(dense, exp0, exp1, exp2, exp3, method=method)
        cho_mo.assign_four_index_transform(cho, exp0, exp1, exp2, exp3, method=method)
        assert cho_mo.is_decoupled
        assert not cho_mo.hermitian_vecs
    else:
        raise ValueError
    assert np.allclose(dense_mo._array, cho_mo.get_dense()._array)
    assert dense_mo.symmetry == cho_mo.symmetry


def test_four_index_transform_8_8_tensordot():
    check_four_index_transform(8, 8, 'tensordot')


def test_four_index_transform_8_8_einsum():
    check_four_index_transform(8, 8, 'einsum')


def test_four_index_transform_4_8_tensordot():
    check_four_index_transform(4, 8, 'tensordot')


def test_four_index_transform_4_8_einsum():
    check_four_index_transform(4, 8, 'einsum')


def test_four_index_transform_2_8_tensordot():
    check_four_index_transform(2, 8, 'tensordot')


def test_four_index_transform_2_8_einsum():
    check_four_index_transform(2, 8, 'einsum')


def test_four_index_transform_1_8_tensordot():
    check_four_index_transform(1, 8, 'tensordot')


def test_four_index_transform_1_8_einsum():
    check_four_index_transform(1, 8, 'einsum')


def test_four_index_transform_8_4_tensordot():
    check_four_index_transform(8, 4, 'tensordot')


def test_four_index_transform_8_4_einsum():
    check_four_index_transform(8, 4, 'einsum')


def test_four_index_transform_4_4_tensordot():
    check_four_index_transform(4, 4, 'tensordot')


def test_four_index_transform_4_4_einsum():
    check_four_index_transform(4, 4, 'einsum')


def test_four_index_transform_2_4_tensordot():
    check_four_index_transform(2, 4, 'tensordot')


def test_four_index_transform_2_4_einsum():
    check_four_index_transform(2, 4, 'einsum')


def test_four_index_transform_1_4_tensordot():
    check_four_index_transform(1, 4, 'tensordot')


def test_four_index_transform_1_4_einsum():
    check_four_index_transform(1, 4, 'einsum')


def test_four_index_transform_8_2_tensordot():
    check_four_index_transform(8, 2, 'tensordot')


def test_four_index_transform_8_2_einsum():
    check_four_index_transform(8, 2, 'einsum')


def test_four_index_transform_4_2_tensordot():
    check_four_index_transform(4, 2, 'tensordot')


def test_four_index_transform_4_2_einsum():
    check_four_index_transform(4, 2, 'einsum')


def test_four_index_transform_2_2_tensordot():
    check_four_index_transform(2, 2, 'tensordot')


def test_four_index_transform_2_2_einsum():
    check_four_index_transform(2, 2, 'einsum')


def test_four_index_transform_1_2_tensordot():
    check_four_index_transform(1, 2, 'tensordot')


def test_four_index_transform_1_2_einsum():
    check_four_index_transform(1, 2, 'einsum')


def test_four_index_transform_8_1_tensordot():
    check_four_index_transform(8, 1, 'tensordot')


def test_four_index_transform_8_1_einsum():
    check_four_index_transform(8, 1, 'einsum')


def test_four_index_transform_4_1_tensordot():
    check_four_index_transform(4, 1, 'tensordot')


def test_four_index_transform_4_1_einsum():
    check_four_index_transform(4, 1, 'einsum')


def test_four_index_transform_2_1_tensordot():
    check_four_index_transform(2, 1, 'tensordot')


def test_four_index_transform_2_1_einsum():
    check_four_index_transform(2, 1, 'einsum')


def test_four_index_transform_1_1_tensordot():
    check_four_index_transform(1, 1, 'tensordot')


def test_four_index_transform_1_1_einsum():
    check_four_index_transform(1, 1, 'einsum')
